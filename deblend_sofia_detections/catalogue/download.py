# -*- coding: future_fstrings -*-
#Functions that look for optical image


from deblend_sofia_detections.deblending.image_manipulation import cut_optical
from deblend_sofia_detections.catalogue.download_classes import NullGate, \
    LW_NedQuery, LW_GaiaQuery
from deblend_sofia_detections.support.errors import DownloadError
from deblend_sofia_detections.support.logging import print_log
from deblend_sofia_detections.support.table_functions import check_table_length
from deblend_sofia_detections.support.support_functions import get_ned_requested_metadata,\
    get_fits_header,write_fits_file
from astropy import units as u

from astroquery.skyview import SkyView
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.table import QTable, Column,vstack,unique
from xml.parsers.expat import ExpatError
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from time import sleep
import urllib
import os
import warnings
import numpy as np
import pickle
import random
import copy
import itertools
from threading import Lock, get_ident
from contextlib import contextmanager




_NULL_GATE = NullGate()
_QUERY_RUN_COUNTER = itertools.count(1)
_GATE_COUNTER_LOCK = Lock()
_ACTIVE_GATE_RUNS = 0
_PEAK_GATE_RUNS = 0

@contextmanager
def _tracked_gate(cfg, service_label, run_id, internet_query_gate):
    global _ACTIVE_GATE_RUNS, _PEAK_GATE_RUNS
    thread_id = get_ident()
    print_log(
        cfg,
        f'{service_label} run {run_id}: waiting for internet gate on thread {thread_id}.',
        case=['verbose'])
    with internet_query_gate:
        with _GATE_COUNTER_LOCK:
            _ACTIVE_GATE_RUNS += 1
            if _ACTIVE_GATE_RUNS > _PEAK_GATE_RUNS:
                _PEAK_GATE_RUNS = _ACTIVE_GATE_RUNS
            active_now = _ACTIVE_GATE_RUNS
            peak_now = _PEAK_GATE_RUNS
        print_log(
            cfg,
            f'{service_label} run {run_id}: acquired internet gate on thread {thread_id}. active={active_now}, peak={peak_now}',
            case=['verbose'])
        try:
            yield
        finally:
            with _GATE_COUNTER_LOCK:
                _ACTIVE_GATE_RUNS -= 1
                active_now = _ACTIVE_GATE_RUNS
            print_log(
                cfg,
                f'{service_label} run {run_id}: released internet gate on thread {thread_id}. active={active_now}',
                case=['verbose'])


def _build_progress_bar(completed, total, width=30):
    if total <= 0:
        total = 1
    fraction = completed / total
    filled = int(round(width * fraction))
    filled = max(0, min(width, filled))
    bar = "#" * filled + "-" * (width - filled)
    return f"[{bar}] {100.0 * fraction:5.1f}%"


    
def build_chunk_grid(cfg, sky_coords, size_in_arcmin, max_chunk_size):
    '''Build a grid of sub-coordinates around the given sky coordinates.

    Parameters
    ----------
    sky_coords : SkyCoord
        The central sky coordinates.
    n_chunks : int
        The number of chunks along each axis.
    chunk_size : Quantity
        The size of each chunk.

    Returns
    -------
    subcoords : list
        A list of sub-coordinates and their corresponding chunk size.
    '''

    subcoords = []
    n_chunks = int(np.ceil((size_in_arcmin/max_chunk_size).decompose().value))
    single_chunk_size = size_in_arcmin / n_chunks
    
    for i in range(n_chunks):
        for j in range(n_chunks):
            ra_offset = ((i - n_chunks/2 + 0.5) * single_chunk_size) / np.cos(sky_coords.dec.to(u.rad).value)
            dec_offset = (j - n_chunks/2 + 0.5) * single_chunk_size
            subcoords.append([SkyCoord(
                ra=sky_coords.ra + ra_offset.to(u.deg),
                dec=sky_coords.dec + dec_offset.to(u.deg),
                frame='fk5'),np.sqrt(2.0*(single_chunk_size/2.)**2)])
    return subcoords

def build_source_chunks(cfg, sources, max_chunk_size):
    '''Build source chunks for the given sources.

    Parameters
    ----------
    cfg : object
        The configuration object.
    sources : list
        A list of source coordinates.
    max_chunk_size : Quantity
        The maximum size of each chunk.

    Returns
    -------
    source_chunks : list
        A list of source chunks.
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")     
        mom0_header = get_fits_header(f'{cfg.sofia.directory}/{cfg.sofia.basename}_mom0.fits')
        mom0_wcs = WCS(mom0_header).celestial
    source_chunks = []
   
    for source in sources: 
        #we need center of the source box
        x, y = 0.5*(source['x_max'] + source['x_min']), 0.5*(source['y_max'] + source['y_min'])
        pixel_size = np.nanmax([(source['x_max'] - x).value
            , (x - source['x_min']).value, 
            (source['y_max'] - y).value, 
            (y - source['y_min']).value])*u.pix
        ra,dec = mom0_wcs.wcs_pix2world(x, y, 0)
        ra_box,dec_box = mom0_wcs.wcs_pix2world(x+pixel_size, y+pixel_size, 0)
        radius = SkyCoord(ra=ra_box*u.deg, dec=dec_box*u.deg).separation(SkyCoord(ra=ra*u.deg, dec=dec*u.deg))
        box_coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='fk5')
        if 2*np.sqrt(0.5*radius**2) > max_chunk_size:           
            tmp_chunks = build_chunk_grid(cfg, box_coords, 2*np.sqrt(0.5*radius**2), max_chunk_size)
            source_chunks.extend(tmp_chunks)
        else:
            source_chunks.append([SkyCoord(
                    ra=ra*u.deg,
                    dec=dec*u.deg,
                    frame='fk5'), radius])  
    return source_chunks

def creating_full_FOV_optical(cfg):
    cfFOV_start =  datetime.now() 
    print_log(cfg, f'Starting full FOV optical image creation at {cfFOV_start}', case=['verbose','screen'])
    SkyView.URL = 'https://skyview.gsfc.nasa.gov/current/cgi/basicform.pl'
    #SkyView.URL = 'https://skyview.gsfc.nasa.gov/current/cgi/query.pl'
   
    print_log(cfg, f'Quering the Sky Survey', case=['verbose'])
    cube_ext = cfg.internal.cube_ext
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mom0_header = get_fits_header(f'{cfg.sofia.directory}/{cfg.sofia.basename}_mom0{cube_ext}')
        mom0_wcs = WCS(mom0_header).celestial

    obj_coords, size_quantity, size_pixels,image_boundaries = get_cutout_region(cfg,
        mom0_header=mom0_header,mom0_wcs=mom0_wcs)
    if not cfg.input.manual_optical_image[0] is None:
        
        print_log(cfg, f'Checking manual input', case=['verbose'])
        for identifier in cfg.input.manual_optical_image:
            if os.path.isfile(identifier):
                print_log(cfg, f'Found manual optical image: {identifier}', case=['verbose'])
            manual_path,manual_file = os.path.split(identifier)
            if manual_path == '':
                manual_path = './'
            cutout = cut_optical(cfg,mom0_header,mom0_wcs,\
                manual_path,
                manual_file)

            cutout_hdr = cutout.wcs.to_header()
            cutout_hdr['COMMENT'] =  f'The original file was  {identifier}'
            write_fits_file(cfg.internal.optical_background,cutout.data,cutout_hdr,overwrite=True)

            return
    print_log(cfg, f'''Obtaining the actual image list with the following parameters:
object coordinates: {obj_coords.to_string('hmsdms')},
radius: {size_quantity},
pixels: {size_pixels},''', case=['verbose'])
    
    if cfg.input.clear_internet_cache:
        SkyView.clear_cache()
    sky_view_list = SkyView.get_image_list(position=obj_coords,
                                radius = size_quantity,
                                coordinates = "J2000",
                                pixels = size_pixels,
                                cache=True,
                                #survey=['WISE 3.4'])
                                #survey=['2MASS-K'])
                                survey=['DSS2 Red'])


    print_log(cfg, f'''Obtained {len(sky_view_list)} images from SkyView
{sky_view_list}
''', case=['verbose'])
    for path in sky_view_list:
        print_log(cfg, f'''Starting the download for {cfg.sofia.original_data_cube}.
This can take a while.''', case=['verbose'])
        filename = path.split("/")[-1]
        print_log(cfg, f'Downloading {filename} from SkyView from {path}')
        urllib.request.urlretrieve(path, f'{cfg.internal.optical_background}')
        if os.path.isfile(f'{cfg.internal.optical_background}'):
            print_log(cfg, "Successfully downloaded the image")
            #os.replace(filename, f'{cfg.internal.ancillary_directory}moment0_full_DSS.fits')
        else:
            print_log(cfg, "Failed to obtain the image from SkyView")
            raise DownloadError(f'''Failed to download the image from SkyView: {path}
Check your internet connection and the SkyView service status.
Note that redownloading the exact same image may fail if it has recently been removed from the SkyView archive.
''')
    cfFOV_end = datetime.now()
    print_log(cfg, f'Finished full FOV optical image creation at {cfFOV_end}', case=['verbose','screen'])
    print_log(cfg, f'Total time taken: {cfFOV_end - cfFOV_start}', case=['verbose','screen'])

def gaia_query_with_retries(coord,size):

    query = f'''SELECT ra,dec,phot_rp_mean_mag FROM gaiadr3.gaia_source 
WHERE ra BETWEEN 180 AND 180.2 AND dec BETWEEN 10.1 AND 10.3 ORDER BY 
phot_rp_mean_mag'''
    job = Gaia.launch_job(query,dump_to_file=False)

'''
def download_gaia_table(cfg,runtime_ctx=None):
    dgf_start = datetime.now()
    print_log(cfg, f'Starting Gaia table download at {dgf_start}', case=['verbose','screen'])
    if runtime_ctx is None:
        internet_query_gate = _NULL_GATE
    else:
        internet_query_gate = runtime_ctx.get('internet_query_gate', _NULL_GATE)

    sky_coords, size_quantity, size_pixels,image_boundaries = get_cutout_region(cfg)
    # we are only running this if the user wants dowloads
    if cfg.internal.gaia_table.lower() == 'none':
        from astroquery.gaia import Gaia
        #Load the gaia table
       
        Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
        #Do not set this to minus one as it can crash
        Gaia.ROW_LIMIT = 5000
        # we want maximum chunks of 0.5 degree
        if cfg.input.gaia_credentials[0] != 'NONE' and cfg.input.gaia_credentials[1] != 'NONE':
            Gaia.login(user=cfg.input.gaia_credentials[0], password=cfg.input.gaia_credentials[1])
        
        gaia_table = Gaia.query_object_async(sky_coords, width=size_quantity*1.2, height=size_quantity*1.2)
        print_log(cfg,f"Found {len(gaia_table)} Gaia sources in the image area. Sorting them"
            ,case=['debug'])
        #Remove galaxy canditates
        gaia_table = gaia_table[gaia_table['in_galaxy_candidates'] == False] 
        print_log(cfg,f"After removing galaxy candidates, {len(gaia_table)} Gaia sources remain.",case=['debug'])
        #remove duplicates
      
        gaia_table = unique(gaia_table, keys=['source_id'])
        print_log(cfg,f"After removing duplicates , {len(gaia_table)} Gaia sources remain.",case=['debug'])
        #gaia_table = gaia_table[gaia_table['in_qso_candidates'] == False]    
        #gaia_table = gaia_table[gaia_table['non_single_star'] == 0] 
        gaia_table.sort('phot_rp_mean_mag')
        if len(gaia_table) > 20000:
            print_log(cfg,"Capping the gaia table at 20000. brightest sources",
                case=['debug'])
            gaia_table = gaia_table[0:20000]
        print_log(cfg,f"Found {len(gaia_table)} Gaia sources in the image area after filtering and sorting. Saving to cache."
            ,case=['verbose'])
      
        with open(f'{cfg.directories.ancillary_directory}/tables/cached_gaia_table.pkl','wb') as tmp:
            pickle.dump(gaia_table,tmp) 
        cfg.internal.gaia_table = f'{cfg.directories.ancillary_directory}/tables/cached_gaia_table.pkl'
    dgf_end = datetime.now()
    print_log(cfg, f'Finished Gaia table download at {dgf_end}', case=['verbose','screen'])
    print_log(cfg, f'Total time taken: {dgf_end - dgf_start}', case=['verbose','screen'])
'''
def download_internet_table(cfg, sources=None, runtime_ctx=None, archive= 'NED'):
    '''
    Download the internet table (e.g., NED) for the given sources.

    Parameters:
    cfg (Config): The configuration object.
    sources (list, optional): List of sources to query. Defaults to None.
    runtime_ctx (dict, optional): Runtime context for managing internet queries. Defaults to None.
    archive (str, optional): Type of internet table to download. Defaults to 'NED'. Option for Now NED and 'Simbad'.

    Returns:
    It writes the downloaded internet table to a cached file and updates the configuration accordingly.
    '''
    download_start = datetime.now()
    print_log(cfg, f'Starting {archive} table download at {download_start}', case=['verbose','screen'])

    if archive.upper() not in ['NED', 'SIMBAD','GAIA']:
        print_log(cfg, f"Unsupported archive type: {archive}", case=['verbose','screen'])
        return
    search = False
    if archive.upper() == 'SIMBAD':
        from astroquery.simbad import Simbad
        query_object = Simbad()  
        if cfg.internal.simbad_table.lower() == 'none':
            search = True

            
        query_function = query_simbad_with_retries
        Simbad.TIMEOUT = 600
        max_chunk_size = 150.0 * u.arcmin
        if cfg.input.clear_internet_cache:
            Simbad.clear_cache()
    elif archive.upper() == 'NED':
        query_object = LW_NedQuery()
        max_chunk_size = 10.0 * u.arcmin
        if cfg.internal.ned_table.lower() == 'none':
            search = True
        query_function = query_ned_with_retries
    elif archive.upper() == 'GAIA':
        query_object = LW_GaiaQuery(credentials= cfg.input.gaia_credentials)
        query_object.login(verbose=cfg.logging.verbose_screen)
        max_chunk_size = 10.0 * u.arcmin
        query_function = query_gaia_with_retries
        if cfg.internal.gaia_table.lower() == 'none':
            search = True

    if not search:
        #check that the cached table is not empty
        with open(f'{cfg.directories.ancillary_directory}/tables/cached_{archive.lower()}_table.pkl','rb') as tmp:
            print(f'Loading cached {archive.lower()} table from {cfg.directories.ancillary_directory}/tables/cached_{archive.lower()}_table.pkl')
            cached_table = pickle.load(tmp)
        try:
            if check_table_length(cached_table) == 0:
                search = True
            del cached_table
        except:
             search = True
        


    # we are only running this if the user wants downloads
    if search:
        # let's avoid mix ups if we  download we first remove
        if os.path.exists(f'{cfg.directories.ancillary_directory}/tables/cached_{archive.lower()}_table.pkl'):
            os.remove(f'{cfg.directories.ancillary_directory}/tables/cached_{archive.lower()}_table.pkl')
        if runtime_ctx is None:
            internet_query_gate = _NULL_GATE
        else:
            internet_query_gate = runtime_ctx.get('internet_query_gate', _NULL_GATE)
        sky_coords, size_in_arcmin, size_pixels,image_boundaries = get_cutout_region(cfg)
        
        print_log(cfg,f"Querying {archive} with the new query sources in the image area, this may take some time...",case=['main'])  
      
        internet_table = None
        query_errors = (ExpatError, ValueError, ConnectionError, TimeoutError, OSError)
        internet_table = run_chunked_query(
            cfg, sky_coords, size_in_arcmin, max_chunk_size,
            query_function, (cfg, query_object, query_errors), archive.upper(),
            internet_query_gate=internet_query_gate,sources=sources)
        

        if internet_table is not None:
            print_log(cfg, f'{archive} query completed with {check_table_length(internet_table)} results', 
                case=['verbose'])
        else:
            print_log(cfg, f'{archive} returned an invalid or temporary response; continuing without {archive}  table.',
                case=['verbose','screen'])
            return
        if 'RA' not in internet_table.colnames or 'DEC' not in internet_table.colnames:
            if 'ra' in internet_table.colnames:
                internet_table.rename_column('ra', 'RA')
            if 'dec' in internet_table.colnames:
                internet_table.rename_column('dec', 'DEC')
        # as astropy is the dumbest project ever they can not be consistant so 
        # we have to correct the units for RA and DEg
        if internet_table['RA'].unit is None:
            internet_table['RA'].unit = u.deg
        if internet_table['DEC'].unit is None:
            internet_table['DEC'].unit = u.deg
        # remove duplicates
        if archive.lower() == 'ned':
            filter_keys = [x for x in ['Object Name', 'RA', 'DEC'] if x in internet_table.colnames]
        elif archive.lower() == 'simbad':
            filter_keys = [x for x in ['main id', 'RA', 'DEC'] if x in internet_table.colnames]
        else:
            filter_keys = ['RA', 'DEC']
        
        if len(filter_keys) > 0:
            internet_table = unique(internet_table, keys=filter_keys)   

       
        #remove the ones that are ouside the cutout region
        tmp =copy.deepcopy(internet_table)
       
        internet_table = tmp[(tmp['RA'] >= image_boundaries['ra_min']) 
            & (tmp['RA'] <= image_boundaries['ra_max']) &
            (tmp['DEC'] >= image_boundaries['dec_min']) & 
            (tmp['DEC'] <= image_boundaries['dec_max'])]

        if check_table_length(internet_table) == 0:
            if image_boundaries['ra_min'] < 0.:
                image_boundaries['ra_min'] = image_boundaries['ra_min'] + 360.0*u.deg
                image_boundaries['ra_max'] = image_boundaries['ra_max'] + 360.0*u.deg
                internet_table = tmp[(tmp['RA'] >= image_boundaries['ra_min']) 
                            & (tmp['RA'] <= image_boundaries['ra_max']) &
                            (tmp['DEC'] >= image_boundaries['dec_min']) & 
                            (tmp['DEC'] <= image_boundaries['dec_max'])]
            
        del tmp
        
        # Astropy is so stupid that it does not provide a QTable from the query
        # so we have to do this as well. 
        if archive.lower() not in ['gaia']:
            result_table = QTable()
            translation_table = get_translation_table(archive)
            requested_columns, requested_dtypes, dummy_units = get_ned_requested_metadata()
            for x in requested_columns:
                if translation_table[x] in internet_table.colnames:
                    tmp_column= internet_table[translation_table[x]]
                    tmp_column[tmp_column.mask] = float('NaN')
                
                    result_table[x] = Column(tmp_column,\
                                        unit=internet_table[translation_table[x]].unit,\
                                        dtype=requested_dtypes[requested_columns.index(x)])
                else:
                    result_table[x] = Column([None for x in range(check_table_length(internet_table))],\
                                        dtype=requested_dtypes[requested_columns.index(x)])
            #select out the galaxies 
            #Simbad has a difference with galaxies wit v and without but ned only slects on type
            objects_to_select_with_v = ['Sy1','BLL','LSB','Bla','BiC','SyG','rG', 'bCG', 'Sy2',
                'SBG', 'LIN', 'QSO', 'H2G', 'EmG', 'AGN', 'G', 'GiP','GiG','GiC','CGG',
                'IG','PaG','GrG','ClG','SCG']
            objects_to_select = ['LSB','BiC','SyG','rG', 'bCG',
                'SBG',  'EmG', 'G', 'GiP','GiG','GiC','CGG',
                'IG','PaG','GrG','SCG']+['G','GPAIR','GTRPL','PofG']
            objects_to_select = set([x.upper() for x in objects_to_select])
            objects_to_select_with_v = set([x.upper() for x in objects_to_select_with_v])
            rows=[]
            for x,v in zip(result_table['Type'], result_table['Velocity']):
                if archive.lower() == 'ned':
                    v = float('NaN')
                if x.upper() in objects_to_select_with_v and not np.isnan(v):
                    rows.append(True)
                elif x.upper() in objects_to_select:
                    rows.append(True)
                else:
                    rows.append(False)

            search_table = result_table[rows]
        else:
            search_table = internet_table
        #select out the galaxies
        if check_table_length(search_table) > 0:
            #Do not cahce empty tables.
            print_log(cfg, f'Caching {archive.lower()} table to {cfg.directories.ancillary_directory}/tables/cached_{archive.lower()}_table.pkl', case=['verbose'])          
            with open(f'{cfg.directories.ancillary_directory}/tables/cached_{archive.lower()}_table.pkl','wb') as tmp:
                 pickle.dump(search_table,tmp) 
            setattr(cfg.internal, f'{archive.lower()}_table', f'{cfg.directories.ancillary_directory}/tables/cached_{archive.lower()}_table.pkl')
        else:
            setattr(cfg.internal, f'{archive.lower()}_table', f'none')
            print_log(cfg, f'Not caching empty {archive.lower()} table', case=['verbose'])
            
   
    download_end = datetime.now()
    print_log(cfg, f'Finished {archive} table download at {download_end}', case=['verbose','screen'])
    print_log(cfg, f'Total time taken: {download_end - download_start}', case=['verbose','screen'])

def get_cutout_region(cfg,mom0_header=None,mom0_wcs=None):
      #First we open the moment header 0 to get the extend of the field 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")     
        if mom0_header is None:
            mom0_header = get_fits_header(f'{cfg.sofia.directory}/{cfg.sofia.basename}_mom0.fits')
        if mom0_wcs is None:
            mom0_wcs = WCS(mom0_header).celestial
    #set the size of the image
    size = np.nanmax([abs(mom0_header['NAXIS1']*mom0_header['CDELT1'])*60.,
                      abs(mom0_header['NAXIS2']*mom0_header['CDELT2'])*60.,])
    size_in_arcmin = u.Quantity(size,u.arcmin)
    #get_image_list seems to guess at the pixel size let's fix it to 3 arcsec 
    beam = mom0_header['BMAJ'] * u.deg
    optical_pixel_scale = beam.to(u.arcsec).value/cfg.general.optical_pixel_scale * u.arcsec
    #If the resolution of our optical images is greater than 5 the deblending becomes hairy 
    if optical_pixel_scale > 4.*u.arcsec:
        optical_pixel_scale = 4. * u.arcsec
    size_pixels = (size_in_arcmin.to(u.arcsec).value/optical_pixel_scale.value).astype(int)
    #obtain the central coordinates
    ra,dec = mom0_wcs.wcs_pix2world(mom0_header['NAXIS1']/2., mom0_header['NAXIS2']/2.,1.)
    obj_coords = SkyCoord(ra= ra* u.degree, dec= dec * u.degree, frame='fk5')

    ra_min,dec_min = mom0_wcs.wcs_pix2world(0, 0, 1)
    ra_max,dec_max = mom0_wcs.wcs_pix2world(mom0_header['NAXIS1'], mom0_header['NAXIS2'], 1)

    if ra_min > ra_max:
        ra_min, ra_max = ra_max, ra_min
    if dec_min > dec_max:
        dec_min, dec_max = dec_max, dec_min
    boundaries = {
        'ra_min': ra_min*u.deg,
        'ra_max': ra_max*u.deg,
        'dec_min': dec_min*u.deg,
        'dec_max': dec_max*u.deg,
    }

    return obj_coords, size_in_arcmin, size_pixels,boundaries


'''   
def get_SIMBAD_translation_table():
    translation_table = {}
    translation_table['Object Name'] = 'main_id'
    translation_table['RA'] = 'ra'
    translation_table['DEC'] = 'dec'
    translation_table['Velocity'] = 'rvz_radvel'
    translation_table['Type'] = 'otype_txt'
    translation_table['Magnitude and Filter'] = 'V'
    translation_table['morph_type'] = 'morph_type'
    translation_table['Distance'] = 'Distance'
    return translation_table
'''
def get_translation_table(archive):
    requested_columns = ['Object Name', 'RA', 'DEC', 'Velocity', 'Type', 'Magnitude and Filter', 'morph_type', 'Distance']
    if archive.lower() == 'simbad':
        named_columns = ['main_id', 'RA', 'DEC', 'rvz_radvel', 'otype_txt', 'V', 'morph_type', 'Distance']
    else:
        named_columns = requested_columns
    translation_table = {}
    for req,named in zip(requested_columns, named_columns):
        translation_table[req] = named
    return translation_table




def make_empty_ned_search_table(include_extra=True):
    requested_columns, requested_dtypes, requested_units = \
        get_ned_requested_metadata(include_extra=include_extra)

    table = QTable(names=requested_columns, dtype=requested_dtypes,
        units=requested_units)
    to_add = ['No object Found', np.nan, np.nan, np.nan,
        'Unknown', None, np.nan]
    if include_extra:
        to_add += [np.nan, np.nan, np.nan]
    table.add_row(to_add)
    return table
   
def print_bar(run_type,completed_chunks, total_chunks, chunk_times, started_at):
    elapsed = datetime.now() - started_at
    average_per_chunk = np.mean(chunk_times)
    eta = average_per_chunk * (total_chunks - completed_chunks)
    progress_bar = _build_progress_bar(completed_chunks, total_chunks)
    print(
        f'\r{run_type} {progress_bar} '
        f'chunk {completed_chunks}/{total_chunks} '
        f'chunk_time={chunk_times[-1]} '
        f'average_chunk_time={average_per_chunk} '
        f'elapsed={elapsed} '
        f'ETA={eta}',
        end='', flush=True)
    
def query_gaia_with_retries(cfg, Gaia, query_errors, coords, chunk_id=None):
    table = None
    for attempt in range(3):
        try:       
            table = Gaia.query_region(coords[0], radius=coords[1])
            break
        except query_errors as e:
            print_log(cfg, f'GAIA query failed on attempt {attempt + 1}/3: {e}', case=['verbose'])
            if attempt < 2:
                sleep(5.)
    return table
    
def query_ned_with_retries(cfg, Ned, query_errors, coords, chunk_id=None):
    table = None
    max_attempts = 4
    for attempt in range(max_attempts):
        try:
            table = Ned.query_region(coords[0],
                radius=coords[1], equinox='J2000.0')
            break
        except query_errors as e:
            run_id = None
            if chunk_id is not None and len(chunk_id) >= 3:
                run_id = chunk_id[2]
            run_str = f' run {run_id}' if run_id is not None else ''
            chunk_str = f' on chunk {chunk_id[0]}/{chunk_id[1]}' if chunk_id is not None else ''
            print_log(cfg,
                f'NED query failed on attempt {attempt + 1}/{max_attempts}{run_str}{chunk_str}: {e}',
                case=['verbose'])
            if attempt < (max_attempts - 1):
                # NED occasionally returns transient malformed responses; backoff helps.
                backoff = min(30.0, 1.1 ** np.sqrt(attempt))
                jitter = random.uniform(0.0, 1.75)
                sleep(backoff + jitter)
        except Exception as e:
            print_log(cfg,
                f'NED query failed With the error: {e}',
                case=['verbose'])
            exit(1)
    return table

def query_simbad_with_retries(cfg, Simbad, query_errors, coords, chunk_id=None):
    table = None
    for attempt in range(3):
        try:
            lamp = Simbad()
            lamp.add_votable_fields('main_id', 'ra', 'dec', 'rvz_radvel', 'otype_txt', 'morph_type', 'V')
            table = lamp.query_region(coords[0], radius=coords[1])
            break
        except query_errors as e:
            print_log(cfg, f'SIMBAD query failed on attempt {attempt + 1}/3: {e}', case=['verbose'])
            if attempt < 2:
                sleep(5.)
    return table

def run_chunked_query(cfg, sky_coord, size_in_arcmin, max_chunk_size, function_to_run,
    function_args, run_type, sources = None,internet_query_gate=_NULL_GATE, new=False):
    print_log(cfg,f'Starting {run_type} internet query. With the new query sources.',
        case=['verbose'])

   
    if size_in_arcmin <= max_chunk_size:
        return function_to_run(*function_args, [sky_coord, size_in_arcmin], chunk_id=(1, 1, None))

    
    chunk_coords,total_chunks = set_chunks(cfg, sky_coord, 
        size_in_arcmin, max_chunk_size, sources=sources)
    chunk_tables = []
    failed_subcoords = []
    chunk_times = []
    non_empty_chunk_count = 0
    empty_chunk_count = 0
    failed_chunk_count = 0
    ncpu = max(1, int(getattr(cfg.general, 'ncpu', 1)))
    mean_size= np.mean([x[1].to(u.arcmin).value for x in chunk_coords])*u.arcmin
    print_log(cfg,
        f'Splitting {run_type} query into {total_chunks} chunks with size {mean_size:.2f}',
        case=['verbose'])
    print_log(cfg, f'Using {ncpu} workers for {run_type} chunk queries.',
        case=['verbose'])

  
    started_at = datetime.now()
    completed_chunks = 0

    if ncpu == 1:
        # Run serially to avoid spawning worker threads for services that are not thread-safe.
        for chunk_no, chunk_coord in enumerate(chunk_coords, start=1):
            chunk_table, chunk_time = run_timed_chunk(
                function_to_run, function_args, chunk_coord,
                chunk_id=(chunk_no, total_chunks))
            if chunk_table is None:
                failed_chunk_count += 1
                failed_subcoords.append((chunk_coord, chunk_no))
            else:
                chunk_times.append(chunk_time)
                if check_table_length(chunk_table) > 0:
                    non_empty_chunk_count += 1
                    chunk_tables.append(chunk_table)
                    chunk_times.append(chunk_time)
                else:
                    empty_chunk_count += 1

            completed_chunks += 1
            print_bar(run_type,completed_chunks, total_chunks, chunk_times, started_at)
    else:
        with ThreadPoolExecutor(max_workers=ncpu) as executor:
            future_to_coord = {
                executor.submit(
                    run_timed_chunk, function_to_run, function_args, chunk_coord,
                    chunk_id=(chunk_no, total_chunks)
                ): (chunk_coord, chunk_no)
                for chunk_no, chunk_coord in enumerate(chunk_coords, start=1)
            }
            for future in as_completed(future_to_coord):
                chunk_coord, chunk_no = future_to_coord[future]
                chunk_table, chunk_time = future.result()
                if chunk_table is None:
                    failed_chunk_count += 1
                    failed_subcoords.append((chunk_coord, chunk_no))
                else:
                    if check_table_length(chunk_table) > 0:
                        non_empty_chunk_count += 1
                        chunk_times.append(chunk_time)
                        chunk_tables.append(chunk_table)
                    else:
                        empty_chunk_count += 1

                completed_chunks += 1
                print_bar(run_type, completed_chunks, total_chunks, chunk_times, started_at)
        # flush the bar      
    print()
    print_log(cfg,
        f'''Summary of table retrieval for {run_type}: 
total={total_chunks}, 
non_empty={non_empty_chunk_count}, 
empty={empty_chunk_count}, 
failed={failed_chunk_count}''',case=['verbose'])

    if failed_chunk_count > 0 and len(failed_subcoords) > 0:
        print_log(cfg,
            f'{service_label}: retrying {failed_chunk_count} failed chunks serially for recovery.',
            case=['verbose','screen'])
        recovered = 0
        still_failed = []
        for sub_coord, chunk_no in failed_subcoords:
            sub_table, _ = run_timed_chunk(
                function_to_run, function_args, sub_coord, 
                chunk_id=(chunk_no, total_chunks))
            if sub_table is None:
                still_failed.append((sub_coord, chunk_no))
            else:
                recovered += 1
                if check_table_length(sub_table) > 0:
                    non_empty_chunk_count += 1
                    chunk_tables.append(sub_table)
                else:
                    empty_chunk_count += 1

        failed_chunk_count = len(still_failed)
        print_log(cfg,
        f'{run_type}: serial recovery recovered {recovered} chunk(s); remaining failed={failed_chunk_count}.',
            case=['verbose'])

    if failed_chunk_count > 0:
        raise RuntimeError(
            f'{run_type} chunked query failed: {failed_chunk_count}/{total_chunks} chunks failed.')

    if len(chunk_tables) == 0:
        return None
    if len(chunk_tables) == 1:
        internet_table = chunk_tables[0]
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            internet_table = vstack(chunk_tables)

    print_log(cfg,
        f'{run_type} chunked query completed with {check_table_length(internet_table)} unfiltered results',
        case=['verbose'])
    return internet_table

def run_timed_chunk(function_to_run, function_args, sub_coords, chunk_id=(1,1,None)):
    chunk_started_at = datetime.now()
    chunk_table = function_to_run(*function_args, sub_coords, chunk_id=chunk_id)
    return chunk_table, (datetime.now() - chunk_started_at)

def set_chunks(cfg, sky_coord, size_in_arcmin,max_chunk_size, sources=None):
    
    if sources is None:
      
        chunk_coords = build_chunk_grid(cfg,sky_coord, size_in_arcmin, max_chunk_size)
       
    else:
        chunk_coords = build_source_chunks(cfg, sources,max_chunk_size)
    total_chunks = len(chunk_coords)
   
    return chunk_coords, total_chunks