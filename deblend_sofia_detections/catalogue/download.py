# -*- coding: future_fstrings -*-
#Functions that look for optical image


from deblend_sofia_detections.deblending.image_manipulation import cut_optical
from deblend_sofia_detections.support.errors import DownloadError
from deblend_sofia_detections.support.logging import print_log
from deblend_sofia_detections.support.table_functions import check_table_length
from deblend_sofia_detections.support.support_functions import get_ned_requested_metadata
from astropy import time, units as u
from astropy.io import fits
from astroquery.skyview import SkyView
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.table import QTable, Column,vstack,unique
from xml.parsers.expat import ExpatError
from concurrent.futures import ThreadPoolExecutor, as_completed

import urllib.request
import os
import warnings
import numpy as np
import pickle
import time


def _format_duration(seconds):
    seconds = int(max(0, round(seconds)))
    mins, secs = divmod(seconds, 60)
    hours, mins = divmod(mins, 60)
    if hours > 0:
        return f"{hours:02d}:{mins:02d}:{secs:02d}"
    return f"{mins:02d}:{secs:02d}"


def _build_progress_bar(completed, total, width=30):
    if total <= 0:
        total = 1
    fraction = completed / total
    filled = int(round(width * fraction))
    filled = max(0, min(width, filled))
    bar = "#" * filled + "-" * (width - filled)
    return f"[{bar}] {100.0 * fraction:5.1f}%"


def _build_chunk_subcoords(sky_coords, n_chunks, chunk_size):
    subcoords = []
    for i in range(n_chunks):
        for j in range(n_chunks):
            ra_offset = ((i - n_chunks/2 + 0.5) * chunk_size) / np.cos(sky_coords.dec.to(u.rad).value)
            dec_offset = (j - n_chunks/2 + 0.5) * chunk_size
            subcoords.append(SkyCoord(
                ra=sky_coords.ra + ra_offset.to(u.deg),
                dec=sky_coords.dec + dec_offset.to(u.deg),
                frame='fk5'))
    return subcoords


def _query_ned_with_retries(cfg, Ned, query_errors, coords, query_size):
    table = None
    for attempt in range(3):
        try:
            table = Ned.query_region(coords,
                radius=np.sqrt(2.0*(query_size/2.)**2), equinox='J2000.0')
            break
        except query_errors as e:
            print_log(cfg, f'NED query failed on attempt {attempt + 1}/3: {e}', case=['verbose','screen'])
            if attempt < 2:
                time.sleep(5.)
    return table


def _query_ned_chunk_timed(cfg, Ned, query_errors, sub_coords, chunk_size):
    chunk_started_at = time.time()
    sub_table = _query_ned_with_retries(cfg, Ned, query_errors, sub_coords, chunk_size)
    return sub_table, (time.time() - chunk_started_at)


def _query_simbad_with_retries(cfg, Simbad, query_errors, coords, query_size):
    table = None
    for attempt in range(3):
        try:
            lamp = Simbad()
            lamp.add_votable_fields('main_id', 'ra', 'dec', 'rvz_radvel', 'otype_txt', 'morph_type', 'V')
            table = lamp.query_region(coords, radius=np.sqrt(2*(query_size/2.)**2))
            break
        except query_errors as e:
            print_log(cfg, f'SIMBAD query failed on attempt {attempt + 1}/3: {e}', case=['verbose','screen'])
            if attempt < 2:
                time.sleep(5.)
    return table


def _query_simbad_chunk_timed(cfg, Simbad, query_errors, sub_coords, chunk_size):
    chunk_started_at = time.time()
    sub_table = _query_simbad_with_retries(cfg, Simbad, query_errors, sub_coords, chunk_size)
    return sub_table, (time.time() - chunk_started_at)


def _query_chunk_timed(worker_fn, worker_args, sub_coords, chunk_size):
    chunk_started_at = time.time()
    sub_table = worker_fn(*worker_args, sub_coords, chunk_size)
    return sub_table, (time.time() - chunk_started_at)


def _run_chunked_query(cfg, sky_coords, size_quantity, max_chunk_size, worker_fn,
    worker_args, service_label):
    if size_quantity <= max_chunk_size:
        return worker_fn(*worker_args, sky_coords, size_quantity)

    n_chunks = int(np.ceil((size_quantity/max_chunk_size).decompose().value))
    chunk_size = size_quantity / n_chunks
    chunk_tables = []
    total_chunks = n_chunks**2
    non_empty_chunk_count = 0
    empty_chunk_count = 0
    failed_chunk_count = 0
    ncpu = max(1, int(getattr(cfg.general, 'ncpu', 1)))

    print_log(cfg,
        f'Splitting {service_label} query into {total_chunks} chunks with size {chunk_size.to(u.arcmin):.2f}',
        case=['verbose','screen'])
    print_log(cfg, f'Using {ncpu} workers for {service_label} chunk queries.', case=['verbose','screen'])

    subcoords = _build_chunk_subcoords(sky_coords, n_chunks, chunk_size)
    started_at = time.time()
    completed_chunks = 0

    with ThreadPoolExecutor(max_workers=ncpu) as executor:
        futures = [
            executor.submit(_query_chunk_timed, worker_fn, worker_args, sub_coord, chunk_size)
            for sub_coord in subcoords
        ]
        for future in as_completed(futures):
            sub_table, chunk_elapsed = future.result()
            if sub_table is None:
                failed_chunk_count += 1
            else:
                if check_table_length(sub_table) > 0:
                    non_empty_chunk_count += 1
                    chunk_tables.append(sub_table)
                else:
                    empty_chunk_count += 1

            completed_chunks += 1
            elapsed = time.time() - started_at
            average_per_chunk = elapsed / completed_chunks
            eta = average_per_chunk * (total_chunks - completed_chunks)
            progress_bar = _build_progress_bar(completed_chunks, total_chunks)
            print(
                f'\r{service_label} {progress_bar} '
                f'chunk {completed_chunks}/{total_chunks} '
                f'chunk_time={chunk_elapsed:5.1f}s '
                f'elapsed={_format_duration(elapsed)} '
                f'ETA={_format_duration(eta)}',
                end='', flush=True)
    print()

    print_log(cfg,
        f'{service_label} chunk summary: total={total_chunks}, non_empty={non_empty_chunk_count}, empty={empty_chunk_count}, failed={failed_chunk_count}',
        case=['verbose','screen'])

    if failed_chunk_count > 0:
        raise RuntimeError(
            f'{service_label} chunked query failed: {failed_chunk_count}/{total_chunks} chunks failed.')

    if len(chunk_tables) == 0:
        return None
    if len(chunk_tables) == 1:
        internet_table = chunk_tables[0]
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            internet_table = vstack(chunk_tables)

    print_log(cfg,
        f'{service_label} chunked query completed with {check_table_length(internet_table)} unfiltered results',
        case=['screen'])
    return internet_table

def creating_full_FOV_optical(cfg):
    SkyView.URL = 'https://skyview.gsfc.nasa.gov/current/cgi/basicform.pl'
    #SkyView.URL = 'https://skyview.gsfc.nasa.gov/current/cgi/query.pl'
   
    print_log(cfg, f'Quering the Sky Survey', case=['verbose'])
    cube_ext = cfg.internal.cube_ext
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mom0_header = fits.getheader(f'{cfg.sofia.directory}/{cfg.sofia.basename}_mom0{cube_ext}')
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
            fits.writeto(cfg.internal.optical_background,cutout.data,cutout_hdr,overwrite=True)

            return
    print_log(cfg, f'''Obtaining the actual image list with the following parameters:
object coordinates: {obj_coords.to_string('hmsdms')},
radius: {size_quantity},
pixels: {size_pixels},''', case=['verbose'])
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



def download_gaia_table(cfg):
    sky_coords, size_quantity, size_pixels,image_boundaries = get_cutout_region(cfg)
    # we are only running this if the user wants dowloads
    if cfg.internal.gaia_table.lower() == 'none':
        from astroquery.gaia import Gaia
        #Load the gaia table
       
        Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
        #Do not set this to minus one as it can crash
        Gaia.ROW_LIMIT = 50000
        # we want maximum chunks of 0.5 degree
        '''
        number_of_chunk = int(np.ceil(size_quantity.to(u.deg).value/1.0))
        total_requests = number_of_chunk**2
        print_log(cfg,f"Querying Gaia for sources in the image area, this may take some time. The query will be split into {total_requests} requests to avoid timeouts."
            ,case=['main','screen'])
        new_size = size_quantity.to(u.deg)/number_of_chunk
        #new_size = 0.1*u.deg
       
        for i in range(number_of_chunk):
            for j in range(number_of_chunk):
                ra_offset = ((i - number_of_chunk/2 + 0.5) * new_size)*np.cos(sky_coords.dec.to(u.rad).value)
                dec_offset = (j - number_of_chunk/2 + 0.5) * new_size
                print_log(cfg,f"Calcuated a ra_offset of {ra_offset} and dec_offset of {dec_offset}",case=['debug'])
                sub_coords = SkyCoord(ra=sky_coords.ra + ra_offset, dec=sky_coords.dec + dec_offset, frame='fk5')
                print_log(cfg,f"Querying Gaia for sources in the sub-region centered at {sub_coords.to_string('hmsdms')} with size {new_size} degrees.",case=['verbose','screen'])
                gaia_table_sub = Gaia.query_object_async(sub_coords, width=new_size*1.5, height=new_size*1.5)
                if i == 0 and j == 0:
                    gaia_table = gaia_table_sub
                else:
                    gaia_table = vstack([gaia_table, gaia_table_sub])
                print_log(cfg,f"Found {len(gaia_table_sub)} Gaia sources in the sub-region. (total so far: {len(gaia_table)})",
                          case=['verbose','screen'])
        '''            
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

def download_ned_table(cfg):
    
    # we are only running this if the user wants dowloads
    if cfg.internal.ned_table.lower() == 'none':
        sky_coords, size_quantity, size_pixels,image_boundaries = get_cutout_region(cfg)
        print_log(cfg,f"Querying NED for sources in the image area, this may take some time...",case=['screen'])  
        from astroquery.ipac.ned import Ned
        Ned.TIMEOUT = 3600
        Ned.clear_cache()
        max_chunk_size = 10.0 * u.arcmin
        internet_table = None
        query_errors = (ExpatError, ValueError, ConnectionError, TimeoutError, OSError)
        internet_table = _run_chunked_query(
            cfg, sky_coords, size_quantity, max_chunk_size,
            _query_ned_with_retries, (cfg, Ned, query_errors), 'NED')

        if internet_table is not None:
            print_log(cfg, f'NED query completed with {check_table_length(internet_table)} results', case=['screen'])
        if internet_table is None:
            print_log(cfg, 'NED returned an invalid or temporary response; continuing without NED counterpart table.',
                case=['verbose'])
            return
       
        # as astropy is the dumbest project ever they can not be consistant so 
        # we have to correct the units for RA and DEg
        internet_table['RA'].unit=u.deg
        internet_table['DEC'].unit=u.deg
        # remove duplicates
        dedupe_keys = [x for x in ['Object Name', 'RA', 'DEC'] if x in internet_table.colnames]
        if len(dedupe_keys) > 0:
            internet_table = unique(internet_table, keys=dedupe_keys)   

       
        #remove the ones that are ouside the cutout region
        internet_table = internet_table[(internet_table['RA'] >= image_boundaries['ra_min']) 
            & (internet_table['RA'] <= image_boundaries['ra_max']) &
            (internet_table['DEC'] >= image_boundaries['dec_min']) & 
            (internet_table['DEC'] <= image_boundaries['dec_max'])]
     
        # Astropy is so stupid that it does not provide a QTable from the query
        # so we have to do this as well. 
        result_table = QTable()
        requested_columns, requested_dtypes, dummy_units = get_ned_requested_metadata()
        for x in requested_columns:
            if x in internet_table.colnames:
                tmp_column= internet_table[x]
                tmp_column[tmp_column.mask] = float('NaN')
            
                result_table[x] = Column(tmp_column,\
                                    unit=internet_table[x].unit,\
                                    dtype=requested_dtypes[requested_columns.index(x)])
            else:
                result_table[x] = Column([None for x in range(check_table_length(internet_table))],\
                                    dtype=requested_dtypes[requested_columns.index(x)])
            
        #select out the galaxies
        objects_to_select = ['G','GPAIR','GTRPL']
        rows = [True if x.upper() in objects_to_select else False for x in result_table['Type']]
        search_table = result_table[rows]
        with open(f'{cfg.directories.ancillary_directory}/tables/cached_ned_table.pkl','wb') as tmp:
            pickle.dump(search_table,tmp) 
        cfg.internal.ned_table = f'{cfg.directories.ancillary_directory}/tables/cached_ned_table.pkl'   

def download_simbad_table(cfg):
    # we are only running this if the user wants dowloads
    if cfg.internal.simbad_table.lower() == 'none':
        sky_coords, size_quantity, size_pixels,image_boundaries = get_cutout_region(cfg)
        print_log(cfg,f"Querying Simbad for sources in the image area, this may take some time...",case=['screen'])  
        from astroquery.simbad import Simbad
        Simbad.TIMEOUT = 3600
        Simbad.clear_cache()
        internet_table = None
        max_chunk_size = 150.0 * u.arcmin
        query_errors = (ExpatError, ValueError, ConnectionError, TimeoutError, OSError)
        internet_table = _run_chunked_query(
            cfg, sky_coords, size_quantity, max_chunk_size,
            _query_simbad_with_retries, (cfg, Simbad, query_errors), 'SIMBAD')

          
        print_log(cfg, f'SIMBAD query completed with originally {check_table_length(internet_table)} results', case=['screen'])  
        if internet_table is None:
            print_log(cfg, 'SIMBAD returned an invalid or temporary response; continuing without a SIMBAD counterpart.',
                case=['verbose'])
            return
        dedupe_keys = [x for x in ['main_id', 'ra', 'dec'] if x in internet_table.colnames]
        if len(dedupe_keys) > 0:
            internet_table = unique(internet_table, keys=dedupe_keys)
        
        # as astropy is the dumbest project ever they can not be consistant so 
        # we have to correct the units for RA and DEg
        internet_table['ra'].unit=u.deg
        internet_table['dec'].unit=u.deg
       
        #remove the ones that are ouside the cutout region
        internet_table = internet_table[(internet_table['ra'] >= image_boundaries['ra_min']) 
            & (internet_table['ra'] <= image_boundaries['ra_max']) &
            (internet_table['dec'] >= image_boundaries['dec_min']) & 
            (internet_table['dec'] <= image_boundaries['dec_max'])]
        # Astropy is so stupid that it does not provide a QTable from the query
        # so we have to do this as well. 
        result_table = QTable()
        translation_table = get_SIMBAD_translation_table()
        requested_columns, requested_dtypes, dummy_units = get_ned_requested_metadata()
  
        for x in requested_columns:
            if translation_table[x] in internet_table.colnames:
                tmp_column= internet_table[translation_table[x]]
                tmp_column[tmp_column.mask] = float('NaN')
                result_table[x] = Column(tmp_column,\
                                    unit=internet_table[translation_table[x]].unit,
                                    dtype=requested_dtypes[requested_columns.index(x)])
            
            else:
                result_table[x] = Column([None for x in range(check_table_length(internet_table))],\
                                    dtype=requested_dtypes[requested_columns.index(x)])
        #select out the galaxies
        objects_to_select_with_v = ['Sy1','BLL','LSB','Bla','BiC','SyG','rG', 'bCG', 'Sy2',
            'SBG', 'LIN', 'QSO', 'H2G', 'EmG', 'AGN', 'G', 'GiP','GiG','GiC','CGG',
            'IG','PaG','GrG','ClG','SCG']
        objects_to_select = ['LSB','BiC','SyG','rG', 'bCG',
            'SBG',  'EmG', 'G', 'GiP','GiG','GiC','CGG',
            'IG','PaG','GrG','SCG']
        objects_to_select = [x.upper() for x in objects_to_select]
        objects_to_select_with_v = [x.upper() for x in objects_to_select_with_v]
        rows=[]
        for x,v in zip(result_table['Type'], result_table['Velocity']):
            if x.upper() in objects_to_select_with_v and not np.isnan(v):
                rows.append(True)
            elif x.upper() in objects_to_select:
                rows.append(True)
            else:
                rows.append(False)

        search_table = result_table[rows]
        print_log(cfg, f'SIMBAD query completed with {check_table_length(search_table)} results after filtering', case=['screen'])
        with open(f'{cfg.directories.ancillary_directory}/tables/cached_simbad_table.pkl','wb') as tmp:
            pickle.dump(search_table,tmp) 
        cfg.internal.simbad_table = f'{cfg.directories.ancillary_directory}/tables/cached_simbad_table.pkl'

def get_cutout_region(cfg,mom0_header=None,mom0_wcs=None):
      #First we open the moment header 0 to get the extend of the field 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")     
        if mom0_header is None:
            mom0_header = fits.getheader(f'{cfg.sofia.directory}/{cfg.sofia.basename}_mom0.fits')
        if mom0_wcs is None:
            mom0_wcs = WCS(mom0_header).celestial
    #set the size of the image
    size = np.nanmax([abs(mom0_header['NAXIS1']*mom0_header['CDELT1'])*60.,
                      abs(mom0_header['NAXIS2']*mom0_header['CDELT2'])*60.,])
    size_quantity= u.Quantity(size,u.arcmin)
    #get_image_list seems to guess at the pixel size let's fix it to 3 arcsec 
    beam = mom0_header['BMAJ'] * u.deg
    optical_pixel_scale = beam.to(u.arcsec).value/cfg.general.optical_pixel_scale * u.arcsec
    #If the resolution of our optical images is greater than 5 the deblending becomes hairy 
    if optical_pixel_scale > 4.*u.arcsec:
        optical_pixel_scale = 4. * u.arcsec
    size_pixels = (size_quantity.to(u.arcsec).value/optical_pixel_scale.value).astype(int)
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

    return obj_coords, size_quantity, size_pixels,boundaries


   
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