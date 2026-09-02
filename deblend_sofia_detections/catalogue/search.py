

from deblend_sofia_detections.support.errors import InputError
from deblend_sofia_detections.support.table_functions import check_table_length,\
    read_manual_table, combine_tables, copy_table_header
from deblend_sofia_detections.support.support_functions import convertRADEC,\
    isquantity, get_nan_for_dtype,calculate_projected_distance,close_variables,\
    get_channel_width,get_ned_requested_metadata,open_fits_file
from deblend_sofia_detections.support.logging import print_log
from deblend_sofia_detections.support.constants import C,rest_HI
from deblend_sofia_detections.support.system_functions import join_path

from astropy.table import QTable,Table,vstack
from datetime import datetime
import astropy.units as u
import copy
import numpy as np
import warnings

import pickle


def create_source_table(source,cfg=None,basename=None,sofia_directory='./'):
    print_log(cfg,f"Processing deblended source with id {source['id'][0]} and name {source['name'][0]} "
                    ,case=['verbose'])
    #First we check if we have previous iteration output
    source = search_counter_part(cfg,source,basename=basename,
        query = 'INTERNET',sofia_directory=sofia_directory,
        insource='sofia')
  
    source = search_counter_part(cfg,source,basename=basename,
        query='Manual',sofia_directory=sofia_directory)
    
    return source    
  
    

def find_counterpart(cfg,source,header_info, sysrange=None, table_source='NED',
                     spectroscopic=True):
    
    if table_source.upper() == 'NED':
        if cfg.internal.ned_table == 'none':
            print_log(cfg, f'No cached NED table found. Please run the download_catalogue function first.', case=['verbose'])
            return QTable()
        with open(f'{cfg.internal.ned_table}', 'rb') as tmp:
            search_table = pickle.load(tmp)
    elif table_source.upper() == 'SIMBAD':
        if cfg.internal.simbad_table == 'none':
            print_log(cfg, f'No cached SIMBAD table found. Please run the download_catalogue function first.', case=['verbose'])
            return QTable()
        with open(f'{cfg.internal.simbad_table}', 'rb') as tmp:
            search_table = pickle.load(tmp)
    elif table_source.upper() == 'MANUAL':
        if cfg.input.manual_input_tables[0] is None:
            return  QTable()
        search_table = read_manual_table(cfg)
       
    else:
        raise InputError(f'Unknown table source {table_source}. Please use NED, SIMBAD or MANUAL.')
    prefix = ''
    for col in source.colnames:
        if col.split('_')[0] == 'sofia':
            prefix = 'sofia_'
            break
            
    coordinates = [source[f'{prefix}ra'],source[f'{prefix}dec']]
    vsys, radius = set_search_radius(cfg,source,header_info,sysrange,
        counterpart_region=cfg.input.counterpart_region)
    
    search_table = sort_on_distance(cfg,search_table, coordinates,
        vsys,header_info=header_info,spectroscopic=spectroscopic)
    print_log(cfg,f'''Searching for {table_source} counterpart for {source[f'{prefix}id']}
Search a radius {radius.to(u.arcsec)} around {", ".join(convertRADEC(cfg,coordinates[0].value,coordinates[1].value))}
The nearest target is {search_table['Object Name'][0]} at a distance of {search_table["Spatial Diff"][0].to(u.arcsec)}
And the velocity difference is {search_table["Velocity Diff"][0].to(u.km/u.s)} to vsys {vsys.to(u.km/u.s)}''',
        case=['verbose'])
    #the closest 200 source should suffice
    search_table = search_table[0:200]
    #let's mask all that are outside the range if we have set a range of velocities
 
    for s_diff,v_diff in zip(search_table['Spatial Diff'],search_table['Velocity Diff']):
        if (np.isnan(v_diff)  or v_diff > sysrange) and spectroscopic:
            search_table.remove_row(np.where(search_table['Spatial Diff'] == s_diff)[0][0]) 
            continue
        #If we want to change this we should change the radius
        if s_diff > radius.to(u.arcsec):
            search_table.remove_row(np.where(search_table['Spatial Diff'] == s_diff)[0][0]) 
   
    # if all are outside the range we return an empty table
    if check_table_length(search_table) == 0: 
       search_table = QTable()
    elif np.isnan(search_table['Velocity Diff'][0]):
        if search_table['Spatial Diff'][0] > 2.*radius.to(u.arcsec):
            search_table = QTable() 
    elif search_table['Spatial Diff'][0] > 2.*radius.to(u.arcsec) or\
       search_table['Velocity Diff'][0] > sysrange:
       search_table = QTable()
   
    return search_table





def search_counter_part(cfg,source,sofia_directory= './',
        basename=None,query ='INTERNET',insource=None):
    '''Look for the optical counterpart of the source'''
    if isinstance(source, (Table, QTable)):
        source = source[0]
    try:
        inid = source['id']
    except:
        inid = source['sofia_id']
   
    print_log(cfg,f'Searching in {query} to find a counterpart for {basename} with id {inid}',
        case=['verbose'])
   
    input_dir = f'{sofia_directory}/{basename}_cubelets'
    
    cube = open_fits_file(f'{input_dir}/{basename}_{inid}_cube.fits',\
        output_verify='warn')
    header_info= {'BMAJ':float(cube[0].header['BMAJ'])*u.deg,
                'pixelsize': float(np.mean([abs(cube[0].header['CDELT1']),
                    abs(cube[0].header['CDELT2'])]))*u.deg,
                'channel_width': get_channel_width(cube[0].header,velocity=True) 
                }
   
    # first try a spectroscopic match
    if query.upper() in ['INTERNET','SIMBAD','NED']:
      
        #spectroscopic_table = find_NED_counterpart(cfg,source, header_info,\
        #    sysrange=150.*u.km/u.s,wide_search=wide_search)
        
        #if ('Object Name' in spectroscopic_table.colnames and
        #    spectroscopic_table['Object Name'][0] == 'No object Found'):
            # something went wrong with the NED query, so lets try cDS
        if query.upper() in ['INTERNET','NED']:
            spectroscopic_table_ned = find_counterpart(cfg,source, header_info,
                sysrange=150.*u.km/u.s, table_source='NED')  
        if query.upper() in ['INTERNET','SIMBAD']:
            spectroscopic_table_simbad = find_counterpart(cfg,source, header_info,
                sysrange=150.*u.km/u.s, table_source='SIMBAD')
        if check_table_length(spectroscopic_table_ned) > 0 and check_table_length(spectroscopic_table_simbad) > 0:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                spectroscopic_table = vstack([spectroscopic_table_ned, spectroscopic_table_simbad])
            vsys, radius = set_search_radius(cfg,source,header_info,150.*u.km/u.s,
                counterpart_region=cfg.input.counterpart_region)
            coord = [source[f'ra'],source[f'dec']]
            spectroscopic_table = sort_on_distance(cfg,spectroscopic_table, 
                coord,vsys,header_info=header_info)
        elif check_table_length(spectroscopic_table_ned) > 0:
            spectroscopic_table = spectroscopic_table_ned
        elif check_table_length(spectroscopic_table_simbad) > 0:
            spectroscopic_table = spectroscopic_table_simbad
        else:
            spectroscopic_table = QTable()    


        search_id = 'INTERNET'
        pref =''
    elif query.upper() == 'MANUAL':
        spectroscopic_table = find_counterpart(cfg,source, header_info,
            sysrange=150.*u.km/u.s,table_source='MANUAL',
            spectroscopic=cfg.input.spectroscopic_manual_counterparts)
        search_id = 'Manual'
        pref = 'sofia_'
        
  
   
    if check_table_length(spectroscopic_table) > 0:
        new_table = spectroscopic_table[0:1]
        final_row = combine_tables(new_table,source,column_indicators=[search_id,insource])
        if np.isnan(new_table['Velocity Diff'][0]):
            final_row[f'{search_id}_spectroscopic'] = False
        else:
            final_row[f'{search_id}_spectroscopic'] = True
           
    else:
        #Producing a same size empty table
        print_log(cfg,f'We found no {search_id} counterpart for {source[pref+"id"]}', 
            case=['verbose'])
        #if query.upper() == 'INTERNET':
        requested_columns, requested_dtypes, requested_units = \
                get_ned_requested_metadata(include_extra=True)
        
        if query.upper() == 'MANUAL':
            if cfg.input.manual_input_tables[0] is None:
                manual_table = QTable(names=['Object Name'],dtype=['str'])
            else:
                manual_table = read_manual_table(cfg,need_velocity=True)
           
            # Add the difference columns 
            add_units = [u.deg, u.km/u.s, u.dimensionless_unscaled]
            for i,col in enumerate(['Spatial Diff','Velocity Diff','Combined Diff']):
                if col not in manual_table.colnames:
                    manual_table.add_column(np.nan,name=col)
                    manual_table[col].unit = add_units[i]
            requested_columns = []
            requested_dtypes = []
            requested_units = []
            for col in manual_table.colnames:
                requested_columns.append(col)
                requested_dtypes.append(manual_table[col].dtype)
                requested_units.append(manual_table[col].unit)
    
        requested_values = []
        for dt in requested_dtypes:
            requested_values.append(get_nan_for_dtype(dt))    
        dummy_table = QTable(names=requested_columns \
            ,dtype=requested_dtypes,units = requested_units)
        dummy_table.add_row(requested_values)
        final_row = combine_tables(dummy_table,source,column_indicators=[search_id,insource])
        final_row[f'{search_id}_spectroscopic'] = False
    
    # We always want to return a table, even if it is empty, so we check if the final row is a table and if not we return the dummy table
    if not isinstance(final_row, (QTable,Table)):
        raise InputError(f'We expected a table but got {type(final_row)}')
    close_variables(spectroscopic_table)
    return final_row



def set_search_radius(cfg,source,header_info,sysrange=None,
        counterpart_region = 'Beam',wide_search = False):
    vsys = None
    pref = ''
    for col in source.colnames:
            if 'sofia_' in col.lower():
                pref = 'sofia_'
    if not sysrange is None:
        vsys = source[pref+'v_sofia'].to(u.km/u.s) #systemic in km/s



    # If we have velocity information we can do a 3 times wider area        
    radius = header_info['BMAJ']/2.

    # but we need to check it doesn't becaome to big for small sources
      

    if counterpart_region.lower() in ['beam']:
        
        pass
        
    elif counterpart_region.lower() in ['3beam']:
        radius = header_info['BMAJ']*3./2.
        
    elif counterpart_region.lower() in ['box']:
        radius = np.nanmax([np.nanmax([source[pref+'x'].value-source[pref+'x_min'].value,\
                         source[pref+'x_max'].value-source[pref+'x'].value,
                         source[pref+'y'].value-source[pref+'y_min'].value,
                         source[pref+'y_max'].value-source[pref+'y'].value
                         ])*header_info['pixelsize'],float(radius.value)])*u.deg
    elif counterpart_region.lower() in ['ellipse']:
        if f'{pref}ell_maj' in source.colnames:
            print_log(cfg,f'''This is the values from sofia
ell_maj = {source[f'{pref}ell_maj'].to(u.deg).value}, radius = {radius.to(u.deg).value}''',case=['debug'])
            radius = np.nanmax([float(source[f'{pref}ell_maj'].to(u.deg).value*0.1),
                                float(radius.to(u.deg).value)])*u.deg
            
    elif counterpart_region.lower() in ['full_ellipse']:
        #for the full elipse we do not allow a wider area with velocity information
        if f'{pref}ell_maj' in source.colnames:
            radius = np.nanmax([float(source[f'{pref}ell_maj'].to(u.deg).value),
                                float(radius.to(u.deg).value)])*u.deg
    else:
        raise InputError(f'We dont know what to do with {counterpart_region} for counterpart_region')
    if wide_search and counterpart_region.lower() not in ['full_ellipse','box']:
        radius *= 5.
    if f'{pref}ell_maj' in source.colnames:
        if radius > source[f'{pref}ell_maj'].to(u.deg):
            print_log(cfg,f'Reducing the search radius to the ellipse major axis {source[f"{pref}ell_maj"]}',
                case=['verbose'])
            radius = float(source[f'{pref}ell_maj'].to(u.deg).value)*u.deg
    
    return vsys, radius



''' Sort our table based on catalogue name'''
def sort_by_name(table):
   # The order we prefer things in 
    preferred_order = [['NGC'], ['UGC','ESO'], ['M']]
    new_table = copy_table_header(  table,)
   
    # add the rows to our table  based on the preferred order
    for names in preferred_order:
        for row in table:
            row_name = row['Object Name']
            for identifier in names:
                if row_name[:len(identifier)].upper() == identifier:
                    if len(identifier) == 1:
                        if not row_name[len(identifier):].isdigit():
                            continue
                    new_table.add_row(row)
    #Attach all final objects that did not have a preffered order
    selected_objects = [x['Object Name'] for x in new_table]
    for row in table:
        if row['Object Name'] not in selected_objects:
            new_table.add_row(row)

    return new_table            

'''Sort the table by distance'''
def sort_on_distance(cfg, table_in, coordinates, vsys, 
        header_info = None, weights = [1.,1.],spectroscopic=True):
    # this stupid table is not ordered so get names, types ra and dec 
    # and sort on distance
    if vsys is None and spectroscopic:
        print_log(cfg,
            f'No systemic velocity found for {table_in["Object Name"][0]}. This is a logical error',
            case=['main'])
        raise InputError(f'No systemic velocity found for {table_in["Object Name"][0]}. This is a logical error')


    only_sort =False
    for x in table_in.colnames:
        if x in ['Spatial Diff']:
            if vsys is None:
                only_sort = True
        elif x in ['Combined Diff']:
            if vsys is not None and not np.isnan(table_in[x][0]):
                only_sort = True
       
    table = copy.deepcopy(table_in)
    if not only_sort:
        print_log(cfg,f'Before sorting: {table_in[0]}',case=['verbose'])
        print_log(cfg,f'This the table type {type(table_in)}',case=['debug'])      
        if np.all(np.array(weights) == 1.) and not header_info is None:
            weights = [header_info['pixelsize'], 
                       header_info['channel_width']]
            
        print_log(cfg,f'Using weights {weights}',case=['debug'])
        if spectroscopic:
            #throw all entries that have velocity nan out of the table, as we cannot use them for the distance calculation
            table = table[~np.isnan(table['Velocity'])]
        # Do not apply the weights to the spatial distance or velocity differences 
        # as they are used a actual physical thresholds 
       
        table['Spatial Diff'] = calculate_projected_distance([table['RA'],table['DEC']],
            coordinates,no_PA = True)
        
   
        if not vsys is None:
            velocities = table['Velocity'].to(u.km/u.s)
           
            if not isquantity(vsys):
                raise InputError(f'vsys is not a quantity but {type(vsys)}')
            
            vsys = vsys.to(u.km/u.s)
           
            table['Velocity Diff'] = np.abs(vsys - velocities)
            table['Combined Diff'] = np.sqrt(
                (table['Velocity Diff'] / weights[1]).decompose()**2 +
                (table['Spatial Diff'] / weights[0]).decompose()**2
            )
        else:
            n = len(table)
            table['Velocity Diff'] = np.full(n, np.nan) * u.km/u.s
            table['Combined Diff'] = np.full(n, np.nan) * u.dimensionless_unscaled
       
    if not spectroscopic:    
        table.sort('Spatial Diff')
    else:
        table.sort('Combined Diff')
    print_log(cfg,f'After sorting:\n {table}',case=['verbose'])
        
    return table
           
