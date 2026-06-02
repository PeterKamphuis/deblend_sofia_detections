from deblend_sofia_detections.support.errors import InputError
from deblend_sofia_detections.support.table_functions import check_table_length,\
    read_manual_table, combine_tables, copy_table_header
from deblend_sofia_detections.support.support_functions import convertRADEC,\
    isquantity, get_nan_for_dtype,calculate_projected_distance,close_variables
from deblend_sofia_detections.support.logging import print_log
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import QTable, Column,Table,vstack


from xml.parsers.expat import ExpatError

import astropy.units as u
import copy
import numpy as np
import time
import warnings


def get_ned_requested_metadata(include_extra=False):
    requested_columns = ['Object Name', 'RA', 'DEC', 'Velocity', 'Type',
        'Magnitude and Filter', 'Distance']
    requested_dtypes = ['U30', float, float, float, object, object, float]
    requested_units = [None, u.deg, u.deg, u.km/u.s, None, None, u.Mpc]

    if include_extra:
        requested_columns += ['Spatial Diff', 'Velocity Diff', 'Combined Diff']
        requested_dtypes += [float, float, float]
        requested_units += [u.deg, u.km/u.s, u.dimensionless_unscaled]

    return requested_columns, requested_dtypes, requested_units


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

def find_internet_counterpart(cfg,source,header_info, sysrange = None,\
        wide_search = False, query = 'ALL'):
    if query.upper() == 'NONE':
        print_log(cfg,'The user does not want to query the internet.',case=['verbose'])
        return make_empty_ned_search_table(), False

    requested_columns, requested_dtypes, dummy_units = get_ned_requested_metadata()
    coordinates = [source['ra'].unmasked,source['dec'].unmasked]
   
    co = SkyCoord(ra=coordinates[0], dec=coordinates[1], frame='fk5')

    # Setup a search table in ned
    vsys, radius = set_search_radius(cfg,source,header_info,sysrange,
        counterpart_region=cfg.general.counterpart_region,wide_search=wide_search)
    
    print_log(cfg, f'''Querying {query} for source {source["name"]} with id {source["id"]}.\n \
Search radius = {radius.to(u.arcmin)}  around {", ".join(convertRADEC(cfg,*[x.value for x in coordinates[:2]]))}''', 
case=['verbose'])
    
    if query.upper() in ['NED','ALL']:
        search_table_ned = find_NED_counterpart(cfg,source,co,radius, 
                                requested_columns, requested_dtypes)
    else:
        search_table_ned = make_empty_ned_search_table(include_extra=False)

    if query.upper() in ['SIMBAD','ALL'] :
        search_table_simbad = find_SIMBAD_counterpart(cfg,source,co,radius, 
            requested_columns, requested_dtypes)
    else:
        search_table_simbad = make_empty_ned_search_table(include_extra=False)
    if check_table_length(search_table_ned) > 0 and check_table_length(search_table_simbad) > 0:
        for col in search_table_ned.colnames:
            print(col)
            print_log(cfg,f'for the column {col} we have {search_table_ned[col].dtype} in ned and {search_table_simbad[col].dtype} in simbad',
                      case=['debug'])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")        
            search_table = vstack([search_table_ned, search_table_simbad])  
    elif check_table_length(search_table_ned) > 0:
        search_table = search_table_ned  
    elif check_table_length(search_table_simbad) > 0:
        search_table = search_table_simbad
    else:
        search_table = make_empty_ned_search_table()
    success = True    
    # if we have set a range of velocities we mask all that that are outside the range
    if not sysrange is None:
        rows =  [True if vsys-sysrange < x < vsys+sysrange else False for x in search_table['Velocity']]
        tmp_table = search_table[rows]
        if check_table_length(tmp_table) == 0:
             print_log(cfg,f'We found no match within the velocity range of {sysrange}, picking the closest object without velocity',
             case=['verbose'])
             search_table['Velocity'] = [float('NaN') if x is None else x for x in search_table['Velocity']]
             success = False
        else:
            search_table = tmp_table
    else:
        # if we do not have a velocity range we have to make sure that we take the None as float
        search_table['Velocity'] = [float('NaN') if x is None else x for x in search_table['Velocity']]
       
    if check_table_length(search_table) > 0:
        if not isquantity(search_table['Velocity']):
            search_table['Velocity'].unit = u.km/u.s
        else:
            if search_table['Velocity'].unit != u.km/u.s:
                try:
                    search_table['Velocity'] = [x.to(u.km/u.s).value for x in search_table['Velocity']]
                    search_table['Velocity'].unit = u.km/u.s
                except:
                    raise InputError(f'the NED counterpart for {source["Name"]} has a weird velocity unit {search_table["Velocity"].unit}')
   
    # first sort on distance, we always do this cause else we do not get the Spatial and combined columns
    #distance to the source not actual distance of the object
    search_table = sort_on_distance(cfg,search_table, coordinates,vsys)
    if len(search_table) > 1:
        # then pick NGC/UGC/ESO/M matches
        search_table = sort_by_name(search_table)

    if check_table_length(search_table) > 0:
        name = search_table['Object Name'][0]
        if ':' in name:
            search_table['Object Name'][0] = name.split(':')[0]
    if search_table['Velocity Diff'].unit != u.km/u.s:
        print_log(cfg,f'Velocity Diff unit is {search_table["Velocity Diff"].unit}')
   
    return search_table,success
    

def find_NED_counterpart(cfg,source,co,radius, requested_columns, requested_dtypes):
    '''
def find_NED_counterpart(cfg,source,header_info, sysrange = None,\
        weights = [1.,1.,1.],wide_search = False ):
    
    requested_columns, dummy_dtypes, dummy_units = get_ned_requested_metadata()
    coordinates = [source['ra'].unmasked,source['dec'].unmasked]
   
    co = SkyCoord(ra=coordinates[0], dec=coordinates[1], frame='fk5')

    # Setup a search table in ned
    vsys, radius = set_search_radius(cfg,source,header_info,sysrange,
        counterpart_region=cfg.general.counterpart_region,wide_search=wide_search)
    
    print_log(cfg, f''Querying NED \
Search a radius {radius.to(u.arcmin)} around {", ".join(convertRADEC(cfg,*[x.value for x in coordinates[:2]]))}'', 
case=['verbose','screen'])
    # get the NED table
    # as astro query is highly unstable only import when we need it
    '''
    from astroquery.ipac.ned import Ned
    Ned.TIMEOUT = 3600
    Ned.clear_cache()
    internet_table = None
    query_errors = (ExpatError, ValueError, ConnectionError, TimeoutError, OSError)


    for attempt in range(3):
        try:
            internet_table = Ned.query_region(co, radius=radius, equinox='J2000.0')
            break
        except query_errors as e:    
            print_log(cfg, f'NED query failed on attempt {attempt + 1}/3: {e}', case=['verbose','screen'])
            if attempt < 2:
                time.sleep(5 * (attempt + 1))
    print_log(cfg, f'NED query completed with {check_table_length(internet_table)} results', case=['screen'])  
    if internet_table is None:
        print_log(cfg, 'NED returned an invalid or temporary response; continuing without a NED counterpart.',
            case=['verbose'])
        return make_empty_ned_search_table()
    
    # as astropy is the dumbest project ever they can not be consistant so 
    # we have to correct the units for RA and DEg
    internet_table['RA'].unit=u.deg
    internet_table['DEC'].unit=u.deg
    # Astropy is so stupid that it does not provide a QTable from the query
    # so we have to do this as well. 
    result_table = QTable()
    
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
    return search_table
    '''
    # if we have set a range of velocities we mask all that that are outside the range
    if not sysrange is None:
        rows =  [True if vsys-sysrange < x < vsys+sysrange else False for x in search_table['Velocity']]
        search_table = search_table[rows]
    else:
        # if we do not have a velocity range we have to make sure that we take the None as float
        search_table['Velocity'] = [float('NaN') if x is None else x for x in search_table['Velocity']]
       
    if check_table_length(search_table) > 0:
        if not isquantity(search_table['Velocity']):
            search_table['Velocity'].unit = u.km/u.s
        else:
            if search_table['Velocity'].unit != u.km/u.s:
                try:
                    search_table['Velocity'] = [x.to(u.km/u.s).value for x in search_table['Velocity']]
                    search_table['Velocity'].unit = u.km/u.s
                except:
                    raise InputError(f'the NED counterpart for {source["Name"]} has a weird unit {search_table["Velocity"].unit}')
   
    # first sort on distance, we always do this cause else we do not get the Spatial and combined columns
    search_table = sort_on_distance(cfg,search_table, coordinates,vsys)
    if len(search_table) > 1:
        # then pick NGC/UGC/ESO/M matches
        search_table = sort_by_name(search_table)

    if check_table_length(search_table) > 0:
        name = search_table['Object Name'][0]
        if ':' in name:
            search_table['Object Name'][0] = name.split(':')[0]
    if search_table['Velocity Diff'].unit != u.km/u.s:
        print_log(cfg,f'Velocity Diff unit is {search_table["Velocity Diff"].unit}')
   
    return search_table
    '''

def find_SIMBAD_counterpart(cfg,source,co,radius, requested_columns,requested_dtypes):
    from astroquery.simbad import Simbad
    Simbad.TIMEOUT = 3600
    Simbad.clear_cache()
    lamp = Simbad()
     
    lamp.add_votable_fields('main_id', 'ra', 'dec', 'rvz_radvel', 'otype_txt','morph_type', 'V')
    internet_table = None
    query_errors = (ExpatError, ValueError, ConnectionError, TimeoutError, OSError)


    for attempt in range(3):
        try:
            internet_table = lamp.query_region(co, radius=radius)
            break
        except query_errors as e:    
            print_log(cfg, f'SIMBAD query failed on attempt {attempt + 1}/3: {e}', case=['verbose','screen'])
            if attempt < 2:
                time.sleep(5 * (attempt + 1))
    print_log(cfg, f'SIMBAD query completed with originally {check_table_length(internet_table)} results', case=['screen'])  
    if internet_table is None:
        print_log(cfg, 'SIMBAD returned an invalid or temporary response; continuing without a SIMBAD counterpart.',
            case=['verbose'])
        return make_empty_ned_search_table()
    
    # as astropy is the dumbest project ever they can not be consistant so 
    # we have to correct the units for RA and DEg
    internet_table['ra'].unit=u.deg
    internet_table['dec'].unit=u.deg
    # Astropy is so stupid that it does not provide a QTable from the query
    # so we have to do this as well. 
    result_table = QTable()
    translation_table = get_SIMBAD_translation_table()
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
    objects_to_select = ['Sy1','BLL','LSB','Bla','BiC','SyG','rG', 'bCG', 'Sy2',
         'SBG', 'LIN', 'QSO', 'H2G', 'EmG', 'AGN', 'G', 'GiP','GiG','GiC','CGG',
         'IG','PaG','GrG','ClG','SCG']
    objects_to_select = [x.upper() for x in objects_to_select]
    rows = [True if x.upper() in objects_to_select else False for x in result_table['Type']]
    search_table = result_table[rows]
    print_log(cfg, f'SIMBAD final query completed with {check_table_length(search_table)} galaxies', case=['screen'])  
    return search_table



def find_manual_counterpart(cfg,source,header_info, sysrange=None):
    if cfg.input.manual_input_tables[0] is None:
        return  QTable()
    manual_table = read_manual_table(cfg)

    coordinates = [source['sofia_ra'],source['sofia_dec']]
    
    vsys, radius = set_search_radius(cfg,source,header_info,sysrange,
        counterpart_region=cfg.general.counterpart_region)
  
    search_table = sort_on_distance(cfg,manual_table, coordinates,vsys)
   
    print_log(cfg,f'''Searching for manual counterpart for {source["sofia_id"]}
Search a radius {radius.to(u.arcsec)} around {", ".join(convertRADEC(cfg,coordinates[0].value,coordinates[1].value))}
The nearest target is {search_table["Name"]} at a distance of {search_table["Spatial Diff"].to(u.arcsec)}
And the velocity difference is {search_table["Velocity Diff"].to(u.km/u.s)} to vsys {vsys.to(u.km/u.s)}''')
   
  
    # if we have set a range of velocities we mask all that that are outside the range
    if search_table['Spatial Diff'][0] > 2.*radius.to(u.arcsec) or\
       search_table['Velocity Diff'][0] > sysrange:
       search_table = QTable()
   
    return search_table
   
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


def search_counter_part(cfg,source,sofia_directory= './',
        basename=None,query ='INTERNET',insource=None,wide_search=False):
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
    cube = fits.open(f'{input_dir}/{basename}_{inid}_cube.fits',\
        output_verify='warn')
    header_info= {'BMAJ':float(cube[0].header['BMAJ'])*u.deg,
                  'pixelsize': float(np.mean([abs(cube[0].header['CDELT1']),\
                                              abs(cube[0].header['CDELT1'])]))\
                                              *u.deg\
                                                 }
   
    # first try a spectroscopic match
    if query.upper() == 'INTERNET':
      
        #spectroscopic_table = find_NED_counterpart(cfg,source, header_info,\
        #    sysrange=150.*u.km/u.s,wide_search=wide_search)
        
        #if ('Object Name' in spectroscopic_table.colnames and
        #    spectroscopic_table['Object Name'][0] == 'No object Found'):
            # something went wrong with the NED query, so lets try cDS
        spectroscopic_table,success = find_internet_counterpart(cfg,source, header_info,\
                sysrange=150.*u.km/u.s,wide_search=wide_search, query=cfg.input.internet_query)
    
        search_id = 'INTERNET'
        pref =''
    elif query.upper() == 'MANUAL':
        spectroscopic_table = find_manual_counterpart(cfg,source, header_info,
            sysrange=150.*u.km/u.s)  
        search_id = 'Manual'
        pref = 'sofia_'
    confirmed = True
   
    if check_table_length(spectroscopic_table) > 0:
        new_table = spectroscopic_table
    else:
       
        print_log(cfg,f'We found no match within the velocity range, picking the closest object without velocity',
            case=['verbose'])
        if query.upper() == 'INTERNET':
            if not success:
                possible_table = spectroscopic_table
            else:
                possible_table,success = find_internet_counterpart(cfg,source, 
                    header_info, query=cfg.input.internet_query)
        if query.upper() == 'MANUAL':
            possible_table = QTable()
        new_table = possible_table
        confirmed = False
        close_variables(possible_table)
    #We don't need the searching rows
    #false_table = True
    #if query.upper() == 'NED':
    #    false_table = False
    if check_table_length(new_table) > 0:
        new_table=new_table[0:1]
        final_row = combine_tables(new_table,source,column_indicators=[search_id,insource])
        final_row[f'{search_id}_spectroscopic'] = confirmed
    else:
       
        print_log(cfg,f'We found no {search_id} counterpart for {source[pref+"id"]}', 
            case=['verbose'])
        if query.upper() == 'INTERNET':
            requested_columns, requested_dtypes, requested_units = \
                get_ned_requested_metadata(include_extra=True)
        elif query.upper() == 'MANUAL':
            manual_table = read_manual_table(cfg,need_velocity=False)
            # Add the difference columns 
            add_units = [u.km/u.s,u.deg, u.km/u.s, u.dimensionless_unscaled]
            for i,col in enumerate(['Velocity','Spatial Diff','Velocity Diff','Combined Diff']):
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
    close_variables(spectroscopic_table, new_table)
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
def sort_on_distance(cfg, table_in, coordinates, vsys, weights = [1.,1.]):
    # this stupid table is not ordered so get names, types ra and dec and sort on distance
    for x in table_in.colnames:
        if x.lower() in ['spatial_diff','velocity_diff', 'combined_diff']:
            raise InputError(f'Column {x} is not allowed in the table, please rename it')     
    
    print_log(cfg,f'Before sorting: {table_in}',case=['verbose'])
    print_log(cfg,f'This the table type {type(table_in)}',case=['debug'])
    table = copy.deepcopy(table_in)
    table['Spatial Diff'] = [calculate_projected_distance([x,\
        y],coordinates,no_PA = True).value for x,y in zip(\
        table['RA'],table['DEC'])]* u.deg
    if not vsys is None:
        velocities = table['Velocity'].to(u.km/u.s)
        if isquantity(velocities):
            velocities = velocities.value
        if isquantity(vsys):
            vsys = vsys.to(u.km/u.s).value
        table['Velocity Diff'] = [float(abs(vsys-z)/weights[1])\
            for z in velocities] * u.km/u.s
        table['Combined Diff']= [float(np.sqrt(x.value**2+y.to(u.arcsec).value**2)) for x,y in\
            zip(table['Velocity Diff'],table['Spatial Diff'])] * u.dimensionless_unscaled
             
        table.sort('Combined Diff')
       
    else:
      
        table['Velocity Diff'] = [float('NaN') for x in table['Spatial Diff']]\
            * u.km/u.s
        table['Combined Diff'] = [float('NaN') for x in table['Spatial Diff']]\
            * u.dimensionless_unscaled
        table.sort('Spatial Diff')
    print_log(cfg,f'After sorting: {table}',case=['verbose'])
        
    return table
           
