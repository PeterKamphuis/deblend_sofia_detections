from deblend_sofia_detections.support.errors import InputError
from deblend_sofia_detections.support.table_functions import check_table_length,\
    read_manual_table, combine_tables, copy_table_header
from deblend_sofia_detections.support.support_functions import convertRADEC,\
    isquantity, get_nan_for_dtype,calculate_projected_distance

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import QTable, Column,Table

from astroquery.ipac.ned import Ned
from xml.parsers.expat import ExpatError

import astropy.units as u
import copy
import numpy as np
import time

Ned.TIMEOUT = 3600


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


def make_empty_ned_search_table():
    requested_columns, requested_dtypes, requested_units = \
        get_ned_requested_metadata(include_extra=True)

    table = QTable(names=requested_columns, dtype=requested_dtypes,
        units=requested_units)
    table.add_row(['No object Found', np.nan, np.nan, np.nan,
        'Unknown', None, np.nan, np.nan, np.nan, np.nan])
    return table


def find_NED_counterpart(cfg,source,header_info, sysrange = None,\
        weights = [1.,1.,1.],wide_search = False ):
    requested_columns, dummy_dtypes, dummy_units = get_ned_requested_metadata()
    coordinates = [source['ra'].unmasked,source['dec'].unmasked]
    #coordinates, radius
    #[source['ra'],source['dec'],source['v_sofia']],\
    #            beam_size[0]/2.
    
    co = SkyCoord(ra=coordinates[0], dec=coordinates[1], frame='fk5')

    # Setup a search table in ned
    vsys, radius = set_search_radius(cfg,source,header_info,sysrange,
        counterpart_region=cfg.general.counterpart_region,wide_search=wide_search)
    if cfg.general.verbose:
        print("Querying NED")
        print(f'Search a radius {radius.to(u.arcmin)} around {", ".join(convertRADEC(*[x.value for x in coordinates[:2]]))}')
    # get the NED table
    internet_table = None
    query_errors = (ExpatError, ValueError, ConnectionError, TimeoutError, OSError)
    for attempt in range(3):
        try:
            internet_table = Ned.query_region(co, radius=radius, equinox='J2000.0')
            break
        except query_errors as e:
            if cfg.general.verbose:
                print(f'NED query failed on attempt {attempt + 1}/3: {e}')
            if attempt < 2:
                time.sleep(5 * (attempt + 1))

    if internet_table is None:
        if cfg.general.verbose:
            print('NED returned an invalid or temporary response; continuing without a NED counterpart.')
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
                                unit=internet_table[x].unit)
        else:
            result_table[x] = Column([None for x in range(check_table_length(internet_table))])
           
  
    #select out the galaxies
    objects_to_select = ['G','GPAIR','GTRPL']
    rows = [True if x.upper() in objects_to_select else False for x in result_table['Type']]
    search_table = result_table[rows]
   
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
    search_table = sort_on_distance(search_table, coordinates,vsys)
    if len(search_table) > 1:
        # then pick NGC/UGC/ESO/M matches
        search_table = sort_by_name(search_table)

    if check_table_length(search_table) > 0:
        name = search_table['Object Name'][0]
        if ':' in name:
            search_table['Object Name'][0] = name.split(':')[0]
    if search_table['Velocity Diff'].unit != u.km/u.s:
        print(f'Velocity Diff unit is {search_table["Velocity Diff"].unit}')
   
    return search_table


def find_manual_counterpart(cfg,source,header_info, sysrange=None):
    if cfg.input.manual_input_tables[0] is None:
        return  QTable()
    manual_table = read_manual_table(cfg)

    coordinates = [source['sofia_ra'],source['sofia_dec']]
    
    vsys, radius = set_search_radius(cfg,source,header_info,sysrange,
        counterpart_region=cfg.general.counterpart_region)
  
    search_table = sort_on_distance(manual_table, coordinates,vsys)
    if cfg.general.verbose:
        print(f'Searching for manual counterpart for {source["sofia_id"]}')
        print(f'Search a radius {radius.to(u.arcsec)} around {", ".join(convertRADEC(coordinates[0].value,coordinates[1].value))}')
        print(f' The nearest target is {search_table["Name"]} at a distance of {search_table["Spatial Diff"].to(u.arcsec)}')
        print(f'And the velocity difference is {search_table["Velocity Diff"].to(u.km/u.s)} to vsys {vsys.to(u.km/u.s)}')
   
  
    # if we have set a range of velocities we mask all that that are outside the range
    if search_table['Spatial Diff'][0] > 2.*radius.to(u.arcsec) or\
       search_table['Velocity Diff'][0] > sysrange:
       search_table = QTable()
   
    return search_table
   



def search_counter_part(cfg,source,sofia_directory= './',
        basename=None,query ='NED',insource=None,wide_search=False):
    '''Look for the optical counterpart of the source'''
    if isinstance(source, (Table, QTable)):
        source = source[0]
    try:
        inid = source['id']
    except:
        inid = source['sofia_id']
    if cfg.general.verbose:
        print(f'Searching in {query} to find a counterpart for {basename} with id {inid}')
   
    input_dir = f'{sofia_directory}/{basename}_cubelets'
    cube = fits.open(f'{input_dir}/{basename}_{inid}_cube.fits',\
        output_verify='warn')
    header_info= {'BMAJ':float(cube[0].header['BMAJ'])*u.deg,
                  'pixelsize': float(np.mean([abs(cube[0].header['CDELT1']),\
                                              abs(cube[0].header['CDELT1'])]))\
                                              *u.deg\
                                                 }
   
    # first try a spectroscopic match
    if query.upper() == 'NED':
      
        spectroscopic_table = find_NED_counterpart(cfg,source, header_info,\
            sysrange=150.*u.km/u.s,wide_search=wide_search)
        
        search_id = 'NED'
        pref =''
    elif query.upper() == 'MANUAL':
        spectroscopic_table = find_manual_counterpart(cfg,source, header_info,
            sysrange=150.*u.km/u.s)  
        search_id = 'Manual'
        pref = 'sofia_'
    confirmed = True
   
    if check_table_length(spectroscopic_table) > 0 and not (
        'Object Name' in spectroscopic_table.colnames and
        spectroscopic_table['Object Name'][0] == 'No object Found'):
        new_table = spectroscopic_table
    else:
        if cfg.general.verbose:
            print(f'We found no match within the velocity range, picking the closest object without velocity')
        if query.upper() == 'NED':
            possible_table = find_NED_counterpart(cfg,source, header_info)
        if query.upper() == 'MANUAL':
            possible_table = QTable()
        new_table = possible_table
        confirmed = False
    #We don't need the searching rows
    #false_table = True
    #if query.upper() == 'NED':
    #    false_table = False
    if check_table_length(new_table) > 0 and not (
        query.upper() == 'NED' and 'Object Name' in new_table.colnames and
        new_table['Object Name'][0] == 'No object Found'):
        new_table=new_table[0:1]
        final_row = combine_tables(new_table,source,column_indicators=[search_id,insource])
        final_row[f'{search_id}_spectroscopic'] = confirmed
    else:
        if cfg.general.verbose:
            print(f'We found no {search_id} counterpart for {source[pref+"id"]}')
        if query.upper() == 'NED':
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
            print(f'This is the values from sofia')
            print(source[f'{pref}ell_maj'].to(u.deg).value,radius.to(u.deg).value)
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
            if cfg.general.verbose:
                print(f'Reducing the search radius to the ellipse major axis {source[f"{pref}ell_maj"]}')
            radius = float(source[f'{pref}ell_maj'].to(u.deg).value)*u.deg
    
    return vsys, radius



''' Sort our table based on catalogue name'''
def sort_by_name(table):
    #print(table)
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
def sort_on_distance(table_in,coordinates,vsys,weights = [1.,1.],print_in=False):
    # this stupid table is not ordered so get names, types ra and dec and sort on distance
    for x in table_in.colnames:
        if x.lower() in ['spatial_diff','velocity_diff', 'combined_diff']:
            raise InputError(f'Column {x} is not allowed in the table, please rename it')     
    if print_in:
        print(f'Before sorting: {table_in}')    
        print(f'This the table type {type(table_in)}')
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
    if print_in:
        print(f'After sorting: {table}')
        
    return table
           
