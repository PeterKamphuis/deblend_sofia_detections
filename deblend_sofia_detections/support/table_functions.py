

from deblend_sofia_detections.support.errors import TableError,InputError
from deblend_sofia_detections.support.support_functions import\
    is_real_unit, isquantity, translate_string_to_unit, convertRADEC
from deblend_sofia_detections.support.logging import print_log
from astropy import units as u
from astropy.table import QTable,Table,Row

import copy
import os
import pickle
import numpy as np
import warnings

def check_table_length(table):
    " Check the length of an astropy table or row"
    if isinstance(table,(QTable)):
        length= len(table)
        if 'Object Name' in table.colnames and length == 1:
            if table['Object Name'][0] == 'No object Found':
                length = 0
    elif isinstance(table,(Table)):
        length= len(table)
    elif isinstance(table,(Row)):
        length = 1
    elif table is None:
        length = 0
    else:
        raise TableError(f'{type(table)} is not a astropy Table or Row')   
    return length

def combine_tables(tableone, tabletwo, column_indicators=[None,None]):
    """
    Combine two tables, ensuring columns match in units and length.

    Parameters:
    - tableone: First table (Astropy Table or Row)
    - tabletwo: Second table (Astropy Table or Row)
    - column_indicators: Optional list of prefixes to differentiate columns
    
    Returns:
    - Combined table with data from both tables.
    """
    
    # Check if tables have the same length 
    
    if check_table_length(tableone) != check_table_length(tabletwo):
        raise ValueError(f"Tables have different lengths ({check_table_length(tableone)}, {check_table_length(tabletwo)}).")
    
    # Initialize lists for final table construction
    input_columns = []
    convert_units = []
    dtypes = []
    inputrows = []
    
    # Loop through both tables
    for i, table in enumerate([tableone, tabletwo]):
        table_units = []
        
        # Handle each column in the table
        for col in table.colnames:
            # Apply column indicators if provided
            colname = f'{column_indicators[i]}_{col}' if not\
                column_indicators[i] is None else col
            if colname in input_columns:
                raise ValueError(f"Column '{colname}' is already in the constructor.")

            input_columns.append(colname)

            # Get the unit of the column if available
            unit = None
            if hasattr(table[col], 'unit'):
                unit = table[col].unit
            else:
                try:
                    if hasattr(table[col][0], 'unit'):
                        unit = table[col][0].unit
                except Exception:
                    unit = None

            # Normalize empty / missing units
            if unit in ['', 'None']:
                unit = None

            # Check for valid units
            if unit is not None and not is_real_unit(unit):
                raise ValueError(f"The unit '{unit}' is not recognized.")

            convert_units.append(unit)
            table_units.append(unit)
            try:
                dtypes.append(table[col].dtype)
            except AttributeError:
                if isinstance(table[col], str):
                    dtypes.append(str)
                else:
                    dtypes.append(float)
                      # Process rows in the table
        for j in range(check_table_length(table)):
            # Handle the row depending on whether it's a Table or Row
            rowin = table[j] if isinstance(table, (Table, QTable)) else table

            newrow = []
            for value, tabunit in zip(rowin, table_units):
                # Convert value to the proper unit if necessary
                if not isquantity(value) and tabunit is not None and not\
                    tabunit == f'str':
                    value *= tabunit
                try:    
                    newrow.append(value.unmasked)
                except AttributeError:
                    # If the value is not a Masked Quantity, just append it
                    newrow.append(value)
            # Add new row to the list of rows
            if i == 0:
                inputrows.append(newrow)
            else:
                inputrows[j] += newrow
    
    # Create the combined table
   
    combined_table = QTable(names=input_columns, units=convert_units, dtype=dtypes)
   
    # Add rows to the combined table
   
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for row in inputrows:
            combined_table.add_row(row)
  
   
    return combined_table

def copy_table_header(input_table):
    '''Copy an astropy table without including any rows'''
    if isinstance(input_table, (QTable, Table)):
        copied_table = copy.deepcopy(input_table[0:1])
        copied_table.remove_row(0)
    elif isinstance(input_table, Row): 
        copied_table = QTable(input_table,copy=True)
        copied_table.remove_row(0)
    else:
        raise TableError(f'Input is not a valid astropy table or row: {type(input_table)}')
   
    return copied_table

def identify_object_names(cfg,table):
    """
    Identify the object name column in a table based on common keywords.
    
    Parameters:
    - table: Astropy Table containing the data.
    
    Returns:
    - Updated table with the identified object name column.
    """
    clear_name_columns = ['name', 'galaxy', 'object', 'NGC', 'IC', 'UGC',
                           'PGC', 'ESO', '2MASX', 'SDSS', 'WISEA']
    check_name_columns = ['id', 'identifier']
    found = False
    for col in table.colnames:
        if col.lower() in clear_name_columns:
            table['Object Name'] = table[col].copy()
            found = True
            break 
    if not found:
        for col in table.colnames:
            if col.lower() in check_name_columns:
                if isinstance(table[col][0], str):
                    try:
                        tmp =float(table[col][0])
                    except ValueError:               
                        table['Object Name'] = table[col].copy()
                        found = True
                        break    
       
    if found:   
        return table
    else:
        raise TableError('No object name column found in the table.')

def identify_velocity_column(cfg,table):
    """
    Identify the velocity column in a table based on common keywords.
    
    Parameters:
    - table: Astropy Table containing the data.
    
    Returns:
    - Updated table with the identified velocity column.
    """
    possible_velocity_columns = ['v_rad', 'v_sofia','cz','v_opt'
        ,'vel','vsys','v_sys','v_hel','v_optical', 'v_helio', 'v_lsr']
    found = False
    for col in table.colnames:
        if col.lower() in possible_velocity_columns:
            if table[col].unit != u.km/u.s:
                try:
                    table[col] = table[col].to(u.km/u.s)
                    table['Velocity'] = table[col].copy()
                    found = True
                    break
                except Exception as e:
                    print_log(cfg, f"Error converting {col} to km/s: {e}", case=['verbose'])
                    pass
            else:
                table['Velocity'] = table[col].copy()
                found = True
                break     
    if found:   
        return table
    else:
        raise TableError('No velocity column found in the table.')
    

def load_table(table_in,fresh_read=False,cfg=None,pickle_output=None):
    table_base,ext = os.path.splitext(table_in)
    table_filename = os.path.basename(os.path.normpath(table_base))
    if pickle_output is not None:
        pickle_file = f'{pickle_output}{table_filename}.pkl'
    else:           
        pickle_file = f'{table_base}.pkl'   
    
    if (not (fresh_read or cfg.input.original_tables) and os.path.exists(pickle_file)): 
        with open(pickle_file,'rb') as tmp:
            table = pickle.load(tmp)
        if type(table) != QTable:
            raise InputError(f'Your pickle file is not a Quantity table, we do not know how to read it.')
    elif ext == '.pkl':
        with open(table_in,'rb') as tmp:
            table = pickle.load(tmp)
    elif ext == '.txt':
        table = read_text_table(cfg,table_in)
        with open(pickle_file,'wb') as tmp:
            pickle.dump(table,tmp)
    elif ext == '.csv':
        table = read_text_table(cfg,table_in,seperator=',')
        with open(pickle_file,'wb') as tmp:
            pickle.dump(table,tmp)
    else:
        raise InputError(f'We do not recognize the extension {ext}, we do not know how to read the file.')
    return table

def read_manual_table(cfg, need_velocity =True):    
    manual_table = None
    for table_in in cfg.input.manual_input_tables:
        if table_in is None:
            manual_table = QTable()
            continue
        if not os.path.isfile(table_in):
            raise InputError(f'Could not find manual input table {table_in}')
          
        print_log(cfg,f'Loading manual input table {table_in}',case=['verbose'])
        manual_table_small = load_table(table_in, cfg=cfg)

        if manual_table is None:
            manual_table = manual_table_small
        else:
            manual_table = combine_tables(manual_table,manual_table_small)  
    
    if 'velocity' not in [x.lower() for x in manual_table.colnames]\
        and need_velocity:
        manual_table = identify_velocity_column(cfg,manual_table) 
    if 'object name' not in [x.lower() for x in manual_table.colnames]:
        manual_table= identify_object_names(cfg,manual_table)
    return manual_table


def read_text_table(cfg,file,seperator=' '):
    sources = None
    with open(file) as tmp:
        lines =tmp.readlines()
    
    for i,line in enumerate(lines):
        print(f"\r Processing the text table: {(i+1)/len(lines)*100.:.1f} % done. ",\
             end =" ",flush = True) 
        
        tmp =line.split(seperator)
     
        if line.strip() in ['','#','%']:
            pass
        elif line[0] in ['#','%']:
            #These two ifs can not be combine as then it will process lines with comments without the seperator
            if len(tmp) > 1:
                tmp =line[1:].split(seperator)
                if tmp[0].strip().lower() in ['name','id','galaxy']:
                    # get the present columns
                    input_columns  = [x.strip() for x in tmp]
                   
                    #determin their location in the line
                    column_locations = []
                    for col in input_columns:
                        column_locations.append(line.find(col)+len(col))
                        columns_triggered = True
                        # Modify certain columns to a unified name
                        if col.lower() in ['name','id','galaxy']:
                            if not 'name' in [x.lower() for x in input_columns] or\
                                col.lower() == 'name':
                                input_columns[input_columns.index(col)] = 'Name'
                            
                        if 'ra' in col.lower():
                            input_columns[input_columns.index(col)] = 'RA'
                        if 'dec' in col.lower():
                            input_columns[input_columns.index(col)] = 'DEC'
                      
                        if ('d25' in col.lower() or  'd_25' in col.lower())\
                            and not ('err' in col.lower() or 'arcsec' in col.lower()):
                            input_columns[input_columns.index(col)] = 'D25'
                        if 'SFR' in col.lower():
                            if not 'SFR' in input_columns:
                                input_columns[input_columns.index(col)] = 'SFR'
                        if 'm' in col.lower() and  '*' in col.lower() :
                            input_columns[input_columns.index(col)] = 'M*' 
                        if 'dist' in col.lower()  or  'D' == col :
                            input_columns[input_columns.index(col)] = 'Distance'
                        if 'sys' in col.lower():
                            input_columns[input_columns.index(col)] = 'Velocity'
                   
                    
                    
                elif columns_triggered:
                    columns_triggered = False
                    
                    input_units  = [x.strip().lower() for x in tmp]
                    convert_units = []
                    for unit in input_units:
                        convert_units.append(translate_string_to_unit(unit))
                    for i,column in enumerate(input_columns):
                        if column.lower() in ['id','name','galaxy','identifier',
                            'object','typ','type','morph','morphology']:
                            convert_units[i] = 'str'

                    dtypes = []
                    for i,unit in enumerate(convert_units):
                        if unit == 'str' or unit == str:
                            dtypes.append(str)
                            convert_units[i] = None
                        elif unit == bool:
                            dtypes.append(bool)
                            convert_units[i] = None
                        elif unit == u.pix:
                            dtypes.append(float)
                        else:
                            dtypes.append(np.float64)
                    # Create the table with the columns 
                    convert_units = [u.dimensionless_unscaled if x in [bool,str,int,'str'] else x for x in convert_units]
                    sources = QTable(names=input_columns,units=convert_units,dtype=dtypes)
                       
        else:
           
            construct_row=[]
            name = tmp[input_columns.index('Name')]
            for i,col in enumerate(input_columns):
                if seperator == ' ':
                    if i == 0:
                        start = 0
                    else:
                        start = column_locations[i-1]
                    end = column_locations[i]
                    value = line[start:end].strip()
                else:
                    value= tmp[i]
                
                if convert_units[i] == 'str':
                    construct_row.append(value.strip())
                elif convert_units[i] in ['hms','dms']:
                
                    if sources[input_columns[i]].unit != u.deg:
                        tmp_column= sources[input_columns[i]]
                        sources[input_columns[i]] = Column(tmp_column.value,\
                                unit=u.deg)
                    if convert_units[i] == 'hms':
                        deg,dummy = convertRADEC(cfg,value,'0d0m0s',invert=True)
                        construct_row.append(deg*u.deg)
                        
                    else:
                      
                        dummy,deg = convertRADEC(cfg,'0h0m0s',value,invert=True)
                        construct_row.append(deg*u.deg)
                else:
                    if dtypes[i] == bool:
                        if value.strip().lower() in ['true','yes','1']:
                            construct_row.append(True)
                        elif value.strip().lower() in ['false','no','0']:
                            construct_row.append(False)
                        else:
                            raise InputError(f'Value |{value}| is not a valid boolean.')
                        #construct_row.append(bool(value))
                      
                    elif dtypes[i] == str:
                        construct_row.append(str(value))
                    elif dtypes[i] == int:
                        construct_row.append(int(value))
                    
                    else:
                        try:
                            construct_row.append(float(value)*convert_units[i])
                        except ValueError:
                            print_log(cfg,
                                f"Value '{value}' could not be converted to float with unit {convert_units[i]}. Setting to NaN.",
                                case = ['debug'])
                            construct_row.append(np.nan * convert_units[i])
            sources.add_row(construct_row)
    print(f"\n")
  
    return sources