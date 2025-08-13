from deblend_sofia_detections.support.errors import TableError,InputError
from deblend_sofia_detections.support.support_functions import\
    is_real_unit, isquantity

from astropy.table import QTable,Table,Row

import copy
import os
import pickle

def check_table_length(table):
    " Check the length of an astropy table or row"
    if isinstance(table,(Table)):
        length= len(table)
    elif isinstance(table,(QTable)):
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
            if hasattr(table[col], 'unit'):
                unit = table[col].unit
            else:
                # Check if the first element has a unit (in case it's a Quantity column)
                if hasattr(table[col][0], 'unit'):
                    unit = table[col][0].unit
                else:
                    unit = None
            
            # Check for valid units
            if not unit is None and not is_real_unit(unit):
                raise ValueError(f"The unit '{unit}' is not recognized.")
            
            convert_units.append(unit)
            table_units.append(unit)
            dtypes.append(table[col].dtype)
        
        # Process rows in the table
        for j in range(check_table_length(table)):
            # Handle the row depending on whether it's a Table or Row
            rowin = table[j] if isinstance(table, Table) else table
            
            newrow = []
            for value, tabunit in zip(rowin, table_units):
                # Convert value to the proper unit if necessary
                
                if not isquantity(value) and not tabunit is None and not\
                    tabunit == f'str':
                    value *=tabunit 
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


def identify_velocity_column(table):
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
                    print(f"Error converting {col} to km/s: {e}")
                    pass
            else:
                table['Velocity'] = table[col].copy()
                found = True
                break     
    if found:   
        return table
    else:
        raise TableError('No velocity column found in the table.')
    


def read_manual_table(cfg, need_velocity =True):    
    manual_table = None
    for table_in in cfg.input.manual_input_tables:
        if table_in is None:
            manual_table = QTable()
            continue
        if not os.path.isfile(table_in):
            raise InputError(f'Could not find manual input table {table_in}')
        if cfg.general.verbose:     
            print(f'Loading manual input table {table_in}')
        with open(f'{table_in}','rb') as tmp:
            try:
                manual_table_small = pickle.load(tmp)
            except Exception as e:
                raise InputError(f''' We are expecting a pickle table to ensure speedy loading. 
Error loading manual input table {table_in}: {e}''')
        if manual_table is None:
            manual_table = manual_table_small
        else:
            manual_table = combine_tables(manual_table,manual_table_small)  
    
    if 'velocity' not in [x.lower() for x in manual_table.colnames]\
        and need_velocity:
        manual_table = identify_velocity_column(manual_table) 
    return manual_table
