
from deblend_sofia_detections import template as templates
from deblend_sofia_detections.support.errors import InputError,SofiaError
from deblend_sofia_detections.support.support_functions import \
    convert_pix_columns_to_arcsec,translate_string_to_unit,get_source_cat_name,\
    get_start_end_locations,convert_pixel_values_to_original
from deblend_sofia_detections.support.system_functions import convert_ps

try:
    from importlib.resources import open_text as pack_open_txt
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    from importlib_resources import open_text as pack_open_txt

from astropy.table import QTable
from astropy.io import votable
from astropy import units as u

import os
import numpy as np
import shutil
import subprocess
import string


def check_parameters(table,variables=None,no_conversion=False):
    '''Check that the parameters in the sofia catalogue are correct and return the velocity column'''
    input_columns = [x.lower() for x in table.colnames]
    if variables is None:
        variables = input_columns
    velocity = None
    possible_velocities = ['v_rad','v_opt','v_app']
    for value in variables:
        trig = False
        if value.lower() in input_columns:
            if value.lower() in possible_velocities and velocity is None:
                velocity = value.lower()   

        elif value.lower() == 'v_sofia':
            for vel in possible_velocities:
                for col in input_columns:
                    if vel.lower() == col and velocity is None:
                        velocity = vel           
            if velocity is None:
                trig = True    
        else:
            trig = True

        if trig:
           raise InputError(f'''SOFIA_CATALOGUE: We cannot find the required column for {value} in the sofia catalogue.
{"":8s}SOFIA_CATALOGUE: This can happen because a) you have tampered with the sofiainput.txt file in the Support directory,
{"":8s}SOFIA_CATALOGUE: b) you are using an updated version of SoFiA2.
''')
  
                #sources[name]['v_sofia'] = sources[name][velocity] 
    #Rename the velocity 
    if not no_conversion:
        for col in table.colnames:
            if velocity == col.lower():
                new_col = col.replace(velocity,'v_sofia')
                table.rename_column(col,new_col)
                if col in variables:
                    variables[variables.index(col)] = new_col
                col = new_col
            if col.lower() not in variables:
                table.remove_column(col)
    return table    
check_parameters.__doc__ =f'''
 NAME:
    check_parameters(Variables,input_columns)

 PURPOSE:
    check wether all  variables are in the sofia catalogue
 
 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Variable = Reguired variable
    input_columns = found columns
    

 OPTIONAL INPUTS:



 OUTPUTS:
    velocity = the velocity found in the catalogue
 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:

'''
'''Run the modified sofia file on the cube'''

def execute_sofia(cfg,run_directory='Sofia_Output',
        sofia_parameter_file='sofia_input.par'):
    indir = os.getcwd()
    os.chdir(f'{run_directory}')
    sfrun = subprocess.Popen(['sofia',sofia_parameter_file], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    sofia_run, sofia_warnings_are_annoying = sfrun.communicate()
    if cfg.general.verbose:
        print(sofia_run.decode("utf-8"))
        print(sofia_warnings_are_annoying.decode("utf-8"))
    if sfrun.returncode == 8:
        with open(f'sofia_output.txt','w') as file:
            file.writelines(sofia_run.decode("utf-8"))
            file.writelines(sofia_warnings_are_annoying.decode("utf-8"))
        os.chdir(indir)    
        return 'No sources found'
    elif sfrun.returncode != 0:
        with open(f'sofia_output.txt','w') as file:
            file.writelines(sofia_run.decode("utf-8"))
            file.writelines(sofia_warnings_are_annoying.decode("utf-8"))
        os.chdir(indir)    
        raise SofiaError(f'Sofia run failed with return code {sfrun.returncode}. Check sofia_output.txt for details.')
     
 
    #Convert the ps files
    all_files_and_directories = os.listdir(f'{run_directory}')
    for file in all_files_and_directories:
        if os.path.splitext(file)[1] in ['.ps','.eps']:
            convert_ps(file) 


    os.chdir(indir)
    return 'Success'



def load_sofia_basename(filename):
    '''Obtain the sofia basename from the input file'''
    input_file = load_sofia_input_file(filename)
    return os.path.basename(os.path.splitext(input_file['input.data'])[0])

def load_sofia_catalogue(filename, variables = None,verbose = False, no_conversion=False):
    '''Read a specified sofia table into a Astropy QTable'''     
    if filename.endswith('.xml'): 
        sources = read_sofia_xml(filename,variables=variables,verbose=verbose)
    else:
        sources = read_sofia_txt(filename,variables=variables,verbose=verbose)

    sources = check_parameters(sources,variables=variables,no_conversion=no_conversion) 
       
    if verbose:
        print(f'We found this table for the file {filename}')
        print(sources)
    return sources

def load_sofia_cat_header(lines):
    '''Read the header of the sofia catalogue file'''
    columns_triggered = False
    row_start = 0
    input_columns = []
    column_locations = []
    for line in lines:
        tmp =line.split()
        # If it is a comment line we increase row start
        if line.startswith('#'):
            row_start += 1
        else:
            # No comment means we have passed the header
            break
            # check if the line is empty or a comment
        if line.strip() == '' or line.strip() == '#':
            pass
        elif tmp[0] == '#' and len(tmp) > 1:
            if tmp[1].strip().lower() in ['name','id']:
                # get the present columns
                input_columns  = [x.strip() for x in tmp[1:]]
               
                #determin their location in the line
                column_locations = []
                for col in input_columns:
                    column_locations.append(line.find(col)+len(col))
                # check that we found all parameters
                
                columns_triggered = True
                
            elif columns_triggered:
                columns_triggered = False
                
                input_units  = [x.strip().lower() for x in tmp[1:]]
                convert_units = []
                for unit in input_units:
                    convert_units.append(translate_string_to_unit(unit))
                for i,column in enumerate(input_columns):
                    if column in ['id','name']:
                        convert_units[i] = 'str'

                dtypes = []
                for unit in convert_units:
                    if unit == 'str':
                        dtypes.append(str)
                    elif unit == u.pix:
                        dtypes.append(float)
                    else:
                        dtypes.append(np.float64)
                
                sources = QTable(names=input_columns,units=convert_units,dtype=dtypes)
            
                    
    return row_start, input_columns, column_locations,convert_units,sources


def load_sofia_input_file(filename='Template'):
    if filename == 'Template':
        with pack_open_txt(templates, 'sofia_template.par') as tmp:
            template = tmp.readlines()
    else:
        with open(filename,'r') as tmp:
            template = tmp.readlines()
    result = {}
    counter = 0
    counter2 = 0
    # Separate the keyword names
    for line in template:
        key = str(line.split('=')[0].strip())
        if key == '':
            result[f"EMPTY{counter}"] = line
            counter += 1
        elif key[0] == '#':
            result[f"HASH{counter2}"] = line
            counter2 += 1
        else:
            result[key] = str(line.split('=')[1].strip())
    return result
load_sofia_input_file.__doc__ ='''
 NAME:
    load_sofia_input_file(filename='Template')

 PURPOSE:
    Read the sofia2 file in Templates into a dictionary

 CATEGORY:
    read_functions

 INPUTS:

 OPTIONAL INPUTS:


 OUTPUTS:
    result = dictionary with the read file

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    split, strip, open

 NOTE:
'''

 
def move_sources(cfg,indir,old_new_ids,originalbasename,basename,original_id,base_dir=''):
    """
    Move the sources to the original directory and update the IDs.
    
    Parameters:
    cfg (Config): The configuration object.
    indir (str): The input directory containing the sources.
    old_new_ids (dict): A dictionary mapping old IDs to new IDs.
    """
    to_move =['_chan.fits','_cube.fits','_mask.fits','_mom0.fits','_mom1.fits',
            '_mom2.fits', '_snr.fits','_spec.txt','_spec_aperture.txt']
    # We need to move the files to the original directory and update the IDs
    for new_id in old_new_ids:
        for g in to_move:
            old_name = f'{indir}{basename}_cubelets/{basename}_{old_new_ids[new_id]}{g}'
            if os.path.exists(old_name):
                new_name = f'{cfg.internal.sofia_directory}/{originalbasename}_cubelets/{originalbasename}_{new_id}{g}'
                #print(f'Moving {old_name} to {new_name}')
                os.rename(old_name, new_name)        
    #And remove the original unsplit source
    original_source_file = f'{cfg.internal.sofia_directory}/{originalbasename}_cubelets/{originalbasename}_{original_id}'
    for g in to_move:
        file_to_remove = f'{original_source_file}{g}'
        #print(f'Trying to remove {file_to_remove}')
        if os.path.exists(file_to_remove):
            #print(f'Removing {file_to_remove}')
            os.remove(file_to_remove)

    #move the watershed output 
    oldname = f'{cfg.internal.sofia_directory}/Watershed_Output/'
    newname = f'{cfg.internal.sofia_directory}/{originalbasename}_cubelets/Watershed_Output_{original_id}/'
    shutil.rmtree(newname) if os.path.exists(newname) else None
    os.rename(oldname, newname)      


def obtain_sofia_id(base_name, cube_name):
    tmp,cube_file = os.path.split(cube_name)
    split_main = cube_file.split(base_name)
    parts = split_main[1].split('_')
    id  = int(parts[1])
    return id,cube_file

def read_sofia_table(cfg,sofia_directory='./',sofia_basename=None,
        no_conversion=False,
        force_text = False):
    '''Locate and read a sofia catolgue into an astropy table '''
    if not os.path.isfile(f'{sofia_directory}/{sofia_basename}_cat.xml')\
        and not os.path.isfile(f'{sofia_directory}/{sofia_basename}_cat.txt'):
        print(f'''No sofia table found in {sofia_directory}/{sofia_basename}_cat.
probably no sources were found or you made an input error. ''')
        return None,sofia_basename,None

    if os.path.isfile(f'{sofia_directory}/{sofia_basename}_cat.xml') and not force_text:
        #If we have an xml we prefer to read that
        table_name = f'{sofia_directory}/{sofia_basename}_cat.xml'
    elif os.path.isfile(f'{sofia_directory}/{sofia_basename}_cat.txt'):
        table_name = f'{sofia_directory}/{sofia_basename}_cat.txt'
    else:
        print(f'''No sofia table found in {sofia_directory} for {cfg.internal.data_cube}.
probably no sources were found.''')
        return None,sofia_basename,None
    req_variables = ['name','f_sum','err_f_sum','id','ell3s_maj',
        'ell3s_min','w20','ra','dec','v_sofia','kin_pa','x','y','z',
        'f_max','ell_maj','ell_min','rms','ell_pa']

    sources = load_sofia_catalogue(table_name,verbose=cfg.general.verbose,
            no_conversion=no_conversion,variables= req_variables) 
    if not no_conversion:
        if cfg.internal.data_cube is None:
            raise InputError('cubename must be specified if no_conversion is True')
        sources = convert_pix_columns_to_arcsec(cfg,sources,f'{cfg.internal.data_directory}/{cfg.internal.data_cube}')

    return sources,sofia_basename,table_name




def read_sofia_xml(filename,variables=None,verbose=False):
    '''Read the sofia xml file into a Astropy QTable'''
  
    if verbose:
        print(f'Reading the sofia catalogue {filename}')
    # Read the xml file
    table = QTable(votable.parse(filename).get_first_table().to_table())
  
    #sources = QTable.read(filename, format='xml')
    actual_colnames = [x for x in table.colnames]
    lower_colnames = [x.lower() for x in actual_colnames]
    
    if 'id' in lower_colnames:
        if not isinstance(table[actual_colnames[lower_colnames.index('id')]][0],(str,np.str_)):
            print(type(table[actual_colnames[lower_colnames.index('id')]][0]))    
            table.add_column([f'{int(x)}' for x in 
                table[actual_colnames[lower_colnames.index('id')]]], name='str_id')
            table.remove_column(actual_colnames[lower_colnames.index('id')])
            table.rename_column('str_id',actual_colnames[lower_colnames.index('id')])
    for col in actual_colnames:
        if isinstance(table[col].unit,(u.UnrecognizedUnit)):
            test = translate_string_to_unit(table[col].unit.to_string())
            table['tmp'] = [x.value for x in table[col]]* test 
            table.remove_column(col)
            table.rename_column('tmp',col)
          
       #if isinstance(unit, str):
        #    table['unit'][table['unit'] == unit] = translate_string_to_unit(unit)       

    #Remove unwanted 
    if verbose:
        print(f'We found the following columns in the sofia catalogue: {actual_colnames}')
    return table


def read_sofia_txt(filename,variables=None,verbose=False):
    with open(filename) as tmp:
        lines = tmp.readlines()
    if variables is None:
        for line in lines:
            tmp= line.split()
            if len(tmp) > 1:
                if tmp[1].strip().lower() in ['name','id']:
                    # get the present columns
                    variables  = [x.strip() for x in tmp[1:]]
                    break    
    row_start, input_columns, column_locations,convert_units,sources = \
        load_sofia_cat_header(lines) 
    if verbose:
        print(f'We found the following columns in the sofia catalogue: {input_columns}')
    for line in lines[row_start:]: 
        if line.strip() == '':
            continue
        #to the table      
        construct_row=[]
        name = get_source_cat_name(line,input_columns,column_locations)
        for i,col in enumerate(input_columns):
            if col == 'name':
                construct_row.append(name)
            else:
                start,end = get_start_end_locations(i,column_locations)
                if col in ['id']:
                    construct_row.append(line[start:end].strip())
                else:
                    #print(f'col: {col} start: {start} end: {end} value: {line[start:end]}')
                    construct_row.append(float(line[start:end].strip())\
                        * convert_units[i])
        sources.add_row(construct_row)
    return sources

def set_sofia(sofia_temp, cube_name, mask, outdir):
    """
    Set the SoFiA parameters for the deblending process.
    
    Parameters:
    sofia_temp (dict): The SoFiA template parameters.
    cube_name (str): The name of the data cube file.
    mask (str): The name of the mask file.
    outdir (str): The directory to save the results.
    
    Returns:
    dict: The updated SoFiA parameters.
    """
    sofia_temp['input.data'] = f'{cube_name}'
    sofia_temp['input.mask'] = f'{mask}'
    sofia_temp['output.directory'] = f'{outdir}/Sofia_Output/' 
    sofia_temp['pipeline.verbose'] = 'true'
    sofia_temp['scfind.enable'] = 'false'
    sofia_temp['linker.enable'] = 'false'
    sofia_temp['reliability.enable'] = 'false'
    sofia_temp['dilation.enable'] = 'false'
    
    return sofia_temp


def update_sofia_catalogue(cfg, cube_name= None,base_name = None, outdir='./',base_dir =''):
    
    #First get the name and id
    if base_name is None:
        base_name = os.path.splitext(os.path.basename(cube_name))[0]
    if cfg.general.verbose:
        print(f'Main name is {cube_name}')
    sofia_id,cube_file_name = obtain_sofia_id(base_name, cube_name) 
    
    #load the original sofia table
 
    sources,sofia_basename,table_name = read_sofia_table(cfg, 
        sofia_directory=cfg.internal.sofia_directory,
        sofia_basename=cfg.internal.sofia_basename,
        no_conversion=True)
    
    #load the split sources
    
    split_sources,split_base_name,split_table_name =  read_sofia_table(cfg, 
        sofia_directory=f'{outdir}/Sofia_Output',
        sofia_basename=os.path.splitext(cube_file_name)[0],no_conversion=True) 
    print(split_sources,outdir,split_base_name,os.path.splitext(cube_file_name)[0])
    # We need to update pixel values to the original cube
  
    split_sources = convert_pixel_values_to_original(split_sources,
                cube_name,f'{cfg.internal.data_directory}/{cfg.internal.data_cube}')
    
    alphabet = list(string.ascii_lowercase)
    split_sources['id'] = split_sources['id'].astype('<U3')
    sources['id'] = sources['id'].astype('<U3')
    for i,row in enumerate(sources):
        if row['id'] == f'{sofia_id}':
            to_remove = i
            break
    if cfg.general.verbose:
        print(f'Removing source {sofia_id} from the sources table at index {to_remove}')
    # Remove the original source from the sources table
    sources.remove_rows(to_remove)
    
    old_new_ids = {}
    for i in range(len(split_sources)):
        if cfg.general.verbose:
            print(f'Splitting source {split_sources["id"][i]} with id {sofia_id}{alphabet[i]}')
        
        old_new_ids[f'{sofia_id}{alphabet[i]}'] = split_sources['id'][i]
        split_sources['id'][i] = f'{sofia_id}{alphabet[i]}'
       
        # Add the split sources to the original sources
        sources.add_row(split_sources[i])
    sources.sort('id')    
  
    #Write the new table to over the old one
    write_sofia_catalogue_xml(sources,f'{cfg.internal.sofia_directory}/{sofia_basename}_cat.xml')
    write_sofia_catalogue_txt(sources,f'{cfg.internal.sofia_directory}/{sofia_basename}_cat.txt')


    # move the split sources to the original directory
    move_sources(cfg,f'{outdir}Sofia_Output/',old_new_ids,sofia_basename
                 ,split_base_name,sofia_id,base_dir=base_dir)



def write_sofia(template,name):
    with open(name,'w') as file:
        for key in template:
            if key[0] == 'E' or key [0] == 'H':
                file.write(template[key])
            else:
                file.write(f"{key} = {template[key]}\n")

write_sofia.__doc__ =f'''
NAME:
sofia
PURPOSE:
write a sofia2 dictionary into file
CATEGORY:
write_functions

INPUTS:
template = sofia template
name = name of the file to write to

OPTIONAL INPUTS:

OUTPUTS:

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
Unspecified

NOTE:
'''


def write_sofia_catalogue_txt(table,filename):
    """
    Update a sofia file format with a table 
    
    Parameters
    ----------
    table : Table
        The table to write.
    filename : str
        The name of the file to write to.
    """
    # Read the current lines from the 
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        tmp= line.split()
        if len(tmp) > 1:
            if tmp[1].strip().lower() in ['name','id']:
                # get the present columns
                variables  = [x.strip() for x in tmp[1:]]
                break     
    
    # obtain the header
    row_start, input_columns, column_locations,convert_units,sources = \
        load_sofia_cat_header(lines) 
    with open(filename, 'w') as f:
        # Write the header
        for line in lines:
            if line.startswith('#'):
                f.write(line)
            else:
                break
        f.write(' \n')
        for row_correct in table: 
          
            line  = ''
            for i,col in enumerate(input_columns):
              
                if col.lower() in [x.lower() for x in row_correct.colnames]:
                    value = row_correct[col]
                    start,end = get_start_end_locations(i, column_locations)
                    #if start == 0:
                    #    start = 1
                    length = int(end - start)-1
                 
                    if isinstance(value, u.Quantity):
                        line += f"{value.value:>{length}.6f}"
                    elif isinstance(value,float):
                        line += f"{value:>{length}.6f}"
                    elif isinstance(value,int):
                        line += f"{value:>{length}d}"
                    elif isinstance(value,str):
                        line += f"{value:>{length}s}"
            f.write(line + '\n')

def write_sofia_catalogue_xml(table,filename):
    """
    Write a sofia catalogue to an xml file.
    
    Parameters
    ----------
    table : Table
        The table to write.
    filename : str
        The name of the file to write to.
    """
    output_votable = votable.from_table(table)
    votable.writeto(output_votable,filename)