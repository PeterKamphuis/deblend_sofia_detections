from deblend_sofia_detections.config.config import defaults
from deblend_sofia_detections.support.errors import InputError
from deblend_sofia_detections.support.logging import print_log,start_new_log
from deblend_sofia_detections.support.system_functions import join_path,create_directory
from deblend_sofia_detections.deblending.sofia_functions import load_sofia_input_file

from omegaconf import OmegaConf,MISSING
from astropy import units as u

import os
import psutil
import shutil
import sys
import copy
import deblend_sofia_detections
import numpy as np
import warnings
try:
    from importlib.resources import files as import_pack_files
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    # For Py<3.9 files is not available
    from importlib_resources import files as import_pack_files
    

def setup_config(argv,single_cube=False):
    if '-v' in argv or '--version' in argv:
        print(f"This is version {deblend_sofia_detections.__version__} of the program.")
        sys.exit()
   
    if '-h' in argv or '--help' in argv or len(argv) == 0:
        print('''
Use package_name in this way:

All config parameters can be set directly from the command line by setting the correct parameters, e.g:
deblend configuration_file=<Your_Input.yml>
''')
        sys.exit()
   

    cfg = OmegaConf.structured(defaults)
    if cfg.general.ncpu == psutil.cpu_count():
        cfg.general.ncpu -= 1
    inputconf = OmegaConf.from_cli(argv)
    cfg_input = OmegaConf.merge(cfg,inputconf)
    
    if cfg_input.print_examples:
        default_name = f'{__name__.split(".")[0]}_default.yml' 
        masked_copy = OmegaConf.masked_copy(cfg,\
                    ['input','sofia','logging','general','directories'])
           
        with open(default_name,'w') as default_write:
            default_write.write(OmegaConf.to_yaml(masked_copy))
        print(f'''We have printed the file {default_name} in {os.getcwd()}.
Exiting {__name__.split(".")[0]}.''')
        my_resources = import_pack_files('deblend_sofia_detections.template')
        data = (my_resources / 'sofia_template.par').read_bytes()
        with open('sofia_template.par','w+b') as default_write:
            default_write.write(data)
        print(f'''We have printed the file  sofia_template.par in {os.getcwd()}.
''')
        sys.exit()
        

    if cfg_input.configuration_file:
        succes = False
        while not succes:
            try:
                yaml_config = OmegaConf.load(cfg_input.configuration_file)
        #merge yml file with defaults
                cfg = OmegaConf.merge(cfg,yaml_config)
                succes = True
            except FileNotFoundError:
                cfg_input.configuration_file = input(f'''
You have provided a config file ({cfg_input.configuration_file}) but it can't be found.
If you want to provide a config file please give the correct name.
Else press CTRL-C to abort.
configuration_file = ''')
    cfg = OmegaConf.merge(cfg,inputconf) 
   
    #open the input parameter file to obtain the data cube and output locations
  
    if single_cube:
        cfg = check_single_cube_input(cfg)
    else:
        cfg = read_parameter_input(cfg)

    
    cfg = directory_check(cfg)
    if cfg.directories.run_directory != os.getcwd():
        os.chdir(cfg.directories.run_directory)
    cfg = background_check(cfg)
    cfg = check_debug_functions(cfg)
    cfg = check_all_required_inputs(cfg)
    cfg = check_internet(cfg)
  
   
    return cfg


def check_all_required_inputs(cfg):
    cfg.input.markers = [x.lower() for x in cfg.input.markers]
    markers = set(cfg.input.markers )
    allowed_markers = ['internet table', 'manual table', 'marker image', 
        'optical image', 'internet image', 'hi peaks','moment0 peaks','cleaned optical image']
    image_count = 0
    for marker in markers:
        tmp = marker.split()
        if marker not in allowed_markers:
            if 'peaks' in marker:
                
                if tmp[1] not in ['hi', 'moment0']:
                    raise InputError(f"Invalid peak type: {tmp[1]}")
                try:
                    num_peaks = int(tmp[0])
                    cfg.internal.maximum_no_peaks = int(np.nanmax([num_peaks, cfg.internal.maximum_no_peaks]))
                except ValueError:
                    raise InputError(f"Invalid format for 'peaks' marker: {marker}")
            else:
                raise InputError(f"Marker '{marker}' is not allowed. Allowed markers are: {allowed_markers}")
        else:
            if len(tmp) == 2:
                if tmp[1] == 'image' and tmp[0] != 'marker':
                    image_count += 1
            if image_count > 1:
                raise InputError(f"Multiple image markers found. Please pick only 1 of 'internet image', 'optical image', 'cleaned optical image'.")
    if 'manual table' in  markers and cfg.input.tables[0] is None:
        print(f"You requested manual markers but have not provided a table.")
        raise InputError(f"You requested manual markers but have not provided a table.")
    
    if 'optical image' in markers and cfg.input.optical_images[0] is None:
        print(f"You requested optical image markers but have not provided an optical image.")
        raise InputError(f"You requested optical image markers but have not provided an optical image.")

    ch_markers = set(['cleaned optical image','optical image'])
    if len(markers.intersection(ch_markers)) > 1:
        print(f"You requested both cleaned optical image and optical image markers. Please choose 1 else we get confused.")
        raise InputError(f"You requested both cleaned optical image and optical image markers. Please choose 1.")
    
    if 'marker image' in markers and cfg.input.marker_images[0] is None:
        print(f"You requested marker image markers but have not provided a marker image.")
        raise InputError(f"You requested marker image markers but have not provided a marker image.")
    ch_markers = set(['cleaned optical image','optical image',
                       'internet image','internet table'])
    if len(markers.intersection(ch_markers)) > 0:
        if cfg.general.counterpart_region == 'No Counterpart':
            print(f"{''.join(['!!!'*10])}\n Warning: Automated markers requested, but counterpart region is set to 'No Counterpart'.\n {''.join(['!!!'*10])}")
           


      
    return cfg
    

def check_debug_functions(cfg):
    cfg.logging.debug_functions = [x.upper() for x in list(cfg.logging.debug_functions)]
    if 'ALL' in cfg.logging.debug_functions:
        cfg.logging.debug_functions = ['ALL']
    if 'NONE' in cfg.logging.debug_functions:
        cfg.logging.debug_functions = ['NONE']

   
    return cfg

def check_internet(cfg):

    download_options = ['gaia','ned','simbad']
    to_download = {}
    for x in cfg.general.internet_repository:
        if x.lower() not in download_options:
            print(f"Warning: {x} is not a recognized download option.")
        else:
            to_download[x.lower()] = True
   
    if not cfg.general.internet:
        print(f"Internet access is disabled in the configuration.")
        for x in to_download:
            to_download[x] = False

    for x in to_download:
        if x == 'gaia':
            if 'cleaned optical image' in cfg.input.markers:
                to_download[x] = False
            continue
        if x.upper() not in cfg.general.internet_repository:
            to_download[x] = False
        if cfg.general.counterpart_region == 'No Counterpart':
            to_download[x] = False

    for x in to_download:
        if not to_download[x]:        
        #check if we already have table
            if os.path.exists(f'{cfg.directories.ancillary_directory}/tables/cached_{x}_table.pkl'):
                setattr(cfg.internal, f'{x}_table', \
                    f'{cfg.directories.ancillary_directory}/tables/cached_{x}_table.pkl')
            else:
                print(f"No cached table found for {x}.")
                setattr(cfg.internal, f'{x}_table', 'none')
                cfg.general.internet_repository.remove(x.upper())
    return cfg
     


def check_single_cube_input(cfg):
    # If we have a parameter file we read it normally but in this it is not required 
    if os.path.isfile(cfg.input.sofia_parameters):
        cfg = read_parameter_input(cfg)
    else:
        if not cfg.sofia.original_data_cube:
            raise InputError(f'''You have not provided a data cube to process.''')
        if os.path.split(cfg.sofia.original_data_cube)[0] == '':
            if cfg.directories.data_directory == '':
                cfg.directories.data_directory = os.getcwd()
        else:
            if cfg.directories.data_directory == '':
                cfg.directories.data_directory = os.path.split(cfg.sofia.original_data_cube)[0]
          
            cfg.sofia.original_data_cube = os.path.split(cfg.sofia.original_data_cube)[1]  
        # check that the cube exists
        if not os.path.isfile(join_path(cfg.directories.data_directory,cfg.sofia.original_data_cube)):
            raise InputError(f'''The data cube {cfg.sofia.original_data_cube} does not exist in the directory {cfg.directories.data_directory}. Please provide a correct data cube and directory.''')   
     
        if cfg.sofia.basename == '':
            cfg.sofia.basename = os.path.splitext(cfg.sofia.original_data_cube)[0].split('_cube')[0]
        # in this case the dta_directory and the sofia_directory have to be hthe same    
        if cfg.sofia.directory == '':
            cfg.sofia.directory = cfg.directories.data_directory
        cube_ext= os.path.splitext(cfg.sofia.original_data_cube)[1]
        if cube_ext.lower() in ['.gz']:
            cube_ext = os.path.splitext(os.path.splitext(cfg.sofia.original_data_cube)[0])[1]+cube_ext

        cfg.internal.cube_ext = cube_ext
        
        if cfg.sofia.original_mask == '':
            cfg.sofia.original_mask = f'{cfg.sofia.basename}_mask{cube_ext}'
        elif os.path.split(cfg.sofia.original_mask)[0] != '':
            # in case the user provides an absolute path
            cfg.sofia.original_mask = os.path.split(cfg.sofia.original_mask)[1]

        if not os.path.isfile(join_path(cfg.sofia.directory,cfg.sofia.original_mask)):
            raise InputError(f'''The mask file {cfg.sofia.original_mask} does not exist in the directory {cfg.sofia.directory}. Please provide a correct mask file.''')
    # As we have no way to reproduce the mask in this case we do not want to modify it.
    file,ext = os.path.splitext(cfg.sofia.original_mask)
    if ext.lower() in ['.gz']:
        file,ext2 = os.path.splitext(file)
        ext = ext2 + ext
    use_mask = f'{file}_deblended{ext}'
    print(f'Using mask file {use_mask}')
    print(f'Original mask file {cfg.sofia.original_mask}')
    print(f'directory {cfg.sofia.directory}')
    print(join_path(cfg.sofia.directory,cfg.sofia.original_mask),join_path(cfg.sofia.directory,use_mask))
    shutil.copyfile(join_path(cfg.sofia.directory,cfg.sofia.original_mask),join_path(cfg.sofia.directory,use_mask))
    cfg.sofia.original_mask = use_mask
    return cfg

def directory_check(cfg):
    dirs = ['data_directory', 'run_directory', 'ancillary_directory', 
        'watershed_directory']
    
    if cfg.sofia.directory[-1] != '/':
        cfg.sofia.directory += '/'
    # make full paths for defaults
    directories_to_check = [cfg.sofia.directory]
    directories_to_create= []   
    if cfg.directories.ancillary_directory == 'ancillary_data':
        cfg.directories.ancillary_directory = join_path(cfg.directories.data_directory,
            cfg.directories.ancillary_directory)
        directories_to_create.append(cfg.directories.ancillary_directory)
    if cfg.directories.watershed_directory == 'Watershed_Output':
        cfg.directories.watershed_directory = join_path(cfg.sofia.directory,
            cfg.directories.watershed_directory)
        directories_to_create.append(cfg.directories.watershed_directory)
   
    cfg.internal.input_log_directory = copy.deepcopy(cfg.logging.log_directory)
    cfg.internal.input_log_file = copy.deepcopy(cfg.logging.log_file)    
    for attr in dirs:
        test_dir = getattr(cfg.directories, attr)
        directories_to_check.append(test_dir)
        if test_dir[-1] != '/':
            test_dir += '/'
            setattr(cfg.directories , attr, test_dir)
    markers = set([x.lower() for x in cfg.input.markers])
    need_internet = False
    if 'internet table' in markers or 'internet image' in markers or 'optical image' in markers:
        need_internet = True
        directories_to_check.append(f'{cfg.directories.ancillary_directory}/tables/')
        directories_to_create.append(f'{cfg.directories.ancillary_directory}/tables/')

    if need_internet and not cfg.general.internet and cfg.general.clear_internet_cache:
        raise InputError("Internet markers requested, but internet access is not available and you want to clear the caches.")
    for test_dir in directories_to_check:
        if not os.path.isdir(test_dir):
            if test_dir in directories_to_create:
                create_directory(test_dir)
            else:
                raise InputError(f"The directory {test_dir} does not exist. Please provide a correct directory.")  
    internet_caches = ['gaia','ned','simbad']
    start_new_log(cfg,basedir=cfg.directories.watershed_directory)
    print_log(cfg, f'''Checked the directories:
        {cfg.directories.data_directory}   
        {cfg.directories.run_directory}
        {cfg.directories.ancillary_directory}
        {cfg.sofia.directory}
    ''')
    for cache in internet_caches:
        if os.path.isfile(f'{cfg.directories.ancillary_directory}/tables/cached_{cache}_table.pkl'):
            if cfg.general.clear_internet_cache:
                os.remove(f'{cfg.directories.ancillary_directory}/tables/cached_{cache}_table.pkl')
            else:
                print_log(cfg, f"Using existing cached {cache} table.", case=['verbose'])
                setattr(cfg.internal, f'{cache}_table', f'{cfg.directories.ancillary_directory}/tables/cached_{cache}_table.pkl')
        else:
            print_log(cfg, f"No cached {cache} table found.", case=['verbose'])
            if not cfg.general.internet and need_internet:
                raise InputError(f"Internet markers requested, but we have no  internet access  and no cached table for {cache}.")
            setattr(cfg.internal, f'{cache}_table', 'none')

  
    
    
    return cfg


def read_parameter_input(cfg):
    
    input_pathname,parameter_file = os.path.split(cfg.input.sofia_parameters)
    cfg.sofia.parameter_file = parameter_file
    if input_pathname == '' or input_pathname[0] != '/':
        input_pathname = join_path(cfg.directories.run_directory,input_pathname)
    cfg.sofia.parameter_path = input_pathname
    if os.path.isfile(f'{cfg.sofia.parameter_path}/{cfg.sofia.parameter_file}') == False:
        print_log(cfg, f"The input parameter file {cfg.sofia.parameter_path}/{cfg.sofia.parameter_file} does not exist.",
                  case=['verbose'])
        cfg.sofia.parameter_path = os.getcwd()
        print_log(cfg, f"Trying {cfg.sofia.parameter_path}/{cfg.sofia.parameter_file} instead.",
                    case=['verbose'])
        
    parameters = load_sofia_input_file(f'{cfg.sofia.parameter_path}/{cfg.sofia.parameter_file}')
    data_path,data_file = os.path.split(parameters['input.data'])
    if cfg.directories.data_directory == '':
        if data_path == '' or data_path[0] != '/':
            cfg.directories.data_directory = join_path(
                input_pathname,data_path)
        else:
            cfg.directories.data_directory = join_path(data_path)
    else:
        warnings.warn(f'''{''.join(['-']*30)} Warning! {''.join(['-']*30)}
Data directory is set to: {cfg.directories.data_directory} by user.
Not applying automatic data directory setting from parameter file.
{''.join(['-']*80)}\n''')
   
    if cfg.sofia.original_data_cube == '':
        cfg.sofia.original_data_cube = data_file
    else:
        warnings.warn(f'''{''.join(['-']*30)} Warning! {''.join(['-']*30)}
Original data cube is set to: {cfg.sofia.original_data_cube} by user.
Not applying automatic original data cube setting from parameter file.
{''.join(['-']*80)}\n''')
    if cfg.sofia.basename == '':
        cfg.sofia.basename = parameters['output.filename']
        if cfg.sofia.basename == '':
            cfg.sofia.basename = os.path.splitext(data_file)[0]
    else:
        warnings.warn(f'''{''.join(['-']*30)} Warning! {''.join(['-']*30)}
Basename is set to: {cfg.sofia.basename} by user.
Not applying automatic basename setting from parameter file.
{''.join(['-']*80)}''')

    if cfg.sofia.directory == '':
        cfg.sofia.directory = join_path( cfg.directories.data_directory 
                                        , parameters['output.directory'])
    else:
        warnings.warn(f'''{''.join(['-']*30)} Warning! {''.join(['-']*30)}
Sofia directory is set to: {cfg.sofia.directory} by user.
Not applying automatic directory setting from parameter file.
{''.join(['-']*80)}\n''')    

    if cfg.sofia.original_mask == '':
        cfg.sofia.original_mask = f'{cfg.sofia.basename}_mask.fits'
    else:
        warnings.warn(f'''{''.join(['-']*30)} Warning! {''.join(['-']*30)}
Original mask is set to: {cfg.sofia.original_mask} by user.
Not applying automatic original mask setting from parameter file.
{''.join(['-']*80)}\n''')
    
    
    return cfg

def background_check(cfg):
    if cfg.input.optical_images[0] is not None:
        cfg.internal.optical_background = f'{cfg.directories.ancillary_directory}/Manual_Optical_Background.fits'
        cfg.internal.cleaned_optical_background = f'{cfg.directories.ancillary_directory}/Cleaned_Manual_Optical_Background.fits'
    else:
        cfg.internal.optical_background = f'{cfg.directories.ancillary_directory}/DSS_Optical_Background.fits'
        cfg.internal.cleaned_optical_background = f'{cfg.directories.ancillary_directory}/Cleaned_DSS_Optical_Background.fits'
    # if we want original images we delete the processed stuf if existing
    if os.path.isfile(cfg.internal.optical_background) and cfg.general.originals:
        if cfg.logging.enable:
            print(f'Deleting existing optical background image: {cfg.internal.optical_background}')
        os.remove(cfg.internal.optical_background)
    
    if cfg.general.originals:
        if cfg.logging.enable:
            print(f'Deleting any existing cleaned optical background image: {cfg.internal.cleaned_optical_background}')
        path,name = os.path.split(cfg.internal.cleaned_optical_background)
        basename = os.path.splitext(name)[0]
        for file in os.listdir(path):
            if file.startswith(basename):
                if cfg.logging.enable:
                    print(f'Deleting existing cleaned optical background image: {join_path(path,file)}')
                os.remove(join_path(path,file))
     
    
    return cfg