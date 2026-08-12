from deblend_sofia_detections.config.config import defaults
from deblend_sofia_detections.support.errors import InputError
from deblend_sofia_detections.support.logging import print_log,start_new_log
from deblend_sofia_detections.support.system_functions import join_path,create_directory
from deblend_sofia_detections.deblending.sofia_functions import load_sofia_input_file
from omegaconf import OmegaConf,MISSING

import os
import psutil
import shutil
import sys
import copy
import deblend_sofia_detections
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

    if '-h' in argv or '--help' in argv:
        print('''
Use package_name in this way:

All config parameters can be set directly from the command line by setting the correct parameters, e.g:
create_package_name def_file=cube.fits error_generator=tirshaker 
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
   
    return cfg

def check_debug_functions(cfg):
    cfg.logging.debug_functions = [x.upper() for x in list(cfg.logging.debug_functions)]
    if 'ALL' in cfg.logging.debug_functions:
        cfg.logging.debug_functions = ['ALL']
    if 'NONE' in cfg.logging.debug_functions:
        cfg.logging.debug_functions = ['NONE']

   
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
        print()
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
        if not os.path.isfile(join_path(cfg.sofia.directory,cfg.sofia.original_mask)):
            raise InputError(f'''The mask file {cfg.sofia.original_mask} does not exist in the directory {cfg.sofia.directory}. Please provide a correct mask file.''')
    # As we have no way to reproduce the mask in this case we do not want to modify it.
    file,ext = os.path.splitext(cfg.sofia.original_mask)
    if ext.lower() in ['.gz']:
        file,ext2 = os.path.splitext(file)
        ext = ext2 + ext
    use_mask = f'{file}_deblended{ext}'
    shutil.copyfile(join_path(cfg.sofia.directory,cfg.sofia.original_mask),join_path(cfg.sofia.directory,use_mask))
    cfg.sofia.original_mask = use_mask
    return cfg

def directory_check(cfg):
    dirs = ['data_directory', 'run_directory', 'ancillary_directory', 'watershed_directory']
    
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
    if cfg.input.internet_query.lower() != 'none':
        directories_to_check.append(f'{cfg.directories.ancillary_directory}/tables/')
        directories_to_create.append(f'{cfg.directories.ancillary_directory}/tables/')
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
            if cfg.input.clear_internet_cache:
                os.remove(f'{cfg.directories.ancillary_directory}/tables/cached_{cache}_table.pkl')
            else:
                print_log(cfg, f"Using existing cached {cache} table.", case=['verbose'])
                setattr(cfg.internal, f'{cache}_table', f'{cfg.directories.ancillary_directory}/tables/cached_{cache}_table.pkl')
        else:
            print_log(cfg, f"No cached {cache} table found.", case=['verbose'])
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
    if data_path == '' or data_path[0] != '/':
        cfg.directories.data_directory = join_path(
            input_pathname,data_path)
    else:
        cfg.directories.data_directory = join_path(data_path)
   
    cfg.sofia.original_data_cube = data_file
    cfg.sofia.basename = parameters['output.filename']
    if cfg.sofia.basename == '':
        cfg.sofia.basename = os.path.splitext(data_file)[0]
    cfg.sofia.directory = join_path( cfg.directories.data_directory 
                                    , parameters['output.directory'])
    
    return cfg

def background_check(cfg):
    if cfg.input.manual_optical_image[0] is not None:
        cfg.internal.optical_background = f'{cfg.directories.ancillary_directory}/Manual_Optical_Background.fits'
        cfg.internal.cleaned_optical_background = f'{cfg.directories.ancillary_directory}/Cleaned_Manual_Optical_Background.fits'
    else:
        cfg.internal.optical_background = f'{cfg.directories.ancillary_directory}/DSS_Optical_Background.fits'
        cfg.internal.cleaned_optical_background = f'{cfg.directories.ancillary_directory}/Cleaned_DSS_Optical_Background.fits'
    # if we want original images we delete the processed stuf if existing
    if os.path.isfile(cfg.internal.optical_background) and cfg.input.original_images:
        if cfg.general.verbose:
            print(f'Deleting existing optical background image: {cfg.internal.optical_background}')
        os.remove(cfg.internal.optical_background)
    
    if cfg.input.original_images:
        if cfg.general.verbose:
            print(f'Deleting any existing cleaned optical background image: {cfg.internal.cleaned_optical_background}')
        path,name = os.path.split(cfg.internal.cleaned_optical_background)
        basename = os.path.splitext(name)[0]
        for file in os.listdir(path):
            if file.startswith(basename):
                if cfg.general.verbose:
                    print(f'Deleting existing cleaned optical background image: {join_path(path,file)}')
                os.remove(join_path(path,file))
     
    
    return cfg