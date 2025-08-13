from deblend_sofia_detections.config.config import defaults
from deblend_sofia_detections.support.errors import InputError
from deblend_sofia_detections.support.system_functions import create_directory
from deblend_sofia_detections.deblending.sofia_functions import load_sofia_input_file
from omegaconf import OmegaConf

import os
import psutil
import sys
import deblend_sofia_detections
try:
    from importlib.resources import files as import_pack_files
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    # For Py<3.9 files is not available
    from importlib_resources import files as import_pack_files
    

def setup_config(argv):
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
                    ['input','general'])
           
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
    cfg = read_parameter_input(cfg)

    cfg = directory_check(cfg)    
    return cfg


def directory_check(cfg):

    for test_dir in [cfg.internal.data_directory, cfg.internal.sofia_directory, cfg.internal.run_directory]:
        if not os.path.isdir(test_dir):
            raise InputError(f'''The directory {test_dir} does not exist.''')
    return cfg


def read_parameter_input(cfg):
    cfg.internal.run_directory = os.getcwd()
    parameters = load_sofia_input_file(cfg.input.sofia_parameters)
    input_pathname,parameter_file = os.path.split(cfg.input.sofia_parameters)
    cfg.internal.sofia_parameter_file = parameter_file
    if input_pathname == '' or input_pathname[0] != '/':
        input_pathname = os.path.join(os.getcwd(),input_pathname)
    data_path,data_file = os.path.split(parameters['input.data'])
    if data_path == '' or data_path[0] != '/':
        cfg.internal.data_directory = f'{input_pathname}{data_path}'
    else:
        cfg.internal.data_directory = data_path
    cfg.internal.ancillary_directory = f'{cfg.internal.data_directory}/ancillary_data'
    cfg.internal.data_cube = data_file
    cfg.internal.sofia_basename = parameters['output.filename']
    if cfg.internal.sofia_basename == '':
        cfg.internal.sofia_basename = os.path.splitext(data_file)[0]
    cfg.internal.sofia_directory = parameters['output.directory']
    if cfg.internal.sofia_directory == '' or cfg.internal.sofia_directory[0] != '/':
        cfg.internal.sofia_directory = os.path.join(input_pathname,cfg.internal.sofia_directory)
    
    return cfg