# OmegaConf setups

from dataclasses import dataclass, field
from omegaconf import OmegaConf,MISSING
from typing import List, Optional

import os
import psutil

@dataclass
class General:
    verbose: bool = True
    try:
        ncpu: int = len(psutil.Process().cpu_affinity())
    except AttributeError:
        ncpu: int = psutil.cpu_count()
    directory: str = os.getcwd()
    multiprocessing: bool = True
    #font_file: str = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"

@dataclass
class defaults:
    print_examples: bool = False
    configuration_file: Optional[str] = None
    general: General = General()



def setup_config(argv):
    cfg = OmegaConf.structured(defaults)
    
    inputconf = OmegaConf.from_cli(argv)
    cfg_input = OmegaConf.merge(cfg,inputconf)
    if cfg_input.print_examples:
        example_name = 'sofia2hdf5.yml'
        with open(example_name,'w') as default_write:
            default_write.write(OmegaConf.to_yaml(cfg))
        print(f'''We have printed the file {example_name} in {os.getcwd()}.
Exiting sofia2hdf5''')
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
    # We only allow the maximum cpus from the command line
    if cfg.general.ncpu == psutil.cpu_count():
        print(f'''The amount of requested cpus ({cfg.general.ncpu}) is equal to the total in the hardware, this is unwise.
Therefore we lower it by one.
If you really want to use the maximum cpus set it from the command line with general.ncpu=''')
        cfg.general.ncpu -= 1
    cfg = OmegaConf.merge(cfg,inputconf) 
    #Further checks go below here

        
    return cfg