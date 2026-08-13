# File for functions related to logging messages to the screen and a log file.
from datetime import datetime
from inspect import stack
from multiprocessing.util import debug

from deblend_sofia_detections.support.system_functions import create_directory, join_path
import os

def linenumber(debug='short'):
    '''get the line number of the print statement in the main.'''
    line = []
    for key in stack():
        if key[1] == 'main.py':
            break
        if key[3] != 'linenumber' and key[3] != 'print_log' and key[3] != '<module>':
            file = key[1].split('/')
            to_add= f"In the function {key[3]} at line {key[2]}"
            if debug == 'long':
                to_add = f"{to_add} in file {file[-1]}."
            else:
                to_add = f"{to_add}."
            line.append(to_add)
    if len(line) > 0:
        if debug == 'long':
            line = ', '.join(line)+f'\n{"":8s}'
        elif debug == 'short':
            line = line[0]+f'\n{"":8s}'
        else:
            line = f'{"":8s}'
    else:
        for key in stack():
            if key[1] == 'main.py':
                line = f"{'('+str(key[2])+')':8s}"
                break
    return line

linenumber.__doc__ =f'''
 NAME:
    linenumber

 PURPOSE:
    get the line number of the print statement in the main. Not sure 
    how well this is currently working.

 CATEGORY:
    log_functions 

 INPUTS:

 OPTIONAL INPUTS:


 OUTPUTS:
    the line number of the print statement

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    If debug = True the full stack of the line print will be given, 
    in principle the first debug message in every function should set 
    this to true and later messages not.
    !!!!Not sure whether currently the linenumber is produced due to 
    the restructuring.
'''

def print_log(cfg,log_statement, case = ['main']):
    if not cfg.logging.enable:
        # If logging is disabled, we run completely silent.
        return

    case_set = set(case)
    debug_case = {'debug_start', 'debug_add', 'debug'}
    has_debug_case = bool(case_set & debug_case)

    debug_functions = [x.lower() for x in cfg.logging.debug_functions]
    no_debug = 'none' in debug_functions
    debug_all = 'all' in debug_functions

    debugging = False
    debug = 'empty'  # Keep indentation spacing when no debug info is added.

    if not no_debug:
        trig = debug_all
        if not trig:
            current_function = None
            for key in stack():
                if key[3] not in ['linenumber', 'print_log', '<module>']:
                    current_function = key[3]
                    break
            if current_function is not None:
                trig = current_function.lower() in debug_functions

        if trig:
            debugging = True
            debug = 'long' if 'debug_start' in case_set else 'short'

    log_statement = f"{linenumber(debug=debug)}{log_statement} \n"

    is_main = 'main' in case_set
    is_screen = 'screen' in case_set
    is_verbose = 'verbose' in case_set

    should_log = (
        cfg.logging.enable_log and (
            is_main
            or is_screen
            or (debugging and has_debug_case)
            or (is_verbose and (cfg.logging.verbose_log or debugging))
        )
    )

    should_screen = (
        is_screen
        or (
            cfg.logging.verbose_screen and (
                is_main
                or (debugging and has_debug_case)
                or is_verbose
            )
        )
    )

    if should_screen:
        print(log_statement)
    if should_log:
        with open(f'{cfg.logging.log_file}','a') as log_file:
            log_file.write(log_statement)
  

print_log.__doc__ =f'''
 NAME:
    print_log
 PURPOSE:
    Print statements to log if existent and screen if Requested
 CATEGORY:
    log_functions 

 INPUTS:
    log_statement = statement to be printed
    Configuration = Standard FAT Configuration

 OPTIONAL INPUTS:


    screen = False
    also print the statement to the screen

 OUTPUTS:
    line in the log or on the screen

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    linenumber, .write

 NOTE:
    If the log is None messages are printed to the screen.
    This is useful for testing functions.
'''

def start_new_log(cfg,basedir = os.getcwd(),source = None):
    if cfg.logging.enable_log and cfg.logging.enable:
        if cfg.internal.input_log_directory[0] != '/':
            cfg.logging.log_directory = join_path(basedir, cfg.internal.input_log_directory)
        cfg.logging.log_file = join_path(cfg.logging.log_directory, cfg.internal.input_log_file)
        if source is not None:
            ext = os.path.splitext(cfg.logging.log_file)[-1]
            cfg.logging.log_file = cfg.logging.log_file.replace(ext,f'_{source}{ext}')  
        # create log directory if it does not exist
        if not os.path.exists(cfg.logging.log_directory):
            create_directory(cfg.logging.log_directory)
       
        if os.path.exists(f'{cfg.logging.log_file}'):
            ext= os.path.splitext(cfg.logging.log_file)[-1]
            os.rename(cfg.logging.log_file, cfg.logging.log_file.replace(ext, f'_previous{ext}'))
        with open(f'{cfg.logging.log_file}','w') as log_file:
            log_file.write(f"Log file for deblend_sofia_detections run on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")