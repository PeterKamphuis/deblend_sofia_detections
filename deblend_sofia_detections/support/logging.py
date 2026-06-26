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
        #if logging is disabled, we run completely silent
        return
    
    debugging = False
    # empty tels line number to just add some spacing in front of the log statement,
    debug= 'empty'
    if cfg.logging.debug:
        trig=False
        if 'ALL' in cfg.logging.debug_functions:
            trig = True
        else:
            # get the function  
            for key in stack():
                if key[3] != 'linenumber' and key[3] != 'print_log' and key[3] != '<module>': 
                    current_function= f"{key[3]}"
                    break
            if current_function.lower() in [x.lower() for x in  cfg.logging.debug_functions]:
                trig=True      
        if trig:
            debugging=True    
            if 'debug_start' in case:
                debug = 'long'
            else:
                debug= 'short'
    log_statement = f"{linenumber(debug=debug)}{log_statement} \n"

    #Now lets check wether we want to print this specifc statement 
    print_screen = False
    print_log = False
    
    if cfg.logging.enable_log: 
        # If we want a log we have to check this specific message
        if 'main' in case or 'screen' in case or\
            (debugging and ('debug_start' in case or 'debug_add' in case or 'debug' in case))\
            or ('verbose' in case and (cfg.logging.verbose_log or debugging)):
            print_log = True
    #if we have a verbose screen we also want to print this specific message to the screen        
    if cfg.logging.verbose_screen:
        # If we want to print to the screen we have to check this specific message
        if 'main' in case or\
            (debugging and ('debug_start' in case or 'debug_add' in case or 'debug' in case))\
            or ('verbose' in case ):
            print_screen = True
    #We always print screen messages to the screen unless logging.enable = false
    if 'screen' in case:
        print_screen = True

    if print_screen:
        print(log_statement)
    if print_log:
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