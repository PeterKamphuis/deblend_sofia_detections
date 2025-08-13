

from deblend_sofia_detections.support.errors import InputError

import os
import time

def create_directory(directory,base_directory):
    split_directory = [x for x in directory.split('/') if x]
    split_directory_clean = [x for x in directory.split('/') if x]
    split_base = [x for x in base_directory.split('/') if x]
    #First remove the base from the directory but only if the first directories are the same
    if split_directory[0] == split_base[0]:
        for dirs,dirs2 in zip(split_base,split_directory):
            if dirs == dirs2:
                split_directory_clean.remove(dirs2)
            else:
                if dirs != split_base[-1]:
                    raise InputError(f"You are not arranging the directory input properly ({directory},{base_directory}).")
    wait = False
    for new_dir in split_directory_clean:
        if not os.path.isdir(f"{base_directory}/{new_dir}"):
            os.mkdir(f"{base_directory}/{new_dir}")
            wait = True
        base_directory = f"{base_directory}/{new_dir}"
    if wait:
        time.sleep(0.1)
create_directory.__doc__ =f'''
 NAME:
    create_directory

 PURPOSE:
    create a directory recursively if it does not exists and strip leading directories when the same fro the base directory and directory to create

 CATEGORY:
    support_functions

 INPUTS:
    directory = string with directory to be created
    base_directory = string with directory that exists and from where to start the check from

 OPTIONAL INPUTS:


 OUTPUTS:

 OPTIONAL OUTPUTS:
    The requested directory is created but only if it does not yet exist

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def convert_ps(image):
    'Convert a ps file to a png file.'
    base = os.path.splitext(image)[0]
    os.system(f'gs -dSAFER -dEPSCrop -r300 -sDEVICE=jpeg -o {base}.png {image}')