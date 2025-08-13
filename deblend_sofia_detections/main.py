# -*- coding: future_fstrings -*-
"""
    A small package to deblend Sofia detections based on the method described in 
    https://ui.adsabs.harvard.edu/abs/2025ApJ...980..157H/abstract

    Most of this is an adapted copy of the accompanying ipython notebook of that paper 
    and as such Qifeng Huang should be considered as a main author of this code

"""

from deblend_sofia_detections.config.functions import setup_config
from deblend_sofia_detections.deblending.deblending import deblend_sofia_detections

import sys
import traceback
import warnings
from multiprocessing import get_context,Manager

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))



def main_with_input(argv):
    cfg = setup_config(argv)
    deblend_sofia_detections(cfg)


def main():
    argv=sys.argv[1:]
    '''Set up the configuration as input by the user'''
    cfg = setup_config(argv)
    deblend_sofia_detections(cfg)
    # for some dumb reason pools have to be called from main
    # !!!!!!!!Starts your Main Here
if __name__ =="__main__":
    main()
