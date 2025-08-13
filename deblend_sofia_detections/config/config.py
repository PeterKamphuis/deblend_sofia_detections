# OmegaConf setups

from dataclasses import dataclass, field
from omegaconf import OmegaConf,MISSING
from typing import List, Optional

import os
import psutil


@dataclass
class Input:
    sofia_parameters: str = 'sofia_input.par' #Give full pathe else expected in the run directory
    manual_input_tables:  List = field(default_factory=lambda: [None])
      # If True we will update the manual input table with the information found online
    manual_optical_image: List = field(default_factory=lambda: [None])

@dataclass
class General:
    verbose: bool = True
    try:
        ncpu: int = len(psutil.Process().cpu_affinity())
    except AttributeError:
        ncpu: int = psutil.cpu_count()
    multiprocessing: bool = True
    optical_pixel_scale: float = 5. # Amount of optical pixels that should cover a beam
    counterpart_region: str = 'Beam' 


@dataclass
class Internal:
    ancillary_directory: str = MISSING
    data_cube: str = MISSING
    data_directory: str = MISSING
    run_directory: str = MISSING
    sofia_catalogue: str = MISSING
    sofia_directory: str = MISSING
    sofia_parameter_file: str = MISSING
    sofia_basename: str = MISSING
  
    #font_file: str = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"

@dataclass
class defaults:
    print_examples: bool = False
    configuration_file: Optional[str] = None
    input: Input = field(default_factory = Input)
    general: General = field(default_factory = General)
    internal: Internal = field(default_factory = Internal)


