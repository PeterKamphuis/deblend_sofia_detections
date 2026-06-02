# OmegaConf setups

from dataclasses import dataclass, field
from omegaconf import OmegaConf,MISSING
from typing import List, Optional

from datetime import datetime
import os
import psutil


@dataclass
class Input:
    sofia_parameters: str = 'sofia_input.par' #Give full pathe else expected in the run directory
    manual_input_tables:  List = field(default_factory=lambda: [None])
    original_tables: bool = False # If True we will use the original tables provided by the user instead looking for pickled ones
    # If the user provides a manual optical image, we will use that as the background for the deblending instead of downloading a DSS image. This can be useful if the user has a better optical image than the DSS one or if they want to use a different wavelength range for the optical image.
    # the images should be in the run_directory the ancillary directory or the data_directory and or be provided with the full path.  
    # At the moment only one image is supported but in the future we can expand this to multiple images and then use the one that best matches the field of view of the data cube.
    manual_optical_image: List = field(default_factory=lambda: [None]) 
    original_images: bool = False # If True we will use the original images provided by the user instead looking for cutouts downloaded from SkyView. This can be useful if the user has better optical images than the DSS ones or if they want to use different wavelength ranges for the optical images. The images should be in the run_directory the ancillary directory or the data_directory and or be provided with the full path.
    # optical and peak deblending not mutually exclusive, if both are True we will use the original images for the deblending and if they are False we will use the cutouts downloaded from SkyView for the deblending.
    use_optical_deblending: bool = True # If True we will use the optical images for the deblending. If False we will not use the optical images for the deblending and only use the cube information. This can be useful if the user does not have good optical images or if they want to test the deblending without using the optical information.
    use_cube_deblending: bool = True # If True we will use the cube for the deblending else the moment 0. This only applies to the optical deblending. If False we will use the moment 0 for the deblending. This can be useful if the user does not have good cube information or if they want to test the deblending without using the cube information.
    
    use_peak_deblending: bool = True # If True we will use the peak deblending for the deblending. If False we will not use the peak deblending and only use the optical information for the deblending. This can be useful if the user does not want to use the peak deblending or if they want to test the deblending without using the peak information.
    maximum_no_peaks: int = 5 # The maximum number of peaks to find in the peak deblending. This can be useful if the user wants to limit the number of peaks to find as it is unlikely that we have more than 5 sources in a cubelet and setting this high can lead to finding spurious peaks.
    compactness: float = 0.0 # The compactness parameter for the watershed algorithm. This can be useful to adjust the compactness of the sources found by the watershed algorithm. A higher value will lead to more compact sources while a lower value will lead to more extended sources. The default value is 0.01 which is a good starting point but can be adjusted based on the specific data and requirements of the user.
    internet_query: str = 'ALL' # The catalogue to query for the internet counterpart search. 
    # This can be useful if the user wants to limit the internet counterpart search to a specific catalogue. 
    # The options are 'NED', 'SIMBAD' and 'ALL'. 
    # If 'ALL' is selected, we will query both NED and SIMBAD for counterparts. 
   
@dataclass
class General:
    try:
        ncpu: int = len(psutil.Process().cpu_affinity())
    except AttributeError:
        ncpu: int = psutil.cpu_count()
    multiprocessing: bool = True
    optical_pixel_scale: float = 5. # Amount of optical pixels that should cover a beam
    counterpart_region: str = 'Ellipse' 
    

@dataclass
class Logging:
    enable: bool = True
    debug: bool = False
    enable_log: bool = False
    verbose_log: bool = True
    verbose_screen: bool = False
    log_directory: str = f'Logs/{datetime.now().strftime("%d-%m-%Y")}'
    log_file: str = 'Log.txt'  # Name of the log file
    debug_functions: List[str] = field(default_factory=lambda: ['ALL']) # List of functions to debug. If 'ALL' is in the list, all functions will be debugged. The function names should be the same as the function names in the code. This can be useful to limit the debugging to specific functions that are of interest to the user.


@dataclass
class Directories:
    ancillary_directory: str = 'ancillary_data'
    data_directory: str = MISSING
    run_directory: str = os.getcwd()
    watershed_directory: str = 'Watershed_Output'
  
@dataclass
class Sofia:
    original_data_cube: str = MISSING
    directory: str = MISSING
    parameter_file: str = MISSING
    parameter_path: str = MISSING
    basename: str = MISSING
    catalogue: str = MISSING

@dataclass
class Internal:
    image_counter: int = 0
    optical_background: str = MISSING
    cleaned_optical_background: str = MISSING
    optical_kernel_fwhm: float = 3.
    sofia: str = 'sofia'
    input_log_directory: str = MISSING
    input_log_file: str = MISSING
    
@dataclass
class defaults:
    print_examples: bool = False
    configuration_file: Optional[str] = None
    input: Input = field(default_factory = Input)
    general: General = field(default_factory = General)
    directories: Directories = field(default_factory = Directories)
    internal: Internal = field(default_factory = Internal)
    sofia: Sofia = field(default_factory = Sofia)
    logging: Logging = field(default_factory = Logging)

