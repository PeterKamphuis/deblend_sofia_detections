# OmegaConf setups

from dataclasses import dataclass, field
from omegaconf import OmegaConf,MISSING
from typing import List, Optional, Set

from datetime import datetime
import os
import psutil


@dataclass
class Input:
    # How to set the markers to use for the deblending. The options are:
    # a) internet table: Set markers based on the downloaded internet table.
    # b) manual table: Set markers based on the manually provided tables
    # c) marker image: Take the provided manual marker image as a direct marker image
    # d) optical image: Take a provided optical image and process it to find markers. 
    # e) cleaned optical image: Use a pre-processed optical image where the background and stars are removed.
    # f) internet image: Same as optical image but the image is downloaded from the specified repository (e.g., DSS) and then processed to find markers.
    # These markers will be incrementally added to the marker image used for the initial deblending. For now only one of the option 'internet image', optical image, or cleaned optical image can be used.
    # If multiple marker types are specified, they will be combined to create the final marker image.
    # The following create a new marker image to start the watershed from, however if the previous markers resulted 
    # in a succesful watershed with multiple sources,the peak markers will retain the marker id's.  
    # g) #no HI peaks: Get markers with the beam width at the highest HI peaks in the cube inside the mask, the number of peak == #no of markers
    # If general.deblend_dimension is set to '2D', the peak markers will be determined in the moment0 map.
    # If general.deblend_dimension is set to '3D', the peak markers will be determined in the full data cube.
    # The peak deblending is a seperate (second if any previous markers were defined) watershed call.
    # To reproduce the original redults from Qifeng+ 2025 One would set input.markers: ['Internet Image', '5 HI peaks']
    # This means, the deblending will use optical markers from the internet image and deblend the 5 HI peaks based on the optical deblend.
 
    markers: List = field(default_factory=lambda: ['Internet Image','5 HI peaks']) # The markers to use for the deblending.

    # The sofia input parameter file. This should be the full path to the sofia input parameter
    # file or else it should be in the run_directory.
    sofia_parameters: str = 'sofia_input.par' #Give full pathe else expected in the run directory
        
    # The manually provided input tables for the deblending. These should be the full path to the tables
    # or else they should be in the run_directory.
    tables:  List = field(default_factory=lambda: [None])
        
    # If the user provides a manual optical image, we will use that as the background for the deblending instead of downloading a DSS image. This can be useful if the user has a better optical image than the DSS one or if they want to use a different wavelength range for the optical image.
    # the images should be in the run_directory the ancillary directory or the data_directory and or be provided with the full path.  
    # At the moment only one image is supported but in the future we can expand this to multiple images and then use the one that best matches the field of view of the data cube.
    optical_images: List = field(default_factory=lambda: [None]) 

    # For more complex arrangements we can provide a direct marker image for the deblending.
    marker_images: List = field(default_factory=lambda: [None])
   
    # Gaia credentials can speed up the internet for querying Gaia catalogues. 
    # The credentials should be provided as a list with the username and password. 
    # If no credentials are provided and the optical image requires masking we will 
    # try to download without logging in
    gaia_credentials: List[str] = field(default_factory=lambda: ['NONE','NONE'])


@dataclass
class Directories:
    # The directory where we stor the the ancillary data, 
    ancillary_directory: str = 'ancillary_data'


    data_directory: str = ''

    # the directory where we run DSD and store the run-specific outputs, such as intermediate files and results.
    run_directory: str = os.getcwd()

    # the directory name that we use within the Sofia_Directory to store our watershed output
    watershed_directory: str = 'Watershed_Output'

@dataclass
class General:
    # General settings for the deblending process.
    # get the number of available CPU cores for multiprocessing
    try:
        ncpu: int = len(psutil.Process().cpu_affinity())
    except AttributeError:
        ncpu: int = psutil.cpu_count()

    # do we want to multiprocess the tasks for faster execution
    multiprocessing: bool = True

    # Amount of optical pixels that should cover a beam
    optical_pixel_scale: float = 5.

    
    # The dimension to use for the deblending. Options are '2D' for two-dimensional images (i.e., moment0)
    #  and '3D' for three-dimensional data cube deblending.
    deblend_dimension: str = '3D' 

    
    # The region to use for finding counterparts. The options are 'Beam', the size of the HI beam,
    # '3Beam', three times the size of the HI beam, 'Box', the rectangular region around the source, 'Ellipse',
    # 10% of the ellipse major axis  as fitted by sofia, or 'Full_Ellipse', the full ellipse as fitted by sofia. The minimum size is the beam.
    # If set to 'No Counterpart' we split all deblended sources without searching for counterparts.
    counterpart_region: str = 'Ellipse' 
    
    # The compactness parameter for the watershed algorithm. 
    # This can be useful to adjust the compactness of the sources found by 
    # the watershed algorithm. A higher value will lead to more compact sources 
    # while a lower value will lead to more extended sources. 
    # The default value is 0.0 which is a good starting point but can be adjusted
    # based on the specific data and requirements of the user.  
    compactness: float = 0.0

    #Limit the number of concurrent internet queries. Set to 0 for no limit.
    #This can be set from outside witt runtime ctx in deblend with_input       
    limit_concurrent_internet: int = 0

    # If set to False manual counterparts will be searched for counterparts without 
    # demanding a spectroscopic match
    spectroscopic_manual_counterparts: bool = True 

    # The catalogue to query for the internet counterpart search. 
    # All listed catalgues are queried if available. The options for now are NED, SIMBAD, GAIA
    # Note that GAIA is used to mask stars in the optical image if 
    # it is not specified the stars will not be masked.
   
    internet_repository: List = field(default_factory=lambda: ['GAIA','NED','SIMBAD'])  

    #the sky survey to download an optincal image from
    # For the options are listed in astroqueries skyview module.
    sky_survey:  List[str] = field(default_factory=lambda: ['DSS2 Red'])  

    # If True we will clear any tables we find in ancillary_data that are from 
    # previous runs of the code. This can be useful if the user wants to 
    # start fresh and not use any previous tables that may have been generated 
    # by previous runs of the code. This will only clear tables that are 
    # generated by this code and not any other tables that may be in the 
    # ancillary_data/tables directory. It will also run the astro query .clear_cache() 
    # method to clear any cached results from previous queries. 
    # When running multiple deblend processes, 
    # this could clear a cache while being used leading to errors.
    clear_internet_cache: bool = False 


    # If originals is True we will use the original images and table
    # provided by the user instead looking for previous cutouts and pickled tables. 
    originals: bool = False

    # Whether the code has internet access. This is used to determine if internet markers can be queried.
    # If not available internet markers cannot be queried.
    internet: bool = True

  

@dataclass
class Logging:
    enable: bool = True 
    enable_log: bool = False
    verbose_log: bool = True
    verbose_screen: bool = False
    log_directory: str = f'Logs/{datetime.now().strftime("%d-%m-%Y")}'
    log_file: str = 'Log.txt'  # Name of the log file
    debug_functions: List[str] = field(default_factory=lambda: ['NONE']) # List of functions to debug. If 'ALL' is in the list, all functions will be debugged. The function names should be the same as the function names in the code. This can be useful to limit the debugging to specific functions that are of interest to the user.
    save_counterpart_table: bool = False # If True we will save the final table of sources with their internet and manual counterparts. This can be useful if the user wants to have a record of the sources that were found and their counterparts. The table will be saved in the ancillary_data/tables directory as matched_table.pkl. If False we will not save the final table of sources with their internet counterparts. This can be useful if the user does not want to save the final table of sources with their internet counterparts. The table will not be saved in the ancillary_data/tables directory as matched_table.pkl.

  
@dataclass
class Sofia:
    original_data_cube: str = ''
    original_mask: str = ''
    directory: str = ''
    parameter_file: str = ''
    parameter_path: str = ''
    basename: str = ''
    catalogue: str = ''

#These are only for internal code book keeping
@dataclass
class Internal:
    image_counter: int = 0
    # The maximum number of HI peaks to find in the peak deblending. 
    # This can be useful if the user wants to limit the number of peaks 
    # to find as it is unlikely that we have more than 5 sources in a cubelet
    #  and setting this high can lead to finding spurious peaks. The user provides this 
    #  as #no HI peaks in input.markers where #no is translated to maximum_no_peaks.
    maximum_no_peaks: int = 0 
    optical_background: str = MISSING
    cleaned_optical_background: str = MISSING
    optical_kernel_fwhm: float = 3.
    sofia: str = 'sofia'
    input_log_directory: str = MISSING
    input_log_file: str = MISSING
    gaia_table: str = 'none'
    ned_table: str = 'none'
    simbad_table: str = 'none'
    cube_ext: str = '.fits'
    weights: List[float] = field(default_factory=lambda: [4., 10.]) #weights in arcsec, km/s
    
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

