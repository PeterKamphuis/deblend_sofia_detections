The  deblend-sofia-detections' yaml file documentation!
======================================

Introduction 
----------------

deblend_sofia_detections is developed to work seamlessly with SoFiA-2 and really only requires the sofia run input .par file. It should then automatically check all detections and download the necessary optical images.
However, to ensure greater control on the deblending process there are optional input parameters that can be specified in the .yml file to customize the behavior of the package.
These parameters can also be specified on the command line and will override the values in the .yml file.
The .yml file will overwrite the defaults and not all parameters need to be listed.

Standalone Keywords
-------------------
*These keyword do not take any prefix*

**print_examples**
 
 *bool, optional, default = False*

  If true the code will print out examples of the .yml file and the .par file and then exit. 

**configuration_file**:
  
 *str, optional, default = deblend_sofia_detections.yml*

  The name of the configuration file to use. If there is no path it is assumed to be in the run directory. 
  logically this parameter only works from the command line and not from the .yml file which it specifies.

Input Keywords
--------------
   
*Specified with input.*

**markers**:

 *list of str, optional, default= ['Internet Image','5 HI peaks']*

  The types of markers to use for the deblending process. 
  deblend_sofia_detections allows the following options to set the markers:
   a) internet table: Set markers based on the downloaded internet table.
   
   b) manual table: Set markers based on the manually provided tables
   
   c) marker image: Take the provided manual marker image as a direct marker image
   
   d) optical image: Take a provided optical image and process it to find markers. 
   
   e) cleaned optical image: Use a pre-processed optical image where the background and stars are removed.
   
   f) internet image: Same as optical image but the image is downloaded from the specified repository (e.g., DSS) and then processed to find markers.
  
  For now only one of the options 'internet image', 'optical image', or 'cleaned optical image' can be spcecified. 
  If multiple marker types are specified, they will be added incrementally to create a single marker image to start the watershed from.
  After a watershed has been performed, with the externally specified markers, the deblend can be refined by starting a new watershed from the peaks of the HI. To ensure coherence with the external marker run, 
  if succesful in finding initial source detections, the peak ID's are merged when they fall inside the same source 
  of the exernal watershed. To use only a watershed starting from the HI peaks do not specify any external peak markers.

  The peak HI watershed run is specified by:

    g) #no HI peaks: Get markers at the highest HI peaks.
    
  In which #no is an integer specifying the number of HI peaks to use as markers. To reproduce the 
  original results from Qifeng+ 2025 one would set input.markers: ['Internet Image', '5 HI peaks'] which would make
  deblend_sofia_detections first use the optical markers from the internet image and then deblend based on the 5 HI peaks.
 

**sofia_parameters**:
 
 *str, optional, default= sofia_input.par*

  The input .par file containing the SoFiA run parameters. If there is no path it is assumed to be in the run directory. 
  If there is a path it should be the full path.

**tables**:

 *list of str, optional, default= []*

  A list of input tables to use for the deblending process. To assist in matching detections to
  optical sources and setting the optical markers.  The data in these tables is taken to be correct without exception.


**optical_images**:

 *list of str, optional, default= [None]*

  A list of optical images to use for the deblending process. If not specified, images will be used downloaded from SkyView. This can be useful if the user has better optical images than the DSS ones or if they want to use different wavelength ranges for the optical images.
  The images should be in the run_directory the ancillary_directory or the data_directory or be provided with the full path.  
  !!!At the moment only one image is supported but in the future we can expand this to multiple images and then use the one that best matche. the field of view of the data cube.!!!


**marker_images**:

 *list of str, optional, default = [None]*

  For more complex arrangements we can provide a direct marker image for the deblending.
  If not specified, no direct marker images will be used.
    

**gaia_credentials**:

 *list of str, optional, default = ['NONE','NONE']*

  The credentials to use for logging into the Gaia archive. By default it is set to ['NONE','NONE'], which means no credentials are used and we do not login.
  Unfortunately, logging in does not appear to speed up the download.
   
General Keywords
--------------

*Specified with general.*

**ncpu**: 

 *int, optional, default = psutil.cpu_count()*

  The number of CPU cores to use for processing.

**multiprocessing**: 

 *bool, optional, default = True*

  Whether to use multiprocessing for processing.

**optical_pixel_scale**: 

 *float, optional, default = 5.0*

  The amount of optical pixels that should cover a beam. If this is more than 4 arcsec
  the optical pixel size will be 4 arcsec. 

**deblend_dimension**:

 *str, optional, default = '3D'*

  The dimension of the file to be deblended. Options are '2D' for two-dimensional images (i.e., moment0)
  and '3D' for three-dimensional data cube deblending.

**counterpart_region**: 

 *str, optional, default = 'Ellipse'*

  The region to use for finding counterparts. The options are 'Beam', the size of the HI beam,
  '3Beam', three times the size of the HI beam, 'Box', the rectangular region around the source, 'Ellipse',
  10% of the ellipse major axis  as fitted by sofia, 'Full_Ellipse', the full ellipse as fitted by sofia, or 
  'No Counterpart' when you simply want to deblend the source without requiring any counterpart. 
  The minimum size is the beam.

**compactness**:

 *float, optional, default = 0.0*
 
  The compactness parameter for the watershed algorithm. This can be useful to adjust the regularity of the sources found by the watershed algorithm. A higher value will lead to more compact sources while a lower value will lead to more extended sources. The default value is 0.0 which is a good starting point but can be adjusted based on the specific data and requirements of the user.
 
 
**limit_concurrent_internet**:

 *int, optional, default = 0*

  Internet queries will be limited to this amount of queries at the same time. This is useful to limit internet traffice when many cores are used.
  When using dsd from other script one can import the main_with_input function and set the limit_concurrent_internet parameter through the runtime_ctx dictionary as {internet_query_gate: Lock()}.


**spectroscopic_manual_counterparts**

 *bool, optional, default = True*

  If set to False the manual input tables will be searched without demanding that the the sources match spectroscopically.

**internet_repository**:

 *List, optional, default = ['GAIA','NED','SIMBAD']*

  The catalogue to query for the internet counterpart search. All listed catalogues are queried if available.
  The options for now are NED, SIMBAD, GAIA. Note that GAIA is used to mask stars in the optical image if 
  it is not specified the stars will not be masked.

**sky_survey**:

 *List, optional, default = ['DSS2 Red']*

  The sky survey to download an optical image from. The options are listed in the astroquery's skyview module. 
  By default, it is set to ['DSS2 Red'].

**clear_internet_cache**:

 *bool, optional, default = False*

  If True we will clear any tables we find in ancillary_data that are from previous runs of the code. This can be useful if the user wants to start fresh and not use any previous tables that may have been generated by previous runs of the code. This will only clear tables that are generated by this code and not any other tables that may be in the ancillary_data/tables directory.
  It will also clear the astroquery caches.


**originals**:

 *bool, optional, default = False*
  
  For purpose of speed deblend_sofia_detections pickles input tables and saves processed optical images. However, this means that if changes are made to input tables they are not propagated to the code.
  If originals=True deblend will re-read the original tables and images provided by the user or redownload them instead of using the pickled and processedones.

**internet**:

 *bool, optional, default = True*

  Whether the code has internet access. This is used to determine if internet markers can be queried.
  If not available internet markers cannot be queried and the code will not be able to download any required data from the internet.
  If this is requested for markers the code will stop else it will execute without downloading any internet data.

Directories
--------------

*Specified with directories.*

**ancillary_directory**:

 *str, optional, default = 'ancillary_data'*

  The directory where the ancillary data will be stored. This includes the optical images downloaded from SkyView and any other ancillary data that is used for the deblending process. By default it is set to 'ancillary_data' in the current working directory but can be changed to any directory.

**data_directory**:

 *str, optional, default = MISSING*

  The directory where the data cubes are stored. This is used to find the data cubes for the deblending process. If not specified it is assumed that the data cubes are in the run_directory.

**run_directory**:
 
 *str, optional, default = os.getcwd()*

  The directory where the deblending process is run. This is used to find the input .par file and any other input files for the deblending process. By default it is set to the current working directory but can be changed to any directory.

**watershed_directory**:  

 *str, optional, default = 'Watershed_Output'*

  The directory where the output of the watershed algorithm will be stored. This includes the deblended masks and any other output files from the watershed algorithm. By default it is set to 'Watershed_Output' in the current working directory but can be changed to any directory.


Sofia Keywords
--------------

*Specified with sofia.*

!!!!!! Normally these would be read from the sofia input parameter file but we can also specify them in the .yml file if we want to override the values in the .par file. This can be useful if we want to run the deblending process on a different data cube than the one used for the SoFiA run or if we want to use a different parameter file for the deblending process.
If we do not have a .par file we can use deblend_cube instead of deblend and specify sofia.original_data_cube=<HI data cube.fits> sofia.original_mask=<blended mask.fits>

**original_data_cube**:

 *str, optional, default = MISSING*

  The original data cube used for the SoFiA run. This is used to find the data cube for the deblending process. If not specified it is assumed that the data cube is in the run_directory or the data_directory.

**original_mask**:

 *str, optional, default = MISSING*

  The original mask used for the SoFiA run. This is used to find the mask for the deblending process. If not specified it is assumed that the mask is in the run_directory or the data_directory.

**directory**:

 *str, optional, default = MISSING*

  The directory where the SoFiA run is located. This is used to find the input .par file and any other input files for the deblending process. If not specified it is assumed that the SoFiA run is in the run_directory.

**parameter_file**:

 *str, optional, default = MISSING*

  The input .par file containing the SoFiA run parameters. This is used to find the input .par file for the deblending process. If not specified it is assumed that the .par file is in the run_directory.

**parameter_path**:

 *str, optional, default = MISSING*

  The full path to the input .par file containing the SoFiA run parameters. This can be used to specify the location of the .par file if it is not in the run_directory.

**basename**: 

 *str, optional, default = MISSING*

  The basename used for the SoFiA run.  

**catalogue**:

 *str, optional, default = MISSING*

  The catalogue file used for the SoFiA run.

Logging Keywords
--------------
These keyword control how the logging is done and where the log files are stored. The log files are stored in the watershed_directory for general logs and in the specified watershed directory for individual sources.
by default but can be changed to any directory.

*Specified with logging.*

**enable**
 
 *bool, optional, default = True*

  If set to True the logging will be enabled. If set to False all logging, including screen messages, will be disabled and no log files will be created.

**enable_log**
 
 *bool, optional, default = False*

  If set to True the log files will be created. If set to False no log files will be created.

**verbose_log**
 
 *bool, optional, default = True*

  If set to True the  logging will be verbose.

**verbose_screen**
 
 *bool, optional, default = False*

  If set to True all messages will also be printed to the screen. If set to False only the log messages will be printed to the log files and only crucial message will be printed to screen 

**log_directory**
 
 *str, optional, default = f'Logs/<date>'*

  The directory where the log files will be stored.

**log_file**
 
 *str, optional, default = 'Log.txt'*

  The name of the log file.

**debug_functions**
 
 *List[str], optional, default = ['NONE']*

  A list of functions for which messages should be printed in debug mode. If 'ALL' is in the list, all functions will. NONE cancels debug messages

**save_counterpart_table**
 
 *bool, optional, default = True*

  If True we will save the final table of sources with their internet and manual counterparts. This can be useful if the user wants to have a record of the sources that were found and their counterparts. The table will be saved in the ancillary_data/tables directory as matched_table.pkl. If False we will not save the final table of sources with their internet counterparts. This can be useful if the user does not want to save the final table of sources with their internet counterparts. The table will not be saved in the ancillary_data/tables directory as matched_table.pkl.
