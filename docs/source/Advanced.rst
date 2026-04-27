The  deblend-sofia-detections' yaml file documentation!
======================================

Introduction 
----------------

deblend-sofia-detections really only requires the sofia run input .par file 
It should then automatically check all detections and download the necessary optical images.
However there are a few more optional input parameters that can be specified in the .yml file to customize the behavior of the package.


@dataclass



  
@dataclass
class defaults:
    print_examples: bool = False
    configuration_file: Optional[str] = None
   
    
    directories: Directories = field(default_factory = Directories)
    internal: Internal = field(default_factory = Internal)
    sofia: Sofia = field(default_factory = Sofia)
Input Keywords
--------------
      
*Specified with input.*

**sofia_parameters**:
 
 *str, optional, default= sofia_input.par*

  The input .par file containing the SOFIA run parameters. If there is no path it is assumed to be in the run directory. 
  If there is a path it should be the full path.

**manual_input_tables**:

   *list of str, optional, default= []*

  A list of input tables to use for the deblending process. To assist in matching detections to
  optical sources and setting the optical markers.  


**original_tables**:

  *bool, optional, default = False*
  
  deblend sofia pickles input table to speed up rerun. However, this mean that if changes are made to text tables they are not propogated to the code.
  If original_tables=True deblend will re-read the original tables provided by the user instead of using the pickled ones.


**manual_optical_image**:

  *list of str, optional, default= []*

  A list of optical images to use for the deblending process. If not specified, images will be used downloaded from SkyView.
  The images should be in the run_directory the ancillary_directory or the data_directory or be provided with the full path.  
  At the moment only one image is supported but in the future we can expand this to multiple images and then use the one that best matche. the field of view of the data cube.
 
  
**original_images**:
  
  *bool, optional, default = False* 
  
  If True deblend will  use the original images provided by the user or redownload a SkyView image instead looking for processed  cutouts downloaded. 
  This can be useful if the user has better optical images than the DSS ones or if they want to use different wavelength ranges for the optical images. The images should be in the run_directory the ancillary directory or the data_directory or be provided with the full path.
    # optical and peak deblending not mutually exclusive, if both are True we will use the original images for the deblending and if they are False we will use the cutouts downloaded from SkyView for the deblending.

**use_optical_deblending**: 
  
  *bool, optional, default = True* 
  
  If True optical markers will be used for the deblending. 
   
**use_cube_deblending**:
  
  *bool, optional, default = True*
  
  If True deblend uses the cube for the deblending otherwise the moment 0. This only applies to the optical deblending. 

**use_peak_deblending**:

  *bool, optional, default = True*
  
  If True deblend use the peak deblending for the deblending. If both use_optical_deblending and use_peak_deblending are True deblend will use the optical_deblended mask to set the values of the peak markers.
  This will ensure that the peaks only belong to a single optical source.

**maximum_no_peaks**:

 *int, optional, default = 5*
 
 The maximum number of peaks to find in the peak deblending. This can be useful if the user wants to limit the number of peaks to find as it is unlikely that we have more than 5 sources in a cubelet and setting this high can lead to finding spurious peaks.

**compactness**:

 *float, optiona, default = 0.0*
 
 The compactness parameter for the watershed algorithm. This can be useful to adjust the regularity of the sources found by the watershed algorithm. A higher value will lead to more compact sources while a lower value will lead to more extended sources. The default value is 0.0 which is a good starting point but can be adjusted based on the specific data and requirements of the user.

General Keywords
--------------

*Specified with general.*

**verbose**: 

  *bool, optional, default = True*

  Whether to print verbose output during processing.


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

**counterpart_region**: 

  *str, optional, default = 'Ellipse'*

  The region to use for finding counterparts. The options are 'Beam', the size of the HI beam,
  '3Beam', three times the size of the HI beam, 'Box', the rectangular region around the source, 'Ellipse',
  10% of the ellipse major axis  as fitted by sofia, or 'Full_Ellipse', the full ellipse as fitted by sofia. The minimum size is the beam.

**debug**:

  *bool, optional, default = False*

  Whether to print debug output during processing. This will print additional information about the processing steps and can be useful for troubleshooting and understanding the processing steps.

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

!!!!!! Normally these would be read from the sofia input parameter file but we can also specify them in the .yml file if we want to override the values in the .par file. This can be useful if we want to run the deblending process on a different data cube than the one used for the SOFIA run or if we want to use a different parameter file for the deblending process.

**original_data_cube**:

  *str, optional, default = MISSING*

  The original data cube used for the SOFIA run. This is used to find the data cube for the deblending process. If not specified it is assumed that the data cube is in the run_directory or the data_directory.

**directory**:

  *str, optional, default = MISSING*

  The directory where the SOFIA run is located. This is used to find the input .par file and any other input files for the deblending process. If not specified it is assumed that the SOFIA run is in the run_directory.

**parameter_file**:

  *str, optional, default = MISSING*

  The input .par file containing the SOFIA run parameters. This is used to find the input .par file for the deblending process. If not specified it is assumed that the .par file is in the run_directory.

**parameter_path**:

  *str, optional, default = MISSING*

  The full path to the input .par file containing the SOFIA run parameters. This can be used to specify the location of the .par file if it is not in the run_directory.

**basename**: 

  *str, optional, default = MISSING*

  The basename used for the SOFIA run.  

**catalogue**:

  *str, optional, default = MISSING*

  The catalogue file used for the SOFIA run.

Internal Keywords
--------------

*Specified with internal.*

**image_counter**:

  *int, optional, default = 0*

  A counter for the number of debug images processed. 

**optical_background**:

  *str, optional, default = MISSING*

  The name of the optical background image used for the deblending process.

**cleaned_optical_background**:

  *str, optional, default = MISSING*

  The name of the cleaned optical background image used for the deblending process.

**optical_kernel_fwhm**:

  *float, optional, default = 3.0*

  The full width at half maximum (FWHM) of the kernel used for smoothing the optical image. 
    