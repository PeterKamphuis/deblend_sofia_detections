.
deblend-sofia-detections' yaml file documentation!
======================================
Introduction 
----------------
deblend-sofia-detections really only requires the sofia run input .par file 
It should then automatically check all detections and download the necessary optical images.
However there are a few more optional input parameters that can be specified in the .yml file to customize the behavior of the package.

Input Keywords
--------------
*Specified with input.*

**sofia_parameters**:
 
 *str, optional, default= sofia_input.par*

  The input .par file containing the SOFIA run parameters. If there is no path it is assume the
  file is in the current working directory. If there is a path it should be the full path.

**manual_input_tables**:

   *list of str, optional, default= []*

  A list of input tables to use for the deblending process. To assist in matching detections to
  optical sources.  

**manual_optical_image**:

  *list of str, optional, default= []*

  A list of optical images to use for the deblending process. If not specified, images will be used downloaded from SkyView.

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

  *str, optional, default = 'Beam'*

  The region to use for finding counterparts. The options are 'Beam', the size of the HI beam,
  '3Beam', three times the size of the HI beam, or 'Box', the rectangular region around the source.
