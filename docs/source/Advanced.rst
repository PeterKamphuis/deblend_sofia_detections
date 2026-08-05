The  deblend-sofia-detections' yaml file documentation!
========================================================

Introduction 
----------------

.. note::

  We thank Peter Kamphuis for creating and openly sharing the original package on
  which this fork is based. This configuration reference covers the committed
  implementation through ``b215ca0``. Automatic DR10 counterparts, targeted
  optical-region deblending, and Gaia-mask diagnostics are included in ``v1.1.0``;
  the position-velocity QA setting is newer than that tag. See
  :doc:`Fork_Differences` for the version-pinned lineage and record the exact
  release or commit used for scientific work.

The minimum deblender configuration identifies the parameter file from the
complete SoFiA run. That file locates the original cube and SoFiA output products;
the separate ``sofia`` executable must also be available. By default the package
checks all catalogue detections and obtains optical images as needed. The YAML
settings below control source selection, counterpart inputs, watershed behavior,
debug products, caching, and optional network services.


The YAML hierarchy follows these top-level configuration groups:

.. code-block:: python

   @dataclass
   class Defaults:
       print_examples: bool = False
       configuration_file: Optional[str] = None
       directories: Directories = field(default_factory=Directories)
       internal: Internal = field(default_factory=Internal)
       sofia: Sofia = field(default_factory=Sofia)

Running ``deblend print_examples=true`` writes all user-facing settings to
``deblend_sofia_detections_default.yml``. The generated starter enables
``auto_query_catalogue``, ``filter_dr10_markers_by_moment0_peaks``, and
``deblend_optical_regions_with_multiple_moment0_peaks``, and disables
``use_peak_deblending``. These four generated-example choices do not change the
compatibility defaults documented below. The generated
``directories.data_directory: ???`` value must be replaced before use.

Input Keywords
--------------
      
*Specified with input.*

**sofia_parameters**:
 
 *str, optional, default= sofia_input.par*

  The input .par file containing the SOFIA run parameters. If there is no path it is assumed to be in the run directory. 
  If there is a path it should be the full path.

**source_ids**:

 *list of str, optional, default = []*

  An optional allowlist of SoFiA catalogue source IDs to deblend. If the list is
  empty, deblend checks every detection in the catalogue. If IDs are supplied,
  only those detections are checked. IDs that are not present in the catalogue
  cause the run to stop with an input error before source processing begins.

  Quote IDs in YAML so they are treated consistently as catalogue identifiers::

    input:
      source_ids:
        - "42"
        - "57"
        - "221"

  This setting limits the detections considered for deblending; it does not make
  a separate SoFiA run. If a selected source is split, the final SoFiA
  re-parameterisation still uses the complete field mask and catalogue.

**manual_input_tables**:

   *list of str, optional, default = [None]*

  Counterpart tables used to match detections and set optical markers. Any
  non-null manual table takes precedence over automatic DR10 querying.

**auto_query_catalogue**:

  *bool, optional, default = False*

  If True and no manual catalogue is supplied, download and cache the relevant
  Legacy Surveys DR10 Tractor rows. A manual catalogue always takes precedence.
  The query footprint comes from the full-field moment-0 celestial WCS and is
  prepared before the per-source loop. Missing input, an invalid TAP response, or
  a network/service failure stops the run unless a compatible cache is available.

**galaxy_types**:

  *list of str, optional, default = [REX, EXP, DEV, SER]*

  The case-insensitive DR10 Tractor morphology allowlist. After this filter, at
  most one row is selected per cyan optical region using the highest finite
  ``flux_g`` value in that region. The resulting association is positional only;
  morphology and brightness do not confirm a common redshift or ownership of H I.

**filter_dr10_markers_by_moment0_peaks**:

  *bool, optional, default = False*

  If True, an automatically selected DR10 row seeds the watershed only when a
  finite positive, beam-scale local maximum in the parent H I moment-0 map maps
  inside the row's exact cyan optical-segmentation label. The elliptical search
  footprint uses ``BMAJ``, ``BMIN``, ``BPA``, and the celestial pixel scale. No
  global flux or S/N threshold is applied. Rejected DR10 rows and their associated
  cyan labels are removed from seeding, while unrelated automatic detections are
  unchanged. Manual catalogues are never filtered, and the option has no effect
  without automatic DR10 querying and optical deblending. This is independent of
  ``use_peak_deblending``. Missing or invalid parent moment-0 celestial WCS or
  synthesized-beam metadata raises a source-level input error.

**deblend_optical_regions_with_multiple_moment0_peaks**:

  *bool, optional, default = False*

  If True, project the beam-scale moment-0 peak map into the raw Photutils cyan
  segmentation and apply multi-threshold optical deblending only to a cyan label
  containing at least two peaks. The pixel membership of labels containing zero
  or one peak is preserved. DR10 type filtering and highest-``flux_g`` selection
  are then repeated independently for each new optical sublabel, followed by the
  moment-0 acceptance filter. This option requires
  ``filter_dr10_markers_by_moment0_peaks: true`` and has no effect when a manual
  catalogue takes precedence or automatic DR10/optical processing is inactive.
  Multiple moment-0 peaks do not by themselves identify multiple galaxies: a
  rotating disc, clumpy emission, noise, residual foreground structure, or beam
  metadata can produce them. The generated sublabels therefore require the same
  channel-map, spectral, position-velocity, and counterpart review as every other
  candidate split.

**optical_deblend_nlevels**:

  *int, optional, default = 32*

  The number of Photutils multi-threshold levels used for each targeted cyan
  region. The value must be at least 2.

**optical_deblend_contrast**:

  *float, optional, default = 0.001*

  The minimum fraction of a targeted region's total flux that a local optical
  peak must contain to become a separate sublabel. Values range from 0 to 1;
  smaller values produce more aggressive splitting.

**optical_deblend_min_pixels**:

  *int, optional, default = 20*

  The minimum number of connected pixels above a Photutils threshold required
  for a targeted optical object to be deblended. The value must be a positive
  integer.

**manual_markers_only**:

  *bool, optional, default = False*

  If True, discard automatically detected optical regions before adding the
  manual catalogue markers. Only manual catalogue positions overlapping the
  parent H I mask can then seed the optical watershed. A manual input table is
  required. Automatic detections remain visible as cyan diagnostics in the QA
  PNG, and the image title identifies the manual-catalogue-only marker mode.
  Yellow catalogue positions outside the parent H I footprint are diagnostic
  only and do not become watershed markers.


**original_tables**:

  *bool, optional, default = False*
  
  Reread manual text tables instead of their cached pickle files. In automatic
  DR10 mode, this also forces a new download instead of reusing the cached CSV and
  query metadata. This setting is unrelated to backups and does not protect SoFiA
  products from being regenerated.


**manual_optical_image**:

  *list of str, optional, default= []*

  A list of optical images to use for the deblending process. If not specified,
  an image will be downloaded from SkyView. Images should be in the
  ``run_directory``, ``ancillary_directory``, or ``data_directory``, or be
  provided with a full path. A manual FITS image must contain a celestial WCS.
  A 2-D image is used directly. RGB, RGBA, and other multi-plane images are
  converted to grayscale by averaging every non-celestial array axis before the
  2-D cutout is made. Verbose mode reports the original shape and collapsed
  axes. Missing files, invalid WCS, and non-overlapping images produce explicit
  input errors. At the moment only one image is supported.
 
  
**original_images**:
  
  *bool, optional, default = False* 
  
  If True deblend will  use the original images provided by the user or redownload a SkyView image instead looking for processed  cutouts downloaded. 
  This can be useful if the user has better optical images than the DSS ones or if they want to use different wavelength ranges for the optical images. The images should be in the run_directory the ancillary directory or the data_directory or be provided with the full path.
  Optical and peak deblending are not mutually exclusive. This refresh setting
  controls image caching; it does not select the watershed route.

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
----------------

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

**continue_on_source_error**:

  *bool, optional, default = True*

  Whether to record an exception from an individual SoFiA source and continue
  with the next selected source. At the end of every run, a detailed report is
  written to ``deblend_failures.log`` in the watershed output directory. Set
  this to False to restore fail-fast behaviour; the partial failure report is
  still written before the exception is raised again.

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

  Whether to print debug output and write additional diagnostic products. Each
  processed optical source receives
  ``debug_products/optical_hi_catalogue_overlay_source_<ID>.png``: a grayscale
  optical cutout with the H I source footprint in purple, faint moment-0
  contours, cyan outlines and crosses at detected optical sources, and yellow
  dots at manual, selected DR10, or accepted NED catalogue positions. DR10 rows
  rejected by the optional moment-0 filter appear as red crosses. The plotted
  catalogue rows are also written to
  ``catalogue_positions_source_<ID>.ecsv``. Plotting failures produce warnings
  and do not stop source processing.
  ``debug_products/optical_hi_components_overlay_source_<ID>.png`` is a
  companion view with separately coloured contours for every raw child from
  the first trial SoFiA parameterisation. Matching coloured circles and labels
  mark the child catalogue RA/Dec centres. These components are captured before
  catalogue-based merging or rejection and are not confirmed galaxies or
  accepted final catalogue sources. This is not a dry-run setting.
  When ``filter_dr10_markers_by_moment0_peaks`` is active, debug output also
  includes the exact parent moment-0 FITS used by the filter, a binary moment-0
  peak-map FITS, and a DR10 ECSV audit table containing identity, type,
  ``flux_g``, optical label, accepted/rejected status, matched peak position and
  value, and rejection reason. The QA overlays show accepted DR10 positions as
  yellow circles and rejected diagnostic-only positions as red crosses.
  When targeted optical deblending is enabled, the exact cyan segmentation maps
  before and after Photutils processing are written as
  ``optical_segmentation_before_targeted_deblend_source_<ID>.fits`` and
  ``optical_segmentation_after_targeted_deblend_source_<ID>.fits``. Their FITS
  headers record raw/final label counts, the number and identities of targeted
  original labels, and the configured Photutils parameters.
  Debug mode also writes ``background_gaia_masked_source_<ID>.fits`` and
  ``gaia_star_mask_source_<ID>.fits``. Their headers record ``GAIA_OK``,
  ``MASKNPIX``, and ``MASKFRAC`` so a failed Gaia query can be distinguished from
  a successful query that found no maskable stars.

**debug_pv_plots**:

  *bool, optional, default = False*

  When True together with ``debug``, write RA-velocity and Dec-velocity versions
  of both the catalogue QA and raw-child component QA plots. This produces up to
  four additional PNGs per source::

    optical_hi_catalogue_pv_ra_velocity_source_<ID>.png
    optical_hi_catalogue_pv_dec_velocity_source_<ID>.png
    optical_hi_components_pv_ra_velocity_source_<ID>.png
    optical_hi_components_pv_dec_velocity_source_<ID>.png

  The RA-velocity projection sums the cubelet over Declination, and the
  Dec-velocity projection sums it over Right Ascension. The grayscale parent-cube
  projection, purple parent-mask footprint, lavender parent H I contours, cyan
  optical positions, yellow/red catalogue positions, and child colours match the
  spatial QA style. Optical and catalogue coordinates are vertical guides because
  DR10 and optical segmentation provide no H I velocity. Raw SoFiA children also
  receive coloured contours and centre points at their measured velocity. For a
  frequency spectral axis, ``RESTFRQ`` or ``RESTFREQ`` is required for the radio-
  velocity conversion. Plotting failures remain warnings and do not stop source
  processing. The setting is opt-in because it can add four rendered figures per
  source in a full-cluster run.

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
-----------------

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
