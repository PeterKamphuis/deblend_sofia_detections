from deblend_sofia_detections.catalogue.download import creating_full_FOV_optical
from deblend_sofia_detections.deblending.image_manipulation import\
    mask_gaia_stars,get_background,split_sources,freq_smooth,subtract_background,\
    mask_source_from_table,add_to_original
from deblend_sofia_detections.deblending.peak_handling import find_peaks
from deblend_sofia_detections.deblending.debug_visualization import \
    catalogue_positions_from_table,write_source_debug_overlay_safely
from deblend_sofia_detections.deblending.sofia_functions import read_sofia_table,\
    obtain_sofia_id, rerun_sofia
from deblend_sofia_detections.support.system_functions import create_directory,join_path
from deblend_sofia_detections.support.support_functions import match_size,\
    close_variables
from deblend_sofia_detections.support.table_functions import read_manual_table
from deblend_sofia_detections.support.errors import InputError
from deblend_sofia_detections.support.failure_reporting import \
    SourceFailure,write_failure_report
from deblend_sofia_detections.support.source_selection import select_source_ids
from astropy.convolution import convolve_fft
from astropy.io import fits
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.wcs import WCS

from skimage.segmentation import watershed
from photutils.segmentation import detect_threshold, detect_sources,  make_2dgaussian_kernel
PROFILING = False  # set to True to enable memory profiling
if PROFILING:
    from memory_profiler import profile
else:
    def profile(stream=None):
        def decorator(func):
            return func
        return decorator

import astropy.units as u
import copy
import numpy as np
import os
import shutil
import traceback
from datetime import datetime
# -*- coding: future_fstrings -*-



def check_source_size(cfg,segments,header):
    # Check the size of the sources in the 2D map
    if 'BMAJ' in header:
        pixel_scale = np.mean(abs(proj_plane_pixel_scales(WCS(header).celestial)))*u.deg  
        beamarea=(np.pi*abs(header['BMAJ']*header['BMIN']))/(4.*np.log(2.))*u.deg
        pix_beam_area = (beamarea/(pixel_scale**2)).value      
    else:
        pix_beam_area = (4.**2)/(4.*np.log(2.))
    for source in np.unique(segments):
        if source == 0:
            continue
        mask = segments == source

        if len(segments.shape) == 3:
            mask = np.sum(mask, axis=0) 
        mask[mask > 0] = 1
        size = np.sum(mask)
        if cfg.general.debug:
            print(f'this is the pixels in the source {source} {size} and the beam area {pix_beam_area}')
        if size < pix_beam_area:
            segments[mask] = 0
              
  
    return segments  

def check_source_surrounded(cfg,mask):
    """Check if the source in the mask is surrounded by another source."""
    results = {}
    sources = np.unique(mask)
    if len(sources)-1 <= 1:
        return mask
    for source in sources:
        if source == 0:
            continue
        # Get the mask for the source
        source_mask = mask == source
        results[f'{source}'] = {'id': source, 'surrounded': False, 
            'mask': source_mask, 'size': np.sum(source_mask),
            'others': {}}
       
        # Check if the source is surrounded by another source
        # loop through the pixels ignoring the edge

    for i,j in zip(*np.nonzero(mask)):
        id = mask[i,j]
        # Check the 8 neighbours
        neighbours = mask[i-1:i+2, j-1:j+2]
                #We are only interested in edge pixels
        if (neighbours == id).all():  # No neighbours
            results[f'{id}']['size'] -= 1
            continue
        if np.all((neighbours == 0) | (neighbours == id)):
            # We don't care about edge to background border 
            continue  # We have a border pixel
        for source in np.unique(neighbours):
            if source == id or source == 0:
                continue
            if f'{source}' not in results[f'{id}']['others']:
                results[f'{id}']['others'][f'{source}'] = 1
            else:
                results[f'{id}']['others'][f'{source}'] += 1.
    for source in results:
        
        #Get the longest border id
        if len(results[source]['others']) == 0:
            continue
        long_border= int(max(results[source]['others'], key=results[source]['others'].get))
        if long_border != 0:
            if results[source]['others'][f'{long_border}'] > results[source]['size']*0.9:
                results[source]['surrounded'] = True
                # If the source is surrounded by another source we set the mask to th id
                if cfg.general.verbose:
                    print(f'''Source {source} is for 95% surrounded by source {long_border}. 
border with {long_border} = {results[source]['others'][f'{long_border}']} pixels, 
total border = {results[source]['size']} pixels. 
So we add it to it.''')
                mask[results[source]['mask']] = long_border    
    return mask



def deblend_on_optical(cfg,data_in,optical_markers_in,outdir='./', optical_header= None,
        base_dir = './',source_id = 'unknown',mask =None):
    """If cube is None we do not deplend on 3D 
    if mom0 is None we do not deblend on 2D
    """
   
    data = copy.deepcopy(data_in)
    optical_markers = copy.deepcopy(optical_markers_in)
    if len(np.unique(optical_markers))-1 <= 1:
        if cfg.general.verbose:
            print(f"We found {len(np.unique(optical_markers))-1} optical markers in the cube {source_id} so we will not use the optical deblending results for the peak deblending.")
        return [False, 1000.]
       
    if len(data[0].data.shape) == 2:
        if cfg.general.verbose:
            print(f"Running 2D deblending based on an optical image for moment 0 \n")
        use_extend = True
        type_ind = 'moment0'
        dimension = 2
    else:
        if cfg.general.verbose:
            print(f"Running 3D deblending based on an optical image for a cube \n")
        use_extend = True
        type_ind = 'cube'
        dimension = 3
    
    if not mask is None:
        if cfg.general.verbose:
            print("Applying the mask to the data")
        data[0].data[mask == 0.] = 0.
    original_deblending = run_watershed(cfg,data[0].data, optical_markers,
        use_extend=use_extend, marker_header=optical_header, data_header=data[0].header,
        outdir=outdir,source_id=source_id)     
    if cfg.general.debug:
        fits.writeto(f"{outdir}debug_products/watershed_uncleaned_mask_{type_ind}_on_optical_source_{source_id}.fits",original_deblending,
            header=optical_header, overwrite=True)  
    if cfg.general.verbose:
        print("Matching the mask to the HI data resolution.")
    new_mask_HI = match_size(data[0].data,original_deblending,max =True)
    if cfg.general.debug:
        fits.writeto(f"{outdir}debug_products/watershed_uncleaned_mask_{type_ind}_in_HI_res.fits",
            new_mask_HI,
            header=data[0].header, overwrite=True)  
    # we need to make sure one source is not continuosly
    # surrounded by the other
   
    # We only want to do this in 2 d as in three it is allowed to be surrounded?
    if cfg.general.verbose:
        print(f"Checking if any source is surrounded by another source in the {dimension}D data.")
    start= datetime.now()
    if dimension == 2:
        new_mask = check_source_surrounded(cfg, new_mask_HI)
    else:
        new_mask = np.zeros_like(new_mask_HI)
        for i in range(new_mask_HI.shape[0]):
            new_mask[i,:,:] = check_source_surrounded(cfg, new_mask_HI[i,:,:])
    end = datetime.now()
    if cfg.general.verbose:
        print(f"Finished checking if any source is surrounded by another source in the {dimension}D data.")
        print(f"Time taken: {end - start}")
    if cfg.general.debug:
        fits.writeto(f"{outdir}debug_products/cleaned_mask_{type_ind}_based_on_optical_markers_source_{source_id}.fits",new_mask,
                header=data[0].header, overwrite=True)  
  
    if cfg.general.verbose:
        print(f"Checking the size of the sources in the {type_ind} and removing sources that are smaller than the beam.")
    new_mask_HI = check_source_size(cfg,new_mask,data[0].header,)
    new_mask_HI = np.array(new_mask_HI,dtype=int)
    hdr = copy.deepcopy(data[0].header)
    if 'BSCALE' in hdr:
        del hdr['BSCALE']
    if 'BZERO' in hdr:
        del hdr['BZERO']
    fits.writeto(f"{outdir}deblended_{type_ind}_mask_based_on_optical.fits",
        new_mask_HI,header=hdr, overwrite=True)

    if len(np.unique(new_mask_HI))-1 <= 1:
        sources = [False, 1000.]
    else:
        sources = [True, len(np.unique(new_mask_HI))-1]
    
    if cfg.general.verbose:
        print(f"Found {sources[1]} sources in the deblended {type_ind} map based on an optical image.")
    close_variables(new_mask,new_mask_HI,original_deblending)
   
    return sources

def deblend_on_peaks(cfg,cube,cube_mask=None,previous_deblend=None,outdir='./',
        source_id = 'unknown'): 
    
    if cfg.general.verbose:
        print(f"Running peak deblending for the cube {source_id}. \n ")
        print("Starting with smoothing the cube in frequency to find the peaks.")
        print(f'We first smooth the cube')
    cube_smooth = freq_smooth(cube[0].data, smooth=4.0)
    if cfg.general.debug:
        fits.writeto(f"{outdir}debug_products/smoothed_cube.fits", cube_smooth,
                header=cube[0].header, overwrite=True)
    wcs = WCS(cube[0].header)
    #Is this always m/s?
    velocity_width = wcs.wcs.cdelt[2]* u.m/u.s
    # We actually do not smooth the mask but let's make a version
    # so we could easily change this
    if not cube_mask[0] is None:
        mask_smooth=cube_mask[0].data
    else:
        mask_smooth = copy.deepcopy(cube[0].data)
        mask_smooth[cube[0].data < 1e-6] = 0.

    if not previous_deblend is None:
        if cfg.general.verbose:
            print("We have a previous deblended mask, so we match it to the current cube size.")
        previous_deblend = match_size(cube[0].data, previous_deblend,max=True)
    
    ## Step 2: peak-3D ##
    if cfg.general.verbose:
        print("Finding peaks in the cube to use as markers for the watershed algorithm.")
    # set a threshold for peak detection    
    threshold = np.zeros_like(cube_smooth) + np.nanstd(cube_smooth[0:2,:,:]) * 3.
    npeaks = cfg.input.maximum_no_peaks  # maximum number of peaks to find as it is unlikely that we have more than 
  
    if 'BMAJ' in cube[0].header:
        pixel_scale = np.mean(abs(proj_plane_pixel_scales(wcs.celestial)))*u.deg  
        pix_fwhm = (cube[0].header['BMAJ']*u.deg/pixel_scale).value
    else:
        #If we have no header we gamble that the fwhm is 3 pixels
        pix_fwhm = 3
    #Make a box size of 2 FWHM
    least_pixels = int(np.abs(((50.*u.km/u.s)/velocity_width.to(u.km/u.s)).value))
    box_size = [least_pixels,int(pix_fwhm),int(pix_fwhm)]  # size of the box to find peaks in
    if cfg.general.verbose:
        print(f"Using a box size of {box_size} to find peaks.")
    peaks,markers3d = find_peaks(cfg,cube_smooth, threshold, box_size=box_size, npeaks=npeaks,      # this is fragile at the moment ...
                mask=mask_smooth,previous_deblend =previous_deblend,outdir=outdir,
                cube_header=cube[0].header,num_processes=cfg.general.ncpu)
    if cfg.general.verbose:
        print(f"Found {len(peaks)} peaks in the cubelet {source_id}.")
   
    res3d = watershed(-cube_smooth, markers3d, mask=np.abs(mask_smooth)>1e-6,
        connectivity=1,compactness=cfg.input.compactness)
    finalhdr = copy.deepcopy(cube[0].header)# Save the results
    if 'BSCALE' in finalhdr:
        del finalhdr['BSCALE']
    if 'BZERO' in finalhdr:
        del finalhdr['BZERO']
    res3d = np.array(res3d, dtype=int)
    final_mask_name = f"{outdir}deblended_cube_mask_based_on_peaks.fits"
    fits.writeto(final_mask_name, res3d,
            header=finalhdr, overwrite=True)

    sources_3D = len(np.unique(res3d))-1
    if cfg.general.verbose:
        print(f"Found {sources_3D} sources in the cubelet {source_id} based on the peaks.")
    if sources_3D <= 1:
        result = [False, sources_3D]
    else:
        result = [True, sources_3D]
    close_variables(cube_smooth,mask_smooth,res3d,markers3d)
    return result

fn = open('profiler_logs/deblend_sofia_detections.log', 'w+') if PROFILING else None
@profile(stream=fn)
def deblend_sofia_detections(cfg):
    """
    Deblend all sources in the given data cube.

    Parameters:
    cfg (Config): The configuration object.
    
    """

    if cfg.general.verbose:
        print(f"Checking the sources in the cube {cfg.sofia.original_data_cube} in the directory {cfg.directories.data_directory}")

    if cfg.general.debug:
        print(f'We use this configuration for the deblending:')
        for key, value in cfg.__dict__.items():
            print(f"{key}: {value}")
    #load the original sofia table
    sources,table_name = read_sofia_table(cfg,
        no_conversion = True)
    source_ids_to_process = select_source_ids(
        sources['id'], cfg.input.source_ids)

    if cfg.general.verbose:
        print(
            f"Processing {len(source_ids_to_process)} of "
            f"{len(sources)} SoFiA detections."
        )
    
    if not os.path.exists(f'{cfg.internal.optical_background}'):
        if cfg.general.verbose:
            print(f"Creating the full FOV optical image for {cfg.sofia.original_data_cube}.")
        creating_full_FOV_optical(cfg)
   
    cubelets_dir = f'{cfg.sofia.directory}/{cfg.sofia.basename}_cubelets/'
    max_source_id_original = np.max([int(x) for x in sources['id']])
    max_source_id = copy.deepcopy(max_source_id_original)
    if cfg.input.use_cube_deblending:
        usemoment=False
        usecube = True
    else:
        usemoment = True
        usecube = False 
    failures = []
    succeeded_count = 0
    failure_report_path = join_path(
        cfg.directories.watershed_directory, 'deblend_failures.log')

    for source_id in source_ids_to_process:
        cube_name = (
            f"{cubelets_dir}{cfg.sofia.basename}_{source_id}_cube.fits")
        try:
            max_source_id = watershed_deblending(cfg,
                            cube_name=cube_name,
                            mask_name=f"{cubelets_dir}{cfg.sofia.basename}_{source_id}_mask.fits",
                            mom0_name=f"{cubelets_dir}{cfg.sofia.basename}_{source_id}_mom0.fits",
                            peak_deblending=cfg.input.use_peak_deblending,
                            optical_deblending= cfg.input.use_optical_deblending,
                            moment0_deblending=usemoment, cube_deblending=usecube,
                            max_source_id=max_source_id)
            succeeded_count += 1
        except Exception as error:
            failure = SourceFailure(
                source_id=str(source_id),
                cube_name=cube_name,
                exception_type=type(error).__name__,
                reason=str(error),
                traceback_text=traceback.format_exc(),
            )
            failures.append(failure)
            print(
                f"Source {source_id} failed with "
                f"{failure.exception_type}: {failure.reason}"
            )

            if cfg.general.continue_on_source_error:
                print("Recording the failure and continuing with the next source.")
            else:
                write_failure_report(
                    failure_report_path,
                    failures,
                    requested_count=len(source_ids_to_process),
                    succeeded_count=succeeded_count,
                )
                print(f"Failure report written to {failure_report_path}")
                raise

    write_failure_report(
        failure_report_path,
        failures,
        requested_count=len(source_ids_to_process),
        succeeded_count=succeeded_count,
    )
    print(
        f"Finished source processing: {succeeded_count} succeeded and "
        f"{len(failures)} failed."
    )
    print(f"Failure report written to {failure_report_path}")
  
    if max_source_id > max_source_id_original:
        rerun_sofia(cfg)
    close_variables(max_source_id,max_source_id_original,sources)

fn = open('profiler_logs/detect_optical_sources.log', 'w+') if PROFILING else None
@profile(stream=fn)
def detect_optical_sources(cfg,mask=None,source_id = 'unknown'):
    """Detect sources in the optical image to use as markers for the watershed algorithm."""
    optical_image = fits.open(cfg.internal.cleaned_optical_background)
    threshold_smooth = detect_threshold(optical_image[0].data, nsigma=3,background= 0.0)
    if cfg.general.verbose:
        print(f"Using a threshold of  {np.mean(threshold_smooth)} for source detection.")  
    npixels = int((np.pi*abs(cfg.internal.optical_kernel_fwhm**2))/(4.*np.log(2.))*3.)
    if cfg.general.verbose:
        print(f"Using npixels = {npixels} for source detection.")
    # our mask has to be the same size as the optical image and True means the pixel should be ignored
    if mask is not None:
        if cfg.general.verbose:
            print("Matching the size of the mask to the optical image.")
        hi_mask = match_size(optical_image[0].data,mask)
        hi_mask[hi_mask < 1e-8] = 0.
        np.ma.make_mask(hi_mask, copy=False)   
        outdir =  join_path(cfg.directories.watershed_directory, f'watershed_source_{source_id}/')
        if cfg.general.debug:
            fits.writeto(f'{outdir}/debug_products/hi_mask_matched_to_optical.fits',
                hi_mask, header=optical_image[0].header, overwrite=True)
        #Detect sources takes a mask where True means the pixel should be ignored 
        #Which is terribly counterintuitive so we reverse the mask
        inv_mask = np.logical_not(hi_mask)
    else:
        inv_mask = None
    
    # try to reduce the value of npixel if only one source is detected
    # As we only know for certain that the target has an optical counterpart we need a single source but
    # not necessarily more but if we have only one there is no point deblending
    segm_deblend = np.zeros(optical_image[0].data.shape)
    while np.max(segm_deblend) < 2 and npixels > 20.:
        segm_deblend = detect_sources(optical_image[0].data, threshold_smooth, npixels=npixels,mask=inv_mask)
        if segm_deblend is None:
            if cfg.general.verbose:
                print("No sources detected in the optical image reducing the size of the pixels.")
            segm_deblend = np.zeros(optical_image[0].data.shape)
        npixels -= 10
      
   
    if np.max(segm_deblend) > 0 and cfg.general.verbose:
        if cfg.general.verbose: 
            print(f"Detected {np.max(segm_deblend)} sources after deblending in the optical image.")
    else:
        if cfg.general.verbose:
            print("No sources detected in the optical image after deblending.")
    # deblending results with background masked
    if np.max(segm_deblend) > 0 and cfg.general.debug:
        fits.writeto(f'{cfg.directories.ancillary_directory}/debug_products/segmentation_map_of_optically_detected_sources_sfs{source_id}_debug_image_{cfg.internal.image_counter}.fits',
            segm_deblend.data, header=optical_image[0].header, overwrite=True)
        cfg.internal.image_counter += 1

    masked_deb = np.ma.masked_array(segm_deblend, np.abs(segm_deblend) < 1e-8)

    return masked_deb,optical_image[0].header,hi_mask

def prepare_background_optical_image(cfg,data,source_id = 'unknown'):
    """ Prepare the background optical image for the deblending on optical sources."""
    #First we set the cleaned name for this source
    path,name = os.path.split(cfg.internal.cleaned_optical_background)
    basename = os.path.splitext(name)[0]
    cfg.internal.cleaned_optical_background = join_path(path,f'{basename}_Source_{source_id}.fits')
    # First check whether our optimal image is already available
    if not os.path.exists(f'{cfg.internal.optical_background}'):
        raise FileNotFoundError(f"The optical background image {cfg.internal.optical_background} does not exist. It should have been created by the function creating_full_FOV_optical.")
    optical_image_header = fits.getheader(f'{cfg.internal.optical_background}')
    wcs_opt = WCS(optical_image_header)
    pixel_scale = np.mean(abs(proj_plane_pixel_scales(wcs_opt)))*u.deg
    fwhm = 5./ pixel_scale.to(u.arcsec).value*2.31
    if fwhm > cfg.internal.optical_kernel_fwhm:
        cfg.internal.optical_kernel_fwhm = float(fwhm)
    if os.path.isfile(f'{cfg.internal.cleaned_optical_background}'):
        if cfg.general.verbose:
            print(f"Found the cleaned optical image at {f'{cfg.internal.cleaned_optical_background}'}, so we will use it for the deblending on optical sources.")
        #lets check if the WCS matches the optical image
        tmp = fits.getheader(f'{cfg.internal.cleaned_optical_background}')
      
        elements = ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2',
                    'CDELT1', 'CDELT2'] #Naxis and CRPIX doesnt have to match as we cut to the HI,'NAXIS1]
        correct_cleaned_image = True
        for elem in elements:
            if elem in ['CTYPE1', 'CTYPE2']:
                if optical_image_header[elem] != tmp[elem]:
                    if cfg.general.verbose:
                        print(f"The WCS of the cleaned optical image does not match the WCS of the original optical image. We will use the original optical image for the deblending on optical sources.")
                        print(f"Element {elem} does not match: {optical_image_header[elem]} != {tmp[elem]}")
                    correct_cleaned_image = False
                    break
            else:
                if np.isclose(optical_image_header[elem],tmp[elem]):
                    if cfg.general.verbose:
                        print(f"The WCS of the cleaned optical image does not match the WCS of the original optical image. We will use the original optical image for the deblending on optical sources.")
                        print(f"Element {elem} does not match: {optical_image_header[elem]} != {tmp[elem]}")
                    correct_cleaned_image = False
                    break
        if not 'GAIA_MASKED' in tmp:
            tmp['GAIA_MASKED'] = False

        if tmp['GAIA_MASKED'] is False:
            correct_cleaned_image = False
            if cfg.general.verbose:
                print(f"The cleaned optical image is not masked for Gaia stars, so we will not use it for the deblending on optical sources.")
           
        # if its the same we do not have to continue to clean it        
        if correct_cleaned_image:
            #If we do not smooth we do not calculated the fwhm

            cfg.internal.cleaned_optical_background = f'{cfg.internal.cleaned_optical_background}'
            return
    else:
        if cfg.general.verbose:
            print(f"Did not find the cleaned optical image at {f'{cfg.internal.cleaned_optical_background}'}, so we will create it for the deblending on optical sources.")
    # If we get here we have to clean our downloaded image    
    if cfg.general.verbose:
        print("Preparing the background optical image for the deblending on optical sources.")
    
    wcs= WCS(data[0].header).celestial
    hi_header = copy.deepcopy(data[0].header)
    bckgrnd,bckgrnd_wcs = get_background(cfg, match_header= hi_header, wcs=wcs)
    # masking the stars in the optical image according to GAIA
    no_gaia_mask = False
    if cfg.general.verbose:
        print("Creating a gaia mask for the optical image.")
    try:    
        gaia_mask = mask_gaia_stars(bckgrnd,bckgrnd_wcs,cfg=cfg) 
    except Exception as e:
        print(f"An error occurred while creating the Gaia mask: {e}")
        print("Continuing without masking Gaia stars.")
        gaia_mask = np.zeros_like(bckgrnd).astype(bool) 
        no_gaia_mask = True
    masked_bckgrnd = copy.deepcopy(bckgrnd)
    masked_bckgrnd_wcs = copy.deepcopy(bckgrnd_wcs)
    if cfg.general.verbose:
        print("Applying a gaia mask for the optical image.")
    masked_bckgrnd[gaia_mask] = np.nan  # Mask Gaia stars in the optical image
    close_variables(gaia_mask,bckgrnd,bckgrnd_wcs)    
    if cfg.general.debug:
        fits.writeto(f'{cfg.directories.ancillary_directory}debug_products/background_gaia_masked_sfs{source_id}_debug_image_{cfg.internal.image_counter}.fits',
            masked_bckgrnd,header=masked_bckgrnd_wcs.to_header(),overwrite=True)
        cfg.internal.image_counter += 1   

    # subtracting the background from the optical image
    if cfg.general.verbose:
        print("Subtracting the background from the optical image.")
    cleaned_optical_image,cleaned_optical_wcs = subtract_background(cfg,
        masked_bckgrnd, masked_bckgrnd_wcs)
    if cfg.general.debug:
        fits.writeto(f'{cfg.directories.ancillary_directory}debug_products/background_gaia_masked_mean_subtracted_sfs{source_id}_debug_image_{cfg.internal.image_counter}.fits',
            cleaned_optical_image,header=cleaned_optical_wcs.to_header(),overwrite=True)
        cfg.internal.image_counter += 1   



    # Finally we want to smooth the image
    if cfg.general.verbose:
        print("Smoothing the optical image with a Gaussian kernel.")
  
    boxin = int(3.*cfg.internal.optical_kernel_fwhm) + 1 if\
        int(3.*cfg.internal.optical_kernel_fwhm) % 2 == 0 else \
        int(3.*cfg.internal.optical_kernel_fwhm)
    ## feel free to adjust the following parameters for better source detection. ##
    #threshold = detect_threshold(optical_image, nsigma=5,background= 0.0)
    if cfg.general.verbose:
        print(f"Using a FWHM of {cfg.internal.optical_kernel_fwhm} pixels for the Gaussian kernel and a box size of {boxin} pixels for smoothing the optical image.")
    kernel = make_2dgaussian_kernel(cfg.internal.optical_kernel_fwhm, size=boxin)
    if cfg.general.verbose:
        print(f"Smoothing the image")
    data_smooth = convolve_fft(cleaned_optical_image, kernel,) # smoothed optical images
    if cfg.general.verbose:
        print(f'writing the final background image for source {source_id} to {cfg.internal.cleaned_optical_background}')
    final_optical_header = copy.deepcopy(cleaned_optical_wcs.to_header())
    if no_gaia_mask:
        final_optical_header['GAIA_MASKED'] = False
    else:
        final_optical_header['GAIA_MASKED'] = True
    fits.writeto(cfg.internal.cleaned_optical_background, data_smooth, 
        header=cleaned_optical_wcs.to_header(), overwrite=True)
    

def load_data(base_name=None,cube_name = None,
              mask_name=None, 
              mom0_name=None):
    if not base_name is None:
        indir,name = os.path.split(base_name)
        if indir == '':
            indir = '.'
    else:
        indir = '.'
        name = 'unknown'
    # load data
    if cube_name is None:
        try:
            cube = fits.open(f"{indir}/{name}_cube.fits")
        except FileNotFoundError:
            cube = None
    else:
        cube = fits.open(f"{cube_name}")
    if mask_name is None:
        try :
            mask = fits.open(f"{indir}{name}_mask.fits")
        except FileNotFoundError:
            mask = None
    else:
        mask = fits.open(f"{mask_name}")

    if not mom0_name is None:
        mom0 = fits.open(f"{mom0_name}")
    else:
        mom0 = None
    return cube, mask, mom0

def obtain_final_mask(cfg,cube_name, results,outdir = './'):
    """ Create a final mask for the deblended sources."""
    # This is not based on the results but on the input
    if cfg.input.use_peak_deblending and results['hi_peaks'][0]:
        final_mask_name = f"{outdir}deblended_cube_mask_based_on_peaks.fits"
    elif cfg.input.use_cube_deblending and cfg.input.use_optical_deblending\
        and results['optical_cube'][0]:
        final_mask_name = f"{outdir}deblended_cube_mask_based_on_optical.fits"
    elif cfg.input.use_optical_deblending and results['optical_moment0'][0]: 
        path,name = os.path.split(cube_name)
        basename = os.path.splitext(name)[0].removesuffix('_cube')
        if cfg.general.verbose:
            print(f'''Using the 2D optical mask for the deblending.
We will check if it is ok and if so we will apply it to the original mask.
looking in {path} for the original mask {basename}_mask.fits''')
       
        original_mask = fits.open(f'{path}/{basename}_mask.fits',
                                  do_not_scale_image_data=True)
        twod_mask = fits.open(f"{outdir}deblended_moment0_mask_based_on_optical.fits"
            ,do_not_scale_image_data=True)
       
        final_mask_name = f"{outdir}deblended_moment0_mask_based_on_optical.fits"
        if twod_mask[0].header['NAXIS'] != 3:
            
           
            if len(np.unique(twod_mask[0].data))-1 <= 1:
                print("The optical 2D mask is not ok, we will not use it for the peak deblending.")
                final_mask_name = None
            else:

                tmp = np.zeros(original_mask[0].data.shape).astype(np.int8)

                for i in range(len(original_mask[0].data[:,0,0])):
                    if np.sum(original_mask[0].data[i,:,:]) > 0:
                        tmp[i,:,:] = twod_mask[0].data * original_mask[0].data[i,:,:]
                original_mask[0].data = tmp
            #Because astropy is sooo arragant and pedantic it could've been written by dutch people
                original_mask[0].scale('int32')

            #print(f"Writing the 2D mask to {original_mask[0].header['BZERO']} {original_mask[0].header['BSCALE']}")
        
                fits.writeto(f"{outdir}moment0_mask_based_on_optical.fits",
                    original_mask[0].data,original_mask[0].header,overwrite=True)
               
       
    else:
        final_mask_name = None
    return final_mask_name

def run_watershed(cfg, data, markers, use_extend=False, marker_header = None, 
        data_header = None, outdir='./',source_id = 'unknown'):
   
    if use_extend:
        if cfg.general.verbose:
            print(f'''Matching the data to the size of  the markers.
This means scaling the data by {data.shape[-1]/markers.shape[-1]} to match the pixel size of the optical image.
''')
        # If we are using the extended optical image, we need to match the size of the
        # if the data is 3D we need etend the markers first 
        if len(data.shape) == 3 and len(markers.shape) == 2:
             new_markers = np.zeros((data.shape[0],markers.shape[0],markers.shape[1])).astype(np.int8)
             for i in range(data.shape[0]):
                new_markers[i] = markers
             markers = new_markers
        data_ext = match_size(markers,data)
        header = marker_header
    else:
        if cfg.general.verbose:
            print(f'''Matching the markers to the size of  the data.
This means scaling the markers by {markers.shape[-1]/data.shape[-1]} to match the pixel size of the data.
''')
        # here match size etends the markers to 3d if required  
        markers_data = np.asarray(markers) 
        if cfg.general.debug:
            fits.writeto(f"{outdir}/debug_products/markers_input_for_dimension_3.fits",
                markers_data, header=marker_header, overwrite=True)     
        markers = match_size(data, markers)
        markers = np.asarray(markers).astype(np.int8)
        data_ext = copy.deepcopy(data)
        header = data_header
    connect = 2
    dimension = int(len(data_ext.shape))
    markers_map = copy.deepcopy(markers)
    if dimension == 3:
        connect += 2
        markers_map = markers_map[0,:,:]
    if cfg.general.debug:
        markers_data = np.asarray(markers) 
        fits.writeto(f"{outdir}/debug_products/markers_used_for_dimension_{dimension}.fits",
            markers_data, header=header, overwrite=True)
    watershed_output_1 = watershed(-data_ext, markers, mask=np.abs(data_ext)>1e-7, 
        connectivity=connect,compactness=cfg.input.compactness)
    if cfg.general.verbose:
        print(f'completed initial watershed and found {np.unique(watershed_output_1)} segments ') 
    if cfg.general.debug:
        fits.writeto(f"{outdir}/debug_products/initial_watershed_mask_dimension_{dimension}.fits",
            watershed_output_1, header=header, overwrite=True)
    # clean markers that grow less than beam  pixels
    pixel_scale = np.mean(abs(proj_plane_pixel_scales(WCS(header))))*u.deg
    if cfg.general.debug:
        print(f"Pixel scale is {pixel_scale}. data header BEAM is {data_header['BMAJ']} {data_header['BMIN']}")
    if 'BMAJ' in data_header:
        beamarea=(np.pi*abs(data_header['BMAJ']*data_header['BMIN']))/(4.*np.log(2.))*u.deg
        min_growth = (beamarea/(pixel_scale**2)).value
    else:
        min_growth = (4.**2)/(4.*np.log(2.))
   
    if dimension == 3:
        # If we have a 3rd axis we demand that they grow a beam in at least 3 channels
        min_growth *=3.
        # Remove markers that do  not grow by at least a beam in the channels that have a mask
    counter = 0
    for i in np.unique(watershed_output_1):
        if i > 0 and np.sum(watershed_output_1 == i) < np.sum(markers_map == i) + min_growth:
            markers[markers==i] = 0 
            if cfg.general.verbose:
                print(f'Removing source {i} that grows less than {min_growth} pixels.')           
            counter += 1

    if cfg.general.verbose:        
        print(f"Removed {counter} sources that grow less than {min_growth} pixels.")
            
    # rerun watershed
    watershed_output_2 = watershed(-data_ext, markers, mask=np.abs(data_ext)>1e-7
            ,connectivity=2,compactness=cfg.input.compactness)
    if cfg.general.debug:
        fits.writeto(f"{outdir}/debug_products/second_watershed_mask_dimension_{dimension}.fits",
            watershed_output_2, header=header, overwrite=True)
    
    if cfg.general.verbose:
        print(f'Finished second watershed')
    # get every segments of the deblend results
    segments = np.zeros_like(markers).astype(np.int8)
    no_source =0
    counter = 1
    for i in np.unique(watershed_output_2):
        if i > 0:
            segments[watershed_output_2 == i] = counter
            counter += 1
            no_source += 1
    if cfg.general.verbose:
        print(f"Found {no_source} sources in run_watershed.\n")
    return segments

def set_optical_markers(cfg,sofia_id, mask, diagnostics=None):
    #If our mask is 3D we sum it to get the 2d mask
    if len(mask.shape) > 2:
        mask = np.nansum(mask, axis=0)
    # detect from the optical image
    detected_optical_markers,detected_optical_markers_header,used_mask = \
        detect_optical_sources(cfg,mask=mask,source_id=sofia_id)
    if diagnostics is not None:
        diagnostics['detected_optical_markers'] = np.ma.asarray(
            detected_optical_markers
        ).filled(0).copy()
        diagnostics['catalogue_positions'] = []
    # and we want to add any source we know to exist
    if not cfg.input.manual_input_tables[0] is None:
        if cfg.general.verbose:
            print("Adding the manual optical source table.")
        source_table= read_manual_table(cfg)
        if diagnostics is not None:
            diagnostics['catalogue_positions'].extend(
                catalogue_positions_from_table(
                    source_table, catalogue="Manual input"
                )
            )
        optical_markers = mask_source_from_table(cfg, detected_optical_markers,
            detected_optical_markers_header, src_table=source_table,
            mask= used_mask)   
    else:
        optical_markers = detected_optical_markers

    if cfg.general.debug:
            fits.writeto(f"{cfg.directories.ancillary_directory}/debug_products/optical_source_markers_sfs{sofia_id}_debug_image_{cfg.internal.image_counter}.fits",
                optical_markers.data,
                header=detected_optical_markers_header, overwrite=True)
    return optical_markers,detected_optical_markers_header 
  
    

def update_original_mask(cfg, original_mask_name=None,final_mask_name=None, id=None):
    if original_mask_name is not None:
        original_mask = fits.open(original_mask_name)
        final_mask = fits.open(final_mask_name)
        original_mask_data = add_to_original(original_mask,final_mask,sofia_id=id)

        # Update the original mask with the new segmentation
     
        print(f'Writing the mask {original_mask_name}')
        original_mask.writeto(original_mask_name, overwrite=True)

fn = open('profiler_logs/watershed_deblending.log', 'w+') if PROFILING else None
@profile(stream=fn)
def watershed_deblending(cfg_in, cube_name = None, 
                         mask_name=None, 
                         mom0_name=None,base_dir= '',
                         moment0_deblending= False,
                         cube_deblending = True,
                         optical_deblending=False,
                         max_source_id = 1,
                         peak_deblending = True,
                         ):
    
    cfg = copy.deepcopy(cfg_in) #Making sure to avoid feedback
    """
    Run the watershed deblending process on the given data cube.
    
    Parameters:
    name (str): Name of the data cube file.
    outdir (str): Directory where the data cube is located.
    """  
    #outdir = f'{cfg.directories.watershed_directory}'
    sofia_id,cube_file_name = obtain_sofia_id(cfg.sofia.basename, cube_name)
    outdir = join_path(cfg.directories.watershed_directory, f'watershed_source_{sofia_id}/')
    

    if not os.path.exists(outdir):
        create_directory(outdir)
      
    else:
        #if the directory already exists we need to clean it because we will 
        #write the results there and we don't want to mix them with previous results
        for file in os.listdir(outdir):
            file_path = os.path.join(outdir, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
    if cfg.general.debug:
        if not os.path.exists(f'{outdir}debug_products'):
            create_directory('debug_products',base_directory=outdir)    
    #We have to clean all the input
    if cfg.general.verbose:
        print(f'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print(f'We start the watershed_deblending for the cube {cube_name}. \n ')
    cube, mask,  mom0 = load_data(cube_name=cube_name, 
        mask_name=mask_name, 
        mom0_name=mom0_name)
    
    results = { 'optical_moment0': [False, 0],
                'optical_cube': [False, 0],
                'hi_peaks': [False, 0]}
    debug_overlay = None
    catalogue_positions = []
    if optical_deblending:
        # If we do deblending based on optical image
        # We first prepare the optical image
        prepare_background_optical_image(cfg,cube,source_id=sofia_id)
        # Then we set markers based on the optical image and our input tables
        # but only within our HI mask 
        marker_diagnostics = {}
        optical_markers,markers_header = set_optical_markers(
            cfg,
            sofia_id,
            mask[0].data,
            diagnostics=marker_diagnostics,
        )
        catalogue_positions.extend(
            marker_diagnostics.get('catalogue_positions', [])
        )
        if cfg.general.debug and mom0 is not None:
            debug_overlay = {
                'optical_image_name': cfg.internal.cleaned_optical_background,
                'moment0_data': mom0[0].data,
                'moment0_header': mom0[0].header,
                'source_mask': mask[0].data,
                'marker_data': marker_diagnostics.get(
                    'detected_optical_markers', optical_markers
                ),
                'catalogue_positions': catalogue_positions,
                'output_name': (
                    f"{outdir}debug_products/"
                    f"optical_hi_catalogue_overlay_source_{sofia_id}.png"
                ),
                'catalogue_output_name': (
                    f"{outdir}debug_products/"
                    f"catalogue_positions_source_{sofia_id}.ecsv"
                ),
                'source_id': sofia_id,
            }
            # Write this once before watershed so a later source-level failure
            # does not remove the most useful optical diagnostic.
            write_source_debug_overlay_safely(**debug_overlay)
        # If we have only one marker or less we do not want to deblend based on the 
        # optical image because it is not useful and can lead to oversegmentation
        if len(np.unique(optical_markers.data)) - 1 <= 1:
            if cfg.general.verbose:
                print(f"We found {len(np.unique(optical_markers.data)) - 1} optical sources in the cube {cube_name} so we will not use the optical deblending results for the peak deblending.")
            # 1000 indicates 1 or less optical sources found so we shouldn't deblend
            # We require optical sources because if we simply deblend on the cube we can split anything
            results['optical_moment0'] = [False, 1000]
            results['optical_cube'] = [False, 1000]
        if moment0_deblending:
            if cfg.general.verbose:
                print(f'Running 2D deblending based on an optical image for mom0 = {mom0_name} \n')
            results['optical_moment0'] = deblend_on_optical(cfg,mom0,optical_markers,
                    outdir=outdir, source_id=sofia_id, optical_header=markers_header)
        if cube_deblending:
            if cfg.general.verbose:
                print(f'Running 3D  deblending based on an optical image for cube {cube_name}. \n')
            cube_smooth = copy.deepcopy(cube)
            cube_smooth[0].data = freq_smooth(cube[0].data, smooth=4.0)
            results['optical_cube'] = deblend_on_optical(cfg,cube_smooth,optical_markers,
                    outdir=outdir, source_id=sofia_id,optical_header=markers_header,
                    mask =mask[0].data)
   
    if peak_deblending:
        if cfg.general.verbose:
            print(f'Running peak deblending for the cube {cube_name}. \n ')

        if results['optical_cube'][0]:
            final_mask_name = f"{outdir}deblended_cube_mask_based_on_optical.fits"
        elif results['optical_moment0'][0]:
            final_mask_name = f"{outdir}deblended_moment0_mask_based_on_optical.fits"
        else:
            print(f"We failed to find any optical sources in the cube {cube_name} so we will not use the optical deblending results for the peak deblending.")
            final_mask_name = None
        
        if not final_mask_name is None:
            if cfg.general.verbose:
                print(f"Using the optical deblending mask {final_mask_name} for the peak deblending")
            
            tmp = fits.open(final_mask_name) 
            previous_deblend = tmp[0].data

        else:
            previous_deblend = None

        results['hi_peaks'] = deblend_on_peaks(cfg,cube, cube_mask=mask, outdir=outdir,
            previous_deblend=previous_deblend,source_id=sofia_id)
  
    final_mask_name = obtain_final_mask(cfg,cube_name, results,outdir=outdir)
  
    if cfg.general.verbose:
        print(f"Final mask name: {final_mask_name}")
        print(f'Which is based on the following deblending results: {results}')
   
   
    if not final_mask_name is None:  
        stil_split = split_sources(
            cfg,
            cube_name,
            final_mask_name,
            outdir=outdir,
            counterpart_positions=catalogue_positions,
            debug_overlay=debug_overlay,
        )  # skip the background
        if debug_overlay is not None:
            # Counterpart matching happens in split_sources. Rewrite the same
            # plot so accepted manual/NED matches appear as yellow points.
            debug_overlay['catalogue_positions'] = catalogue_positions
            write_source_debug_overlay_safely(**debug_overlay)
        # If we are running as a larger sofia run update the catalogue and cubelets
        if stil_split:
            max_source_id =create_final_mask(cfg, input_mask_name=final_mask_name,
                              max_source_id=max_source_id,
                              sofia_source_id=sofia_id,
                              outdir=outdir)
            update_original_mask(cfg, original_mask_name=
                f'{cfg.sofia.directory}/{cfg.sofia.basename}_mask.fits',
                final_mask_name = f'{outdir}utilized_mask.fits',id = sofia_id)
           
    if cfg.general.verbose:
        print(f"Finished the watershed deblending for the cube {cube_name}. \n ")
        print(f'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print(f'We found the following results for the deblending of the cube {cube_name}: {results}')   
    #if simple_copy:
    #    if cfg.general.verbose:
    #        print(f'''No final mask was created, skipping the source splitting. 
#Copying Watershed directory to Watershed_Output_{sofia_id} in {cfg.sofia.basename}_cubelets''')
    #    newname = f'{cfg.sofia.sofia_directory}/{cfg.internal.basename}_cubelets/Watershed_Output_{sofia_id}/'
    #    shutil.rmtree(newname) if os.path.exists(newname) else None
    #    os.rename(outdir, newname)
    return max_source_id
#@profile

def create_final_mask(cfg, input_mask_name=None, max_source_id=1, sofia_source_id=1,
                      outdir='./'):
    """ Create the final mask for the deblended sources."""
    if input_mask_name is not None:
        input_mask = fits.open(input_mask_name)
        ids_in_mask = [int(x) for x in np.unique(input_mask[0].data) if int(x) != 0]
       
        print(f'Found {len(ids_in_mask)} unique IDs in the input mask.')
        translate_ids_to_all = {}
        if int(sofia_source_id) in ids_in_mask:
            translate_ids_to_all[f'{int(sofia_source_id)}'] = int(sofia_source_id)
            skip_id = int(sofia_source_id)
        else:
            translate_ids_to_all[f'{np.min(ids_in_mask)}'] = int(sofia_source_id)
            skip_id = int(np.min(ids_in_mask))
        print(f'Skipping ID: {skip_id}')
       
        for id in ids_in_mask:
            if id != skip_id:
                max_source_id += 1
                translate_ids_to_all[f'{id}'] = max_source_id
        print(f' translation: {translate_ids_to_all}')

        for id in translate_ids_to_all:
            input_mask[0].data[input_mask[0].data == int(id)] = translate_ids_to_all[id]
        input_mask[0].data = input_mask[0].data.astype(np.int32)
        input_mask.writeto(f"{outdir}utilized_mask.fits",overwrite=True)

    return max_source_id
