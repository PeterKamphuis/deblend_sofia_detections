

from deblend_sofia_detections.catalogue.download import download_full_FOV_optical
from deblend_sofia_detections.deblending.image_manipulation import\
    mask_gaia_stars,get_background,split_sources,freq_smooth,subtract_background,\
    mask_source_from_table
from deblend_sofia_detections.deblending.peak_handling import find_peaks
from deblend_sofia_detections.deblending.sofia_functions import read_sofia_table,\
    obtain_sofia_id,update_sofia_catalogue
from deblend_sofia_detections.support.system_functions import create_directory
from deblend_sofia_detections.support.support_functions import match_size

from astropy.convolution import convolve
from astropy.io import fits
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.wcs import WCS

from skimage.segmentation import watershed
from photutils.segmentation import detect_threshold, detect_sources, deblend_sources, make_2dgaussian_kernel


import astropy.units as u
import copy
import numpy as np
import os
import shutil

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

        if len(segments.shape) == 2:
            size = np.sum(mask)
            print(f'this is the pixels in the source {source} {size} and the beam area {pix_beam_area}')
            if size < pix_beam_area:
                segments[mask] = 0
        elif len(segments.shape) == 3:
            size = np.sum(mask, axis=(1,2))
            for i in range(len(size)):
                if size[i] < pix_beam_area:
                    print(size[i],pix_beam_area)
                    exit()
  
    return segments  

def check_source_surrounded(cfg,mask):
    """Check if the source in the mask is surrounded by another source."""
    results = {}
    sources = np.unique(mask)
    for source in sources:
        if source == 0:
            continue
        # Get the mask for the source
        source_mask = mask == source
        results[f'{source}'] = {'id': source, 'surrounded': False, 
            'mask': source_mask, 'size': source_mask.sum(),
            'others': {}}
        # Check if the source is surrounded by another source
        # loop through the pixels ignoring the edge
   
    for i in range(1,mask.shape[0]-1):
        for j in range(1,mask.shape[1]-1):
            if mask[i,j] == 0:
                continue
            else:
                id = mask[i,j]
                
                # Check the 8 neighbours
                neighbours = mask[i-1:i+2, j-1:j+2]
                #We are only interested in edge pixels
                if (neighbours == id).all():  # No neighbours
                    results[f'{id}']['size'] -= 1
                for source in sources:
                    if source == id:
                        continue
                    if (neighbours == source).any():  # All neighbours are the same source
                        if f'{source}' not in results[f'{id}']['others']:
                            results[f'{id}']['others'][f'{source}'] = 0
                        results[f'{id}']['others'][f'{source}'] += np.sum(neighbours == source)/9.
    for source in results:
        #Get the longest border id
        long_border= int(max(results[source]['others'], key=results[source]['others'].get))
        if long_border != 0:
            if results[source]['others'][f'{long_border}'] > results[source]['size']/9.:
                results[source]['surrounded'] = True
                # If the source is surrounded by another source we set the mask to th id
                if cfg.general.verbose:
                    print(f"Source {source} is surrounded by source {long_border}. So we add it to it.")
                mask[results[source]['mask']] = long_border    
    return mask

def deblend_moment0(data,optical_image,source_map=None, use_extend = False,
            two_dimensional =False,cfg=None,outdir='./',hdr=None):
    """ Deblends the moment 0 map of HI data using watershed segmentation.
    Parameters:
    data (np.ndarray): The HI data cube or moment 0 map.
    optical_image: optical image
    source_map (np.ndarray): optical source map 
    use_extend (bool): If true the HI is sampled on the pixel size of optical
    two_dimensional (bool): If true the data is treated as a 2D image
    ratio (int): The ratio of the pixel size of the optical map to the HI data
    """
   
       
    # find the central regions of galaxies and set them as markers
    if two_dimensional:
        markers_opt = np.zeros_like(optical_image).astype(np.int8)
        for i in np.unique(source_map):
            #print(f"Processing segment {i} of the optical source map.")
            if i == 0 or np.ma.is_masked(i) or i is None: 
                continue
            thresh = 1.
            segment = (source_map == i)
            mean, std = np.nanmean(optical_image[segment]), np.nanstd(optical_image[segment])
            if np.isnan(mean) or np.isnan(std):
                print(f"Warning: mean or std is NaN for segment {i}. Skipping this segment.")
                pass 
            else:
              
                markers_opt[np.logical_and(optical_image > mean + std, segment)] = i
    else:
        markers_opt = copy.deepcopy(source_map)
    fits.writeto(f"{outdir}optical_markers.fits",markers_opt,
                header=hdr, overwrite=True) 
   
    new_source_map = run_watershed(data, markers_opt,use_extend=use_extend
        ,cfg=cfg)  

    return np.array(new_source_map,dtype=float)



def deblend_on_optical(optical_image,cube = None, img=None,optical_name=None,
                       outdir='./',base_dir = './', cfg=None,source_id = 'unknown',
                       source_table=None, mom0=None,cube_mask =None):
    """If cube is None we do not deplend on 3D 
    if mom0 is None we do not deblend on 2D
    """
    if cube is not None:
        wcs     = WCS(cube[0].header).celestial
        hi_header = copy.deepcopy(cube[0].header)
        hi_mask = match_size(optical_image[0].data,np.nansum(cube[0].data,axis=0))
    else:
        wcs = WCS(mom0[0].header).celestial
        hi_header = copy.deepcopy(mom0[0].header)
        hi_mask = match_size(optical_image[0].data,mom0[0].data)

   
    #Obtain the background image we are using
    
    wcs_opt = WCS(optical_image[0].header).celestial
    bckgrnd,bckgrnd_wcs = get_background(img=img,wcs_opt=wcs_opt,
                            optical_name=optical_name,
                            match_header= hi_header,
                            wcs=wcs)
    bckgrnd = bckgrnd.astype(np.float32)
    
   
    #____________________ 2. Optical sources
    # First we mask the stars in the optical image according to GAIA
    if cfg.general.verbose:
        print("Creating a gaia mask for the optical image.")
    mask_stars = True
    if os.path.isfile(f'{cfg.internal.ancillary_directory}/masked_background_{source_id}.fits'):
        tmp = fits.open(f'{cfg.internal.ancillary_directory}/masked_background_{source_id}.fits')
        #Ensure that this is the same WCS as the optical image
        masked_bckgrnd_wcs = WCS(tmp[0].header).celestial

        bckgrnd_header = bckgrnd_wcs.to_header()
        elements = ['CTYPE1', 'CTYPE2', 'CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2',
                    'CDELT1', 'CDELT2']
       
        mask_stars = False
       
        for el in elements:
            if tmp[0].header[el] != bckgrnd_header[el]:
                print(f"Element {el} does not match: {tmp[0].header[el]} != {bckgrnd_header[el]}")
                mask_stars = True
                if cfg.general.verbose:
                    print(f"The WCS of the masked background {masked_bckgrnd_wcs} does not match the WCS of the background {bckgrnd_wcs}.")
                break
        if not mask_stars:
            masked_bckgrnd = tmp[0].data
       

    if mask_stars:  
        gaia_mask = mask_gaia_stars(bckgrnd,bckgrnd_wcs,cfg=cfg) 
        masked_bckgrnd = copy.deepcopy(bckgrnd)
        if cfg.general.verbose:
            print("Applying a gaia mask for the optical image.")
        masked_bckgrnd[gaia_mask] = np.nan  # Mask Gaia stars in the optical image
        # subtract the background
        #masked_bckgrnd = subtract_background(masked_bckgrnd, bckgrnd_wcs)

        fits.writeto(f'{cfg.internal.ancillary_directory}/masked_background_{source_id}.fits',masked_bckgrnd
                    ,header=bckgrnd_wcs.to_header(),
                    overwrite=True)   
        del gaia_mask     
    # We can now use the masked data to create a mask for the HI cube
    #if cfg.general.verbose:
    #    print("Plotting the masked image with mom 0 overlay.")    
    #It seems that parts of galaxies are blanked, let's see if this becomes a problem 
    if cfg.general.verbose:
        print("Looking for optical sources in the background image.")
    # Create a mask based on the hi map so we only use sources within the HI detection
    if cube is not None:
        hi_mask = match_size(masked_bckgrnd,np.nansum(cube[0].data,axis=0))
    else:
        hi_mask = match_size(masked_bckgrnd,mom0[0].data)
    hi_mask[hi_mask < 1e-8] = 0.
    np.ma.make_mask(hi_mask, copy=False)    
   
    source_markers_2d = detect_optical_sources(masked_bckgrnd, bckgrnd_wcs,
        cfg=cfg,base_dir= base_dir,mask= hi_mask
        ,outdir=outdir,hdr=optical_image[0].header)
    fits.writeto(f"{outdir}detected_sources.fits",source_markers_2d.data,
                header=bckgrnd_wcs.to_header(), overwrite=True)
   
    # If we have an input table of sources lets add them
    if not source_table is None:
        if cfg.general.verbose:
            print("Adding the manual optical source table.")
        source_markers_2d  = mask_source_from_table(masked_bckgrnd, bckgrnd_wcs,
                                            src_table=source_table,
                                            mask= source_markers_2d )
    if cfg.general.verbose:
        print(f'We found {len(np.unique(source_markers_2d))-1} sources in the optical image.')
   
    if len(np.unique(source_markers_2d))-1 <= 1:
        if cfg.general.verbose:
            print("Only one or no source found in the optical. Skipping the deblending.")
        return [False,1000]
    
    # from here on 3D differs from  2D deblending
    sources_2D = [False,None]
   
    if not mom0 is None:
        masked_mom0 = np.ma.masked_array(mom0[0].data, np.abs(mom0[0].data) < 1e-8)
        segments_2d_map = deblend_moment0(masked_mom0,masked_bckgrnd,
                                        source_map= source_markers_2d , 
                                        use_extend=True, two_dimensional=True,
                                        cfg= cfg,outdir=outdir,hdr = optical_image[0].header,)
        segments_2d_map= np.array(segments_2d_map, dtype=int)
        fits.writeto(f"{outdir}optical_2D_watershed_all.fits",segments_2d_map,
                    header=optical_image[0].header, overwrite=True)  
        # we need to make sure one source is not continuosly
        # surrounded by the other
        new_mask = check_source_surrounded(cfg, segments_2d_map)
        
       
        fits.writeto(f"{outdir}optical_2D_watershed_masks_original.fits",new_mask,
                    header=optical_image[0].header, overwrite=True)  
        
        segments_2d_map_HI = match_size(mom0[0].data,segments_2d_map,max =True)
        segments_2d_map_HI = check_source_size(cfg,segments_2d_map_HI,mom0[0].header,)
       
        fits.writeto(f"{outdir}optical_2D_watershed_masks.fits",segments_2d_map_HI,
                    header=mom0[0].header, overwrite=True)
      
        if len(np.unique(segments_2d_map))-1 <= 1:
            sources_2D = [False, 1000.]
        else:
            sources_2D = [True, len(np.unique(segments_2d_map_HI))-1]
        '''
        plot_source_image(masked_bckgrnd,bckgrnd_wcs,masked_mom0,
                          outdir=outdir,name='mom0_on_stars')
        plot_source_image(segments_2d_map,bckgrnd_wcs,masked_mom0,
                          outdir=outdir,name='mom0_on_2Dsegment',cmap='Set1')
        '''
        if cfg.general.verbose:
            print(f"Found {sources_2D} sources in the 2D deblending process.")
        del segments_2d_map
        del segments_2d_map_HI
        del masked_mom0

    sources_3D = [False, None]
    

    # Normally we only do the 3D deblending if we have more than one source    
    if not cube is None:
        #Smooth the cube in frequecncy
        cube_smooth = freq_smooth(cube[0].data, smooth=4.0)
        # We actually do not smooth the mask but let's make a version
        # so we could easily change this
        if not cube_mask[0] is None:
            mask_smooth=cube_mask[0].data
        else:
            mask_smooth = copy.deepcopy(cube[0].data)
            mask_smooth[cube[0].data < 1e-6] = 0.

        # The 3D markers are initially the same as the 2D markers in every channel
        #Copy the markers to the 3D sized cube first in 2D and then put them in every channel
        markersin_HI_resolution = match_size(cube[0].data[0,:,:], source_markers_2d,max=True)       
        markers3d = np.zeros(cube_smooth.shape).astype(np.int8)
        mask_channels = 0.
        for freq in range(len(cube_smooth)):
            markers3d[freq] =  markersin_HI_resolution
            #Lets check the channels that have a mask
            if np.nansum(mask_smooth[freq]) > 0:
                mask_channels += 1
        #Then apply the watershed algorithm
        res3d0 = watershed(-cube_smooth, markers3d, mask=np.abs(mask_smooth)>1e-6,
                   connectivity=4)
        # The segment need to grow by at least a beam in the channels that have a mask
        pixel_scale = np.mean(abs(proj_plane_pixel_scales(wcs)))*u.deg
        if 'BMAJ' in cube[0].header:
           beamarea=(np.pi*abs(cube[0].header['BMAJ']*cube[0].header['BMIN']))/(4.*np.log(2.))*u.deg
           pix_beam_area = (beamarea/(pixel_scale**2)).value
           pix = pix_beam_area
           #pix = (cube[0].header['BMAJ']*u.deg/pixel_scale).value
        else:
           pix = (4.**2)/(4.*np.log(2.))
        min_growth = pix * mask_channels # Minimum growth in pixels
        # Remove markers that do  not grow by at least a beam in the channels that have a mask
        counter = 1
        for i in np.unique(res3d0):
            if i > 0 and np.sum(res3d0 == i) < np.sum(markers3d == i) + min_growth:
                markers3d[res3d0==i] = 0 
                if cfg.general.verbose:
                    print(f'Removing source {i} that grows less than {min_growth} pixels.') 
            else:
                markers3d[res3d0==i] = counter
                counter += 1
        # rerun watershed
        res3d0 = watershed(-cube_smooth, markers3d, mask=np.abs(mask_smooth)>1e-6,
                   connectivity=2)
        # Save the segments
        fits.writeto(f"{outdir}optical_3D_watershed_mask.fits",res3d0,
                    header=cube[0].header, overwrite=True)
        if len(np.unique(res3d0))-1 <= 1:
            sources_3D = [False, len(np.unique(res3d0))-1]
        else:
            sources_3D = [True, len(np.unique(res3d0))-1]

        del res3d0
        del markers3d
        del markersin_HI_resolution  
        del mask_smooth
        del cube_smooth  
    del source_markers_2d
    del masked_bckgrnd
    del bckgrnd
    del hi_mask
    del bckgrnd_wcs
    del wcs_opt
    del wcs
    if not sources_2D[1] is None and not sources_3D[1] is None:
        if sources_2D[1] == sources_3D[1] and sources_2D[0]:  
            return sources_2D
        else:
            return [False,sources_3D[1]]
    elif not sources_2D[1] is None:
        return sources_2D
    elif not sources_3D[1] is None:
        return sources_3D
    else:
        return [False,0]


def deblend_on_peaks(cfg,cube,cube_mask=None,previous_deblend=None,outdir='./',
        source_id = 'unknown'): 
    cube_smooth = freq_smooth(cube[0].data, smooth=4.0)
    fits.writeto(f"{outdir}smoothed_cube.fits", cube_smooth,
                header=cube[0].header, overwrite=True)
    wcs = WCS(cube[0].header)
    velocity_width = wcs.wcs.cdelt[2] * u.m/u.s
    # We actually do not smooth the mask but let's make a version
    # so we could easily change this
    if not cube_mask[0] is None:
        mask_smooth=cube_mask[0].data
    else:
        mask_smooth = copy.deepcopy(cube[0].data)
        mask_smooth[cube[0].data < 1e-6] = 0.
    if not previous_deblend is None:
        if len(previous_deblend.shape) == 2:
            previous_deblend = match_size(cube[0].data[0,:,:], previous_deblend,max=True)
            tmp = np.zeros(cube_smooth.shape).astype(np.int8)
            for i in range(len(cube_smooth)):
                tmp[i] = previous_deblend
            previous_deblend = tmp
        else:
            previous_deblend = match_size(cube[0].data, previous_deblend,max=True)
    
    ## Step 2: peak-3D ##
    # set a threshold for peak detection    
    threshold = np.zeros_like(cube_smooth) + np.nanstd(cube_smooth[0:2,:,:]) * 3.
    npeaks = 5  # maximum number of peaks to find as it is unlikely that we have more than 
    # 5 sources in a cubelet we do not want to set this high

  
    if 'BMAJ' in cube[0].header:
        pixel_scale = np.mean(abs(proj_plane_pixel_scales(wcs.celestial)))*u.deg  
        #print(pixel_scale,cube[0].header['BMAJ'],cube[0].header['CDELT1'])
        #exit()
        pix_fwhm = (cube[0].header['BMAJ']*u.deg/pixel_scale).value
    else:
        #If we have no header we gamble that the fehm is 3 pixels
        pix_fwhm = 3
    #Make a box size of 2 FWHM
    least_pixels = int(((50.*u.km/u.s)/velocity_width.to(u.km/u.s)).value)
    box_size = [least_pixels,int(pix_fwhm),int(pix_fwhm)]  # size of the box to find peaks in

    peaks,markers3d = find_peaks(cfg,cube_smooth, threshold, box_size=box_size, npeaks=npeaks,      # this is fragile at the moment ...
                mask=mask_smooth,previous_deblend =previous_deblend,outdir=outdir,
                cube_header=cube[0].header,num_processes=cfg.general.ncpu)
  
  
    res3d = watershed(-cube_smooth, markers3d, mask=np.abs(mask_smooth)>1e-6)
    finalhdr = copy.deepcopy(cube[0].header)# Save the results

    res3d = np.array(res3d, dtype=int)
    final_mask_name = f"{outdir}peak_watershed_masks.fits"
    fits.writeto(final_mask_name, res3d,
            header=finalhdr, overwrite=True)

    sources_3D = len(np.unique(res3d))-1
    if sources_3D <= 1:
        result = [False, sources_3D]
    else:
        result = [True, sources_3D]
    del cube_smooth
    del mask_smooth
    del res3d
    del markers3d
    
    return result

def deblend_sofia_detections(cfg):
    """
    Deblend all sources in the given data cube.

    Parameters:
    cfg (Config): The configuration object.
    
    """

    if cfg.general.verbose:
        print(f"Checking the sources in the cube {cfg.internal.data_cube} in the directory {cfg.internal.data_directory}")

    
    #load the original sofia table
    sources,sofia_basename,table_name = read_sofia_table(cfg,
        sofia_directory=cfg.internal.sofia_directory,
        sofia_basename=cfg.internal.sofia_basename,
        no_conversion = True)

    if not os.path.exists( f'{cfg.internal.data_directory}/ancillary_data/moment0_full_DSS.fits'):
        if cfg.general.verbose:
            print(f"Downloading the full FOV optical image for {cfg.internal.data_cube}.")
        download_full_FOV_optical(cfg)
    cubelets_dir = f'{cfg.internal.sofia_directory}/{sofia_basename}_cubelets/'

    for id in sources['id']:
        id = int(id)
        watershed_deblending(cfg,
                        cube_name=f"{cubelets_dir}{sofia_basename}_{id}_cube.fits",
                        mask_name=f"{cubelets_dir}{sofia_basename}_{id}_mask.fits",
                        optical_name=f"{cfg.internal.data_directory}/ancillary_data/moment0_full_DSS.fits",
                        mom0_name=f"{cubelets_dir}{sofia_basename}_{id}_mom0.fits",
                        peak_deblending=True,optical_deblending= True,
                        main_name = cfg.internal.data_cube,
                        two_dim_deblending=True,
                        cfg_in=cfg)
        
def detect_optical_sources(optical_image, optical_wcs
        ,cfg=None,base_dir='',mask=None,outdir='./',hdr=None,subtract = True):
    
    if subtract:
        optical_image = subtract_background(optical_image, optical_wcs)
        fits.writeto(f'{outdir}background_subtracted_optical_image.fits',
                optical_image, header=optical_wcs.to_header(), overwrite=True)
    # source detection and deblending using photutils. 
    pixel_scale = np.mean(abs(proj_plane_pixel_scales(optical_wcs)))*u.deg
    fwhm = 7./pixel_scale.to(u.arcsec).value
    if fwhm < 3:
        fwhm = 3.0
    boxin = int(3.*fwhm) + 1 if int(3.*fwhm) % 2 == 0 else int(3.*fwhm)
    
    ## feel free to adjust the following parameters for better source detection. ##
    #threshold = detect_threshold(optical_image, nsigma=5,background= 0.0)
    kernel = make_2dgaussian_kernel(fwhm, size=boxin)
    data_smooth = convolve(optical_image, kernel,) # smoothed optical images
    
    #data_smooth = optical_image
    threshold_smooth = detect_threshold(data_smooth, nsigma=3,background= 0.0)
    if cfg.general.verbose:
        print(f"Using a threshold of  {np.mean(threshold_smooth)} for source detection.")
    fits.writeto(f'{outdir}/smoothed_optical_image.fits', 
                data_smooth, header=optical_wcs.to_header(), overwrite=True)
    npixels = int((np.pi*abs(fwhm**2))/(4.*np.log(2.))*3.)
    if cfg.general.verbose:
        print(f"Using npixels={npixels} for source detection.")
    #Detect sources takes a mask where True means the pixel should be ignored 
    #Which is terribly counterintuitive so we reverse the mask
    inv_mask = np.logical_not(mask) if mask is not None else None
    segm = detect_sources(data_smooth, threshold_smooth, npixels=npixels,mask=inv_mask)
    if segm is None:
        print("No sources detected in the optical image. Exiting.")

        segm_deblend = np.zeros(data_smooth.shape)
    #elif np.max(segm) == 1:
    #    if cfg.general.verbose:
    #        print(f"Detected {np.unique(segm)-1} sources in the optical image. Retrying with smooth threshold")
    #    segm = detect_sources(data_smooth, threshold_smooth, npixels=npixels,mask=inv_mask)
    #    if segm is None:
    #        print("No sources detected in the optical image. Exiting.")
    #        segm_deblend = np.zeros(data_smooth.shape)
    #    else:
    #        segm_deblend = segm
    else:
        segm_deblend = segm
    if np.sum(segm_deblend) != 0:
        fits.writeto(f'{outdir}/first_segmentation.fits',
            segm_deblend.data, header=optical_wcs.to_header(), overwrite=True)

    #    if cfg.general.verbose:
    #        print(f"Detected {np.unique(segm)-1} sources before deblending in the optical image.")
    #    segm_deblend = deblend_sources(data_smooth, segm, npixels=npixels, nlevels=16, contrast=0.1)
    #    if cfg.general.verbose:
    #        print(f"Detected {np.unique(segm_deblend)-1} sources in the optical image.")
    # try to reduce the value of npixel if only one source is detected
    # As we only know for certain that the target has an optical counterpart we need a single source but
    # not necessarily more
    while np.max(segm_deblend) < 1 and npixels > 20.:
        if cfg.general.verbose:
            print(f"Detected only {np.max(segm_deblend)} sources, reducing npixels from {npixels} to {npixels-10}.")
        npixels -= 10
        segm = detect_sources(data_smooth, threshold_smooth, npixels=npixels,mask=inv_mask)
        if segm is None:
            print("No sources detected in the optical image. Exiting.")
            segm_deblend = np.zeros(data_smooth.shape)
        #else:
        #    segm_deblend = deblend_sources(data_smooth, segm, npixels=npixels, nlevels=16, contrast=0.1)

    # deblending results with background masked
    masked_deb = np.ma.masked_array(segm_deblend, np.abs(segm_deblend) < 1e-8)

    return masked_deb


    

def load_data(name=None, indir=None,cube_name = None,
              mask_name=None, optical_name=None,
              mom0_name=None):
    
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
    if optical_name is None:
        try:
            
            opti = f"{indir}{name}_optical.fits"
        except FileNotFoundError:
            print(f"Optical file {name}_optical.fits not found in {indir}.")
            opti = None
    else:
        basename = os.path.splitext(os.path.basename(optical_name))[0]
        opt_path = os.path.dirname(optical_name)
      
        if os.path.isfile(f'{opt_path}/{basename}_bg_subtracted.fits'):
            opti = fits.open(f'{opt_path}/{basename}_bg_subtracted.fits')
        else:
            opti = fits.open(optical_name)
            subtr = subtract_background(opti[0].data, WCS(opti[0].header))
            fits.writeto(f'{opt_path}/{basename}_bg_subtracted.fits', subtr, opti[0].header, overwrite=True)
            opti[0].data= subtr
    if not mom0_name is None:
        mom0 = fits.open(f"{mom0_name}")
    else:
        mom0 = None
    return cube, mask, opti, mom0

def obtain_final_mask(cfg,cube_name, results,outdir = './'):
    """ Create a final mask for the deblended sources."""
    if results['hi_peaks'][0]:
        final_mask_name = f"{outdir}peak_watershed_masks.fits"
    elif results['optical_3d'][0]:
        final_mask_name = f"{outdir}optical_3D_watershed_masks.fits"
    elif results['optical_2d'][0]: 
        path,name = os.path.split(cube_name)
        basename = os.path.splitext(name)[0]
        original_mask = fits.open(f'{path}/{basename.split("_cube")[0]}_mask.fits',
                                  do_not_scale_image_data=True)
        twod_mask = fits.open(f"{outdir}optical_2D_watershed_masks.fits"
            ,do_not_scale_image_data=True)
       
        final_mask_name = f"{outdir}optical_2D_watershed_masks.fits"
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
        
                fits.writeto(f"{outdir}optical_2D_watershed_masks.fits",
                    original_mask[0].data,original_mask[0].header,overwrite=True)
               
       
    else:
        final_mask_name = None
    return final_mask_name

def run_watershed(mom0_ext, markers2d,use_extend=False, cfg = None):
    ratio = abs(markers2d.shape[-1] / mom0_ext.shape[-1])
    if use_extend:
        if cfg.general.verbose:
            print(f"Using the extended optical image with ratio {ratio} to match the moment-0 map.")
        # If we are using the extended optical image, we need to match the size of the
        mom0_ext = match_size(markers2d,mom0_ext)
    else:
        if cfg.general.verbose:
            print(f"Using the extended optical image with ratio {ratio} to match the moment-0 map.")
    
        markers2d = match_size(mom0_ext, markers2d)
   
    #print(f'Using these two array  {mom0_ext.shape} and {markers2d.shape}')
    # initial watershed (res: result)
    
    res2d = watershed(-mom0_ext, markers2d, mask=np.abs(mom0_ext)>1e-6, connectivity=2)
    if cfg.general.verbose:
        print(f'completed initial watershed and found {np.unique(res2d)} segments ') 
    # clean markers that grow less than min_growth pixels
    if len(mom0_ext.shape) == 3:
        min_growth = 5.*mom0_ext.shape[0] * ratio**2
    else:
        min_growth = 5. * ratio**2
    no_source = 0
    for i in np.unique(res2d):
        if i > 0 and np.sum(res2d == i) < np.sum(markers2d == i) + min_growth:
            markers2d[res2d==i] = 0
            no_source += 1
    if cfg.general.verbose:        
        print(f"Removed {no_source} sources that grow less than {min_growth} pixels.")
            
    # rerun watershed
    res2d = watershed(-mom0_ext, markers2d, mask=np.abs(mom0_ext)>1e-6,connectivity=2)
    if cfg.general.verbose:
        print(f'Finished second watershed')
    # get every segments of the deblend results
    segments = np.zeros_like(markers2d).astype(np.int8)
    no_source =0
    counter = 1
    for i in np.unique(res2d):
        if i > 0:
            segments[res2d == i] = counter
            counter += 1
            no_source += 1
    if cfg.general.verbose:
        print(f"Found {no_source} sources in run_watershed.\n")
    return segments


def watershed_deblending(cfg, name = None, cube_name = None, 
                         mask_name=None, optical_name=None, 
                         mom0_name=None,base_dir= '',
                         source_table=None, two_dim_deblending= False,
                         optical_deblending=False, cfg_in = None,
                         peak_deblending = True,
                         main_name=None):
    
    cfg = copy.deepcopy(cfg_in) #Making sure to avoid feedback
    """
    Run the watershed deblending process on the given data cube.
    
    Parameters:
    name (str): Name of the data cube file.
    outdir (str): Directory where the data cube is located.
    """
    outdir = f'{cfg.internal.sofia_directory}/Watershed_Output/' 
    if not os.path.exists(outdir):
        create_directory('Watershed_Output',cfg.internal.sofia_directory)
    else:
        if os.path.isfile(f'{outdir}/sofia_input.par'):
            os.remove(f'{outdir}/sofia_input.par')
        if os.path.isfile(f'{outdir}/optical_2d_watershed_masks.fits'):
            os.remove(f'{outdir}/optical_2d_watershed_masks.fits')
        if os.path.isfile(f'{outdir}/optical_3d_watershed_masks.fits'):
            os.remove(f'{outdir}/optical_3d_watershed_masks.fits')
        if os.path.isfile(f'{outdir}/peak_watershed_masks.fits'):
            os.remove(f'{outdir}/peak_watershed_masks.fits')
        
        

   
    #We have to clean all the input

    sofia_id,cube_file_name = obtain_sofia_id(cfg.internal.sofia_basename, cube_name) 
    cube, mask, opti, mom0 = load_data(name=name, indir=outdir,
                                            cube_name=cube_name, 
                                            mask_name=mask_name, 
                                            optical_name=optical_name, 
                                            mom0_name=mom0_name)
    
    results = { 'optical_2d': [False, 0],
                'optical_3d': [False, 0],
                'hi_peaks': [False, 0]}
    if optical_deblending:
        if two_dim_deblending:
            results['optical_2d'] = deblend_on_optical(opti, optical_name=optical_name,
                       outdir=outdir,base_dir =base_dir, cfg=cfg,source_id=sofia_id,
                       source_table=source_table, mom0=mom0,cube_mask =mask)
        else:
            results['optical_3d'] = deblend_on_optical(opti, cube=cube, optical_name=optical_name,
                       outdir=outdir,base_dir =base_dir, cfg=cfg,source_id=sofia_id,
                       source_table=source_table, mom0=mom0,cube_mask =mask)
   
  
    if results['optical_2d'][1] == 1000 or results['optical_3d'][1] == 1000:
            # 1000 indicates 1 or less optical sources found so we shouldn't deblend
            # We require optical sources because if we simply deblend on the cube we can split anything
            final_mask_name = None
    else:

        if peak_deblending:
           
            if results['optical_2d'][0]:
                final_mask_name = f"{outdir}optical_2D_watershed_masks.fits"
            elif results['optical_3d'][0]:
                final_mask_name = f"{outdir}optical_3D_watershed_masks.fits"
            else:
                final_mask_name = None
            
            if not final_mask_name is None:
                if cfg.general.verbose:
                    print(f"Using the optical deblending mask {final_mask_name} for the peak deblending")
               
                tmp = fits.open(final_mask_name) 
                previous_deblend = tmp[0].data

            else:
                previous_deblend = None

            results['hi_peaks'] = deblend_on_peaks(cfg,cube, cube_mask=mask, outdir=outdir,
                                                previous_deblend=previous_deblend)

        final_mask_name = obtain_final_mask(cfg,cube_name, results,outdir=outdir)
    
    if cfg.general.verbose:
        print(f"Final mask name: {final_mask_name}")
        print(f'Which is based on the following deblending results: {results}')
    simple_copy = False
    if not final_mask_name is None:  
        os.system(f'cp {final_mask_name} {outdir}utilized_mask.fits')
       
        stil_split = split_sources(cfg,cube_name, final_mask_name, outdir=outdir)  # skip the background
        # If we are running as a larger sofia run update the catalogue and cubelets
        if not cfg is None and stil_split:
            update_sofia_catalogue(cfg,cube_name=cube_name,
                base_name=cfg.internal.sofia_basename,
                outdir=outdir,base_dir=base_dir,) 
        else:
            simple_copy = True
    else:
        simple_copy = True
    if simple_copy:
        if cfg.general.verbose:
            print(f'''No final mask was created, skipping the source splitting. 
Copying Watershed directory to Watershed_Output_{sofia_id} in {cfg.internal.sofia_basename}_cubelets''')
        newname = f'{cfg.internal.sofia_directory}/{cfg.internal.sofia_basename}_cubelets/Watershed_Output_{sofia_id}/'
        shutil.rmtree(newname) if os.path.exists(newname) else None
        os.rename(outdir, newname)
#@profile


#