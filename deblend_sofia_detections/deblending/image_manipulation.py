
from deblend_sofia_detections.catalogue.search import \
    search_counter_part

from deblend_sofia_detections.deblending.sofia_functions import \
    load_sofia_input_file,set_sofia,write_sofia,read_sofia_table,\
    execute_sofia,closest_sofia_source
from deblend_sofia_detections.support.system_functions import \
    create_directory
from deblend_sofia_detections.support.support_functions import \
    close_variables,get_channel_width
from deblend_sofia_detections.support.table_functions import check_table_length
from deblend_sofia_detections.support.logging import print_log


from astropy.convolution import convolve,Gaussian1DKernel
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.nddata.utils import NoOverlapError
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales


from photutils.background import Background2D # Background2D is used for background subtraction

from deblend_sofia_detections.support.profiling import profile


import astropy.units as u
import copy
import numpy as np
import os
import pickle
import shutil



def add_to_original(original_data, cut_data,sofia_id = 1):
    original_wcs = WCS(original_data[0].header)
    cut_wcs = WCS(cut_data[0].header)
    cut_origin = cut_wcs.wcs_pix2world(0,0,0,1)
    original_coord = original_wcs.wcs_world2pix(*cut_origin,1.)
    expanded_new = np.zeros_like(original_data[0].data)
    expanded_new[int(original_coord[2]):int(original_coord[2])+cut_data[0].data.shape[0],
              int(original_coord[1]):int(original_coord[1])+cut_data[0].data.shape[1],
              int(original_coord[0]):int(original_coord[0])+cut_data[0].data.shape[2]] = cut_data[0].data
    original_data[0].data[original_data[0].data == int(sofia_id)] = 0
    original_data[0].data[expanded_new > 0] = expanded_new[expanded_new > 0]
    return original_data



def cut_optical(cfg,hdr_over,wcs,imdir,image):
    '''Cut out the optical image'''
    #load a smaller part from a larger fits image
    optical_image=fits.open(f'{imdir}/{image}',verify_output='ignore')
   
    try:
        hdr = optical_image[0].header
        data = optical_image[0].data
    except:
        try:
            hdr = optical_image.header
            data = optical_image.data
        except:
            return None
    opt_wcs= WCS(hdr)
    sizecut = [hdr_over['NAXIS1'], hdr_over['NAXIS2']]
    centralpix = [hdr_over['NAXIS1']/ 2., hdr_over['NAXIS2']/ 2.]
    rascr, decscr = wcs.wcs_pix2world(*centralpix,1.)
    obj_coords = SkyCoord(ra= rascr* u.degree, dec=decscr * u.degree, frame='fk5')
    size = u.Quantity((sizecut[1]* 3600 * abs(hdr_over['CDELT2']),\
                       sizecut[0]* 3600 * abs(hdr_over['CDELT1'])), u.arcsec)
    try:
        optical_cutout = Cutout2D(data, obj_coords, size, wcs=opt_wcs)
    except NoOverlapError:
        print_log(cfg,f'No overlap between the optical image and the SOFIA image',case=['verbose'])
        optical_cutout = None
    optical_image.close()
    return optical_cutout



def freq_smooth(cube, bin_size=0, smooth=0):
    '''
    bin or smooth the data cube along the frequency axis.
    '''
    if bin_size:    # bin
        shape = cube.shape
        cube_new = np.zeros((shape[0]//bin_size, shape[1], shape[2]))
        for freq in range(cube.shape[0]//bin_size):
            cube_new[freq] = np.sum(cube[freq*bin_size : (freq+1)*bin_size], axis=0)
    elif smooth:    # smooth
        cube_new = np.zeros_like(cube)
        kernel = Gaussian1DKernel(smooth)
        for i in range(cube.shape[1]):
            for j in range(cube.shape[2]):
                cube_new[:,i,j] = convolve(cube[:,i,j], kernel)
        close_variables(kernel)
    else: cube_new = cube  # do nothing
    return cube_new


def get_background(cfg, match_header=None, wcs=None):   
    optdir,optfile = os.path.split(cfg.internal.optical_background)
    optical = cut_optical(cfg,match_header,wcs,optdir,optfile)
    bckgrnd_wcs = optical.wcs
    bckgrnd = optical.data
    bckgrnd = bckgrnd.astype(np.float32)
    return bckgrnd, bckgrnd_wcs


def mask_gaia_stars(cfg,optical_image, optical_wcs):
    """
    Masks Gaia stars in an astronomical FITS image.
    optical_image: The cutout image containing optical data.
    optical_wcs: The WCS of the optical image.
    radius_arcsec: The radius in arcseconds to mask around each star.
    We are always using the pickled and dowloaded table for the whole image.
    If it doesn't exist we never run this function.
    """
    star_mask = np.zeros_like(optical_image.data).astype(bool) 
    
    # Astro queries are woefully unstable so we only import them when we need them
    
  
    # Run astroquery.
    # This may take some time for large images. 
    # In such case, you can save this table and reload it next time.
    h,w = optical_image.data.shape
    coords = optical_wcs.pixel_to_world(h/2.-0.5, w/2.-0.5)
    pixel_scale = np.mean(proj_plane_pixel_scales(optical_wcs))*u.deg
    radius_pixels = 5./ pixel_scale.to(u.arcsec).value  # Convert arcsec to pixels
   
    print_log(cfg,f''' We find a pixel scale of {pixel_scale.to(u.arcsec)}
Which means we use a basic masking radius of {radius_pixels} pixels.''',
        case=['verbose'])
    # We have already matched the optical image to the size of our HI detection so we can just search that area
    # Query Gaia catalog
 
    with open(f'{cfg.internal.gaia_table}','rb') as tmp:
            gaia_table = pickle.load(tmp)
   
    if len(gaia_table) == 0:
        print_log(cfg,"No Gaia sources found in the full image area. Returning an empty mask.",case=['verbose'])
        return star_mask, True   
    else:
        print_log(cfg,f"Found {len(gaia_table)} Gaia sources in the full image area.",case=['verbose'])
    # generate star masks

    ra_range = [optical_wcs.wcs_pix2world(0, 0, 1)[0], optical_wcs.wcs_pix2world(w, h, 1)[0]]
    ra_range = [min(ra_range), max(ra_range)]  # Ensure ra_range is in increasing order
    dec_range = [optical_wcs.wcs_pix2world(0, 0, 1)[1], optical_wcs.wcs_pix2world(w, h, 1)[1]]
   
    # exclude stars outside the image bounds
    gaia_table = gaia_table[(gaia_table["ra"] >= ra_range[0]-0.1) & (gaia_table["ra"] <= ra_range[1]+0.1) &
                            (gaia_table["dec"] >= dec_range[0]-0.1) & (gaia_table["dec"] <= dec_range[1]+0.1)]
    gaia_table.sort('phot_rp_mean_mag')  # Sort by brightness (smaller magnitude is brighter)
    gaia_table = gaia_table[0:int(len(gaia_table)*.5)]  # Limit to the upper half of the brightness distribution
    star_coords = SkyCoord(ra=gaia_table["ra"], dec=gaia_table["dec"],
                            unit=(u.deg, u.deg), frame='fk5')
   
    x, y = optical_wcs.world_to_pixel(star_coords)
    # these are magnitude so the smaller they are the brighter the star is
    individual_radius = (np.median(gaia_table["phot_rp_mean_mag"])/
        gaia_table["phot_rp_mean_mag"])**3 * radius_pixels  # Example scaling factor for radius
    individual_radius[individual_radius <radius_pixels] = radius_pixels

    # Mask stars
    print_log(cfg,f"Masking {len(x)} stars in the optical image.",case=['verbose','screen'])
    yy, xx = np.indices(optical_image.data.shape)
    
    # Memory-efficient chunked approach
    x_arr = np.array(x)
    y_arr = np.array(y)
    r_arr = np.array(individual_radius)
    
    # Filter stars that are within the image bounds + maximum radius
    max_radius = np.max(r_arr)
    height, width = optical_image.data.shape
    
    # Pre-filter stars that could possibly affect the image
    #valid_mask = ((x_arr >= -max_radius) & (x_arr < width + max_radius) & 
    #              (y_arr >= -max_radius) & (y_arr < height + max_radius))
    
    valid_mask = ((x_arr >= 0) & (x_arr < width) & 
                  (y_arr >= 0) & (y_arr < height))

    x_valid = x_arr[valid_mask]
    y_valid = y_arr[valid_mask]
    r_valid = r_arr[valid_mask]
    


    print_log(cfg,f"Processing {len(x_valid)} valid stars out of {len(x_arr)} total stars",case=['verbose'])
    
    # Process stars in chunks to avoid memory issues
    chunk_size = min(35, len(x_valid))  # Adjust based on available memory
    
    for i in range(0, len(x_valid), chunk_size):
        if cfg.logging.enable:
            print(f"\r Processing chunks  {(i/len(x_valid))*100.:.1f} % done. ",\
                  end=" ", flush=True)
            
        end_idx = min(i + chunk_size, len(x_valid))
        
        # Get chunk of stars
        x_chunk = x_valid[i:end_idx]
        y_chunk = y_valid[i:end_idx]
        r_chunk = r_valid[i:end_idx]
        
        # Vectorized computation for this chunk
        x_stars = x_chunk[:, np.newaxis, np.newaxis]
        y_stars = y_chunk[:, np.newaxis, np.newaxis]
        r_stars = r_chunk[:, np.newaxis, np.newaxis]
        
        # Compute distances for this chunk
        distances_sq = (xx[np.newaxis, :, :] - x_stars)**2 + (yy[np.newaxis, :, :] - y_stars)**2
        chunk_masks = distances_sq < r_stars**2
        
        # Update the star mask with OR operation
        star_mask |= np.any(chunk_masks, axis=0)
 
    # Mask the stars in the optical image
   
    print_log(cfg,f'''\r Processing chunks 100.0 % done.\n 
Created the star mask to the optical image.''',case=['screen'])
    return star_mask, False

def mask_source_from_table(cfg,optical_markers,optical_header,mask=None, 
        src_table = None):
    optical_wcs = WCS(optical_header)
    if mask is None:
        masked_deb = np.full_like(optical_markers, 0)
    else: 
        masked_deb = copy.deepcopy(mask)
 
    if src_table is None:

        print_log(cfg,"No source table provided. Not adding to the mask"
            ,case=['verbose'])
        return optical_markers
    # input source table (e.g., SGA2020)
    #src_table = Table(names=['ra', 'dec', 'PA', 'sma', 'e'],
    #                data=np.array([[158.9368, -28.7691, 107.8, 23.8, 0.36], 
    #                                [158.9026, -28.7686, 154.7, 24.7, 0.65]]))
   
    pixel_scale = proj_plane_pixel_scales(optical_wcs)[0] * u.deg
    seg_start = np.max(optical_markers) if np.any(optical_markers) else 1
    if 'PA' not in src_table.colnames:
            
        print_log(cfg,"No PA column found in the source table. Adding a default PA column with 0 degrees."
            ,case=['debug'])
        src_table['PA'] = np.full(len(src_table), 0.) * u.deg
  
    maj_sizes = ["sma",'major_axis','maj_ang_size','maj_angsize']
    min_sizes = ["smb",'minor_axis','min_ang_size','min_angsize','e','ellipticity']
    source_counter = 1
    # Pre-filter rows with NaN coordinates/PA and batch-convert to pixel coords
    all_indices = np.array(range(check_table_length(src_table)))
    valid = all_indices[~(np.isnan(src_table["RA"][all_indices]) |
                          np.isnan(src_table["DEC"][all_indices]) |
                          np.isnan(src_table["PA"][all_indices]))]
    if len(valid) > 0:
        all_coords = SkyCoord(ra=src_table["RA"][valid], dec=src_table["DEC"][valid], unit='deg')
        all_xcens, all_ycens = optical_wcs.world_to_pixel(all_coords)
        in_bounds = ((all_xcens >= 0) & (all_xcens <= masked_deb.shape[1]) &
                     (all_ycens >= 0) & (all_ycens <= masked_deb.shape[0]) &
                     ~np.isnan(all_xcens) & ~np.isnan(all_ycens))
        valid = valid[in_bounds]
        all_xcens = all_xcens[in_bounds]
        all_ycens = all_ycens[in_bounds]
    else:
        all_xcens = all_ycens = np.array([])
    for loop_idx, i in enumerate(valid):
        xcen = all_xcens[loop_idx]
        ycen = all_ycens[loop_idx]

        sma = 10.* pixel_scale.to(u.arcsec)
        for size in maj_sizes:
            if size in src_table.colnames:
                if not np.isnan(src_table[size][i]):
                    sma = src_table[size][i].to(u.arcsec)
                    break
        smb = float('NaN')
        for size in min_sizes:
            if size in src_table.colnames:
                if not np.isnan(src_table[size][i]):
                    if size in ['e','ellipticity']:
                        smb = sma * (1-src_table[size][i])
                    else:
                        smb = src_table[size][i].to(u.arcsec)
                    break
        if np.isnan(smb):
            print_log(cfg,f"No valid minor axis size found for source {src_table['RA'][i], src_table['DEC'][i]}. Defaulting to a circle with radius =  {sma}."
                ,case=['debug'])
            smb = sma   
        print_log(cfg,f"Processing source {src_table['RA'][i], src_table['DEC'][i]} with PA {src_table['PA'][i]}, sma {sma} and smb {smb}. Pixel coordinates {xcen, ycen}"
                ,case=['debug'])    
        #if we have a source at the location we remove it for maximum control
        if optical_markers[int(ycen),int(xcen)] != 0:
            id = optical_markers[int(ycen),int(xcen)]
            optical_markers[optical_markers == id] = 0
            
        sma_pix = abs(sma.to(u.arcsec).value / pixel_scale.to(u.arcsec).value)
        smb_pix = abs(smb.to(u.arcsec).value / pixel_scale.to(u.arcsec).value)
        yy, xx = np.indices(masked_deb.shape)
        theta = (src_table["PA"][i].to(u.rad)).value+np.radians(90)  # Convert PA to radians and add 90 degrees to align with the major axis
        dx = xx - xcen
        dy = yy - ycen

        xrot = dx * np.cos(theta) + dy * np.sin(theta)
        yrot = -dx * np.sin(theta) + dy * np.cos(theta)

        ellipse = np.isfinite(xrot) & np.isfinite(yrot) & \
            ((xrot / sma_pix) ** 2 + (yrot / smb_pix) ** 2 <= 1)
     
        optical_markers[ellipse] = source_counter + seg_start
        source_counter += 1
    optical_markers[mask == 0] = 0
    return optical_markers


def split_sources(cfg_in,cube_name, mask, 
        outdir='./', catalogue = False):
    """
    Split the sources in the deblended 3D data cube.
    
    Parameters:
    cube (astropy.io.fits.HDUList): The data cube to split.
    res3d (np.ndarray): The 3D segmentation map.
    dir (str): The directory to save the results.
    catalogue (bool): If True, save the source catalogue.
    """

    
    cfg = copy.deepcopy(cfg_in) #Making sure to avoid feedback
    path,name = os.path.split(cube_name)
    basename = os.path.splitext(name)[0]
    sofia_temp = load_sofia_input_file()
    
    if not os.path.exists(f'{outdir}/Sofia_Output/'):
        create_directory('Sofia_Output',base_directory=outdir)

   
    shutil.copy2(mask,f"{outdir}/Sofia_Output/tmp_mask.fits")
    sofia_temp= set_sofia(sofia_temp, cube_name, f"{outdir}/Sofia_Output/tmp_mask.fits",outdir) 


    write_sofia(sofia_temp,f'{outdir}/Sofia_Output/sofia_input.par')
    #Run Sofia
    matched = False
    counter = 0
    maskhdr= fits.getheader(f"{outdir}/Sofia_Output/tmp_mask.fits",verify_output='ignore')
    header_info = {'pixelsize': float(np.mean([abs(maskhdr['CDELT1']),
                    abs(maskhdr['CDELT2'])]))*u.deg,
                       'channel_width': get_channel_width(maskhdr)}
    close_variables(maskhdr)
    while not matched:
        # Run sofia
        print_log(cfg,f"Running SoFiA on {cube_name} with mask {mask} in {outdir}/Sofia_Output/sofia_input.par",
            case=['verbose','screen'])
        execute_sofia(cfg,run_directory=f'{outdir}/Sofia_Output/')
        #read the ouput table
       
        print_log(cfg,f"Reading the SoFiA output table from {outdir} the cube {name}",
            case=['verbose'])
        split_sources,table_name =  read_sofia_table(cfg, 
            sofia_directory=f'{outdir}/Sofia_Output/',sofia_basename=basename,
            no_conversion=False) 

        print_log(cfg,f"Read the SoFiA output table from {outdir} the cube {name} with {len(split_sources)} sources.",
            case=['verbose','screen'])
        if split_sources is None:
            raise ValueError(f"SoFiA did not produce an output table for {cube_name}. Please check the SoFiA output for errors.")
        id = []
        replace_id = []
        present_id = [int(x) for x in split_sources['id']]  
        counterparts = {}  
        watername= os.path.splitext(os.path.split(table_name)[-1])[0].split('_cat')[0]
        shutil.copy2(f"{outdir}/Sofia_Output/{watername}_mask.fits",f"{outdir}/Sofia_Output/tmp_mask.fits")
        for source in split_sources:
            print_log(cfg,f"Processing deblended source with id {source['id']} and name {source['name']} "
                ,case=['verbose','screen'])
            #First we check if we have previous iteration output

            source = search_counter_part(cfg,source,basename=watername,
                query = 'INTERNET',sofia_directory=f'{outdir}/Sofia_Output/',
                insource='sofia')
            source = search_counter_part(cfg,source,basename=watername,
                    query='Manual',sofia_directory=f'{outdir}/Sofia_Output/') 
         
            if source['Manual_spectroscopic'] and not source['Manual_Name'][0]\
                in [x for x in counterparts]:
                source['Name'] =  source['Manual_Name'][0]
            elif source['INTERNET_spectroscopic'] and not \
                source['INTERNET_Object Name'][0] in [x for x in counterparts]:
                source['Name'] =  source['INTERNET_Object Name'][0]  
            else:
                source['Name'] =  source['sofia_name'][0]
            source_row = source[0]
            if source_row['Name'] == source_row['sofia_name']:
                print_log(cfg,f"SPLIT_SOURCE: Source id {source_row['sofia_id']} with name {source_row['Name']} has no counterpart in the catalogue. Replacement needed."
                    ,case=['verbose','screen']) 
            else:
                counterparts[source_row['Name']] = source_row['sofia_id']
                print_log(cfg,f"SPLIT_SOURCE: Source id {source_row['sofia_id']} with name {source_row['Name']} has a counterpart in the catalogue. No replacement needed."
                    ,case=['verbose','screen'])
         
          
            if source_row['Name'] == source_row['sofia_name']:
                rep = closest_sofia_source(cfg,source_row['sofia_id'],split_sources,
                    header_info = header_info)
                replace_id.append([int(source_row['sofia_id']),int(rep)])
                if source_row['sofia_id'] == 4:
                    print(source_row,rep)
                    exit()
                 
            else:
                id.append(source_row['sofia_id'])
            
        print_log(cfg,f"Found {len(id)} sources with a counterpart in the catalogue.", case=['verbose','screen'])
        counter += 1
        if counter > 50:
            print_log(cfg,f"Warning: More than 50 matching counterparts for {name}.", case=['verbose'])
            matched = True
        
        maskin= fits.open(f"{outdir}/Sofia_Output/tmp_mask.fits",verify_output='ignore')
       
        #maskin= fits.open(mask)
        print_log(cfg,f"Found {np.unique(maskin[0].data).size-1} sources in the mask. Found {len(id)} sources with a counterpart in the catalogue."
            , case=['verbose','screen'])
            
        if len(id) == len(split_sources):
            print_log(cfg,f"The id and split_sources lengths match. No further deblending needed.", case=['verbose', 'screen'])
            matched = True
        elif np.unique(maskin[0].data).size-1 == 1:
            print_log(cfg,f"Only one source found in the mask {mask}. No deblending needed.", case=['verbose', 'screen'])
            matched = True
        else:
            print_log(cfg,f'''The mask has the following source {np.unique(maskin[0].data)}
the counterparts map {counterparts}''', case=['verbose','screen'])
            for pair in replace_id:
                print_log(cfg,f"Replacing source {pair[0]} with {pair[1]} in the mask.", case=['verbose','screen'])
                maskin[0].data[maskin[0].data == pair[0]] = pair[1]
            fits.writeto(f'{outdir}/Sofia_Output/tmp_mask.fits',maskin[0].data,maskin[0].header,
                 overwrite=True,output_verify='ignore')
     
    shutil.copy2(f'{outdir}/Sofia_Output/tmp_mask.fits',f'{outdir}/final_mask.fits')    
    
    if np.unique(maskin[0].data).size-1 == 1 or counter  > 50:
        ret_val= False
    else:
        ret_val = True
    close_variables(maskin,split_sources)
    return ret_val


@profile('profiler_logs/subtract_background.log')
def subtract_background(cfg,image,wcs):
    """
    Subtracts the background from an image using a 2D background estimation.
    
    Parameters:
    image (np.ndarray): The input image data.
    wcs (astropy.wcs.WCS): The WCS of the image.
    
    Returns:
    np.ndarray: The image with the background subtracted.
    """
  
    boxin = int(3.*cfg.internal.optical_kernel_fwhm) + 1 if\
        int(3.*cfg.internal.optical_kernel_fwhm) % 2 == 0\
        else int(3.*cfg.internal.optical_kernel_fwhm)
    
    box_size = [boxin, boxin]  # box size for background estimation
    background = Background2D(image, box_size)
    new_image = image - background.background
    new_wcs = copy.deepcopy(wcs)
    close_variables(image,background,wcs)
    return new_image,new_wcs

     