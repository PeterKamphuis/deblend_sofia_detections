from deblend_sofia_detections.catalogue.search import \
    search_counter_part

from deblend_sofia_detections.deblending.sofia_functions import \
    load_sofia_input_file,set_sofia,write_sofia,read_sofia_table,\
    execute_sofia
from deblend_sofia_detections.deblending.debug_visualization import \
    matched_counterpart_positions,trial_hi_components_from_table,\
    write_hi_component_debug_overlay_safely
from deblend_sofia_detections.support.system_functions import \
    create_directory
from deblend_sofia_detections.support.support_functions import \
    close_variables
from deblend_sofia_detections.support.optical_image import \
    celestial_numpy_axes,collapse_optical_data
from deblend_sofia_detections.support.errors import InputError


from astropy.convolution import convolve,Gaussian1DKernel
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.nddata.utils import NoOverlapError
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astroquery.gaia import Gaia

from photutils.background import Background2D # Background2D is used for background subtraction

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


def _celestial_numpy_axes(full_wcs, data_ndim):
    """Map the WCS celestial pixel axes to NumPy array-axis indices."""
    try:
        return celestial_numpy_axes(
            full_wcs.get_axis_types(),
            full_wcs.axis_correlation_matrix,
            full_wcs.pixel_n_dim,
            data_ndim,
        )
    except (AttributeError, IndexError, TypeError, ValueError):
        return None


def cut_optical(hdr_over,wcs,dir,image,debug=False,verbose=False):
    """Cut a 2-D celestial image from a 2-D or multi-plane optical FITS."""
    optical_path = os.path.join(dir, image)
    try:
        # Scaled integer FITS images (BSCALE/BZERO) cannot always be accessed
        # through Astropy's memory mapping, so load them normally here.
        optical_image = fits.open(optical_path, memmap=False)
    except (OSError, ValueError) as error:
        raise InputError(
            f"Could not open optical FITS image {optical_path}: {error}"
        ) from error

    try:
        image_hdu = next(
            (
                hdu
                for hdu in optical_image
                if hdu.data is not None and np.ndim(hdu.data) >= 2
            ),
            None,
        )
        if image_hdu is None:
            raise InputError(
                f"Optical FITS image {optical_path} contains no image HDU "
                "with at least two dimensions."
            )
        hdr = image_hdu.header.copy()
        data = np.asanyarray(image_hdu.data)
        try:
            full_wcs = WCS(hdr)
            opt_wcs = full_wcs.celestial
        except Exception as error:
            raise InputError(
                f"Optical FITS image {optical_path} has an invalid WCS: "
                f"{error}"
            ) from error
        if not full_wcs.has_celestial:
            raise InputError(
                f"Optical FITS image {optical_path} has no celestial WCS."
            )
        original_shape = data.shape
        celestial_axes = _celestial_numpy_axes(full_wcs, data.ndim)
        try:
            data, collapsed_axes = collapse_optical_data(
                data, celestial_numpy_axes=celestial_axes
            )
        except ValueError as error:
            raise InputError(
                f"Could not convert optical FITS image {optical_path} with "
                f"shape {original_shape} into a 2-D celestial image: {error}"
            ) from error
        if collapsed_axes and verbose:
            print(
                f"Optical image {optical_path} has shape {original_shape}. "
                f"Averaging non-celestial array axis/axes "
                f"{list(collapsed_axes)} produced the 2-D image "
                f"{data.shape}."
            )
    except Exception:
        optical_image.close()
        raise

    sizecut = [hdr_over['NAXIS1'], hdr_over['NAXIS2']]
    centralpix = [hdr_over['NAXIS1']/ 2., hdr_over['NAXIS2']/ 2.]
    rascr, decscr = wcs.wcs_pix2world(*centralpix,1.)
    obj_coords = SkyCoord(ra= rascr* u.degree, dec=decscr * u.degree, frame='fk5')
    size = u.Quantity((sizecut[1]* 3600 * abs(hdr_over['CDELT2']),\
                       sizecut[0]* 3600 * abs(hdr_over['CDELT1'])), u.arcsec)
    try:
        optical_cutout = Cutout2D(
            data, obj_coords, size, wcs=opt_wcs, copy=True
        )
    except NoOverlapError as error:
        raise InputError(
            f"Optical image {optical_path} does not overlap the SoFiA field."
        ) from error
    except Exception as error:
        raise InputError(
            f"Could not cut the 2-D optical image {optical_path} with shape "
            f"{data.shape}: {error}"
        ) from error
    finally:
        optical_image.close()
    return optical_cutout

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
    optical = cut_optical(
        match_header,
        wcs,
        optdir,
        optfile,
        debug=cfg.general.debug,
        verbose=cfg.general.verbose,
    )
    bckgrnd_wcs = optical.wcs
    bckgrnd = optical.data
    bckgrnd = bckgrnd.astype(np.float32)
    return bckgrnd, bckgrnd_wcs


def mask_gaia_stars(optical_image, optical_wcs, 
                    gaia_table=None, cfg= None,):
    """
    Masks Gaia stars in an astronomical FITS image.
    optical_image: The cutout image containing optical data.
    optical_wcs: The WCS of the optical image.
    radius_arcsec: The radius in arcseconds to mask around each star.
    gaia_table: Optional pre-loaded Gaia table. If None, it will be queried.
    """
    #Load the gaia table
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    #Do not set this to minus one as it can crash
    #
    Gaia.ROW_LIMIT = 50000
    #Gaia.clear_cache()
    # Run astroquery.
    # This may take some time for large images. 
    # In such case, you can save this table and reload it next time.
    h,w = optical_image.data.shape
    coords = optical_wcs.pixel_to_world(h/2.-0.5, w/2.-0.5)
    pixel_scale = np.mean(proj_plane_pixel_scales(optical_wcs))*u.deg
    
  
    radius_pixels = 5./ pixel_scale.to(u.arcsec).value  # Convert arcsec to pixels
    if cfg.general.verbose:
        print(f''' We find a pixel scale of {pixel_scale.to(u.arcsec)}
Which means we use a basic masking radius of {radius_pixels} pixels.''')
    # We have already matched the optical image to the size of our HI detection so we can just search that area
    # Query Gaia catalog
    if gaia_table is None:
        gaia_table = Gaia.query_object_async(coords, width=w*pixel_scale*1.5, height=h*pixel_scale*1.5)
    else:
        gaia_table = Table.read(gaia_table)
    if cfg.general.verbose:
        print(f"Found {len(gaia_table)} Gaia sources in the image area. Sorting them")
    #Remove galaxy canditates
    gaia_table = gaia_table[gaia_table['in_galaxy_candidates'] == False] 
    #gaia_table = gaia_table[gaia_table['in_qso_candidates'] == False]    
    #gaia_table = gaia_table[gaia_table['non_single_star'] == 0] 
    gaia_table.sort('phot_rp_mean_mag')
    gaia_table = gaia_table[0:int(len(gaia_table)*0.5)]
    if len(gaia_table) > 5000:
        if cfg.general.debug:
            print("Capping the gaia table at 5000.")
        gaia_table = gaia_table[0:5000]
  
       # generate star masks
    star_mask  = np.zeros_like(optical_image.data).astype(bool)
    if len(gaia_table) == 0:
        if cfg.general.verbose:
            print("No Gaia sources found in the image area. Returning an empty mask.")
        return star_mask   
    
    star_coords = SkyCoord(ra=gaia_table["ra"], dec=gaia_table["dec"],
                            unit=(u.deg, u.deg), frame='icrs')
    
    x, y = optical_wcs.world_to_pixel(star_coords)
    # these are magnitude so the smaller they are the brighter the star is
    individual_radius = (  np.median(gaia_table["phot_rp_mean_mag"])/
        gaia_table["phot_rp_mean_mag"])**3 * radius_pixels  # Example scaling factor for radius
    individual_radius[individual_radius <radius_pixels] = radius_pixels

    # Mask stars
    if cfg.general.verbose:
        print("Masking stars in the optical image.")
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
    if cfg.general.verbose:
        print(f"Processing {len(x_valid)} valid stars out of {len(x_arr)} total stars")
    
    # Process stars in chunks to avoid memory issues
    chunk_size = min(35, len(x_valid))  # Adjust based on available memory
    
    for i in range(0, len(x_valid), chunk_size):
        if cfg.general.verbose:
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
    if cfg.general.verbose:
        print(f"\r Processing chunks 100.0 % done.\n ")
        print("Created the star mask to the optical image.")
    return star_mask

def mask_source_from_table(cfg,optical_markers,optical_header,mask=None, 
        src_table = None):
    optical_wcs = WCS(optical_header)
    if mask is None:
        masked_deb = np.full_like(optical_markers, 0)
    else: 
        masked_deb = copy.deepcopy(mask)
 
    if src_table is None:
        if cfg.general.verbose:
            print("No source table provided. Not adding to the mask")
        return optical_markers
    # input source table (e.g., SGA2020)
    #src_table = Table(names=['ra', 'dec', 'PA', 'sma', 'e'],
    #                data=np.array([[158.9368, -28.7691, 107.8, 23.8, 0.36], 
    #                                [158.9026, -28.7686, 154.7, 24.7, 0.65]]))

    pixel_scale = proj_plane_pixel_scales(optical_wcs)[0] * u.deg
    seg_start = np.max(optical_markers) if np.any(optical_markers) else 1
    if 'PA' not in src_table.colnames:
        if cfg.general.debug:
            print("No PA column found in the source table. Adding a default PA column with 0 degrees.")
        src_table['PA'] = np.full(len(src_table), 0.) * u.deg
  
    maj_sizes = ["sma",'major_axis','maj_ang_size']
    min_sizes = ["smb",'minor_axis','min_ang_size','e','ellipticity']
    for i in range(len(src_table)):
        if any([np.isnan(src_table["RA"][i]), np.isnan(src_table["DEC"][i]),
                np.isnan(src_table["PA"][i])]):
            continue
        
        sma = 10 * pixel_scale.to(u.arcsec)
        for size in maj_sizes:
            if size in src_table.colnames:
                if not np.isnan(src_table[size][i]):
                    sma = src_table[size][i]
                    break
        smb = float('NaN')
        for size in min_sizes:
                if size in src_table.colnames:
                    if not np.isnan(src_table[size][i]):
                        if size in ['e','ellipticity']:
                            smb = sma * (1-src_table[size][i])
                        else:
                            smb = src_table[size][i]
                        break
        if np.isnan(smb):
            if cfg.general.debug:
                print(f"No valid minor axis size found for source {src_table['RA'][i], src_table['DEC'][i]}. Using a default value of {sma} arcsec.")
            smb = sma   
        gal_coord = SkyCoord(ra=src_table["RA"][i], dec=src_table["DEC"][i], unit='deg')
        xcen, ycen = optical_wcs.world_to_pixel(gal_coord)
        if cfg.general.debug:
            print(f"Processing source {src_table['RA'][i], src_table['DEC'][i]} with PA {src_table['PA'][i]}, sma {sma} and smb {smb}. Pixel coordinates {xcen, ycen}")
        if xcen > masked_deb.shape[1] or ycen > masked_deb.shape[0] or xcen < 0\
            or ycen < 0 or np.isnan(xcen) or np.isnan(ycen):
            if cfg.general.debug:
                print(f"Source {src_table['RA'][i], src_table['DEC'][i]} is outside the image bounds. Skipping.")
            continue
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

        optical_markers[ellipse] = i + 1 + seg_start
    optical_markers[mask == 0] = 0
    return optical_markers


def split_sources(cfg_in,cube_name, mask, 
        outdir='./', catalogue = False, counterpart_positions=None,
        automatic_counterpart_table=None, debug_overlay=None):
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



    sofia_temp= set_sofia(sofia_temp, cube_name, mask,outdir) 


    write_sofia(sofia_temp,f'{outdir}/Sofia_Output/sofia_input.par')
    #Run Sofia
    matched = False
    counter = 0
    while not matched:
        # Rune sofia
        execute_sofia(cfg,run_directory=f'{outdir}/Sofia_Output/')
        #read the ouput table
        if cfg.general.verbose:
            print(f"Reading the SoFiA output table from {outdir} the cube {name}")
        split_sources,table_name =  read_sofia_table(cfg, 
            sofia_directory=f'{outdir}/Sofia_Output/',sofia_basename=basename,
            no_conversion=False) 
        if split_sources is None:
            raise ValueError(f"SoFiA did not produce an output table for {cube_name}. Please check the SoFiA output for errors.")
        id = []
        replace_id = []
      
        present_id = [int(x) for x in split_sources['id']]
       
        watername= os.path.splitext(os.path.split(table_name)[-1])[0].split('_cat')[0]
        if counter == 0 and debug_overlay is not None:
            components = trial_hi_components_from_table(
                split_sources,
                cubelet_directory=(
                    f"{outdir}/Sofia_Output/{watername}_cubelets"
                ),
                basename=watername,
            )
            if components:
                component_overlay = {
                    key: debug_overlay[key]
                    for key in (
                        "optical_image_name",
                        "moment0_data",
                        "moment0_header",
                        "source_mask",
                        "marker_data",
                        "catalogue_positions",
                        "marker_mode",
                        "source_id",
                    )
                }
                component_overlay.update(
                    {
                        "components": components,
                        "output_name": (
                            f"{outdir}debug_products/"
                            "optical_hi_components_overlay_source_"
                            f"{debug_overlay['source_id']}.png"
                        ),
                    }
                )
                write_hi_component_debug_overlay_safely(
                    **component_overlay
                )
        for source in split_sources:
            if cfg.general.verbose:
                print(f"Processing source with id {source['id']} and name {source['name']}")
            source = search_counter_part(cfg,source,basename=watername,
                query = 'NED',sofia_directory=f'{outdir}/Sofia_Output/',
                insource='sofia')

            source = search_counter_part(
                cfg,
                source,
                basename=watername,
                query='DR10',
                sofia_directory=f'{outdir}/Sofia_Output/',
                automatic_catalogue=automatic_counterpart_table,
            )
         
            source = search_counter_part(cfg,source,basename=watername,
                    query='Manual',sofia_directory=f'{outdir}/Sofia_Output/')
            if counterpart_positions is not None:
                counterpart_positions.extend(
                    matched_counterpart_positions(source)
                )
            
          
            if source['Manual_spectroscopic']:
                source['Name'] =  source['Manual_Name'][0]
            elif source['DR10_counterpart']:
                source['Name'] = source['DR10_Name'][0]
            elif source['NED_spectroscopic']:
                source['Name'] =  source['NED_Object Name'][0]  
            else:
                source['Name'] =  source['sofia_name'][0]
            source_row = source[0]
            if cfg.general.verbose:
                print(f"Processing source {source_row['Name']} with id {source_row['sofia_id']}")

            if source_row['Name'] == source_row['sofia_name']:
                if len(id) == 0:
                    rep = np.min(present_id)
                    if int(rep) == int(source_row['sofia_id']):
                        rep = np.max(present_id)
                else:
                    rep = id[-1]
                replace_id.append([int(source_row['sofia_id']),int(rep)])
            else:
                id.append(source_row['sofia_id'])
        if cfg.general.verbose:
            print(f"Found {len(id)} sources with a counterpart in the catalogue.")
        counter += 1
        if counter > 50:
            if cfg.general.verbose:
                print(f"Warning: More than 50 matching counterparts for {name}.")
            matched = True
        
        maskin= fits.open(f"{outdir}/Sofia_Output/{watername}_mask.fits")
        #maskin= fits.open(mask)
        if cfg.general.verbose:
            print(f"Found {np.unique(maskin[0].data).size-1} sources in the mask. Found {len(id)} sources with a counterpart in the catalogue.")
      
        if len(id) == len(split_sources):
            matched = True
        elif np.unique(maskin[0].data).size-1 == 1:
            if cfg.general.verbose:
                print(f"Only one source found in the mask {mask}. No deblending needed.")
            matched = True
        else:
            
            for pair in replace_id:
                if cfg.general.verbose:
                    print(f"Replacing source {pair[0]} with {pair[1]} in the mask.")
                maskin[0].data[maskin[0].data == pair[0]] = pair[1]
            maskin.writeto(f'{outdir}final_mask.fits', overwrite=True)
        
        
    
    if np.unique(maskin[0].data).size-1 == 1 or counter  > 50:
        ret_val= False
    else:
        ret_val = True
    close_variables(maskin,split_sources)
    return ret_val

fn = open('profiler_logs/subtract_background.log', 'w+') if PROFILING else None
@profile(stream=fn)
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
