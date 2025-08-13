from deblend_sofia_detections.catalogue.search import \
    search_counter_part
from deblend_sofia_detections.config.functions import setup_config
from deblend_sofia_detections.deblending.sofia_functions import \
    load_sofia_input_file,set_sofia,write_sofia,read_sofia_table,\
    execute_sofia
from deblend_sofia_detections.support.system_functions import \
    create_directory

from astropy.convolution import convolve,Gaussian1DKernel
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.nddata.utils import NoOverlapError
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astroquery.gaia import Gaia

from photutils.aperture import EllipticalAperture
from photutils.background import Background2D # Background2D is used for background subtraction

import astropy.units as u

import copy
import numpy as np
import os

def cut_optical(hdr_over,wcs,dir,image):
    '''Cut out the optical image'''
    #load a smaller part from a larger fits image
    optical_image=fits.open(f'{dir}/{image}')
   
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
    sizecut = np.max([hdr_over['NAXIS1'], hdr_over['NAXIS2']])
    centralpix = [hdr_over['NAXIS1']/ 2., hdr_over['NAXIS2']/ 2.]
    rascr, decscr = wcs.wcs_pix2world(*centralpix,1.)
    obj_coords = SkyCoord(ra= rascr* u.degree, dec=decscr * u.degree, frame='fk5')
    size = u.Quantity((sizecut* 3600 * abs(hdr_over['CDELT2']),\
                       sizecut* 3600 * abs(hdr_over['CDELT2'])), u.arcsec)
    try:
        optical_cutout = Cutout2D(data, obj_coords, size, wcs=opt_wcs)
    except NoOverlapError:
        print(f'No overlap between the optical image and the SOFIA image')
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
        del kernel
    else: cube_new = cube  # do nothing
    return cube_new


def get_background(img=None, wcs_opt=None, optical_name=None, match_header=None,
                   wcs=None):   
    if img is not None:
        bckgrnd = img
        bckgrnd_wcs = wcs_opt
    else:
        optdir,optfile = os.path.split(optical_name)
      
        optical = cut_optical(match_header,wcs,\
                optdir,optfile)
        bckgrnd_wcs = optical.wcs
        bckgrnd = optical.data
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
    Gaia.ROW_LIMIT = -1
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
        gaia_table = Gaia.query_object_async(coords, width=w*pixel_scale, height=h*pixel_scale)
    else:
        gaia_table = Table.read(gaia_table)
   
    #Remove galaxy canditates
    gaia_table = gaia_table[gaia_table['in_galaxy_candidates'] == False] 
    #gaia_table = gaia_table[gaia_table['in_qso_candidates'] == False]    
    #gaia_table = gaia_table[gaia_table['non_single_star'] == 0] 
    gaia_table.sort('phot_rp_mean_mag')
 
    gaia_table = gaia_table[0:int(len(gaia_table)*0.5)]
    if cfg.general.verbose:
        print(f"Found {len(gaia_table)} Gaia sources in the image area.")
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
    valid_mask = ((x_arr >= -max_radius) & (x_arr < width + max_radius) & 
                  (y_arr >= -max_radius) & (y_arr < height + max_radius))
    
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

def mask_source_from_table(optical_image,optical_wcs,mask=None, src_table = None):
    if mask is None:
        masked_deb = np.full_like(optical_image, 0)
    else:
        masked_deb = copy.deepcopy(mask)

    if src_table is None:
        print("No source table provided. Not adding to the mask")
        return masked_deb
    
    # input source table (e.g., SGA2020)
    #src_table = Table(names=['ra', 'dec', 'PA', 'sma', 'e'],
    #                data=np.array([[158.9368, -28.7691, 107.8, 23.8, 0.36], 
    #                                [158.9026, -28.7686, 154.7, 24.7, 0.65]]))

    pixel_scale = proj_plane_pixel_scales(optical_image.wcs)[0] * u.deg
    seg_start = np.max(masked_deb) if np.any(masked_deb) else 0
    for i in range(len(src_table)):
        gal_coord = SkyCoord(ra=src_table["ra"][i], dec=src_table["dec"][i], unit='deg')
        xcen, ycen = optical_wcs.world_to_pixel(gal_coord)
        aper = EllipticalAperture((xcen, ycen), 
                                a=src_table["sma"][i] / pixel_scale, 
                                b=src_table["sma"][i] * (1-src_table["e"][i]) / pixel_scale, 
                                theta=(90+src_table["PA"][i])*u.deg)
        segment = aper.to_mask(method='center')
        masked_deb[int(ycen)-segment.shape[0]//2+1:int(ycen)+1-segment.shape[0]//2+segment.shape[0],
                int(xcen)-segment.shape[1]//2+1:int(xcen)+1-segment.shape[1]//2+segment.shape[1]]\
                = segment.data * (i + 1+ seg_start)    
    return masked_deb


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
        create_directory('Sofia_Output',outdir)

     

    sofia_temp= set_sofia(sofia_temp, cube_name, mask,outdir) 


    write_sofia(sofia_temp,f'{outdir}/Sofia_Output/sofia_input.par')
    #Run Sofia
    matched = False
    while not matched:
        # Rune sofia
        execute_sofia(cfg,run_directory=f'{outdir}/Sofia_Output/')
        #read the ouput table
        if cfg.general.verbose:
            print(f"Reading the SoFiA output table from {outdir} the cube {name}")
        split_sources,split_base_name,split_table_name =  read_sofia_table(cfg, 
            sofia_directory=f'{outdir}/Sofia_Output/',sofia_basename=basename,
            no_conversion=False) 
        id = []
        replace_id = []
        present_id = [int(x) for x in split_sources['id']]
        counter = 0
        for source in split_sources:
           
            watername= f'{split_base_name}'
        
            source = search_counter_part(cfg,source,basename=watername,
                query = 'NED',sofia_directory=f'{outdir}/Sofia_Output/',
                insource='sofia')
          
            source = search_counter_part(cfg,source,basename=watername,
                    query='Manual',sofia_directory=f'{outdir}/Sofia_Output/')
            if source['Manual_spectroscopic']:
                source['Name'] =  source['Manual_Name']
            elif source['NED_spectroscopic']:
                source['Name'] =  source['NED_Object Name']  
            else:
                source['Name'] =  source['sofia_name']
            if cfg.general.verbose:
                print(f"Processing source {source['Name'][0]} with id {source['sofia_id'][0]}")

            if source['Name'] == source['sofia_name']:
                if len(id) == 0:
                    rep = np.min(present_id)
                    if int(rep) == int(source['sofia_id']):
                        rep = np.max(present_id)
                else:
                    rep = id[-1]
                replace_id.append([int(source['sofia_id']),int(rep)])
            else:
                id.append(source['sofia_id'])
        if cfg.general.verbose:
            print(f"Found {len(id)} sources with a counterpart in the catalogue.")
            print(id,replace_id)
        counter += 1
        if counter > 50:
            print(f"Warning: More than 50 times matching counterpart for {name}.")
            matched = True
        maskin= fits.open(mask)
        if len(id) == len(split_sources):
            matched = True
        elif np.unique(maskin[0].data).size-1 == 1:
            if cfg.general.verbose:
                print(f"Only one source found in the mask {mask}. No deblending needed.")
            matched = True
        else:
            
            for pair in replace_id:
                maskin[0].data[maskin[0].data == pair[0]] = pair[1]
            maskin.writeto(mask, overwrite=True)
        
    
    if np.unique(maskin[0].data).size-1 == 1 or counter  > 50:
        del maskin
        return False
    else:
        del maskin
        return True
    
def subtract_background(image,wcs):
    """
    Subtracts the background from an image using a 2D background estimation.
    
    Parameters:
    image (np.ndarray): The input image data.
    wcs (astropy.wcs.WCS): The WCS of the image.
    
    Returns:
    np.ndarray: The image with the background subtracted.
    """
    pixel_scale = np.mean(abs(proj_plane_pixel_scales(wcs)))*u.deg
    fwhm = 5./ pixel_scale.to(u.arcsec).value*2.31
    boxin = int(3.*fwhm) + 1 if int(3.*fwhm) % 2 == 0 else int(3.*fwhm)
    
    box_size = [boxin, boxin]  # box size for background estimation
    background = Background2D(image, box_size)
    
    return image - background.background    

     