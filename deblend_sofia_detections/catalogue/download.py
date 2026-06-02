# -*- coding: future_fstrings -*-
#Functions that look for optical image

from deblend_sofia_detections.deblending.image_manipulation import cut_optical
from deblend_sofia_detections.support.errors import DownloadError
from deblend_sofia_detections.support.logging import print_log
from astropy import units as u
from astropy.io import fits
from astroquery.skyview import SkyView
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import urllib.request
import os
import warnings
import numpy as np



def creating_full_FOV_optical(cfg):
    SkyView.URL = 'https://skyview.gsfc.nasa.gov/current/cgi/basicform.pl'
    #SkyView.URL = 'https://skyview.gsfc.nasa.gov/current/cgi/query.pl'
   
    print_log(cfg, f'Quering the Sky Survey', case=['verbose'])
    #First we open the moment header 0 to get the extend of the field 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")     
        mom0_header = fits.getheader(f'{cfg.sofia.directory}/{cfg.sofia.basename}_mom0.fits')
        mom0_wcs = WCS(mom0_header).celestial
    #set the size of the image
    size = np.nanmax([abs(mom0_header['NAXIS1']*mom0_header['CDELT1'])*60.,
                      abs(mom0_header['NAXIS2']*mom0_header['CDELT2'])*60.,])
    size_quantity= u.Quantity(size,u.arcmin)
    #get_image_list seems to guess at the pixel size let's fix it to 3 arcsec 
    beam = mom0_header['BMAJ'] * u.deg
    optical_pixel_scale = beam.to(u.arcsec).value/cfg.general.optical_pixel_scale * u.arcsec
    #If the resolution of our optical images is greater than 5 the deblending becomes hairy 
    if optical_pixel_scale > 4.*u.arcsec:
        optical_pixel_scale = 4. * u.arcsec
    size_pixels = (size_quantity.to(u.arcsec).value/optical_pixel_scale.value).astype(int)
    #obtain the central coordinates
    ra,dec = mom0_wcs.wcs_pix2world(mom0_header['NAXIS1']/2., mom0_header['NAXIS2']/2.,1.)
    obj_coords = SkyCoord(ra= ra* u.degree, dec= dec * u.degree, frame='fk5')

  
    if not cfg.input.manual_optical_image[0] is None:
        
        print_log(cfg, f'Checking manual input', case=['verbose'])
        for identifier in cfg.input.manual_optical_image:
            if os.path.isfile(identifier):
                print_log(cfg, f'Found manual optical image: {identifier}', case=['verbose'])
            manual_path,manual_file = os.path.split(identifier)
            if manual_path == '':
                manual_path = './'
            cutout = cut_optical(mom0_header,mom0_wcs,\
                manual_path,
                manual_file)

            cutout_hdr = cutout.wcs.to_header()
            cutout_hdr['COMMENT'] =  f'The original file was  {identifier}'
            fits.writeto(cfg.internal.optical_background,cutout.data,cutout_hdr,overwrite=True)

            return


    i
    print_log(cfg, f'''Obtaining the actual image list with the following parameters:
object coordinates: {obj_coords.to_string('hmsdms')},
radius: {size_quantity},
pixels: {size_pixels},''', case=['verbose'])
    SkyView.clear_cache()
    sky_view_list = SkyView.get_image_list(position=obj_coords,
                                   radius = size_quantity,
                                   coordinates = "J2000",
                                   pixels = size_pixels,
                                   cache=True,
                                   #survey=['WISE 3.4'])
                                   #survey=['2MASS-K'])
                                   survey=['DSS2 Red'])

 
    print_log(cfg, f'''Obtained {len(sky_view_list)} images from SkyView
{sky_view_list}
              ''', case=['verbose'])
    for path in sky_view_list:
        print_log(cfg, f'''Starting the download for {cfg.sofia.original_data_cube}.
This can take a while.''', case=['verbose'])
        filename = path.split("/")[-1]
        print_log(cfg, f'Downloading {filename} from SkyView from {path}')
        urllib.request.urlretrieve(path, f'{cfg.internal.optical_background}')
        if os.path.isfile(f'{cfg.internal.optical_background}'):
            print_log(cfg, "Successfully downloaded the image")
            #os.replace(filename, f'{cfg.internal.ancillary_directory}moment0_full_DSS.fits')
        else:
            print_log(cfg, "Failed to obtain the image from SkyView")
            raise DownloadError(f'''Failed to download the image from SkyView: {path}
Check your internet connection and the SkyView service status.
Note that redownloading the exact same image may fail if it has recently been removed from the SkyView archive.
''')
