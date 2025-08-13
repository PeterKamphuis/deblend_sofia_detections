

from deblend_sofia_detections import report_version

from astropy.table import QTable
from astropy.io import fits

from multiprocessing import get_context

import numpy as np

def filter_peaks(maxima, border_width = None,npeaks=np.inf,
        previous_deblend=None):
      

    z_peaks = []
    y_peaks = []
    x_peaks = []
    mask_values = []
    peak_values = []
    for peaks in maxima:
        if previous_deblend is not None:
            if not np.isnan(previous_deblend[peaks[0],peaks[1],peaks[2]]):
            # if we have a previous deblend we use the mean of the previous deblend
            # in the vicinity of the peak to determine the current source
                # this is to avoid that we have a lot of sources with only one pixel
                # which is not what we want
                # We use the mean of the previous deblend in the vicinity of the peak
                current = previous_deblend[peaks[0],peaks[1],peaks[2]]
                print(f"Using previous deblend value {current} for peak {peaks}")
                if current in mask_values:
                    current = 0.
            else:
                current = 0.
        else:
            if len(mask_values) == 0:
                current = 1
            else:
                current = np.max(mask_values)+1

        if current != 0.:
            mask_values.append(current)
            peak_values.append(peaks[3])
            z_peaks.append(peaks[0])
            y_peaks.append(peaks[1])
            x_peaks.append(peaks[2])
        if len(mask_values) >= npeaks:
            break

   
    meta = {'version': report_version(),}
    colnames = ['z_peak', 'y_peak', 'x_peak', 'peak_value', 'mask_values']
    coldata = [z_peaks, y_peaks, x_peaks, peak_values,mask_values]
    table = QTable(coldata, names=colnames, meta=meta)
  
    return table




def find_peaks(cfg,data, threshold, box_size=[3,3,3], mask=None,
               border_width=None, npeaks=np.inf, previous_deblend=None,
               num_processes=6,outdir='./',cube_header=None):
    local_maxima = []
    shape = data.shape
    print(f'Data shape {shape}')
    #maskd the data if a mask is present
    if mask is not None:
        data[mask < 0.5] = float('NaN')
    #box_size = [z,y,x]
    boxside = np.array(np.ceil(np.array(box_size) / 2.0),dtype=int)
    box = np.array([[boxside[0]-1, boxside[0]], [boxside[1]-1,boxside[1]], 
                    [boxside[2]-1,boxside[2]]],dtype=int)
   
    with get_context("spawn").Pool(processes=num_processes) as pool:
        local_maxima = pool.starmap(is_local_maxima, [(data, threshold,box, z, y, x)
            for z in range(0, shape[0]) for y in range(0, shape[1]) 
            for x in range(0, shape[2])])
    local_maxima = np.array(local_maxima).reshape(shape[0], shape[1], shape[2])    
    #print(local_maxima)
    coords = np.where(local_maxima == True)
   
    local_maxima_coords = [[int(z), int(y), int(x), float(data[z, y, x])] 
            for z, y, x in zip(*coords)]
    test_cube = np.zeros(shape)
    local_maxima_coords = sorted(local_maxima_coords, key=lambda x: x[3], reverse=True)
    
    peak_table = filter_peaks(local_maxima_coords, border_width=border_width,
        previous_deblend=previous_deblend,npeaks=npeaks)


  #Make a new markers map based on the peaks
    markers3d = np.zeros(shape).astype(np.int8)

    for peak in peak_table:
        z,y,x = peak["z_peak"], peak["y_peak"], peak["x_peak"]
        # place markers in the vicinity of the flux peaks
        # We  use mask_smooth here because it provides the correct mask previous mask number 
        # for peak, i.e. peaks are grouped with the old source detection. 
        if cfg.general.verbose:
            print(f"Placing markers for peak {peak['mask_values']} at {x},{y},{z} ")

        markers3d[z-box[0,0]:z+box[0,1],y-box[1,0]:y+box[1,1],
            x-box[2,0]:x+box[2,1]] = peak["mask_values"]

        # you can add some markers manually
        # (manual markers here ...)
    fits.writeto(f"{outdir}peak3d_markers.fits",markers3d,
                header=cube_header, overwrite=True)

    return peak_table,markers3d
      


def is_local_maxima(arr,threshold,box,z, y, x):
    if arr[z, y, x] <= threshold[z, y, x]:
        return False

    if np.isnan(arr[z, y, x]):
        return False

    subarray = arr[z-box[0,0]:z+box[0,1], y-box[1,0]:y+box[1,1], x-box[2,0]:x+box[2,1]]

    if np.isnan(subarray).all():
        return False
    elif np.nanmax(subarray) == arr[z, y, x]:
        return True
    else:
        return False