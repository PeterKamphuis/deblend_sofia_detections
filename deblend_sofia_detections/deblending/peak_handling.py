from tabnanny import verbose

from deblend_sofia_detections.support.logging import print_log
from deblend_sofia_detections.support.support_functions import write_fits_file
from astropy.table import QTable


from multiprocessing import get_context
from itertools import islice
import numpy as np
import datetime

from scipy.ndimage import maximum_filter

def filter_peaks(cfg,maxima,npeaks=np.inf,previous_deblend=None):
    n = min(len(maxima), int(npeaks))
    if previous_deblend is not None:
        prev = previous_deblend
        used= set()
        keep = []
        current = np.empty(len(maxima), dtype=float)
        current[:] = 0.0
        
        for i, (z, y, x, _) in enumerate(maxima):
            if np.isnan(prev[int(z), int(y), int(x)]):
                current = 0
            else:
                current = int(prev[int(z), int(y), int(x)])
                if current in used:
                    current = 0

            if current != 0:
                used.add(current)
                keep.append(i)

            if len(used) >= n:
                break
        
        maxima = maxima[keep]
        mask_values = np.array(
            [int(prev[int(z), int(y), int(x)]) for z, y, x, _ in maxima],
            dtype=int
            )
        
    else:
        # assign mask values in order
     
        maxima = maxima[:n]
        mask_values = np.arange(1, n + 1, dtype=int)

    peak_table = QTable(
        {
            "z_peak": maxima[:, 0].astype(int),
            "y_peak": maxima[:, 1].astype(int),
            "x_peak": maxima[:, 2].astype(int),
            "peak_value": maxima[:, 3],
            "mask_values": mask_values,
        }
    )  
    '''
    print(f'We found {peak_table} peaks in the cubelet')
    z_peaks = []
    y_peaks = []
    x_peaks = []
    mask_values = []
    peak_values = []
    for peaks in maxima:
        if previous_deblend is not None:
            if not np.isnan(previous_deblend[int(peaks[0]),int(peaks[1]),int(peaks[2])]):
            # if we have a previous deblend we use the mean of the previous deblend
            # in the vicinity of the peak to determine the current source
                # this is to avoid that we have a lot of sources with only one pixel
                # which is not what we want
                # We use the mean of the previous deblend in the vicinity of the peak
                current = previous_deblend[int(peaks[0]),int(peaks[1]),int(peaks[2])]
               
                print_log(cfg,f"Using previous deblend value {current} for peak {peaks}"
                    ,case=['verbose'])
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
            z_peaks.append(int(peaks[0]))
            y_peaks.append(int(peaks[1]))
            x_peaks.append(int(peaks[2]))
        if len(mask_values) >= npeaks:
            break

   
    meta = {'version': report_version(),}
    colnames = ['z_peak', 'y_peak', 'x_peak', 'peak_value', 'mask_values']
    coldata = [z_peaks, y_peaks, x_peaks, peak_values,mask_values]
    table = QTable(coldata, names=colnames, meta=meta)
    print(f'We found {table} the old table')
    exit()
    '''
    return peak_table

def chunked_iterable(iterable, size):
    """Yield successive chunks from iterable of given size."""
    it = iter(iterable)
    while True:
        chunk = list(islice(it, size))
        if not chunk:
            break
        yield chunk


def find_peaks(cfg,data, threshold, box_size=[3,3,3], mask=None,
               npeaks=np.inf, previous_deblend=None,
               num_processes=6,outdir='./',cube_header=None):
    
    shape = data.shape
    
    print_log(cfg,f'Data shape {shape}',case=['debug'])
    #maskd the data if a mask is present
    if mask is not None:
        data[mask < 0.5] = np.nan
    #box_size = [z,y,x]
    boxside = np.array(np.ceil(np.array(box_size) / 2.0),dtype=int)
    box = np.array([[boxside[0]-1, boxside[0]], [boxside[1]-1,boxside[1]], 
                    [boxside[2]-1,boxside[2]]],dtype=int)
   
    print_log(cfg,f'doing local maxima detection with {num_processes} processes',case=['verbose'])
    print_log(cfg,f'box size {box_size}, box {box}',case=['verbose'])
 
    start_time = datetime.datetime.now()
    peak_table_array = find_peaks_fast(data, threshold, box_size)
    finish_time = datetime.datetime.now()
    
    print_log(cfg,f'We found our peaks in {finish_time - start_time}',case=['verbose'])
    print_log(cfg,f'We now sort the peaks and filter them',case=['verbose'])
    # sorted is an intrinsic python funcion    
    #local_maxima_coords = sorted(local_maxima_coords, key=lambda x: x[3], reverse=True)
    peak_table = filter_peaks(cfg, peak_table_array, previous_deblend=previous_deblend
        ,npeaks=npeaks)
  
  
    if 'FIND_PEAKS' in cfg.logging.debug_functions or 'ALL' in cfg.logging.debug_functions:
        #write the peak table to debug directory
        peak_table.write(f'{cfg.logging.log_directory}/peak_table.ecsv', overwrite=True)

    #Make a new markers map based on the peaks
    markers3d = np.zeros(shape).astype(np.int8)

    for peak in peak_table:
        z,y,x = peak["z_peak"], peak["y_peak"], peak["x_peak"]
        # place markers in the vicinity of the flux peaks
        # We  use mask_smooth here because it provides the correct mask previous mask number 
        # for peak, i.e. peaks are grouped with the old source detection. 
       
        print_log(cfg,f"Placing markers for peak {peak['mask_values']} at {x},{y},{z} "
            ,case=['verbose'])

        markers3d[z-box[0,0]:z+box[0,1],y-box[1,0]:y+box[1,1],
            x-box[2,0]:x+box[2,1]] = int(peak["mask_values"])

        # you can add some markers manually
        # (manual markers here ...)
   
    print_log(cfg,f"Saving peak markers to {cfg.logging.log_directory}/peak3d_markers.fits",case=['debug'])
    if 'FIND_PEAKS' in cfg.logging.debug_functions or 'ALL' in cfg.logging.debug_functions:    
        write_fits_file(f"{cfg.logging.log_directory}/peak3d_markers_{cfg.internal.image_counter}.fits",markers3d,
            header=cube_header, overwrite=True)
        cfg.internal.image_counter += 1

    return peak_table,markers3d
      
def find_peaks_fast(data, threshold, box_size):
    # Apply maximum filter
    local_max = (data == maximum_filter(data, size=box_size, mode='constant'))
    # Apply threshold
    detected_peaks = local_max & (data > threshold)
   
    z, y, x = np.nonzero(detected_peaks)
    values = data[z, y, x]
    order = np.argsort(values)[::-1]
    peaks = np.column_stack([z[order], y[order], x[order], values[order]]).astype(float)

    #local_maxima_coords = list(zip(z, y, x, values))
    return peaks
    
      
def is_local_maxima(arr_box,threshold,box):
    if arr_box[box[0,0],box[1,0],box[2,0]] <= threshold:
        return False
    if np.isnan(arr_box[box[0,0],box[1,0],box[2,0]]):
        return False
    if np.isnan(arr_box).all():
        return False
    elif np.nanmax(arr_box) == arr_box[box[0,0],box[1,0],box[2,0]]:
        return True
    else:
        return False



'''
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
'''