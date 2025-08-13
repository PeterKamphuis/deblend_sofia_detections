from deblend_sofia_detections.support.errors import InputError,UnitError,\
    RegriddingError


from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from scipy.ndimage import map_coordinates

import astropy.units as u
import copy
import numpy as np
import os
import re



def calculate_projected_distance(coord1,coord2,no_PA=False): 
    '''Calculate the projected distance between two coordinates'''
    sk1 = SkyCoord(*coord1)
    sk2 = SkyCoord(*coord2)
    separation = sk1.separation(sk2).degree*u.deg
    if no_PA:
        return separation
   
    projected_PA=(np.degrees(np.arcsin(np.radians(((coord1[1].to(u.deg)-coord2[1].to(u.deg))\
                                    /separation.to(u.deg)).value)))+90.)*u.deg
    
    if coord1[0]-coord2[0] > 0.:
        projected_PA=360.*u.deg-projected_PA
  
    return separation,projected_PA.to(u.deg)


def convert_pix_columns_to_arcsec(cfg,table,file):
    #first open the cube
    hdr  = fits.getheader(file)
    pixsize= np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])]*u.deg/u.pix)
    for col in table.colnames:
        if table[col].unit == u.pix and not col[0] in ['x','y','z'] and \
            not col[-1] in ['x','y','z'] :
            table[col] = table[col]*pixsize
    return table

def convert_pixel_values_to_original(intable, cube_file_name,original_cube):
    original_hdr  = fits.getheader(original_cube)
    hdr = fits.getheader(cube_file_name)
    original_wcs = WCS(original_hdr)
    wcs_in = WCS(hdr)
    to_convert = ['','_min','_max']
    for conv_set in to_convert:
        for r in range(len(intable)):
            original_values = [intable[f'x{conv_set}'][r],
                               intable[f'y{conv_set}'][r],
                               intable[f'z{conv_set}'][r]]
            original_coord = wcs_in.wcs_pix2world(*original_values,1.)
            new_values = original_wcs.wcs_world2pix(*original_coord,1.)
            intable[f'x{conv_set}'][r] = new_values[0]*u.pix
            intable[f'y{conv_set}'][r] = new_values[1]*u.pix
            intable[f'z{conv_set}'][r] = new_values[2]*u.pix
    return intable

def convertRADEC(RAin,DECin,invert=False, colon=False, verbose=False):
    if verbose:
        print(f'''CONVERTRADEC: Starting conversion from the following input.
{'':8s}RA = {RAin}
{'':8s}DEC = {DECin}
''')
    RA = copy.deepcopy(RAin)
    DEC = copy.deepcopy(DECin)
    if not invert:
        try:
            _ = (e for e in RA)
        except TypeError:
            RA= [RA]
            DEC =[DEC]
        for i in range(len(RA)):
            xpos=RA
            ypos=DEC
            xposh=int(np.floor((xpos[i]/360.)*24.))
            xposm=int(np.floor((((xpos[i]/360.)*24.)-xposh)*60.))
            xposs=(((((xpos[i]/360.)*24.)-xposh)*60.)-xposm)*60
            yposh=int(np.floor(np.absolute(ypos[i]*1.)))
            yposm=int(np.floor((((np.absolute(ypos[i]*1.))-yposh)*60.)))
            yposs=(((((np.absolute(ypos[i]*1.))-yposh)*60.)-yposm)*60)
            sign=ypos[i]/np.absolute(ypos[i])
            if colon:
                RA[i]="{}:{}:{:2.2f}".format(xposh,xposm,xposs)
                DEC[i]="{}:{}:{:2.2f}".format(yposh,yposm,yposs)
            else:
                RA[i]="{}h{}m{:2.2f}".format(xposh,xposm,xposs)
                DEC[i]="{}d{}m{:2.2f}".format(yposh,yposm,yposs)
            if sign < 0.: DEC[i]='-'+DEC[i]
        if len(RA) == 1:
            RA = str(RA[0])
            DEC = str(DEC[0])
    else:
        if isinstance(RA,str):
            RA=[RA]
            DEC=[DEC]

        xpos=RA
        ypos=DEC

        for i in range(len(RA)):
            # first we split the numbers out
            tmp = re.split(r"[a-z,:]+",xpos[i])
            RA[i]=(float(tmp[0])+((float(tmp[1])+(float(tmp[2])/60.))/60.))*15.
            tmp = re.split(r"[a-z,:'\"]+",ypos[i])
            if float(tmp[0]) != 0.:
                DEC[i]=float(np.absolute(float(tmp[0]))+((float(tmp[1])+\
                    (float(tmp[2])/60.))/60.))*float(tmp[0])/np.absolute(float(tmp[0]))
            else:
                DEC[i] = float(np.absolute(float(tmp[0])) + ((float(tmp[1])\
                             + (float(tmp[2]) / 60.)) / 60.))
                if tmp[0][0] == '-':
                    DEC[i] = float(DEC[i])*-1.
        if len(RA) == 1:
            RA= float(RA[0])
            DEC = float(DEC[0])
        else:
            RA =np.array(RA,dtype=float)
            DEC = np.array(DEC,dtype=float)
    return RA,DEC

convertRADEC.__doc__ =f'''
 NAME:
    convertRADEC

 PURPOSE:
    convert the RA and DEC in degre to a string with the hour angle

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    RAin = RA to be converted
    DECin = DEC to be converted

 OPTIONAL INPUTS:


    invert=False
    if true input is hour angle string to be converted to degree

    colon=False
    hour angle separotor is : instead of hms

 OUTPUTS:
    converted RA, DEC as string list (hour angles) or numpy float array (degree)

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_nan_for_dtype(dtype):
    """Return appropriate NaN value for given dtype"""
    dtype = np.dtype(dtype)
    
    if dtype.kind == 'f':  # floating point
        return np.nan
    elif dtype.kind in ['i', 'u']:  # integer types
        return np.iinfo(dtype).min  # or use a sentinel value like -999
    elif dtype.kind == 'b':  # boolean
        return False  # or None
    elif dtype.kind in ['U', 'S']:  # string types
        return 'NaN'  # empty string or 'NaN'
    elif dtype.kind == 'O':  # object
        return None
    elif dtype.kind == 'M':  # datetime
        return np.datetime64('NaT')
    elif dtype.kind == 'm':  # timedelta
        return np.timedelta64('NaT')
    else:
        return None

def get_source_cat_name(line,input_columns,column_locations):
    '''Read out the source name from the catalogue input line'''
    if input_columns.index('name')-1 < 0.:
        start = 0
    else:
        start = column_locations[input_columns.index('name')-1]
    end = column_locations[input_columns.index('name')]
    name = line[start:end].strip()
    name = name.strip('"')
    name = '_'.join(name.split())
    return name

def get_start_end_locations(i,column_locations):
    if i == 0:
        start = 0
    else:
        start = column_locations[i-1]
    end = column_locations[i]
    return start,end

def is_real_unit(unit):
    # None is what we consider a valid unit for strings and such
    #If you want to check for a quanity use isQuantity
    if unit is None:
        return True
    try:
        u.Unit(unit)  # Try to create a unit from the string
        return True
    except ValueError:
        return False
    
def isiterable(variable):
    '''Check whether variable is iterable'''
    #First check it is not a string as those are iterable
    if isinstance(variable,str):
        return False
    try:
        iter(variable)
    except TypeError:
        return False

    return True
isiterable.__doc__ =f'''
 NAME:
    isiterable

 PURPOSE:
    Check whether variable is iterable

 CATEGORY:
    support_functions

 INPUTS:
    variable = variable to check

 OPTIONAL INPUTS:

 OUTPUTS:
    True if iterable False if not

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def isquantity(value):
    verdict= True
    if not isinstance(value,u.quantity.Quantity):
        if isiterable(value):
        #if it is an iterable we make sure it it is an numpy array
            if not isinstance(value,np.ndarray) and not value is None:
                value = quantity_array(value)
            else:
                verdict = False
        else:
            verdict = False
       
    return verdict

def match_size(matcharray, inarray,max=False):
    """
    Match the size of inarray to matcharray by regridding.
    Parameters:
    matcharray (np.ndarray): The array whose shape we want to match.
    inarray (np.ndarray): The array we want to regrid.
    Returns:
    np.ndarray: The regridded array with the same shape as matcharray.
    """

    if matcharray.shape == inarray.shape:
        return inarray
    else:
        if len(matcharray.shape) == 3 and len(inarray.shape) == 2:
            New_Shape = matcharray.shape[1:3]
        else:
            New_Shape = matcharray.shape
       
        return regrid_array(inarray, New_Shape,max= max)



def read_unit_part(input_string,transform_in):
    '''
    Read a part of the input string and return the unit, transformation and power.
    '''

    proc_transform = None
    proc_power = None
        # If we have a backslash we need to split the string
    tmp = input_string.split(transform_in)
    proc_unit = tmp[0].strip()
    #if we still have a trans form we pass
    if len(proc_unit.split('/')) > 1 or len(proc_unit.split('*')) > 1:
        #print(f'WARNING: The unit {proc_unit} contains a transformation {transform_in} that is not supported. Skipping this part.')
        proc_unit = None
        remainder = transform_in
        pass
    else:
        # If we have a power we need to split the string
        remainder = input_string.removeprefix(f'{proc_unit}').strip()
        if len(remainder) > 0:
            remainder = remainder.removeprefix(transform_in)
       
        if '^' in proc_unit:
            proc_unit,proc_power = proc_unit.split('^')
            proc_power = float(proc_power)
        else:
            proc_power = 1
        proc_transform = transform_in
        #remainder = tmp[1] if len(tmp) > 1 else ''
    
    return proc_unit,proc_transform,proc_power,remainder


def regrid_array(oldarray, Out_Shape,max = False):
    oldshape = np.array(oldarray.shape)
    newshape = np.array(Out_Shape, dtype=float)
    ratios = oldshape/newshape
        # calculate new dims
    nslices = [ slice(0,j) for j in list(newshape) ]
    #make a list with new coord
    new_coordinates = np.mgrid[nslices]
    #scale the new coordinates
    for i in range(len(ratios)):
        new_coordinates[i] *= ratios[i]
    #create our regridded array
    if max:
        order = 0
    else:
        order = 1
    newarray = map_coordinates(oldarray, new_coordinates,order=order)
    if any([x != y for x,y in zip(newarray.shape,newshape)]):
        raise RegriddingError(f'''Something went wrong when regridding.
This newarray {newarray.shape} and requested {newshape}''')
    return newarray
regrid_array.__doc__ =f'''
 NAME:
    regridder
 PURPOSE:
    Regrid an array into a new shape through the ndimage module
 CATEGORY:
    fits_functions

 INPUTS:
    oldarray = the larger array
    newshape = the new shape that is requested

 OPTIONAL INPUTS:

 OUTPUTS:
    newarray = regridded array

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    scipy.ndimage.map_coordinates, np.array, np.mgrid

 NOTE:
'''

def translate_string_to_unit(input,invert=False):
    '''
    translation_dict = {'ARCSEC': u.arcsec,
                        'ARCMIN': u.arcmin,
                        'DEGREE': u.degree,
                        'DEG': u.degree,
                        'DMS': 'dms',
                        'HMS': 'hms',
                        'MPC': u.Mpc,
                        'KPC': u.kpc,
                        'PC': u.pc,
                        'KM/S': u.km/u.s,
                        'M/S': u.m/u.s,
                        'M_SOLAR': u.Msun, 
                        'M_SUN': u.Msun,
                        'L_SOLAR': u.Lsun,
                        'L_SOLAR/PC^2': u.Lsun/u.pc**2,
                        'M_SOLAR/PC^2': u.Msun/u.pc**2,
                        'L_SOLAR/PC^3': u.Lsun/u.pc**3,
                        'M_SOLAR/PC^3': u.Msun/u.pc**3, 
                        'M_SOLAR/YR': u.Msun/u.yr,
                        'MAG/ARCSEC^2': u.mag/u.arcsec**2,
                        'PIX': u.pix,
                        'JY': u.Jy,
                        'JY/BEAM': u.Jy/u.beam,
                        'JY/BEAM*KM/S': u.Jy/u.beam*u.km/u.s,
                        'JY/BEAM*M/S': u.Jy/u.beam*u.m/u.s,
                        'JY*KM/S': u.Jy*u.km/u.s,
                        'JY*M/S': u.Jy*u.m/u.s,
                        'YR': u.yr,
                        '-':u.dimensionless_unscaled,
                        '':u.dimensionless_unscaled,
                        'SomethingIsWrong': None,
                        'UNKOWN': None,}
    '''
    translation_dict = {'ARCSEC': u.arcsec,
                        'ARCMIN': u.arcmin,
                        'DEGREE': u.degree,
                        'DEG': u.degree,
                        'DMS': 'dms',
                        'HMS': 'hms',
                        'MPC': u.Mpc,
                        'KPC': u.kpc,
                        'PC': u.pc,
                        'KM': u.km,
                        'M': u.m,
                        'S': u.s,
                        'M_SOLAR': u.Msun, 
                        'M_SUN': u.Msun,
                        'L_SOLAR': u.Lsun,
                        'L_SUN': u.Lsun,
                        'MAG': u.mag,
                        'PIX': u.pix,
                        'JY': u.Jy,
                        'MJY': u.mJy, #This has degeneracy with MJY
                        'BEAM': u.beam,
                        'YR': u.yr,
                        'LOG': u.dex,
                        'LOG10': u.dex,
                        '-':u.dimensionless_unscaled,
                        '':u.dimensionless_unscaled,
                        ' ':u.dimensionless_unscaled,
                        'SomethingIsWrong': None,
                        'STR': str,
                        'BOOL': bool,
                        'INT': int,
                        'UNKNOWN': u.dimensionless_unscaled,}
    output =False
    # First we need to seperate the main units and transformations
    if invert:
        output = unit_to_string(input,translation_dict)   
    else:
        output = string_to_unit(input,translation_dict)
      
             
    if output is False:
        raise UnitError(f'The unit {input} is not recognized for a valid translation.')
    else:
        return output
    
def unit_to_string(input,translation_dict):
    ''' Convert an astropy unit to a string.
    The input should be an astropy unit, e.g. u.Jy/u.beam*u.km/u.s**2.
    The function will return the corresponding string representation of the unit.
    If the input is not a valid unit, it will raise a UnitError
    '''
    
    if f'{input}' == '':
        return ''
    if input is None:
        return ''
    try:
        input.decompose()   
    except u.UnitTypeError:
        raise InputError(f'The input {input} is not a valid astropy unit.')
    add_log = False
   
    if isinstance(input, u.DexUnit):
        add_log = True
        input = input.physical_unit  
   
    units = input.bases
    powers = input.powers
    
    output = list(translation_dict.keys())[list(translation_dict.values()).index(units[0])]
    if powers[0] != 1:
        if powers[0] < 0:
            output = f'1/{output}^{abs(powers[0])}'
        else:
            output = f'{output}^{powers[0]}'
    if len(units) > 1:
        for unit,power in zip(units[1:],powers[1:]):
            if unit in list(translation_dict.values()):
                unit_str = list(translation_dict.keys())[list(translation_dict.values()).index(unit)]
            else:
                raise UnitError(f'The unit {unit} in {input} is not recognized for a valid translation.')
            if power < 0:
                if abs (power) > 1:
                    output += f'/{unit_str}^{abs(power)}'
                else:
                    output += f'/{unit_str}'
            elif power > 1:
                output += f'*{unit_str}^{power}'
            else:
                output += f'*{unit_str}'
    if add_log: 
        output = f'LOG({output})'
    return output
   

def string_to_unit(input,translation_dict):
    ''' Convert a string to an astropy unit.
    The input string should be in the form of a unit, e.g. 'JY/BEAM*KM/S**2'.
    The function will parse the string and return the corresponding astropy unit.
    If the input is not a valid unit, it will raise a UnitError

    For now this does not handle composite log units like 'LOG10(JY/BEAM*KM/S**2)'.
    '''
    unit_in = copy.deepcopy(input)
    
    if 'log' in unit_in.lower() or 'lg' in unit_in.lower():
        tmp = re.split(r'\(|\)', unit_in) 
        if len(tmp) != 3:
            raise UnitError(f'The unit {unit_in} is not a valid log unit. It should be in the form of LOG10(JY/BEAM*KM/S**2).')
        for i in range(len(tmp)):
            if tmp[i].strip() == 'log' or tmp[i].strip() == 'lg':
                conv_unit = process_composite_string(tmp[i+1],translation_dict) 
            output = u.dex(conv_unit)
    else:
        output = process_composite_string(unit_in,translation_dict)        # If we do not have a log unit we can just process the string    
      
    return output

def process_composite_string(input_string,translation_dict):
    input_proc = input_string.strip().upper() 
    #First make sure that powers are ^ not **
    input_proc = re.sub(r'\*\*',r'^',input_proc)
    transforms = []
    powers = []
    units = []
    if len(input_proc) == 0:
        units = ['']
        transforms = ['']
        powers = [1]
    else:
        while (len(input_proc) > 0):
            for tr in [r'/', r'*']:
                # another check as the transform loop is in the while
                if len(input_proc) > 0:
                    proc_unit,proc_transform,proc_power,remainder =\
                        read_unit_part(input_proc,tr)        
                    if proc_unit is None:
                        pass
                    else:
                        units.append(proc_unit)
                        transforms.append(proc_transform)
                        powers.append(proc_power)
                        input_proc = remainder.strip()
    # convert the units to the astropy units    
    convert_units =[]   
   
    for unit_in in units:       
        if unit_in in list(translation_dict.keys()):
            convert_units.append(translation_dict[unit_in])
        else:
            raise UnitError(f'The unit {unit_in} in {input} is not recognized for a valid translation.')
    #Combine the units in the correct manner
    output = convert_units[0]
    if len(convert_units) > 1:
        for unit_in,comb,power in zip(convert_units[1:],transforms,powers[1:]):
            if power > 1:
                unit_in = unit_in**power
            if comb == '*':
                output *= unit_in
            elif comb == '/':
                output /= unit_in
    return output

def quantity_array(list,unit):
    #Because astropy is coded by nincompoops Units to not convert into numpy arrays well.
    #It seems impossible to convert a list of Quantities into a quantity  with a list or np array
    #This means we have to pull some ticks when using numpy functions because they don't accept lists of Quantities
    # Major design flaw in astropy unit and one think these nincompoops could incorporate a function like this 
    #Convert a list of quantities into quantity with a numpy array
    return np.array([x.to(unit).value for x in list],dtype=float)*unit 
