
import os
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d

def get_extrema_row(ABBA):
    
    max_row = -np.inf ; min_row = np.inf
    avgs = np.zeros(ABBA.size)
    for i, row in enumerate(ABBA):
        avgs[i] = np.average(row) 
        if avgs[i] == 0: continue
        max_row = max(max_row, avgs[i])
        min_row = min(min_row, avgs[i])

    max_i = np.where(avgs==max_row)[0][0]
    min_i = np.where(avgs==min_row)[0][0]
    
    return max_i, min_i

def pixels_to_wavelength(file_header):
 
    naxis1 = file_header['naxis1']
    crpix1 = file_header['crpix1']
    crval1 = file_header['crval1']
    cdelt1 = file_header['cdelt1']
    
    wavelengths = crval1 + (np.arange(1, naxis1 + 1) - crpix1) * cdelt1
    return wavelengths
    
def combine_ABBA(ABBA, max_i, min_i):
    sp1 = np.sum(ABBA[max_i-5:max_i+5, :], axis=0)
    sp2 = -np.sum(ABBA[min_i-5:min_i+5, :], axis=0)

    return sp1 + sp2

def smooth(y, radius=10):
    smoothed_data = []
    for i, _ in enumerate(y):
        if i <= radius//2 or i >= len(y)-radius//2: continue
        smoothed_data += [np.median(y[i-radius//2:i+radius//2])]
        
    
    return np.array([np.nan]*(radius//2 + 1) + smoothed_data + [np.nan]*(radius//2))

def get_spectrum(ABBA_fname, grism, SN_position=None):
    
    ABBA_file = fits.open(ABBA_fname)
    ABBA = ABBA_file[1].data
    ABBA_header = ABBA_file[1].header
    ABBA_file.close()
    
    if SN_position is None:
        max_i, min_i = get_extrema_row(ABBA)
    else:
        min_i, max_i = SN_position
    
    wavelengths = pixels_to_wavelength(ABBA_header)
    spectrum = combine_ABBA(ABBA, max_i, min_i)
    exptime = ABBA_header['exptime']
    spectrum = spectrum / (2 * exptime)
    
    return spectrum, wavelengths, ABBA_header

def get_atmospheric_spectrum(wavelengths):

    tell_path = os.path.join(os.environ['EMIR_PIPE'], 'telluric.txt')
    telluric_flux = [] ; telluric_wave = []
    with open(tell_path, 'r') as f:
        for line in f:
            spl = line.split()
            telluric_flux += [float(spl[1])]
            telluric_wave += [float(spl[0])*10]

    atmos_interpolation = interp1d(telluric_wave, telluric_flux, kind="linear")
    return atmos_interpolation(wavelengths)



