
import os
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d

from matplotlib import pyplot as plt

EXTRACTION_HALF_HEIGHT = 5
BACKGROUND_MARGIN = 5

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
    try:
        cdelt1 = file_header['cdelt1']
    except KeyError:
        cdelt1 = file_header['cd1_1']

    wavelengths = crval1 + (np.arange(1, naxis1 + 1) - crpix1) * cdelt1
    return wavelengths
    
def combine_ABBA(ABBA, max_i, min_i, half_height=EXTRACTION_HALF_HEIGHT):
    sp1 = np.sum(ABBA[max_i-half_height:max_i+half_height, :], axis=0)
    sp2 = -np.sum(ABBA[min_i-half_height:min_i+half_height, :], axis=0)

    return sp1 + sp2


def get_extraction_rows(ABBA, SN_position=None):
    if SN_position is None:
        return get_extrema_row(ABBA)
    if type(SN_position) == tuple:
        min_i, max_i = SN_position
        return max_i, min_i

    new_im = np.empty_like(ABBA)
    new_im[:] = np.nan
    new_im[SN_position-60:SN_position+60,:] = ABBA[SN_position-60:SN_position+60,:]
    return get_extrema_row(new_im)


def estimate_spectrum_error(ABBA, max_i, min_i,
                            half_height=EXTRACTION_HALF_HEIGHT,
                            background_margin=BACKGROUND_MARGIN):
    nrows, ncols = ABBA.shape
    background_mask = np.ones(nrows, dtype=bool)

    for center in (max_i, min_i):
        lo = max(0, center - half_height - background_margin)
        hi = min(nrows, center + half_height + background_margin)
        background_mask[lo:hi] = False

    background = ABBA[background_mask, :]
    if background.size == 0:
        sigma_pix = np.full(ncols, np.nan)
    else:
        median = np.nanmedian(background, axis=0)
        mad = np.nanmedian(np.abs(background - median), axis=0)
        sigma_pix = 1.4826 * mad

        fallback = np.nanstd(background, axis=0, ddof=1)
        invalid = ~np.isfinite(sigma_pix) | (sigma_pix == 0)
        sigma_pix[invalid] = fallback[invalid]

    n_extract_rows = 4 * half_height
    return np.sqrt(n_extract_rows) * sigma_pix

def smooth(y, radius=10):
    smoothed_data = []
    for i, _ in enumerate(y):
        if i <= radius//2 or i >= len(y)-radius//2: continue
        smoothed_data += [np.median(y[i-radius//2:i+radius//2])]
        
    
    return np.array([np.nan]*(radius//2 + 1) + smoothed_data + [np.nan]*(radius//2))

def fits_to_data(fname):
    file = fits.open(fname)
    try:
        ydata = file[1].data
        header = file[1].header
        file.close()
    except IndexError:
        ydata = file[0].data
        header = file[0].header
        file.close()

    return pixels_to_wavelength(header), ydata

def get_spectrum(ABBA_fname, grism, SN_position=None, debug=False, return_error=False):
    
    ABBA_file = fits.open(ABBA_fname)
    try:
    	ABBA = ABBA_file[1].data
    	ABBA_header = ABBA_file[1].header
    	ABBA_file.close()
    except IndexError:
        ABBA = ABBA_file[0].data
        ABBA_header = ABBA_file[0].header
        ABBA_file.close()

    max_i, min_i = get_extraction_rows(ABBA, SN_position=SN_position)
    
    wavelengths = pixels_to_wavelength(ABBA_header)

    if debug:
        fig, axs = plt.subplots(2, figsize=(20,5))

        r = 20
        axs[0].imshow(ABBA)
        axs[0].set_ylim(max_i-r, max_i+r)
        axs[1].imshow(ABBA)
        axs[1].set_ylim(min_i-r, min_i+r)
        plt.suptitle("%s raw data"%grism)
        plt.show()

    spectrum = combine_ABBA(ABBA, max_i, min_i)
    spectrum_error = estimate_spectrum_error(ABBA, max_i, min_i)
    exptime = ABBA_header['exptime']
    spectrum = spectrum / (2 * exptime)
    spectrum_error = spectrum_error / (2 * exptime)

    if return_error:
        return spectrum, spectrum_error, wavelengths, ABBA_header
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

def get_header(o, idx):
	
    f = [f for f in os.listdir(o.abba_dir) if f[0] !='.'][0]
    fname = o.abba_dir + '/' + f
    with fits.open(fname) as f:
        h = f[idx].header

    return h

def get_hydrogen_lines(widths=None):
    paschen = [18750, 12820, 10940, 10050, 9546,]
    bracket = [40510,26250,21660,19440,18170,14580,]
    etc = [15550, 15730, 15900, 16100, 16400, 16800, 17350,
		   20000, 20100, 20500, 20700, 20950, 21050, 21150]
    all_lines = sorted(paschen + bracket + etc)
    all_lines[15] = 0 ; all_lines[16] = 0

    all_lines[2] -= 30
    all_lines[6] -= 10
    all_lines[7] -= 15
    all_lines[8] += 20
    all_lines[9] += 20
    all_lines[10] += 20
    all_lines[11] += 15

    widths = widths if widths is not None else [9, 15, 15, 10, 3, 8, 9,
            									7, 7, 7, 5, 5, 7, 1, 
            									1, 9, 4, 4, 3, 3, 3,
            									3, 7, 1, 3]


    return all_lines, widths


def get_std_magnitudes(std_obs):
    """Query Simbad for the 2MASS J, H, K magnitudes of the standard star.

    The star name is taken from the OBJECT FITS keyword of the standard
    observation (e.g. 'HIP16897_SN2024xdq' → 'HIP16897').

    Parameters
    ----------
    std_obs : Observation
        The standard-star Observation object.

    Returns
    -------
    magnitudes : dict
        {'J': float, 'H': float, 'K': float}
    """
    h = get_header(std_obs, 1)
    raw_name = h['OBJECT'].split('_')[0]          # 'HIP16897_SN2024xdq' → 'HIP16897'
    simbad_id = urllib.parse.quote(raw_name)       # URL-safe

    url = ('https://simbad.u-strasbg.fr/simbad/sim-id'
           '?output.format=votable'
           '&Ident=' + simbad_id +
           '&output.params=flux(J),flux(H),flux(K)')

    with urllib.request.urlopen(url, timeout=20) as resp:
        tree = ET.parse(resp)

    ns = {'vo': 'http://www.ivoa.net/xml/VOTable/v1.2'}
    root = tree.getroot()
    tr = root.find('.//vo:TABLEDATA/vo:TR', ns)
    if tr is None:
        raise ValueError("Simbad returned no results for '%s'" % raw_name)

    j_mag, h_mag, k_mag = [float(td.text) for td in tr.findall('vo:TD', ns)]
    magnitudes = {'J': j_mag, 'H': h_mag, 'K': k_mag}
    print("2MASS magnitudes for %s  →  J=%.3f  H=%.3f  K=%.3f"
          % (raw_name, j_mag, h_mag, k_mag))
    return magnitudes
