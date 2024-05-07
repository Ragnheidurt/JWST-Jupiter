from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import get_pkg_data_filename
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.signal import savgol_filter
from scipy import signal
import sys
sys.path.append('C:/Users/dansf/OneDrive/Documents/KTH/Thesis/Code/JWST_extrasolar_aurora/Functions')
from average_errors import average_errors
from short_interval import short_interval
from brightest_spot import brightest_spot
from scale_model import scale_model
from model_data import model_shortening

#nrs1 = get_pkg_data_filename('Data_BD\MAST_2024-04-30T08_43_33.520Z\MAST_2024-04-30T08_43_33.520Z\JWST\jw02124001001_03102_00001\jw02124001001_03102_00001_nrs1_x1d.fits')
#nrs2 = get_pkg_data_filename('Data_BD\MAST_2024-04-30T08_43_33.520Z\MAST_2024-04-30T08_43_33.520Z\JWST\jw02124001001_03102_00001\jw02124001001_03102_00001_nrs2_x1d.fits')
#nrs1 = get_pkg_data_filename('Data_BD\MAST_2024-04-30T08_55_19.229Z\MAST_2024-04-30T08_55_19.229Z\JWST\jw02124001001_03102_00002\jw02124001001_03102_00002_nrs1_x1d.fits')
#nrs2 = get_pkg_data_filename('Data_BD\MAST_2024-04-30T08_55_19.229Z\MAST_2024-04-30T08_55_19.229Z\JWST\jw02124001001_03102_00002\jw02124001001_03102_00002_nrs2_x1d.fits')
#nrs1 = get_pkg_data_filename('Data_BD\MAST_2024-04-30T08_55_19.229Z\MAST_2024-04-30T08_55_19.229Z\JWST\jw02124001001_03102_00003\jw02124001001_03102_00003_nrs1_x1d.fits')
#nrs2 = get_pkg_data_filename('Data_BD\MAST_2024-04-30T08_55_19.229Z\MAST_2024-04-30T08_55_19.229Z\JWST\jw02124001001_03102_00003\jw02124001001_03102_00003_nrs2_x1d.fits')
nrs1 = get_pkg_data_filename('Data_BD\MAST_2024-04-30T09_08_48.992Z\MAST_2024-04-30T09_08_48.992Z\JWST\jw02124051001_03102_00001\jw02124051001_03102_00001_nrs1_x1d.fits')
nrs2 = get_pkg_data_filename('Data_BD\MAST_2024-04-30T09_08_48.992Z\MAST_2024-04-30T09_08_48.992Z\JWST\jw02124051001_03102_00001\jw02124051001_03102_00001_nrs2_x1d.fits')

# Get the data
nrs1_data = fits.getdata(nrs1,ext=1)
nrs2_data = fits.getdata(nrs2,ext=1)

# Get the header of the data
nrs1_header = fits.getheader(nrs1,ext=1)
nrs2_header = fits.getheader(nrs2,ext=1)
print(nrs1_header)

#print(nrs1_data.dtype)

# Extract wavelength and flux from the data
nrs1_wavelength = nrs1_data['WAVELENGTH']
nrs1_flux = nrs1_data['FLUX']
nrs1_sb = nrs1_data['SURF_BRIGHT']
nrs1_flux_error = nrs1_data['FLUX_ERROR']
nrs1_flux_poisson = nrs1_data['FLUX_VAR_POISSON']

nrs2_wavelength = nrs2_data['WAVELENGTH']
nrs2_flux = nrs2_data['FLUX']
nrs2_sb = nrs2_data['SURF_BRIGHT']
nrs2_flux_error = nrs2_data['FLUX_ERROR']

# Focus on a specific interval
'''
low_limit = 3.992
up_limit = 3.996
'''
low_limit = 3.25
up_limit = 3.45

nrs1_lengd = 0
for i in range(nrs1_wavelength.shape[0]):
    if(nrs1_wavelength[i]>=low_limit and nrs1_wavelength[i]<= up_limit):
        nrs1_lengd = nrs1_lengd +1


nrs2_lengd = 0
for i in range(nrs2_wavelength.shape[0]):
    if(nrs2_wavelength[i]>=low_limit and nrs2_wavelength[i]<= up_limit):
        nrs2_lengd = nrs2_lengd +1

nrs1_close_wave = np.zeros(nrs1_lengd,dtype='float32')
nrs1_close_value = np.zeros(nrs1_lengd,dtype='float32')
count = 0
for i in range(nrs1_wavelength.shape[0]):
    if nrs1_wavelength[i] >= low_limit and nrs1_wavelength[i] <= up_limit:
        nrs1_close_wave[count] = nrs1_wavelength[i]
        nrs1_close_value[count] = nrs1_flux[i]
        count = count+1


nrs2_close_wave = np.zeros(nrs2_lengd,dtype='float32')
nrs2_close_value = np.zeros(nrs2_lengd,dtype='float32')
count = 0
for i in range(nrs2_wavelength.shape[0]):
    if nrs2_wavelength[i] >= low_limit and nrs2_wavelength[i] <= up_limit:
        nrs2_close_wave[count] = nrs2_wavelength[i]
        nrs2_close_value[count] = nrs2_flux[i]
        count = count+1


plt.figure(1)
plt.plot(nrs1_wavelength,nrs1_flux)
plt.plot(nrs2_wavelength,nrs2_flux)
#plt.fill_between(nrs1_wavelength,nrs1_flux-nrs1_flux_error, nrs1_flux+nrs1_flux_error,alpha = 0.5, label='error')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Brightness')

plt.show()