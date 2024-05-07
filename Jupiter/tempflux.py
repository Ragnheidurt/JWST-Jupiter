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

path = 'C:/Users/dansf/OneDrive/Documents/KTH/Thesis/Code/JWST-Jupiter/Jupiter/Data/'
files = ['200-200.txt', '200-600.txt', '200-1000.txt', '600-200.txt', '600-600.txt', '600-1000.txt', '1000-200.txt', '1000-600.txt', '1000-1000.txt']
flux = np.empty((len(files)), dtype=object)
wave = np.empty((len(files)), dtype=object)

nrs1_jupiter = get_pkg_data_filename('Sure/MAST_2024-03-04T13_01_54.243Z/MAST_2024-03-04T13_01_54.243Z/JWST/jw01373003001_03105_00001/jw01373003001_03105_00001_nrs1_s3d.fits')
nrs2_jupiter = get_pkg_data_filename('Sure/MAST_2024-03-04T13_01_54.243Z/MAST_2024-03-04T13_01_54.243Z/JWST/jw01373003001_03105_00001/jw01373003001_03105_00001_nrs2_s3d.fits')

# Extract the data from the files
nrs1_data = fits.getdata(nrs1_jupiter,ext=1)
nrs2_data = fits.getdata(nrs2_jupiter,ext=1)
nrs1_error = fits.getdata(nrs1_jupiter,ext=2)
nrs2_error = fits.getdata(nrs2_jupiter,ext=2)


# Extract header for both data sets
nrs1_header = fits.getheader(nrs1_jupiter,ext=1)
nrs2_header = fits.getheader(nrs2_jupiter,ext=1)
error_header = fits.getheader(nrs1_jupiter,ext=2)

# Get the start and end wavelengths for both sensors
nrs1_wavestart = nrs1_header['WAVSTART']
nrs1_wavend = nrs1_header['WAVEND']
nrs2_wavestart = nrs2_header['WAVSTART']
nrs2_wavend = nrs2_header['WAVEND']

# Calculate the interval and size of wavelength data
nrs1_dist = nrs1_wavend-nrs1_wavestart
nrs1_interval = nrs1_dist/nrs1_data.shape[0]
nrs2_dist = nrs2_wavend-nrs2_wavestart
nrs2_interval = nrs2_dist/nrs2_data.shape[0]

# Create array with wavelength
nrs1_wave = np.zeros(nrs1_data.shape[0])
nrs2_wave = np.zeros(nrs2_data.shape[0])


for i in range(nrs1_data.shape[0]):
    if i != 0: 
        nrs1_wave[i] = nrs1_wave[i-1] + nrs1_interval
    else:
        nrs1_wave[i] = nrs1_wavestart

for i in range(nrs2_data.shape[0]):
    if i != 0:
        nrs2_wave[i] = nrs2_wave[i-1] + nrs2_interval
    else:
        nrs2_wave[i] = nrs2_wavestart

nrs1_avg, nrs1_error, index1 = brightest_spot(nrs1_data, 0, nrs1_error, 1500)
nrs2_avg, nrs2_error, index2 = brightest_spot(nrs2_data,0,nrs2_error,2500)
nrs1_avg = nrs1_data[:,26,30]
nrs2_avg = nrs2_data[:,26,30]

for k in range(nrs1_avg.shape[0]):
    if nrs1_avg[k] < -100 and k!=0:
        nrs1_avg[k] = nrs1_avg[k-1]

for k in range(nrs2_avg.shape[0]):
    if nrs2_avg[k] < -100 and k!=0:
        nrs2_avg[k] = nrs2_avg[k-1]


for i in range(len(files)):
    filepath = path + files[i]
    with open(filepath,'r') as file:
        lines = file.readlines()

    wavenumber = []
    absorption_intensity1 = []
    vacuum_wavelength = []
    column_density = []
    stuff = []

    for line in lines:
        data = line.split()

        wavenumber.append(float(data[0]))
        absorption_intensity1.append(float(data[1]))
        vacuum_wavelength.append(float(data[2]))
        column_density.append(float(data[3]))
        stuff.append(data[4])

    wavenumber = np.array(wavenumber)
    absorption_intensity1 = np.array(absorption_intensity1)
    vacuum_wavelength = np.array(vacuum_wavelength)
    column_density = np.array(column_density)
    stuff = np.array(stuff)
    vacuum_wavelength = vacuum_wavelength*1e-10
    low_limit = 4.5e-06
    up_limit = 5e-06
    model_wave_short,model_flux_short = model_shortening(vacuum_wavelength, low_limit, up_limit, column_density)
    scaled_model = scale_model(column_density,nrs1_avg, nrs2_avg, vacuum_wavelength,nrs2_wave)

    #flux[i] = model_flux_short
    #wave[i] = model_wave_short
    flux[i] = scaled_model
    wave[i] = vacuum_wavelength


plt.figure(0)
plt.plot(wave[0],flux[0],alpha=0.9, linestyle='--')
plt.plot(wave[1],flux[1],alpha=0.8, linestyle='--')
plt.plot(wave[2],flux[2],alpha=0.7, linestyle='--')
plt.plot(wave[3],flux[3],alpha=0.6, linestyle='--')
plt.plot(wave[4],flux[4],alpha=0.5, linestyle='--')
plt.plot(wave[5],flux[5],alpha=0.4, linestyle='--')
plt.plot(wave[6],flux[6],alpha=0.3, linestyle='--')
plt.plot(wave[7],flux[7],alpha=0.2, linestyle='--')
plt.plot(wave[8],flux[8],alpha=0.1, linestyle='--')

fig, axs = plt.subplots(3, 3)

# Plot data on each subplot
count = 0
for i in range(3):
    for j in range(3):
        index = i * 3 + j  # Calculate index for accessing data
        axs[i, j].plot(wave[count],flux[count], linestyle='--')  # Plot data
        axs[i,j].plot(nrs1_wave,nrs1_avg, linestyle='--', alpha=0.5)
        axs[i,j].plot(nrs2_wave,nrs2_avg, linestyle='--', alpha=0.5)
        #axs[i, j].set_title(f"Plot {index + 1}")  # Set subplot title
        axs[i,j].set_title(f"Rot-Vib {files[index]}")
        count = count + 1

# Adjust layout to prevent overlap
plt.tight_layout()

plt.figure(2)
plt.plot(nrs1_wave,nrs1_avg)
plt.plot(nrs2_wave,nrs2_avg)
plt.plot(wave[7],flux[7],linestyle='--', alpha=0.5, color='red')



plt.show()
