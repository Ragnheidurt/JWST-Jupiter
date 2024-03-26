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

    flux[i] = model_flux_short
    wave[i] = model_wave_short
    #flux[i] = column_density
    #wave[i] = vacuum_wavelength


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
        axs[i, j].plot(wave[count],flux[count])  # Plot data
        axs[i, j].set_title(f"Plot {index + 1}")  # Set subplot title
        count = count + 1

# Adjust layout to prevent overlap
plt.tight_layout()

plt.show()
