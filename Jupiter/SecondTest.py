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


# Read in files from both sensors
path = "Sure/MAST_2024-03-04T13_01_54.243Z/MAST_2024-03-04T13_01_54.243Z/JWST/"
files = [["jw01373003001_0310b_00001/jw01373003001_0310b_00001_nrs1_s3d.fits", "jw01373003001_0310b_00001/jw01373003001_0310b_00001_nrs2_s3d.fits"],
         ["jw01373003001_0310b_00002/jw01373003001_0310b_00002_nrs1_s3d.fits", "jw01373003001_0310b_00002/jw01373003001_0310b_00002_nrs2_s3d.fits"],
         ["jw01373003001_03101_00002/jw01373003001_03101_00002_nrs1_s3d.fits", "jw01373003001_03101_00002/jw01373003001_03101_00002_nrs2_s3d.fits"],
         ["jw01373003001_03103_00001/jw01373003001_03103_00001_nrs1_s3d.fits", "jw01373003001_03103_00001/jw01373003001_03103_00001_nrs2_s3d.fits"],
         ["jw01373003001_03103_00002/jw01373003001_03103_00002_nrs1_s3d.fits", "jw01373003001_03103_00002/jw01373003001_03103_00002_nrs2_s3d.fits"],
         ["jw01373003001_03105_00001/jw01373003001_03105_00001_nrs1_s3d.fits", "jw01373003001_03105_00001/jw01373003001_03105_00001_nrs2_s3d.fits"],
         ["jw01373003001_03105_00002/jw01373003001_03105_00002_nrs1_s3d.fits", "jw01373003001_03105_00002/jw01373003001_03105_00002_nrs2_s3d.fits"],
         ["jw01373003001_03107_00001/jw01373003001_03107_00001_nrs1_s3d.fits", "jw01373003001_03107_00001/jw01373003001_03107_00001_nrs2_s3d.fits"],
         ["jw01373003001_03107_00002/jw01373003001_03107_00002_nrs1_s3d.fits", "jw01373003001_03107_00002/jw01373003001_03107_00002_nrs2_s3d.fits"],
         ["jw01373003001_03109_00001/jw01373003001_03109_00001_nrs1_s3d.fits", "jw01373003001_03109_00001/jw01373003001_03109_00001_nrs2_s3d.fits"],
         ["jw01373003001_03109_00002/jw01373003001_03109_00002_nrs1_s3d.fits", "jw01373003001_03109_00002/jw01373003001_03109_00002_nrs2_s3d.fits"]]

nrs = []
nrs_data = np.empty((len(files), 2), dtype=object)
nrs_header = np.empty((len(files),2), dtype=object)
nrs_wavstart = np.empty((len(files),2), dtype=float)
nrs_wavend = np.empty((len(files),2), dtype=float)
nrs_interval = np.empty((len(files),2), dtype=float)
nrs_error = np.empty((len(files),2),dtype=object)

for i in range(len(files)):
    row = []
    filename1 = path+files[i][0]
    filename2 = path+files[i][1]
    row.append(get_pkg_data_filename(filename1))
    row.append(get_pkg_data_filename(filename2))
    nrs.append(row)


for i in range(len(files)):
    nrs_data[i,0] = fits.getdata(nrs[i][0])
    nrs_data[i,1] = fits.getdata(nrs[i][1])
    nrs_header[i,0] = fits.getheader(nrs[i][0],ext=1)
    nrs_header[i,1] = fits.getheader(nrs[i][1],ext=1)
    nrs_error[i,0] = fits.getdata(nrs[i][0],ext = 2)
    nrs_error[i,1] = fits.getdata(nrs[i][1],ext = 2)
    nrs_wavstart[i,0] = nrs_header[i,0]['WAVSTART']
    nrs_wavstart[i,1] = nrs_header[i,1]['WAVSTART']
    nrs_wavend[i,0] = nrs_header[i,0]['WAVEND']
    nrs_wavend[i,1] = nrs_header[i,1]['WAVEND']
    dist1 = nrs_wavend[i,0]-nrs_wavstart[i,0]
    dist2 = nrs_wavend[i,1]-nrs_wavstart[i,1]
    nrs_interval[i,0] = dist1/nrs_data[i,0].shape[0]
    nrs_interval[i,1] = dist2/nrs_data[i,1].shape[0]

nrs_wave = np.empty((len(files),2), dtype=object)

for i in range(len(files)):
    wave = np.empty(nrs_data[i,0].shape[0])
    for j in range(nrs_data[i,0].shape[0]):
        if j != 0:
            wave[j] = wave[j-1]+nrs_interval[i,0]
        else:
            wave[j] = nrs_wavstart[i,0]
    nrs_wave[i,0] = wave
    wave = np.empty(nrs_data[i,1].shape[0])
    for j in range(nrs_data[i,1].shape[0]):
        if j != 0:
            wave[j] = wave[j-1]+nrs_interval[i,1]
        else:
            wave[j] = nrs_wavstart[i,1]
    nrs_wave[i,1] = wave
    #divisors = [i for i in range(1, nrs_data[i,0].shape[0] + 1) if nrs_data[i,0].shape[0] % i == 0]
    # Assuming nrs_data is a 2D NumPy array
    for k in range(nrs_data.shape[0]):
        divisors = [j for j in range(1, nrs_data[i, 0].shape[0] + 1) if nrs_data[i, 0].shape[0] % j == 0]
    # Do something with divisors


    error = nrs_error[i,0][:,26,30]
    nrs1_wave_average, nrs1_average, nrs1_error_average = average_errors(nrs_wave[i,0],nrs_data[i,0][:,26,30],error,divisors[4])
    plt.figure(i)
    plt.plot(nrs1_wave_average,nrs1_average)
    plt.fill_between(nrs1_wave_average,nrs1_average-nrs1_error_average, nrs1_average+nrs1_error_average,alpha = 0.5, color='green',label='error')

    #plt.plot(nrs_wave[i,0],nrs_data[i,0][:,26,30])
    #plt.plot(nrs_wave[i,0],nrs_error[i,0][:,26,30])
    #plt.plot(nrs_wave[i,1],nrs_data[i,1][:,26,30])

plt.show()
        
low_limit = 3.983e-06
up_limit = 3.997e-06

for i in range(len(files)):
    start_index1 = np.where(nrs_wave[i,0] >= low_limit)[0][0]
    end_index1 = np.where(nrs_wave[i,0] <= up_limit)[0][-1]
    start_index2 = np.where(nrs_wave[i,1] >= low_limit)[0][0]
    end_index2 = np.where(nrs_wave[i,1] <= up_limit)[0][-1]
    interval_wavelengths1 = nrs_wave[i,0][start_index1:end_index1 + 1]
    interval_intensities1 = nrs_data[i,0][start_index1:end_index1 + 1]
    interval_wavelengths2 = nrs_wave[i,1][start_index2:end_index2 + 1]
    interval_intensities2 = nrs_data[i,1][start_index2:end_index2 + 1]
    plt.figure(i)
    plt.plot(interval_wavelengths1,interval_intensities1[:,26,30])
    plt.plot(interval_wavelengths2,interval_intensities2[:,26,30])
    #plt.axvline(x=3953e-9, color='r', linestyle='--', label='Vertical Line')
    plt.axvline(x=3985.5e-9, color='r', linestyle='--', label='Vertical Line')

plt.show()


for i in range(len(files)):
    max1 = np.argmax(nrs_data[i,0], axis=0)
    max2 = np.argmax(nrs_data[i,1], axis=0)
    print(max1)
    fig, axs = plt.subplots(1,2)
    plt.figure(i)
    axs[0].imshow(nrs_data[i,0][max1[0],:,:])
    axs[1].imshow(nrs_data[i,1][max2[0],:,:])

plt.show()

fig, axs = plt.subplots(2, 2)


plt.figure(1)
axs[0,0].imshow(nrs_data[0,0][1500,:,:])
axs[0,1].imshow(nrs_data[0,1][2000,:,:])
axs[1,0].imshow(nrs_data[1,0][1000,:,:])
axs[1,1].imshow(nrs_data[1,1][2000,:,:])


plt.show()