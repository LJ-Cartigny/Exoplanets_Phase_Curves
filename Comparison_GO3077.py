#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# June 2025
# Comparison of the GO 3077 observations with the simulations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import h5py
from astropy.time import Time

mpl.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 20,
    'axes.titlesize': 20,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    # 'figure.figsize': (16, 9),
    # 'lines.linewidth': 2,
    # 'grid.alpha': 0.5,
    # 'grid.linestyle': '--',
})

def load_data(filename, key):
    """
    Load the data from the HDF5 file.

    :param filename: The path to the HDF5 file.
    :type filename: str

    :param key: The key to access the data in the HDF5 file.
    :type key: str

    :return: The data loaded from the file.
    :rtype: np.ndarray
    """
    
    with h5py.File(filename, 'r') as f:
        # print("Available keys in the file:", list(f.keys()))
        data = f[key][:]
    
    return data

def main():

    planets = "bcdefgh"
    filter = "F1500W"
    model = "sphinx"
    unit = "mJy"

    # Load the simulated data

    program_ID, visit, t_start, t_end, filter_obs, flux_obs, err_obs = np.loadtxt("JWST_Obs_times.txt", delimiter=',', skiprows=2, unpack=True,dtype=str)

    flux_obs = flux_obs.astype(float)
    err_obs = err_obs.astype(float)

    t_start = Time(t_start, format='isot', scale='tdb')
    t_end = Time(t_end, format='isot', scale='tdb')

    t_start.format = 'jd'
    t_end.format = 'jd'

    t_start = t_start.jd
    t_end = t_end.jd

    t_start -= 2450000
    t_end -= 2450000

    i = 0

    while program_ID[i] != 'GO_3077':
        i += 1

    t0_visit1 = t_start[i]
    t0_visit2 = t_start[i+1]


    visit1_time_sphinx, visit1_flux_sphinx = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_"+filter+"_"+model+"_"+unit+"_"+str(t0_visit1)+".txt",unpack=True, skiprows=1, delimiter=',')
    visit2_time_sphinx, visit2_flux_sphinx = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_"+filter+"_"+model+"_"+unit+"_"+str(t0_visit2)+".txt",unpack=True, skiprows=1, delimiter=',')

    visit1_time_sphinx -= 0.5
    visit2_time_sphinx -= 0.5

    offset_sphinx = -0.19

    model = 'phoenix'
    visit1_time_phoenix, visit1_flux_phoenix = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_"+filter+"_"+model+"_"+unit+"_"+str(t0_visit1)+".txt",unpack=True, skiprows=1, delimiter=',')
    visit2_time_phoenix, visit2_flux_phoenix = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_"+filter+"_"+model+"_"+unit+"_"+str(t0_visit2)+".txt",unpack=True, skiprows=1, delimiter=',')

    visit1_time_phoenix -= 0.5
    visit2_time_phoenix -= 0.5

    offset_phoenix = -0.225


    # Load the GO 3077 observations

    # visit= 'Visit1'

    # visit1_flux_obs = load_data("Data_GO_3077/S3_miri_photometry_T1bc_ap10_bg15_45_SpecData_"+visit+".h5", 'aplev')
    # visit1_flux_err_obs = load_data("Data_GO_3077/S3_miri_photometry_T1bc_ap10_bg15_45_SpecData_"+visit+".h5", 'aperr')

    # visit1_time_obs = load_data("Data_GO_3077/S3_miri_photometry_T1bc_ap10_bg15_45_SpecData_"+visit+".h5", 'time') - 50000

    # visit= 'Visit2'
    # visit2_flux_obs = load_data("Data_GO_3077/S3_miri_photometry_T1bc_ap10_bg15_45_SpecData_"+visit+".h5", 'aplev')
    # visit2_flux_err_obs = load_data("Data_GO_3077/S3_miri_photometry_T1bc_ap10_bg15_45_SpecData_"+visit+".h5", 'aperr')
    # visit2_time_obs = load_data("Data_GO_3077/S3_miri_photometry_T1bc_ap10_bg15_45_SpecData_"+visit+".h5", 'time') -50000

    visit= 'Visit12'

    visit12_flux_obs = load_data("Data_GO_3077/S3_miri_photometry_T1bc_ap12_bg20_45_SpecData_"+visit+".h5", 'aplev')
    visit12_flux_err_obs = load_data("Data_GO_3077/S3_miri_photometry_T1bc_ap12_bg20_45_SpecData_"+visit+".h5", 'aperr')
    visit12_time_obs = load_data("Data_GO_3077/S3_miri_photometry_T1bc_ap12_bg20_45_SpecData_"+visit+".h5", 'time') -50000

    offset = 0


    # Plot the data

    plt.figure(figsize=(16, 9))

    # plt.errorbar(visit1_time_obs, visit1_flux_obs, yerr=visit1_flux_err_obs, fmt='o', label='Visit 1 observation', color='blue', markersize=2, capsize=3)
    # plt.errorbar(visit2_time_obs, visit2_flux_obs, yerr=visit2_flux_err_obs, fmt='o', label='Visit 2 observation', color='orange', markersize=2, capsize=3)
    plt.errorbar(visit12_time_obs, visit12_flux_obs+offset, yerr=visit12_flux_err_obs, fmt='o', label='Visit 1 & 2 observations', color='green', markersize=2, capsize=3)

    plt.plot(visit1_time_sphinx, visit1_flux_sphinx+offset_sphinx, label='Visit 1 simulation (corrected SPHINX - 0.19 mJy)', color='blue', linestyle=':',linewidth=2,zorder=10)
    plt.plot(visit2_time_sphinx, visit2_flux_sphinx+offset_sphinx, label='Visit 2 simulation (corrected SPHINX - 0.19 mJy)', color='orange', linestyle=':',linewidth=2,zorder=10)
    plt.plot(visit1_time_phoenix, visit1_flux_phoenix+offset_phoenix, label='Visit 1 simulation (corrected PHOENIX - 0.225 mJy)', color='blue', linestyle='-',linewidth=2,zorder=10)
    plt.plot(visit2_time_phoenix, visit2_flux_phoenix+offset_phoenix, label='Visit 2 simulation (corrected PHOENIX - 0.225 mJy)', color='orange', linestyle='-',linewidth=2,zorder=10)

    plt.title('Phase curves comparison for the GO 3077 program')
    plt.xlabel(r'Time ($BMJD_{TBD}$)')
    plt.ylabel('Flux (mJy)')
    plt.grid()
    plt.legend()
    plt.xlim(np.min(visit12_time_obs), np.max(visit12_time_obs))
    plt.ylim(2.33,2.41)
    # plt.ylim(2.54,2.61)
    ax = plt.gca()
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

    # plt.savefig('GO_3077_curve_comparison_v2_zoom.png', bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()