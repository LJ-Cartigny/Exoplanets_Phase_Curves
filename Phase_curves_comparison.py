#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# March 2025
# Phase curves comparisons

import numpy as np
import matplotlib.pyplot as plt
from TRAPPIST1_parameters import *



def main():
    simu_file_b = "Phase_curve_v1_output/phase_curve_b_bis.txt"
    # simu_phase_b, simu_flux_b = np.loadtxt(simu_file_b, unpack=True)
    # simu_phase_b=simu_phase_b/P_b+omega_b/(2*np.pi)-0.75

    simu_time_b, simu_flux_b = np.loadtxt(simu_file_b, unpack=True)

    ref_file_b = "References_b_and_c/lc_michael_bestfit_model_b.dat"
    ref_phase_b, ref_flux_b = np.loadtxt(ref_file_b, usecols=(0,1), unpack=True)
    ref_time_b = (ref_phase_b-omega_b/(2*np.pi)+1.75)*P_b
    ref_flux_b = (ref_flux_b-1)*1e6


    simu_file_c = "Phase_curve_v1_output/phase_curve_c_bis.txt"
    # simu_phase_c, simu_flux_c = np.loadtxt(simu_file_c, unpack=True)
    # simu_phase_c=simu_phase_c/P_c+omega_c/(2*np.pi)-0.75

    simu_time_c, simu_flux_c = np.loadtxt(simu_file_c, unpack=True)

    ref_file_c = "References_b_and_c/lc_michael_bestfit_model_c.dat"
    ref_phase_c, ref_flux_c = np.loadtxt(ref_file_c, usecols=(0,1), unpack=True)
    ref_time_c = (ref_phase_c-omega_c/(2*np.pi)+1.75)*P_c
    ref_flux_c = (ref_flux_c-1)*1e6



    plt.figure()
    plt.subplot(121)
    # plt.plot(simu_phase_b-1, simu_flux_b, label="Simulation")
    # plt.plot(ref_phase_b, ref_flux_b, label="Reference")
    plt.plot(simu_time_b, simu_flux_b, label="Simulation")
    plt.plot(ref_time_b, ref_flux_b, label="Reference")
    # plt.xlabel("Phase")
    plt.xlabel("Time (days)")
    plt.ylabel("$F_{planet}/F_{star}$ (ppm)")
    # plt.xlim(np.min(ref_phase_b), np.max(ref_phase_b))
    plt.xlim(np.min(ref_time_b), np.max(ref_time_b))
    plt.title("Comparison of phase curves of TRAPPIST-1 b")
    plt.ylim(0,)
    plt.legend()
    plt.grid()

    plt.subplot(122)
    # plt.plot(simu_phase_c-1, simu_flux_c, label="Simulation")
    # plt.plot(ref_phase_c, ref_flux_c, label="Reference")
    plt.plot(simu_time_c, simu_flux_c, label="Simulation")
    plt.plot(ref_time_c, ref_flux_c, label="Reference")
    # plt.xlabel("Phase")
    plt.xlabel("Time (days)")
    plt.ylabel("$F_{planet}/F_{star}$ (ppm)")
    # plt.xlim(np.min(ref_phase_c), np.max(ref_phase_c))
    plt.xlim(np.min(ref_time_c), np.max(ref_time_c))
    plt.title("Comparison of phase curves of TRAPPIST-1 c")
    plt.ylim(0,)
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()