#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# March 2026
# Detect quick flux variations that could be interesting for follow-up observations

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from pathlib import Path
# racine du projet (Exoplanets_Phase_Curves/)
Code_files_DIR = Path(__file__).resolve().parents[1] / "Code_files"

from Code_files.Phase_curve_TTV import phase_TTV
from Code_files.Flux_wavelength import flux_ratio_miri, planet_equilibirium_temperature
from Code_files.TRAPPIST1_parameters import *

def detect_flux_variations(t0, nb_days, planets, filter, t_lapse, min_var, nb_points=10000, save_output = True, output_path = "Flux_variations", show_plot=False):
    """
    Detect quick flux variations in the total phase curve that could be interesting for follow-up observations.

    :param t0: Initial time (in BJD_TBD - 2450000)
    :type t0: float

    :param nb_days: Number of days to consider
    :type nb_days: int

    :param planets: Planets to consider (e.g., "defgh")
    :type planets: str

    :param filter: MIRI filter to use for the phase curve (F1500W or F1280W)
    :type filter: str

    :param t_lapse: Targeted lapse of time during which to detect variations (in hours)
    :type t_lapse: float

    :param min_var: Minimum flux variation threshold (in ppm)
    :type min_var: float

    :param nb_points: Number of points to use for the phase curve calculation (default: 10000)
    :type nb_points: int

    :param save_output: Whether to save the output in a .txt file and plot in a .png file (default: True)
    :type save_output: bool

    :param output_path: The path to the output file (default: "Flux_variations")
    :param output_path: str

    :param show_plot: Whether to show the plot (default: True)
    :param show_plot: bool

    :return: array containing detected flux variations with their corresponding times
    :rtype: numpy.ndarray
    """

    if filter not in ["F1500W", "F1280W"]:
        raise ValueError("Invalid filter. Choose either 'F1500W' or 'F1280W'.")
    
    if nb_days/nb_points > t_lapse/24:
        print(f"Warning: The time resolution ({nb_days/nb_points} days) is too low to detect variations within the specified time lapse ({t_lapse/24} days). Consider increasing nb_points or decreasing t_lapse.")
        # print(f"Warning: The number of points ({nb_points}) is too low to detect variations within the specified time lapse ({t_lapse} hours) over the total time span ({nb_days} days). Consider increasing nb_points or decreasing t_lapse.")
        raise ValueError("The number of points is too low to detect variations within the specified time lapse. Increase nb_points or decrease t_lapse.")
    
    t_end = t0 + nb_days
    
    t = np.linspace(t0, t_end, nb_points)

    # For TRAPPIST-1 b

    if 'b' in planets:

        P_b_TTV, transit_peaks_b = np.loadtxt(Code_files_DIR/"Files_TTV/TTV_b.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)
        r_b = a_b

        phase_b_TTV, t_b_TTV = phase_TTV(P_b_TTV,t0,t_end,transit_peaks_b,nb_points)

        T_b = planet_equilibirium_temperature(T_eff_star, R_star, a_b, redistribution=0)

        flux_ratio_b = flux_ratio_miri(filter, R_b, R_star, T_b) 

        phase_curve_b_TTV = flux_ratio_b * phase_b_TTV

    
    # For TRAPPIST-1 c

    if 'c' in planets:

        P_c_TTV, transit_peaks_c = np.loadtxt(Code_files_DIR/"Files_TTV/TTV_c.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)
        r_c = a_c

        phase_c_TTV, t_c_TTV = phase_TTV(P_c_TTV,t0,t_end,transit_peaks_c,nb_points)

        T_c = planet_equilibirium_temperature(T_eff_star, R_star, a_c, redistribution=0)

        flux_ratio_c = flux_ratio_miri(filter, R_c, R_star, T_c) 

        phase_curve_c_TTV = flux_ratio_c * phase_c_TTV

    
    # For TRAPPIST-1 d

    if 'd' in planets:

        P_d_TTV, transit_peaks_d = np.loadtxt(Code_files_DIR/"Files_TTV/TTV_d.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)
        r_d = a_d

        phase_d_TTV, t_d_TTV = phase_TTV(P_d_TTV,t0,t_end,transit_peaks_d,nb_points)

        T_d = planet_equilibirium_temperature(T_eff_star, R_star, a_d, redistribution=0)

        flux_ratio_d = flux_ratio_miri(filter, R_d, R_star, T_d) 

        phase_curve_d_TTV = flux_ratio_d * phase_d_TTV

    
    # For TRAPPIST-1 e

    if 'e' in planets:

        P_e_TTV, transit_peaks_e = np.loadtxt(Code_files_DIR/"Files_TTV/TTV_e.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)
        r_e = a_e

        phase_e_TTV, t_e_TTV = phase_TTV(P_e_TTV,t0,t_end,transit_peaks_e,nb_points)

        T_e = planet_equilibirium_temperature(T_eff_star, R_star, a_e, redistribution=0)

        flux_ratio_e = flux_ratio_miri(filter, R_e, R_star, T_e) 

        phase_curve_e_TTV = flux_ratio_e * phase_e_TTV

    
    # For TRAPPIST-1 f

    if 'f' in planets:

        P_f_TTV, transit_peaks_f = np.loadtxt(Code_files_DIR/"Files_TTV/TTV_f.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)
        r_f = a_f

        phase_f_TTV, t_f_TTV = phase_TTV(P_f_TTV,t0,t_end,transit_peaks_f,nb_points)

        T_f = planet_equilibirium_temperature(T_eff_star, R_star, a_f, redistribution=0)

        flux_ratio_f = flux_ratio_miri(filter, R_f, R_star, T_f) 

        phase_curve_f_TTV = flux_ratio_f * phase_f_TTV

    # For TRAPPIST-1 g

    if 'g' in planets:

        P_g_TTV, transit_peaks_g = np.loadtxt(Code_files_DIR/"Files_TTV/TTV_g.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)
        r_g = a_g

        phase_g_TTV, t_g_TTV = phase_TTV(P_g_TTV,t0,t_end,transit_peaks_g,nb_points)

        T_g = planet_equilibirium_temperature(T_eff_star, R_star, a_g, redistribution=0)

        flux_ratio_g = flux_ratio_miri(filter, R_g, R_star, T_g) 

        phase_curve_g_TTV = flux_ratio_g * phase_g_TTV

    
    # For TRAPPIST-1 h

    if 'h' in planets:

        P_h_TTV, transit_peaks_h = np.loadtxt(Code_files_DIR/"Files_TTV/TTV_h.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)
        r_h = a_h

        phase_h_TTV, t_h_TTV = phase_TTV(P_h_TTV,t0,t_end,transit_peaks_h,nb_points)

        T_h = planet_equilibirium_temperature(T_eff_star, R_star, a_h, redistribution=0)

        flux_ratio_h = flux_ratio_miri(filter, R_h, R_star, T_h) 

        phase_curve_h_TTV = flux_ratio_h * phase_h_TTV


    # Calculate the total phase curve by summing the contributions of the selected planets

    total_phase_curve = np.zeros(nb_points)

    if 'b' in planets:
        total_phase_curve += phase_curve_b_TTV
    if 'c' in planets:
        total_phase_curve += phase_curve_c_TTV
    if 'd' in planets:
        total_phase_curve += phase_curve_d_TTV
    if 'e' in planets:
        total_phase_curve += phase_curve_e_TTV
    if 'f' in planets:
        total_phase_curve += phase_curve_f_TTV
    if 'g' in planets:
        total_phase_curve += phase_curve_g_TTV
    if 'h' in planets:
        total_phase_curve += phase_curve_h_TTV



    # Looking for the extrema


    variations = []
    
    maxima, max_prop = find_peaks(total_phase_curve)
    # print(maxima)
    # print(f"Number of maxima: {len(maxima)}")
    # print(max_prop)
    minima, min_prop = find_peaks(-total_phase_curve)
    # print(minima)
    # print(f"Number of minima: {len(minima)}")

    if maxima[0]>minima[0]: # If the first extrmum is a minimum
        for i in range(len(maxima)):
            if total_phase_curve[maxima[i]] - total_phase_curve[minima[i]] > min_var: # Comparing with previous minimum
                variations.append((t[minima[i]],total_phase_curve[maxima[i]]-total_phase_curve[minima[i]],(t[maxima[i]]-t[minima[i]])*24))
            if i < len(minima)-1 and total_phase_curve[maxima[i]] - total_phase_curve[minima[i+1]] > min_var: # Comparing with folowing minimum
                variations.append((t[maxima[i]],total_phase_curve[maxima[i]] - total_phase_curve[minima[i+1]],(t[minima[i+1]]-t[maxima[i]])*24))
    else:
        for i in range(len(maxima)): # If the first extremum is a maximum
            if i != 0 and total_phase_curve[maxima[i]] - total_phase_curve[minima[i-1]] > min_var: # Comparing with previous minimum
                variations.append((t[minima[i-1]],total_phase_curve[maxima[i]]-total_phase_curve[minima[i]],(t[maxima[i]]-t[minima[i-1]])*24))
            # if i >= len(minima):
            #     break
            if i < len(minima)-1 and total_phase_curve[maxima[i]] - total_phase_curve[minima[i]] > min_var: # Comparing with following minimum
                variations.append((t[maxima[i]],total_phase_curve[maxima[i]]-total_phase_curve[minima[i]],(t[minima[i]]-t[maxima[i]])*24))

    if save_output:
        np.savetxt(f"{output_path}.txt", variations, header = "t0 (BJD_TBD - 2450000), Variation (ppm), Duration (hours)", fmt="%f, %f, %f")
        
    variations = np.array(variations).T


    plt.figure(figsize=(16,10))
    plt.plot(t, total_phase_curve, label="Total phase curve")
    plt.scatter(variations[0]+variations[2]/48,variations[1], marker="x", color='red', label="Quick variation")
    plt.title("Quick variations detected in total phase curve")
    plt.xlabel("Time")
    plt.ylabel("Flux (ppm)")
    plt.grid()
    plt.legend()
    if save_output:
        plt.savefig(f"{output_path}.png", bbox_inches="tight")
    if show_plot:
        plt.show()

    return variations


def main():
    # Example usage
    t0 = 11192.5000000  # Initial time in BJD_TBD - 2450000 (June 1st, 2026)
    nb_days = 500   # Number of days to consider
    planets = "defgh"  # Planets to consider
    filter = "F1500W"  # MIRI filter to use
    t_lapse = 10  # Lapse of time during which to detect variations (in hours)
    min_var = 400  # Minimum flux variation threshold (in ppm)
    nb_points = 100000  # Number of points to use for the phase curve calculation

    save_output = False

    output_path = Code_files_DIR/f"Flux_variations_output/Variations_{int(t0)}_{nb_days}d_{min_var}ppm"

    show_plot = False

    variations = detect_flux_variations(t0, nb_days, planets, filter, t_lapse, min_var, nb_points, save_output=save_output, output_path=output_path, show_plot=show_plot)
    # print(variations.shape)
    print("Detected flux variations:")
    for time, variation, duration in variations.T:
        print(f"Time: {time}, Variation: {variation} ppm, Duration: {duration} hours")

    print(f"Total number of detected variations: {variations.shape[1]}")

if __name__ == "__main__":
    main()