#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# March 2025
# Exoplanet phase curve with transit timing variations

import numpy as np
import matplotlib.pyplot as plt
from TRAPPIST1_parameters import *
from Phase_curve_v1 import star_planet_separation, flux_star, flux_planet, luminosity_planet_dayside, phase_curve
from Transits import eclipse, eclipse_impact_parameter
from Orbital_motion import compute_true_anomaly

def phase_TTV(P_TTV,t0,t_end,transit_peaks,nb_points):
    """
    Computes the phase of the planet taking the TTVs into account for a circular orbit and the adjusted time array

    :param P_TTV: the modified orbital periods of the planet due to the TTVs (in days)
    :type P_TTV: numpy.ndarray

    :param t0: the initial time (in days)
    :type t0: float

    :param t_end: the final time (in days)
    :type t_end: float

    :param transit_peaks: the peaks of the transits (in days)
    :type transit_peaks: numpy.ndarray

    :param nb_points: the number of points for the phase curve
    :type nb_points: int

    :return: phases_TTV, t
    :rtype: numpy.ndarray, numpy.ndarray
    """

    phases_TTV = np.zeros(nb_points)
    i = 0
    while transit_peaks[i] < t0:
        i += 1
    t_first_peak = transit_peaks[i+1]
    while t_end > transit_peaks[i]:
        i += 1
    t_last_peak = transit_peaks[i]
    t = np.linspace(t_first_peak, t_last_peak, nb_points)
    
    t_nearest_peak = np.zeros(nb_points)
    j = 1
    for i in range(nb_points):
        if t[i] < transit_peaks[j]:
            t_nearest_peak[i] = transit_peaks[j-1]
        else:
            t_nearest_peak[i] = transit_peaks[j]
            j += 1
        phases_TTV[i] = np.sin((t[i]+t_nearest_peak[i])/P_TTV[j-1]*2*np.pi)/2 + 0.5
        

    # for i in range(nb_points):
    #     for j in range(len(transit_peaks)-1):
    #         if t[i] >= transit_peaks[j] and t[i] < transit_peaks[j+1]:
    #             phases_TTV[i] = (t[i]-transit_peaks[j])/P_TTV[j]
    #     #phases_TTV[i] = (t[i]-t0)/P_TTV[i]

    return phases_TTV, t


def main():
    
    t0 = 7260 # starting time in BJD_TBD - 2450000
    nb_days = 20 # simulation duration in days
    t_end = t0 + nb_days

    nb_points = 10000

    t = np.linspace(t0, t_end, nb_points)
    

    # For TRAPPIST-1 b

    P_b_TTV, transit_peaks_b = np.loadtxt("Files_TTV/TTV_b.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)
    #print(P_b_TTV)
    
    # nu_b = compute_true_anomaly(t, P_b, t0, omega_b)
    # r_b = star_planet_separation(a_b,e_b,nu_b)
    r_b = a_b # We assume a circular orbit as there are too many errors for an elliptical one
    flux_star_b = flux_star(L_star,r_b)
    flux_b = flux_planet(flux_star_b)
    L_b = luminosity_planet_dayside(flux_b,R_b)
    b_b = eclipse_impact_parameter(a_b,i_b,e_b,R_star,omega_b)
    
    phase_b_TTV, t_b_TTV = phase_TTV(P_b_TTV,t0,t_end,transit_peaks_b,nb_points)
    print(phase_b_TTV)
    eclipse_b = eclipse(P_b,a_b,R_star,R_b,i_b,np.arccos(phase_b_TTV)/(2*np.pi),e_b,omega_b,b_b)
    
    phase_curve_b_TTV = phase_curve(L_star,L_b,R_star,R_b,phase_b_TTV,eclipse_b)

    plt.figure()
    plt.plot(t_b_TTV,phase_curve_b_TTV,label="b")
    plt.xlabel("Time ($BJD_{TBD} - 2450000$)")
    plt.ylabel("$L_{planet}/L_{star}$ (ppm)")
    plt.title("Phase curves of planets of TRAPPIST-1 as black bodies with TTVs")
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    main()