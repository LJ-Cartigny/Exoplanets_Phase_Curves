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
    Computes the phase of the planet taking into account the modification of the period due to TTVs starting from the nearest transit peak from t0

    :param P_TTV: the modified orbital periods of the planet due to the TTVs (in days)
    :type P_TTV: numpy.ndarray

    :param t0: the initial time (in BJD_TBD - 2450000)
    :type t0: float

    :param t_end: the final time (in BJD_TBD - 2450000)
    :type t_end: float

    :param transit_peaks: the peaks of the transits (in BJD_TBD - 2450000)
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

    if np.abs(t0 - transit_peaks[i]) < np.abs(t0 - transit_peaks[i-1]):
        t_first_transit = transit_peaks[i]
        if type(P_TTV) == np.ndarray:
            P = P_TTV[i]
        else:
            P = P_TTV
        # P = P_TTV[i]
    else:
        t_first_transit = transit_peaks[i-1]
        if type(P_TTV) == np.ndarray:
            P = P_TTV[i]
        else:
            P = P_TTV
        # P = P_TTV[i-1]
    
    t = np.linspace(t0, t_end, nb_points)

    for i in range(nb_points):
        phases_TTV[i] = np.sin((t[i]-t_first_transit)/P * 2*np.pi - np.pi/2)/2 + 0.5

    return phases_TTV, t


def phase_curve_simulation(t0, nb_days, nb_points=10000, planets='bcdefgh', Keplerian=False, total=True, plot=True, save_plot=False, save_txt=False):
    """
    Simulates the phase curves of the planets of TRAPPIST-1 for a given number of days starting from t0 taking into account the modified periods due to TTVs.
    We assume circular orbits as otherwise the code does not manage to solve the Kepler equation to compute the true anomaly due to the modified periods.

    :param t0: the initial time (in BJD_TBD - 2450000)
    :type t0: float

    :param nb_days: the number of days to simulate
    :type nb_days: int

    :param nb_points: the number of points for the phase curves (default: 10000)
    :type nb_points: int

    :param planets: the planets to simulate (default: 'bcdefgh')
    :type planets: str

    :param Keplerian: whether to use the Keplerian periods or not (default: False)
    :type Keplerian: bool

    :param total: whether to plot the total phase curve or not (default: True)
    :type total: bool

    :param plot: whether to plot the phase curves or not (default: True)
    :type plot: bool

    :param save_plot: whether to save the plot or not (default: False)
    :type save_plot: bool

    :param save_txt: whether to save the phase curves as txt files or not (default: False)
    :type save_txt: bool

    :return: None
    """
    
    t_end = t0 + nb_days

    t = np.linspace(t0, t_end, nb_points)
    

    # For TRAPPIST-1 b

    if 'b' in planets:

        P_b_TTV, transit_peaks_b = np.loadtxt("Files_TTV/TTV_b.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)
        
        if Keplerian:
            P_b_TTV = P_b
        
        r_b = a_b
        flux_star_b = flux_star(L_star,r_b)
        flux_b = flux_planet(flux_star_b)
        L_b = luminosity_planet_dayside(flux_b,R_b)
        b_b = eclipse_impact_parameter(a_b,i_b,e_b,R_star,omega_b)
        
        phase_b_TTV, t_b_TTV = phase_TTV(P_b_TTV,t0,t_end,transit_peaks_b,nb_points)
        eclipse_b = eclipse(P_b,a_b,R_star,R_b,i_b,np.arccos(phase_b_TTV)/(2*np.pi),e_b,omega_b,b_b)
        
        phase_curve_b_TTV = phase_curve(L_star,L_b,R_star,R_b,phase_b_TTV,eclipse_b)

        if save_txt:
            np.savetxt("Phase_curve_TTV_output/phase_curve_b_TTV_"+str(t0)+".txt", np.column_stack((t_b_TTV, phase_curve_b_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), L_b/L_star (ppm)', comments='')


    # For TRAPPIST-1 c

    if 'c' in planets:

        P_c_TTV, transit_peaks_c = np.loadtxt("Files_TTV/TTV_c.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)

        if Keplerian:
            P_c_TTV = P_c

        r_c = a_c
        flux_star_c = flux_star(L_star,r_c)
        flux_c = flux_planet(flux_star_c)
        L_c = luminosity_planet_dayside(flux_c,R_c)
        b_c = eclipse_impact_parameter(a_c,i_c,e_c,R_star,omega_c)

        phase_c_TTV , t_c_TTV = phase_TTV(P_c_TTV,t0,t_end,transit_peaks_c,nb_points)
        eclipse_c = eclipse(P_c,a_c,R_star,R_c,i_c,np.arccos(phase_c_TTV)/(2*np.pi),e_c,omega_c,b_c)

        phase_curve_c_TTV = phase_curve(L_star,L_c,R_star,R_c,phase_c_TTV,eclipse_c)

        if save_txt:
            np.savetxt("Phase_curve_TTV_output/phase_curve_c_TTV_"+str(t0)+".txt", np.column_stack((t_c_TTV, phase_curve_c_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), L_c/L_star (ppm)', comments='')


    # For TRAPPIST-1 d
    
    if 'd' in planets:

        P_d_TTV, transit_peaks_d = np.loadtxt("Files_TTV/TTV_d.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)

        if Keplerian:
            P_d_TTV = P_d

        r_d = a_d
        flux_star_d = flux_star(L_star,r_d)
        flux_d = flux_planet(flux_star_d)
        L_d = luminosity_planet_dayside(flux_d,R_d)
        b_d = eclipse_impact_parameter(a_d,i_d,e_d,R_star,omega_d)

        phase_d_TTV , t_d_TTV = phase_TTV(P_d_TTV,t0,t_end,transit_peaks_d,nb_points)
        eclipse_d = eclipse(P_d,a_d,R_star,R_d,i_d,np.arccos(phase_d_TTV)/(2*np.pi),e_d,omega_d,b_d)

        phase_curve_d_TTV = phase_curve(L_star,L_d,R_star,R_d,phase_d_TTV,eclipse_d)

        if save_txt:
            np.savetxt("Phase_curve_TTV_output/phase_curve_d_TTV_"+str(t0)+".txt", np.column_stack((t_d_TTV, phase_curve_d_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), L_d/L_star (ppm)', comments='')


    # For TRAPPIST-1 e

    if 'e' in planets:

        P_e_TTV, transit_peaks_e = np.loadtxt("Files_TTV/TTV_e.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)

        if Keplerian:
            P_e_TTV = P_e

        r_e = a_e
        flux_star_e = flux_star(L_star,r_e)
        flux_e = flux_planet(flux_star_e)
        L_e = luminosity_planet_dayside(flux_e,R_e)
        b_e = eclipse_impact_parameter(a_e,i_e,e_e,R_star,omega_e)

        phase_e_TTV , t_e_TTV = phase_TTV(P_e_TTV,t0,t_end,transit_peaks_e,nb_points)
        eclipse_e = eclipse(P_e,a_e,R_star,R_e,i_e,np.arccos(phase_e_TTV)/(2*np.pi),e_e,omega_e,b_e)

        phase_curve_e_TTV = phase_curve(L_star,L_e,R_star,R_e,phase_e_TTV,eclipse_e)

        if save_txt:
            np.savetxt("Phase_curve_TTV_output/phase_curve_e_TTV_"+str(t0)+".txt", np.column_stack((t_e_TTV, phase_curve_e_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), L_e/L_star (ppm)', comments='')


    # For TRAPPIST-1 f

    if 'f' in planets:

        P_f_TTV, transit_peaks_f = np.loadtxt("Files_TTV/TTV_f.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)

        if Keplerian:
            P_f_TTV = P_f

        r_f = a_f
        flux_star_f = flux_star(L_star,r_f)
        flux_f = flux_planet(flux_star_f)
        L_f = luminosity_planet_dayside(flux_f,R_f)
        b_f = eclipse_impact_parameter(a_f,i_f,e_f,R_star,omega_f)

        phase_f_TTV , t_f_TTV = phase_TTV(P_f_TTV,t0,t_end,transit_peaks_f,nb_points)
        eclipse_f = eclipse(P_f,a_f,R_star,R_f,i_f,np.arccos(phase_f_TTV)/(2*np.pi),e_f,omega_f,b_f)

        phase_curve_f_TTV = phase_curve(L_star,L_f,R_star,R_f,phase_f_TTV,eclipse_f)

        if save_txt:
            np.savetxt("Phase_curve_TTV_output/phase_curve_f_TTV_"+str(t0)+".txt", np.column_stack((t_f_TTV, phase_curve_f_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), L_f/L_star (ppm)', comments='')


    # For TRAPPIST-1 g

    if 'g' in planets:

        P_g_TTV, transit_peaks_g = np.loadtxt("Files_TTV/TTV_g.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)

        if Keplerian:
            P_g_TTV = P_g

        r_g = a_g
        flux_star_g = flux_star(L_star,r_g)
        flux_g = flux_planet(flux_star_g)
        L_g = luminosity_planet_dayside(flux_g,R_g)
        b_g = eclipse_impact_parameter(a_g,i_g,e_g,R_star,omega_g)

        phase_g_TTV , t_g_TTV = phase_TTV(P_g_TTV,t0,t_end,transit_peaks_g,nb_points)
        eclipse_g = eclipse(P_g,a_g,R_star,R_g,i_g,np.arccos(phase_g_TTV)/(2*np.pi),e_g,omega_g,b_g)

        phase_curve_g_TTV = phase_curve(L_star,L_g,R_star,R_g,phase_g_TTV,eclipse_g)

        if save_txt:
            np.savetxt("Phase_curve_TTV_output/phase_curve_g_TTV_"+str(t0)+".txt", np.column_stack((t_g_TTV, phase_curve_g_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), L_g/L_star (ppm)', comments='')


    # For TRAPPIST-1 h

    if 'h' in planets:

        P_h_TTV, transit_peaks_h = np.loadtxt("Files_TTV/TTV_h.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)

        if Keplerian:
            P_h_TTV = P_h

        r_h = a_h
        flux_star_h = flux_star(L_star,r_h)
        flux_h = flux_planet(flux_star_h)
        L_h = luminosity_planet_dayside(flux_h,R_h)
        b_h = eclipse_impact_parameter(a_h,i_h,e_h,R_star,omega_h)

        phase_h_TTV , t_h_TTV = phase_TTV(P_h_TTV,t0,t_end,transit_peaks_h,nb_points)
        eclipse_h = eclipse(P_h,a_h,R_star,R_h,i_h,np.arccos(phase_h_TTV)/(2*np.pi),e_h,omega_h,b_h)

        phase_curve_h_TTV = phase_curve(L_star,L_h,R_star,R_h,phase_h_TTV,eclipse_h)

        if save_txt:
            np.savetxt("Phase_curve_TTV_output/phase_curve_h_TTV_"+str(t0)+".txt", np.column_stack((t_h_TTV, phase_curve_h_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), L_h/L_star (ppm)', comments='')


    # Total signal

    if total:
        phase_curve_total = np.zeros(nb_points)
        if 'b' in planets:
            phase_curve_total += phase_curve_b_TTV
        if 'c' in planets:  
            phase_curve_total += phase_curve_c_TTV
        if 'd' in planets:
            phase_curve_total += phase_curve_d_TTV
        if 'e' in planets:
            phase_curve_total += phase_curve_e_TTV
        if 'f' in planets:
            phase_curve_total += phase_curve_f_TTV
        if 'g' in planets:
            phase_curve_total += phase_curve_g_TTV
        if 'h' in planets:
            phase_curve_total += phase_curve_h_TTV

    if save_txt:
        np.savetxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_TTV_"+str(t0)+".txt", np.column_stack((t, phase_curve_total)), delimiter=',', header='Time (BJD_TBD - 2450000), L_total/L_star (ppm)', comments='')


    # Plotting the phase curves

    if plot:
        plt.figure(figsize=(16,9))
        if 'b' in planets:
            plt.plot(t_b_TTV,phase_curve_b_TTV,label="b")
        if 'c' in planets:
            plt.plot(t_c_TTV,phase_curve_c_TTV,label="c")
        if 'd' in planets:
            plt.plot(t_d_TTV,phase_curve_d_TTV,label="d")
        if 'e' in planets:
            plt.plot(t_e_TTV,phase_curve_e_TTV,label="e")
        if 'f' in planets:
            plt.plot(t_f_TTV,phase_curve_f_TTV,label="f")
        if 'g' in planets:
            plt.plot(t_g_TTV,phase_curve_g_TTV,label="g")
        if 'h' in planets:
            plt.plot(t_h_TTV,phase_curve_h_TTV,label="h")
        if total:
            plt.plot(t,phase_curve_total,'--',label="Total")
        plt.xlabel(r"Time ($BJD_{TBD} - 2450000$)")
        plt.ylabel(r"$L_{planet}/L_{star}$ (ppm)")
        plt.title("Phase curves of planets of TRAPPIST-1 as black bodies with TTVs")
        plt.legend()
        plt.grid()

        if save_plot:
            plt.savefig("Phase_curve_TTV_plots/phase_curve_TTV_"+planets+"_"+str(t0)+".png")

        # Show the plot
        plt.show()


def main():

    t0 = 9785 # starting time in BJD_TBD - 2450000
    nb_days = 30 # simulation duration in days

    # Call the function to simulate the phase curves
    phase_curve_simulation(t0, nb_days)



if __name__ == "__main__":
    main()