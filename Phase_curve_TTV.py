#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# March 2025
# Exoplanet phase curve with transit timing variations

import numpy as np
import matplotlib.pyplot as plt
from TRAPPIST1_parameters import *
from Phase_curve_v1 import star_planet_separation, flux_star, flux_planet, luminosity_planet_dayside, phase_curve
from Flux_wavelength import flux_ratio_miri, planet_equilibirium_temperature, flux_planet_miri, integrate_flux_model_mJy
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


def phase_curve_simulation(t0, nb_days, nb_points=10000, planets='bcdefgh', redistribution=0, filter=None, model='sphinx', unit='ppm', Keplerian=False, total=True, plot=True, save_plot=False, save_txt=False):
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

    :param redistribution: the redistribution efficiency between the day side and night side (default: 0)
    :type redistribution: float

    :param filter: the filter to use (default: None). If None, the bolometric fluxes, expressed in ppm, are used relatively to the stellar flux with the planets considered as bare rocks.
    :type filter: str or None

    :param model: the model to use for the stellar flux (default: 'sphinx'). If 'sphinx', the flux is computed using the SPHINX model. If 'phoenix', the flux is computed using the PHOENIX model.
    :type model: str

    :param unit: the unit of the phase curve (default: 'ppm'). If 'ppm', the fluxes of the planets will be computed relatively to the stellar flux in ppm. If 'mJy', the planetary fluxes will be computed in absolute value in "mJy". Will be set automatically to 'mJy' if the model is 'phoenix'.
    :type unit: str

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

    if model == 'phoenix':
        unit = 'mJy'

    if unit == 'mJy':
        flux_star_mJy = integrate_flux_model_mJy(filter,model=model) # Compute the flux of the star in mJy
    
    if redistribution != 0 and redistribution != 1:
        raise ValueError("For now the code only supports a redistribution effciency equal to 0 or 1.")
    

    # For TRAPPIST-1 b

    if 'b' in planets:

        P_b_TTV, transit_peaks_b = np.loadtxt("Files_TTV/TTV_b.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)
        
        if Keplerian:
            P_b_TTV = P_b
        
        r_b = a_b
        
        b_b = eclipse_impact_parameter(a_b,i_b,e_b,R_star,omega_b)

        phase_b_TTV, t_b_TTV = phase_TTV(P_b_TTV,t0,t_end,transit_peaks_b,nb_points)
        eclipse_b = eclipse(P_b,a_b,R_star,R_b,i_b,np.arccos(phase_b_TTV)/(2*np.pi),e_b,omega_b,b_b)

        if filter==None:
            flux_star_b = flux_star(L_star,r_b)
            flux_b = flux_planet(flux_star_b)
            L_b = luminosity_planet_dayside(flux_b,R_b)
            
            phase_curve_b_TTV = phase_curve(L_star,L_b,R_star,R_b,phase_b_TTV,eclipse_b)
        
        else:
            T_b = planet_equilibirium_temperature(T_eff_star, R_star, a_b, redistribution=redistribution)

            if unit == 'ppm': # Compute the flux relatively to the star in ppm
                flux_ratio_b = flux_ratio_miri(filter, R_b, R_star, T_b) 
                if redistribution == 0:
                    phase_curve_b_TTV = flux_ratio_b * phase_b_TTV * (-1*eclipse_b+1)
                else:
                    phase_curve_b_TTV = flux_ratio_b * (-1*eclipse_b+1)
            
            elif unit == 'mJy': # Compute the flux in absolute value in mJy
                flux_ratio_b = flux_ratio_miri(filter, R_b, R_star, T_b)*1e-6
                if redistribution == 0:
                    phase_curve_b_TTV = (flux_ratio_b * phase_b_TTV * (-1*eclipse_b+1) +1) * flux_star_mJy
                else:
                    phase_curve_b_TTV = (flux_ratio_b * (-1*eclipse_b+1) +1) * flux_star_mJy

            else:
                raise ValueError("The unit must be 'ppm' or 'mJy'.")

            # if redistribution == 0:
            #     phase_curve_b_TTV = flux_ratio_b * phase_b_TTV * (-1*eclipse_b+1)
            # else:
            #     phase_curve_b_TTV = flux_ratio_b * (-1*eclipse_b+1)
                

        if save_txt:
            if filter==None:
                np.savetxt("Phase_curve_TTV_output/phase_curve_b_"+"_bolometric_"+str(t0)+".txt", np.column_stack((t_b_TTV, phase_curve_b_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), F_b/F_star (ppm)', comments='')
            else:
                if unit == "ppm":
                    header_flux = "F_b/F_star (ppm)"
                else:
                    header_flux = "F_star + F_b (mJy)"
                if redistribution==0:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_b_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_b_TTV, phase_curve_b_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')
                else:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_b_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_b_TTV, phase_curve_b_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')


    # For TRAPPIST-1 c

    if 'c' in planets:

        P_c_TTV, transit_peaks_c = np.loadtxt("Files_TTV/TTV_c.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)

        if Keplerian:
            P_c_TTV = P_c

        r_c = a_c
        
        b_c = eclipse_impact_parameter(a_c,i_c,e_c,R_star,omega_c)

        phase_c_TTV , t_c_TTV = phase_TTV(P_c_TTV,t0,t_end,transit_peaks_c,nb_points)
        eclipse_c = eclipse(P_c,a_c,R_star,R_c,i_c,np.arccos(phase_c_TTV)/(2*np.pi),e_c,omega_c,b_c)

        if filter==None:
            flux_star_c = flux_star(L_star,r_c)
            flux_c = flux_planet(flux_star_c)
            L_c = luminosity_planet_dayside(flux_c,R_c)
            
            phase_curve_c_TTV = phase_curve(L_star,L_c,R_star,R_c,phase_c_TTV,eclipse_c)

        else:
            T_c = planet_equilibirium_temperature(T_eff_star, R_star, a_c, redistribution=redistribution)
            
            if unit == 'ppm': # Compute the flux relatively to the star in ppm
                flux_ratio_c = flux_ratio_miri(filter, R_c, R_star, T_c)
                if redistribution == 0:
                    phase_curve_c_TTV = flux_ratio_c * phase_c_TTV * (-1*eclipse_c+1)
                else:
                    phase_curve_c_TTV = flux_ratio_c * (-1*eclipse_c+1)
            
            elif unit == 'mJy': # Compute the flux in absolute value in mJy
                flux_ratio_c = flux_ratio_miri(filter, R_c, R_star, T_c)*1e-6
                if redistribution == 0:
                    phase_curve_c_TTV = (flux_ratio_c * phase_c_TTV * (-1*eclipse_c+1) + 1) * flux_star_mJy 
                else:
                    phase_curve_c_TTV = (flux_ratio_c * (-1*eclipse_c+1) + 1) * flux_star_mJy

            else:
                raise ValueError("The unit must be 'ppm' or 'mJy'.")


        if save_txt:
            if filter==None:
                np.savetxt("Phase_curve_TTV_output/phase_curve_c_"+"_bolometric_"+str(t0)+".txt", np.column_stack((t_c_TTV, phase_curve_c_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), F_c/F_star (ppm)', comments='')
            else:
                if unit == "ppm":
                    header_flux = "F_c/F_star (ppm)"
                else:
                    header_flux = "F_star + F_c (mJy)"
                if redistribution==0:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_c_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_c_TTV, phase_curve_c_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')
                else:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_c_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_c_TTV, phase_curve_c_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')
      

    # For TRAPPIST-1 d
    
    if 'd' in planets:

        P_d_TTV, transit_peaks_d = np.loadtxt("Files_TTV/TTV_d.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)

        if Keplerian:
            P_d_TTV = P_d

        r_d = a_d
        
        b_d = eclipse_impact_parameter(a_d,i_d,e_d,R_star,omega_d)

        phase_d_TTV , t_d_TTV = phase_TTV(P_d_TTV,t0,t_end,transit_peaks_d,nb_points)
        eclipse_d = eclipse(P_d,a_d,R_star,R_d,i_d,np.arccos(phase_d_TTV)/(2*np.pi),e_d,omega_d,b_d)

        if filter==None:
            flux_star_d = flux_star(L_star,r_d)
            flux_d = flux_planet(flux_star_d)
            L_d = luminosity_planet_dayside(flux_d,R_d)
            
            phase_curve_d_TTV = phase_curve(L_star,L_d,R_star,R_d,phase_d_TTV,eclipse_d)

        else:
            T_d = planet_equilibirium_temperature(T_eff_star, R_star, a_d, redistribution=redistribution)

            if unit == 'ppm': # Compute the flux relatively to the star in ppm
                flux_ratio_d = flux_ratio_miri(filter, R_d, R_star, T_d)
                if redistribution == 0:
                    phase_curve_d_TTV = flux_ratio_d * phase_d_TTV * (-1*eclipse_d+1)
                else:
                    phase_curve_d_TTV = flux_ratio_d * (-1*eclipse_d+1)
            
            elif unit == 'mJy': # Compute the flux in absolute value in mJy
                flux_ratio_d = flux_ratio_miri(filter, R_d, R_star, T_d)*1e-6
                if redistribution == 0:
                    phase_curve_d_TTV = (flux_ratio_d * phase_d_TTV * (-1*eclipse_d+1) + 1) * flux_star_mJy
                else:
                    phase_curve_d_TTV = (flux_ratio_d * (-1*eclipse_d+1) + 1) * flux_star_mJy
            
            else:
                raise ValueError("The unit must be 'ppm' or 'mJy'.")


        if save_txt:
            if filter==None:
                np.savetxt("Phase_curve_TTV_output/phase_curve_d_"+"_bolometric_"+str(t0)+".txt", np.column_stack((t_d_TTV, phase_curve_d_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), F_d/F_star (ppm)', comments='')
            else:
                if unit == "ppm":
                    header_flux = "F_d/F_star (ppm)"
                else:
                    header_flux = "F_star + F_d (mJy)"
                if redistribution==0:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_d_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_d_TTV, phase_curve_d_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')
                else:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_d_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_d_TTV, phase_curve_d_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')
     

    # For TRAPPIST-1 e

    if 'e' in planets:

        P_e_TTV, transit_peaks_e = np.loadtxt("Files_TTV/TTV_e.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)

        if Keplerian:
            P_e_TTV = P_e

        r_e = a_e
        
        b_e = eclipse_impact_parameter(a_e,i_e,e_e,R_star,omega_e)

        phase_e_TTV , t_e_TTV = phase_TTV(P_e_TTV,t0,t_end,transit_peaks_e,nb_points)
        eclipse_e = eclipse(P_e,a_e,R_star,R_e,i_e,np.arccos(phase_e_TTV)/(2*np.pi),e_e,omega_e,b_e)

        if filter==None:
            flux_star_e = flux_star(L_star,r_e)
            flux_e = flux_planet(flux_star_e)
            L_e = luminosity_planet_dayside(flux_e,R_e)
            
            phase_curve_e_TTV = phase_curve(L_star,L_e,R_star,R_e,phase_e_TTV,eclipse_e)

        else:
            T_e = planet_equilibirium_temperature(T_eff_star, R_star, a_e, redistribution=redistribution)

            if unit == 'ppm': # Compute the flux relatively to the star in ppm
                flux_ratio_e = flux_ratio_miri(filter, R_e, R_star, T_e) 
                if redistribution == 0:
                    phase_curve_e_TTV = flux_ratio_e * phase_e_TTV * (-1*eclipse_e+1)
                else:
                    phase_curve_e_TTV = flux_ratio_e * (-1*eclipse_e+1)
            
            elif unit == 'mJy': # Compute the flux in absolute value in mJy
                flux_ratio_e = flux_ratio_miri(filter, R_e, R_star, T_e)*1e-6
                if redistribution == 0:
                    phase_curve_e_TTV = (flux_ratio_e * phase_e_TTV * (-1*eclipse_e+1) + 1) * flux_star_mJy 
                else:
                    phase_curve_e_TTV = (flux_ratio_e * (-1*eclipse_e+1) + 1) * flux_star_mJy

            else:
                raise ValueError("The unit must be 'ppm' or 'mJy'.")


        if save_txt:
            if filter==None:
                np.savetxt("Phase_curve_TTV_output/phase_curve_e_"+"_bolometric_"+str(t0)+".txt", np.column_stack((t_e_TTV, phase_curve_e_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), F_e/F_star (ppm)', comments='')
            else:
                if unit == "ppm":
                    header_flux = "F_e/F_star (ppm)"
                else:
                    header_flux = "F_star + F_e (mJy)"
                if redistribution==0:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_e_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_e_TTV, phase_curve_e_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')
                else:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_e_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_e_TTV, phase_curve_e_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')


    # For TRAPPIST-1 f

    if 'f' in planets:

        P_f_TTV, transit_peaks_f = np.loadtxt("Files_TTV/TTV_f.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)

        if Keplerian:
            P_f_TTV = P_f

        r_f = a_f
        
        b_f = eclipse_impact_parameter(a_f,i_f,e_f,R_star,omega_f)

        phase_f_TTV , t_f_TTV = phase_TTV(P_f_TTV,t0,t_end,transit_peaks_f,nb_points)
        eclipse_f = eclipse(P_f,a_f,R_star,R_f,i_f,np.arccos(phase_f_TTV)/(2*np.pi),e_f,omega_f,b_f)

        if filter==None:
            flux_star_f = flux_star(L_star,r_f)
            flux_f = flux_planet(flux_star_f)
            L_f = luminosity_planet_dayside(flux_f,R_f)
            
            phase_curve_f_TTV = phase_curve(L_star,L_f,R_star,R_f,phase_f_TTV,eclipse_f)

        else:
            T_f = planet_equilibirium_temperature(T_eff_star, R_star, a_f, redistribution=redistribution)

            if unit == 'ppm': # Compute the flux relatively to the star in ppm
                flux_ratio_f = flux_ratio_miri(filter, R_f, R_star, T_f)
                if redistribution == 0:
                    phase_curve_f_TTV = flux_ratio_f * phase_f_TTV * (-1*eclipse_f+1)
                else:
                    phase_curve_f_TTV = flux_ratio_f * (-1*eclipse_f+1)

            elif unit == 'mJy': # Compute the flux in absolute value in mJy
                flux_ratio_f = flux_ratio_miri(filter, R_f, R_star, T_f)*1e-6
                if redistribution == 0:
                    phase_curve_f_TTV = (flux_ratio_f * phase_f_TTV * (-1*eclipse_f+1) + 1) * flux_star_mJy
                else:
                    phase_curve_f_TTV = (flux_ratio_f * (-1*eclipse_f+1) + 1) * flux_star_mJy

            else:
                raise ValueError("The unit must be 'ppm' or 'mJy'.")

            

        if save_txt:
            if filter==None:
                np.savetxt("Phase_curve_TTV_output/phase_curve_f_"+"_bolometric_"+str(t0)+".txt", np.column_stack((t_f_TTV, phase_curve_f_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), F_f/F_star (ppm)', comments='')
            else:
                if unit == "ppm":
                    header_flux = "F_f/F_star (ppm)"
                else:
                    header_flux = "F_star + F_f (mJy)"
                if redistribution==0:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_f_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_f_TTV, phase_curve_f_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')
                else:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_f_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_f_TTV, phase_curve_f_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')


    # For TRAPPIST-1 g

    if 'g' in planets:

        P_g_TTV, transit_peaks_g = np.loadtxt("Files_TTV/TTV_g.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)

        if Keplerian:
            P_g_TTV = P_g

        r_g = a_g
        
        b_g = eclipse_impact_parameter(a_g,i_g,e_g,R_star,omega_g)

        phase_g_TTV , t_g_TTV = phase_TTV(P_g_TTV,t0,t_end,transit_peaks_g,nb_points)
        eclipse_g = eclipse(P_g,a_g,R_star,R_g,i_g,np.arccos(phase_g_TTV)/(2*np.pi),e_g,omega_g,b_g)

        if filter==None:
            flux_star_g = flux_star(L_star,r_g)
            flux_g = flux_planet(flux_star_g)
            L_g = luminosity_planet_dayside(flux_g,R_g)
            
            phase_curve_g_TTV = phase_curve(L_star,L_g,R_star,R_g,phase_g_TTV,eclipse_g)

        else:
            T_g = planet_equilibirium_temperature(T_eff_star, R_star, a_g, redistribution=redistribution)

            if unit == 'ppm': # Compute the flux relatively to the star in ppm
                flux_ratio_g = flux_ratio_miri(filter, R_g, R_star, T_g)
                if redistribution == 0:
                    phase_curve_g_TTV = flux_ratio_g * phase_g_TTV * (-1*eclipse_g+1)
                else:
                    phase_curve_g_TTV = flux_ratio_g * (-1*eclipse_g+1)

            elif unit == 'mJy': # Compute the flux in absolute value in mJy
                flux_ratio_g = flux_ratio_miri(filter, R_g, R_star, T_g)*1e-6
                if redistribution == 0:
                    phase_curve_g_TTV = (flux_ratio_g * phase_g_TTV * (-1*eclipse_g+1) + 1) * flux_star_mJy
                else:
                    phase_curve_g_TTV = (flux_ratio_g * (-1*eclipse_g+1) + 1) * flux_star_mJy

            else:
                raise ValueError("The unit must be 'ppm' or 'mJy'.")


        if save_txt:
                if filter==None:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_g_"+"_bolometric_"+str(t0)+".txt", np.column_stack((t_g_TTV, phase_curve_g_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), F_g/F_star (ppm)', comments='')
                else:
                    if unit == "ppm":
                        header_flux = "F_g/F_star (ppm)"
                    else:
                        header_flux = "F_star + F_g (mJy)"
                    if redistribution==0:
                        np.savetxt("Phase_curve_TTV_output/phase_curve_g_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_g_TTV, phase_curve_g_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')
                    else:
                        np.savetxt("Phase_curve_TTV_output/phase_curve_g_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_g_TTV, phase_curve_g_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')


    # For TRAPPIST-1 h

    if 'h' in planets:

        P_h_TTV, transit_peaks_h = np.loadtxt("Files_TTV/TTV_h.txt", delimiter=',',skiprows=1,usecols=(1,2),unpack=True)

        if Keplerian:
            P_h_TTV = P_h

        r_h = a_h
        
        b_h = eclipse_impact_parameter(a_h,i_h,e_h,R_star,omega_h)

        phase_h_TTV , t_h_TTV = phase_TTV(P_h_TTV,t0,t_end,transit_peaks_h,nb_points)
        eclipse_h = eclipse(P_h,a_h,R_star,R_h,i_h,np.arccos(phase_h_TTV)/(2*np.pi),e_h,omega_h,b_h)

        if filter==None:
            flux_star_h = flux_star(L_star,r_h)
            flux_h = flux_planet(flux_star_h)
            L_h = luminosity_planet_dayside(flux_h,R_h)
            
            phase_curve_h_TTV = phase_curve(L_star,L_h,R_star,R_h,phase_h_TTV,eclipse_h)

        else:
            T_h = planet_equilibirium_temperature(T_eff_star, R_star, a_h, redistribution=redistribution)

            if unit == 'ppm': # Compute the flux relatively to the star in ppm
                flux_ratio_h = flux_ratio_miri(filter, R_h, R_star, T_h)
                if redistribution == 0:
                    phase_curve_h_TTV = flux_ratio_h * phase_h_TTV * (-1*eclipse_h+1)
                else:
                    phase_curve_h_TTV = flux_ratio_h * (-1*eclipse_h+1)

            elif unit == 'mJy': # Compute the flux in absolute value in mJy
                flux_ratio_h = flux_ratio_miri(filter, R_h, R_star, T_h)*1e-6
                if redistribution == 0:
                    phase_curve_h_TTV = (flux_ratio_h * phase_h_TTV * (-1*eclipse_h+1) + 1) * flux_star_mJy
                else:
                    phase_curve_h_TTV = (flux_ratio_h * (-1*eclipse_h+1) + 1) * flux_star_mJy

            else:
                raise ValueError("The unit must be 'ppm' or 'mJy'.")


        if save_txt:
            if filter==None:
                np.savetxt("Phase_curve_TTV_output/phase_curve_h_"+"_bolometric_"+str(t0)+".txt", np.column_stack((t_h_TTV, phase_curve_h_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), F_h/L_star (ppm)', comments='')
            else:
                if unit == "ppm":
                    header_flux = "F_h/F_star (ppm)"
                else:
                    header_flux = "F_star + F_h (mJy)"
                if redistribution==0:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_h_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_h_TTV, phase_curve_h_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')
                else:
                    np.savetxt("Phase_curve_TTV_output/phase_curve_h_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t_h_TTV, phase_curve_h_TTV)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')
    

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
        
        if unit == 'mJy':
            if len(planets)>1:
             phase_curve_total -= flux_star_mJy*(len(planets)-1)

    if save_txt:
        if filter==None:
            np.savetxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_bolometric_"+str(t0)+".txt", np.column_stack((t, phase_curve_total)), delimiter=',', header='Time (BJD_TBD - 2450000), F_total/F_star (ppm)', comments='')
        else:
            if unit == "ppm":
                header_flux = "F_planets/F_star (ppm)"
            else:
                header_flux = "F_star + F_planets (mJy)"
            if redistribution==0:
                np.savetxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t, phase_curve_total)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')
            else:
                np.savetxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", np.column_stack((t, phase_curve_total)), delimiter=',', header='Time (BJD_TBD - 2450000), '+header_flux, comments='')
        

    # Plotting the phase curves

    if plot or save_plot:
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
            plt.plot(t,phase_curve_total,'--',color='grey',label="Total")
        plt.xlabel(r"Time ($BJD_{TBD} - 2450000$)")
        if unit == 'ppm':
            plt.ylabel(r"$F_{planet}/F_{star}$ (ppm)")
        elif unit == 'mJy':
            plt.ylabel(r"$F_{star} + F_{planet}$ (mJy)")
        else:
            raise ValueError("The unit must be 'ppm' or 'mJy'.")

        if filter==None:
            plt.title("Phase curves of planets of TRAPPIST-1 as bare rocks with bolometric fluxes")
        else:
            if redistribution==0:
                plt.title("Phase curves of planets of TRAPPIST-1 as bare rocks with MIRI "+filter+" filter using the "+model+" model")
            else:
                plt.title("Phase curves of planets of TRAPPIST-1 with atmospheres with MIRI "+filter+" filter using the "+model+" model")
        plt.legend()
        plt.grid()

        if save_plot:
            if filter==None:
                plt.savefig("Phase_curve_TTV_plots/phase_curve_TTV_"+planets+"_bolometric_"+str(t0)+".png", bbox_inches='tight')
            else:
                if redistribution==0:
                    plt.savefig("Phase_curve_TTV_plots/phase_curve_TTV_"+planets+"_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".png", bbox_inches='tight')
                else:
                    plt.savefig("Phase_curve_TTV_plots/phase_curve_TTV_"+planets+"_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".png", bbox_inches='tight')

        if plot:
            plt.show()
        else:
            plt.close()


def main():

    t0 = 9875 # starting time in BJD_TBD - 2450000
    nb_days = 30 # simulation duration in days

    # Call the function to simulate the phase curves
    phase_curve_simulation(t0, nb_days, planets='defgh', redistribution=1, filter="F1500W", Keplerian=True, total=True, plot=True, save_plot=False, save_txt=False)



if __name__ == "__main__":
    main()