#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# January 2025
# Phase curves

import numpy as np
import matplotlib.pyplot as plt
from Orbital_motion import compute_true_anomaly

def phase_angle(omega, nu, i):
    """
    Determines the phase angle of a planet from its orbital parameters (in rad).

    :param omega: the argument of pericentre (in rad)
    :type omega: float

    :param nu: the true anomaly (in rad)
    :type nu: float

    :param i: the inclination (in rad)
    :type i: float

    :return: alpha
    :rtype: float
    """

    alpha = np.arccos(np.sin(omega+nu)*np.sin(i))
    return alpha

def phase_function(alpha):
    """
    Determines the phase function of a Lambert sphere.

    :param alpha: the phase angle (in rad)
    :type alpha: float

    :return: g
    :rtype: float
    """

    g = (np.sin(alpha)+(np.pi-alpha)*np.cos(alpha))/np.pi
    return g

def star_planet_separation(a,e,nu):
    """
    Determines the distance between a planet and its star using its orbital parameters.

    :param a: the semimajor axis (in m)
    :type a: float

    :param e: the eccentricity
    :type e: float

    :param nu: the true anomaly (in rad)
    :type nu: float

    :return: r
    :rtype: float
    """

    r = (a*(1-e**2))/(1+e*np.cos(nu))
    return r

def flux_star(L,d):
    """
    Determines the flux received from a star (in W/m^2) at a distance d.

    :param L: the star luminosity (in W)
    :type L: float

    :param d: the distance (in m)
    :type d: float

    :return: F
    :rtype: float
    """

    F = L/(4*np.pi*d**2)
    return F

def flux_planet(F_star):
    """
    Determines the flux reemitted by a planet (in W/m^2) from the one it receives from its star considering the planet is a black body.

    :param F_star: the flux received by the planet from its star (in W/m^2)
    :type F_star: float

    :return: F_planet
    :rtype: float
    """

    F_planet = F_star/4
    return F_planet

def surface(R):
    """
    Determines the surface of a sphere of radius R.

    :param R: the radius (in m)
    :type R: float

    :return: S
    :rtype: float
    """

    S = 4*np.pi*R**2
    return S

def luminosity_planet_dayside(F_planet,R_planet):
    """
    Determines the luminosity of the dayside of a planet from the flux it reemits and its radius.

    :param F_planet: the flux reemitted by the planet's dayside (in W/m^2)
    :type F_planet: float

    :param R_planet: the planet radius (in m)
    :type R_planet: float

    :return: L_planet
    :rtype: float
    """

    L_planet = F_planet * surface(R_planet)/2
    return L_planet

def phase_curve(L_star, L_planet, phase_planet):
    """
    Determines the phase curve of a planet from its luminosity, its star's luminosity and its phase function expressed as the ratio between the planet and star's luminosities in ppm.

    :param L_star: the star luminosity (in W)
    :type L_star: float

    :param L_planet: the planet luminosity (in W)
    :type L_planet: float

    :param phase_planet: the phase function of the planet
    :type phase_planet: float

    :return: curve
    :rtype: float
    """

    curve = L_planet/L_star*phase_planet*10**6
    return curve


def main():
    # Constants
    L_Sun = 3.83E26 # Solar luminosity in W (from Wikipedia)
    M_Sun = 2E30 # Solar mass in kg (from Wikipedia)
    R_Sun = 6.96E8 # Solar radius in m (from Wikipedia)
    R_Earth = 6378E3 # Earth radius in m (from Wikipedia)

    t = np.linspace(0,30,10000) # time in days


    # For TRAPPIST-1 (using NASA Exoplanet Archive)
    T_eff_star = 2520 # in K (from Ducrot et al. 2020)
    R_star = 0.1234 * R_Sun # in m (from Ducrot et al. 2020)
    M_star = 0.0898 * M_Sun # in kg (from Ducrot et al. 2020)
    L_star = 10**(-3.26)*L_Sun # in W (from Ducrot et al. 2020)

    # For TRAPPIST-1 b (using NASA Exoplanet Archive)
    a_b = 20.13 * R_star # in m (from Ducrot et al. 2020)
    P_b = 1.51088432 # in days (from Ducrot et al. 2020)
    i_b = 89.28 * np.pi/180 # in rad (from Ducrot et al. 2020)
    omega_b = 336.86 * np.pi/180 # in rad (from Grimm et al. 2018)
    e_b = 0.00622 # (from Grimm et al. 2018)

    R_b = 1.017 * R_Earth # in m (from Grimm et al. 2018)

    nu_b = compute_true_anomaly(0,e_b,P_b,t)

    alpha_b = phase_angle(omega_b,nu_b,i_b)
    phase_b = phase_function(alpha_b)

    r_b = star_planet_separation(a_b,e_b,nu_b)

    flux_star_b = flux_star(L_star,r_b)
    flux_b = flux_planet(flux_star_b)
    L_b = luminosity_planet_dayside(flux_b,R_b)

    phase_curve_b = phase_curve(L_star,L_b,phase_b)


    #For TRAPPIST-1 c (using NASA Exoplanet Archive)

    a_c = 27.57 * R_star # in m (from Ducrot et al. 2020)
    P_c = 2.42179346 # in days (from Ducrot et al. 2020)
    i_c = np.radians(89.47) # in rad (from Ducrot et al. 2020)
    omega_c = np.radians(282.45) # in rad (from Grimm et al. 2018)
    e_c = 0.00654 # (from Grimm et al. 2018)

    R_c = 1.095 * R_Earth # in m (from Grimm et al. 2018)

    nu_c = compute_true_anomaly(0,e_c,P_c,t)
    alpha_c = phase_angle(omega_c,nu_c,i_c)
    phase_c = phase_function(alpha_c)

    r_c = star_planet_separation(a_c,e_c,nu_c)

    flux_star_c = flux_star(L_star,r_c)
    flux_c = flux_planet(flux_star_c)
    L_c = luminosity_planet_dayside(flux_c,R_c)

    phase_curve_c = phase_curve(L_star,L_c,phase_c)

    
    #For TRAPPIST-1 d (using NASA Exoplanet Archive)

    a_d = 38.85 * R_star # in m (from Ducrot et al. 2020)
    P_d = 4.04978035 # in days (from Ducrot et al. 2020)
    i_d = np.radians(89.65) # in rad (from Ducrot et al. 2020)
    omega_d = np.radians(-8.73) # in rad (from Grimm et al. 2018)
    e_d = 0.00837 # (from Grimm et al. 2018)

    R_d = 0.784 * R_Earth # in m (from Grimm et al. 2018)

    nu_d = compute_true_anomaly(0,e_d,P_d,t)
    alpha_d = phase_angle(omega_d,nu_d,i_d)
    phase_d = phase_function(alpha_d)

    r_d = star_planet_separation(a_d,e_d,nu_d)

    flux_star_d = flux_star(L_star,r_d)
    flux_d = flux_planet(flux_star_d)
    L_d = luminosity_planet_dayside(flux_d,R_d)

    phase_curve_d = phase_curve(L_star,L_d,phase_d)


    #For TRAPPIST-1 e (using NASA Exoplanet Archive)

    a_e = 51.0 * R_star # in m (from Ducrot et al. 2020)
    P_e = 6.09956479 # in days (from Ducrot et al. 2020)
    i_e = np.radians(89.663) # in rad (from Ducrot et al. 2020)
    omega_e = np.radians(108.37) # in rad (from Grimm et al. 2018)
    e_e = 0.000510 # (from Grimm et al. 2018)

    R_e = 0.910 * R_Earth # in m (from Grimm et al. 2018)

    nu_e = compute_true_anomaly(0,e_e,P_e,t)
    alpha_e = phase_angle(omega_e,nu_e,i_e)
    phase_e = phase_function(alpha_e)

    r_e = star_planet_separation(a_e,e_e,nu_e)

    flux_star_e = flux_star(L_star,r_e)
    flux_e = flux_planet(flux_star_e)
    L_e = luminosity_planet_dayside(flux_e,R_e)

    phase_curve_e = phase_curve(L_star,L_e,phase_e)


    #For TRAPPIST-1 f (using NASA Exoplanet Archive)

    a_f = 67.10 * R_star # in m (from Ducrot et al. 2020)
    P_f = 9.20659399 # in days (from Ducrot et al. 2020)
    i_f = np.radians(89.666) # in rad (from Ducrot et al. 2020)
    omega_f = np.radians(368.81) # in rad (from Grimm et al. 2018)
    e_f = 0.01007 # (from Grimm et al. 2018)

    R_f = 1.046 * R_Earth # in m (from Grimm et al. 2018)

    nu_f = compute_true_anomaly(0,e_f,P_f,t)
    alpha_f = phase_angle(omega_f,nu_f,i_f)
    phase_f = phase_function(alpha_f)

    r_f = star_planet_separation(a_f,e_f,nu_f)

    flux_star_f = flux_star(L_star,r_f)
    flux_f = flux_planet(flux_star_f)
    L_f = luminosity_planet_dayside(flux_f,R_f)

    phase_curve_f = phase_curve(L_star,L_f,phase_f)


    #For TRAPPIST-1 g (using NASA Exoplanet Archive)

    a_g = 81.7 * R_star # in m (from Ducrot et al. 2020)
    P_g = 12.35355570 # in days (from Ducrot et al. 2020)
    i_g = np.radians(89.698) # in rad (from Ducrot et al. 2020)
    omega_g = np.radians(191.34) # in rad (from Grimm et al. 2018)
    e_g = 0.00208 # (from Grimm et al. 2018)

    R_g = 1.148 * R_Earth # in m (from Grimm et al. 2018)

    nu_g = compute_true_anomaly(0,e_g,P_g,t)
    alpha_g = phase_angle(omega_g,nu_g,i_g)
    phase_g = phase_function(alpha_g)

    r_g = star_planet_separation(a_g,e_g,nu_g)

    flux_star_g = flux_star(L_star,r_g)
    flux_g = flux_planet(flux_star_g)
    L_g = luminosity_planet_dayside(flux_g,R_g)

    phase_curve_g = phase_curve(L_star,L_g,phase_g)


    #For TRAPPIST-1 h (using NASA Exoplanet Archive)

    a_h = 107.9 * R_star # in m (from Ducrot et al. 2020)
    P_h = 18.76727450 # in days (from Ducrot et al. 2020)
    i_h = np.radians(89.763) # in rad (from Ducrot et al. 2020)
    omega_h = np.radians(338.92) # in rad (from Grimm et al. 2018)
    e_h = 0.00567 # (from Grimm et al

    R_h = 0.773 * R_Earth # in m (from Grimm et al. 2018)

    nu_h = compute_true_anomaly(0,e_h,P_h,t)
    alpha_h = phase_angle(omega_h,nu_h,i_h)
    phase_h = phase_function(alpha_h)
    
    r_h = star_planet_separation(a_h,e_h,nu_h)

    flux_star_h = flux_star(L_star,r_h)
    flux_h = flux_planet(flux_star_h)
    L_h = luminosity_planet_dayside(flux_h,R_h)

    phase_curve_h = phase_curve(L_star,L_h,phase_h)


    # Plot

    plt.figure()
    plt.plot(t,phase_curve_b,label="b")
    plt.plot(t,phase_curve_c,label="c")
    plt.plot(t,phase_curve_d,label="d")
    plt.plot(t,phase_curve_e,label="e")
    plt.plot(t,phase_curve_f,label="f")
    plt.plot(t,phase_curve_g,label="g")
    plt.plot(t,phase_curve_h,label="h")
    plt.xlabel("Time (days)")
    plt.ylabel("$L_{planet}/L_{star}$ (ppm)")
    plt.title("Phase curves of planets of TRAPPIST-1 as black bodies")
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()