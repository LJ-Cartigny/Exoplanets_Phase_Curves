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


def main():
    # Constants
    L_Sun = 3.83E26 # Solar luminosity in W (from Wikipedia)
    M_Sun = 2E30 # Solar mass in kg (from Wikipedia)
    R_Sun = 6.96E8 # Solar radius in m (from Wikipedia)
    R_Earth = 6378E3 # Earth radius in m (from Wikipedia)

    t = np.linspace(0,10,1000) # time in days


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

    #nu_b = np.linspace(0,8*np.pi,10000) # in rad
    nu_b = compute_true_anomaly(0,e_b,P_b,t)
    #print(nu_b.shape)

    alpha_b = phase_angle(omega_b,nu_b,i_b)
    #print(nu_b.shape)
    phase_b = phase_function(alpha_b)
    #print(phase_b.shape)

    r_b = star_planet_separation(a_b,e_b,nu_b)

    flux_star_b = flux_star(L_star,r_b)
    flux_b = flux_planet(flux_star_b)
    L_b = luminosity_planet_dayside(flux_b,R_b)



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

    # print(L_b/L_star)
    # print(L_c/L_star)

    # Plot

    plt.figure()
    #plt.plot(nu_b, flux_b/flux_star_b*phase_b)
    #plt.plot(t, flux_b/flux_star_b*phase_b, label="b")
    #plt.plot(t, flux_c/flux_star_c*phase_c, label="c")
    plt.plot(t,L_b/L_star * phase_b * 10**6,label="b")
    plt.plot(t,L_c/L_star * phase_c * 10**6,label="c")
    #plt.xlabel("True anomaly $\\nu$ (rad)")
    plt.xlabel("Time (days)")
    plt.ylabel("$L_{planet}/L_{star}$ (ppm)")
    plt.title("Phase curves of TRAPPIST-1 b and c as black bodies")
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()