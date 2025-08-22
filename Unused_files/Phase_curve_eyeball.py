#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# February 2025
# Phase curves with an eyeball planet

import numpy as np
import matplotlib.pyplot as plt
from Orbital_motion import compute_true_anomaly
from Phase_curve_v1 import phase_angle, phase_function, star_planet_separation, surface_sphere, phase_curve, flux_star, flux_planet
from TRAPPIST1_parameters import *


# def flux_star_zenith(L,d,theta):
#     """
#     Determines the flux received from a star (in W/m^2) at a distance d with a zenith angle theta.

#     :param L: the star luminosity (in W)
#     :type L: float

#     :param d: the distance (in m)
#     :type d: float

#     :param theta: the zenith angle (in rad)
#     :type theta: float

#     :return: F
#     :rtype: float
#     """

#     F = flux_star(L,d)*np.cos(theta)
#     return F

def flux_planet_latitude(F_star,lat):
    """
    Determines the flux reemitted by a planet (in W/m^2) at a given latitude.

    :param F_star: the flux received from the star (in W/m^2)
    :type F_star: float

    :param lat: the latitude (in rad)
    :type lat: float

    :return: F
    :rtype: float
    """

    F = F_star*np.cos(lat)
    return F


def fraction_ring(lat_in,phase):
    """
    Determines the fraction of a ring of the dayside of an eyeball planet visible at a given phase.

    :param lat_in: the latitude of the inner radius of the ring (in rad)
    :type lat_in: float

    :param phase: the phase of the planet (in rad)
    :type phase: float

    :return: f
    :rtype: float
    """

    f = []
    for theta in phase:
        if np.abs(np.tan(lat_in) * np.tan(theta))<1:
            f.append(1/np.pi * np.arccos(-np.tan(lat_in) * np.tan(theta)))
        else:
            f.append(0)

    # if np.abs(np.tan(lat_in) * np.tan(phase))<1:
    #     f = 1/np.pi * np.arccos(-np.tan(lat_in) * np.tan(phase))
    # else:
    #     f = 0
    
    return f


def fraction_disk(lat_out,phase):
    """
    Determines the fraction of the central disk of the dayside of an eyeball planet visible at a given phase.

    :param lat_out: the latitude of the radius of the disk (in rad)
    :type lat_out: float

    :param phase: the phase of the planet (in rad)
    :type phase: float

    :return: f
    :rtype: float
    """
    f = []
    for theta in phase:
        if lat_out<theta:
            f.append(1)
        elif np.abs(np.tan(lat_out) * np.tan(theta))<1:
            f.append(1/np.pi * np.arccos(-np.tan(lat_out) * np.tan(theta)))
        else:
            f.append(0)

    # if lat_out<phase:
    #     f = 1
    # elif np.abs(np.tan(lat_out) * np.tan(phase))<1:
    #     f = 1/np.pi * np.arccos(-np.tan(lat_out) * np.tan(phase))
    # else:
    #     f = 0
    
    return np.array(f)


def area_ring(r_in,r_out):
    """
    Determines the area of a ring between two radii.

    :param r_in: the inner radius
    :type r_in: float

    :param r_out: the outer radius
    :type r_out: float

    :return: A
    :rtype: float
    """

    A = np.pi*(r_out**2-r_in**2)
    return A


def area_disk(r):
    """
    Determines the area of a disk of radius r.

    :param r: the radius
    :type r: float

    :return: A
    :rtype: float
    """

    A = np.pi*r**2
    return A


def luminosity_planet_eyeball(F_planet,R_planet,n,phase):
    """
    Determines the luminosity of an eyeball planet with its dayside divided in a central disk and n rings at a given phase.

    :param F_planet: the flux reemitted by the planet (in W/m^2)
    :type F_planet: float

    :param R_planet: the radius of the planet (in m)
    :type R_planet: float

    :param n: the number of rings
    :type n: int

    :param phase: the phase of the planet (in rad)
    :type phase: float

    :return: L_planet
    :rtype: float
    """

    lat = np.linspace(0,np.pi/2,n+1)

    r_disk = R_planet*np.sin(lat[1])
    A_disk = area_disk(r_disk)
    f_disk = fraction_disk(lat[1],phase)

    L_disk = F_planet[0]*A_disk*f_disk

    L_rings = 0
    for i in range(2,n):
        r_in = R_planet*np.sin(lat[i])
        r_out = R_planet*np.sin(lat[i+1])
        A_ring = area_ring(r_in,r_out)
        f_ring = fraction_ring(lat[i],phase)
        L_rings += F_planet*A_ring*f_ring
    
    L_planet = L_disk + L_rings

    return L_planet


def main():
    
    t_end = 30 # simulation duration in days
    nb_points = 1000 # number of points in the time and phase arrays

    t = np.linspace(0,t_end,nb_points) # time array in days

    #nb_rings = 100

    #latitudes = np.linspace(0,np.pi/2,nb_rings+1)


    # For TRAPPIST-1 b
    
    nu_b = compute_true_anomaly(0,e_b,P_b,t)

    alpha_b = phase_angle(omega_b,nu_b,i_b)
    phase_b = phase_function(alpha_b)

    r_b = star_planet_separation(a_b,e_b,nu_b)

    flux_star_b = flux_star(L_star,r_b)
    flux_b = flux_planet_latitude(flux_star_b,phase_b)
    L_b = luminosity_planet_eyeball(flux_b,R_b,nb_points,phase_b)

    phase_curve_b = L_b/L_star*10**6



    # Plot

    plt.figure()
    plt.plot(t,phase_b,label="b")
    plt.show()

    plt.figure()
    plt.plot(t,phase_curve_b,label="b")
    plt.xlabel("Time (days)")
    plt.ylabel("$L_{planet}/L_{star}$ (ppm)")
    plt.title("Phase curves of planets of TRAPPIST-1 as eyeball black bodies")
    plt.legend()
    plt.grid()
    plt.show()



if __name__ == "__main__":
    main()

