#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# February 2025
# Exoplanet transit

import numpy as np
import matplotlib.pyplot as plt
from TRAPPIST1_parameters import *
#from Phase_curve_v1 import phase_angle, phase_function, star_planet_separation, surface_sphere, phase_curve, flux_star, flux_planet

def transit_depth(R_planet, R_star):
    """
    Determines the depth of an exoplanet transit.

    :param R_planet: the radius of the planet (in m)
    :type Rp: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :return: delta_F
    :rtype: float
    """

    delta_F = (R_planet/R_star)**2
    return delta_F

def transit_impact_parameter(a, i, e, R_star, omega):
    """
    Determines the impact parameter of an exoplanet transit.

    :param a: the semimajor axis (in m)
    :type a: float

    :param i: the inclination (in rad)
    :type i: float

    :param e: the eccentricity
    :type e: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :param omega: the argument of pericentre (in rad)
    :type omega: float

    :return: b
    :rtype: float
    """

    b = a/R_star * np.cos(i) * (1-e**2)/(1+e*np.sin(omega))
    return b

def eclipse_impact_parameter(a, i, e, R_star, omega):
    """
    Determines the impact parameter of an exoplanet eclipse.

    :param a: the semimajor axis (in m)
    :type a: float

    :param i: the inclination (in rad)
    :type i: float

    :param e: the eccentricity
    :type e: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :param omega: the argument of pericentre (in rad)
    :type omega: float

    :return: b
    :rtype: float
    """

    b = a/R_star * np.cos(i) * (1-e**2)/(1-e*np.sin(omega))
    return b


def total_transit_duration(P, a, R_star, R_planet, i, e, omega, b):
    """
    Determines the total duration of an exoplanet transit (in s).

    :param P: the orbital period (in s)
    :type P: float

    :param a: the semimajor axis (in m)
    :type a: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :param R_planet: the radius of the planet (in m)
    :type R_planet: float

    :param i: the inclination (in rad)
    :type i: float

    :param e: the eccentricity
    :type e: float

    :param omega: the argument of pericentre (in rad)
    :type omega: float

    :param b: the transit impact parameter
    :type b: float

    :return: t_total
    :rtype: float
    """

    t_total = P/np.pi * np.arcsin(R_star/a * np.sqrt(((1+R_planet/R_star)**2 - (a/R_star * np.cos(i))**2) / (1-np.cos(i)**2)))

    #t_total = P/np.pi * np.arcsin(np.sqrt((1+R_planet/R_star)**2 - b**2) / (1-np.cos(i)**2)) * np.sqrt(1-e**2)/(1+e*np.sin(omega))

    return t_total


def flat_transit_duration(P, a, R_star, R_planet, i, e, omega, b):
    """
    Determines the flat duration of an exoplanet transit (in s).

    :param P: the orbital period (in s)
    :type P: float

    :param a: the semimajor axis (in m)
    :type a: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :param R_planet: the radius of the planet (in m)
    :type R_planet: float

    :param i: the inclination (in rad)
    :type i: float

    :param e: the eccentricity
    :type e: float

    :param omega: the argument of pericentre (in rad)
    :type omega: float

    :param b: the transit impact parameter
    :type b: float

    :return: t_flat
    :rtype: float
    """

    t_flat = P/np.pi * np.arcsin(np.sin(total_transit_duration(P, a, R_star, R_planet, i, e, omega, b)*np.pi/P) * np.sqrt(((1-R_planet/R_star)**2-(a/R_star*np.cos(i))**2)/((1+R_planet/R_star)**2-(a/R_star*np.cos(i))**2)))
    
    #t_flat = P/np.pi * np.arcsin(np.sqrt((1-R_planet/R_star)**2 - b**2) / (1-np.cos(i)**2)) * np.sqrt(1-e**2)/(1+e*np.sin(omega))

    return t_flat


def total_eclipse_duration(P, a, R_star, R_planet, i, e, omega, b):
    """
    Determines the total duration of an exoplanet eclipse (in s).

    :param P: the orbital period (in s)
    :type P: float

    :param a: the semimajor axis (in m)
    :type a: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :param R_planet: the radius of the planet (in m)
    :type R_planet: float

    :param i: the inclination (in rad)
    :type i: float

    :param e: the eccentricity
    :type e: float

    :param omega: the argument of pericentre (in rad)
    :type omega: float

    :param b: the eclipse impact parameter
    :type b: float

    :return: t_total
    :rtype: float
    """

    t_total = P/np.pi * np.arcsin(np.sqrt((1+R_planet/R_star)**2 - b**2) / (1-np.cos(i)**2)) * np.sqrt(1-e**2)/(1-e*np.sin(omega))

    return t_total


def flat_eclipse_duration(P, a, R_star, R_planet, i, e, omega, b):
    """
    Determines the flat duration of an exoplanet eclipse (in s).

    :param P: the orbital period (in s)
    :type P: float

    :param a: the semimajor axis (in m)
    :type a: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :param R_planet: the radius of the planet (in m)
    :type R_planet: float

    :param i: the inclination (in rad)
    :type i: float

    :param e: the eccentricity
    :type e: float

    :param omega: the argument of pericentre (in rad)
    :type omega: float

    :param b: the eclipse impact parameter
    :type b: float

    :return: t_flat
    :rtype: float
    """

    t_flat = P/np.pi * np.arcsin(np.sqrt((1-R_planet/R_star)**2 - b**2) / (1-np.cos(i)**2)) * np.sqrt(1-e**2)/(1-e*np.sin(omega))

    return t_flat


def eclipse_phase(P, a, R_star, R_planet, i, e, omega, b):
    """
    Determines the phases of an exoplanet for which its secondary eclipse starts and ends (centered at 0 or 1).

    :param P: the orbital period (in s)
    :type P: float

    :param a: the semimajor axis (in m)
    :type a: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :param R_planet: the radius of the planet (in m)
    :type R_planet: float

    :param i: the inclination (in rad)
    :type i: float

    :param e: the eccentricity
    :type e: float

    :param omega: the argument of pericentre (in rad)
    :type omega: float

    :param b: the impact parameter
    :type b: float

    :return: phase_eclipse_start, phase_eclipse_end
    :rtype: float
    """

    #t_eclipse = total_transit_duration(P, a, R_star, R_planet, i)#+flat_transit_duration(P, a, R_star, R_planet, i))/2

    # t_eclipse = total_eclipse_duration(P, a, R_star, R_planet, i, e, omega, b)

    t_eclipse = total_transit_duration(P, a, R_star, R_planet, i, e, omega, b)

    phase_eclipse_start = 1-t_eclipse/(2*P)
    phase_eclipse_end = t_eclipse/(2*P)

    return phase_eclipse_start, phase_eclipse_end


def eclipse(P, a, R_star, R_planet, i, phase, e, omega, b):
    """
    Determines if an exoplanet is in eclipse or not at a given phase.

    :param P: the orbital period (in s)
    :type P: float

    :param a: the semimajor axis (in m)
    :type a: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :param R_planet: the radius of the planet (in m)
    :type R_planet: float

    :param i: the inclination (in rad)
    :type i: float

    :param phase: the phase of the exoplanet (in rad)
    :type phase: float

    :param e: the eccentricity
    :type e: float

    :param omega: the argument of pericentre (in rad)
    :type omega: float

    :param b: the impact parameter
    :type b: float

    :return: in_eclipse
    :rtype: bool
    """

    phase_eclipse_start, phase_eclipse_end = eclipse_phase(P, a, R_star, R_planet, i, e, omega, b)

    # in_eclipse = (phase_eclipse_start < phase-np.trunc(phase)) + (phase-np.trunc(phase) < phase_eclipse_end)

    in_eclipse = (phase_eclipse_start < phase-np.trunc(phase)) + (phase-np.trunc(phase) < phase_eclipse_end)

    return in_eclipse


def main():
    print("TRAPPIST-1b")
    # phase_b = np.array([0, 0.25, 0.5, 0.75, 1])
    phase_b = np.linspace(0,1,10001,endpoint=True)
    b_b = eclipse_impact_parameter(a_b, i_b, e_b, R_star, omega_b)
    print("Impact parameter:",b_b)
    print("Transit depth:",transit_depth(R_b, R_star)*100, '%')
    # print(np.sqrt((1+R_b/R_star)**2 - b_b**2)/(1-np.cos(i_b)**2))
    # print("Total eclipse duration:",total_eclipse_duration(P_b*24, a_b, R_star, R_b, i_b, e_b, omega_b, b_b), 'h')
    # print("Flat eclipse duration:",flat_eclipse_duration(P_b*24, a_b, R_star, R_b, i_b, e_b, omega_b, b_b), 'h')
    print("Total transit duration:",total_transit_duration(P_b*24, a_b, R_star, R_b, i_b, e_b, omega_b, b_b), 'h')
    print("Flat transit duration:",flat_transit_duration(P_b*24, a_b, R_star, R_b, i_b, e_b, omega_b, b_b), 'h')
    print(eclipse_phase(P_b*24*3600, a_b, R_star, R_b, i_b, e_b, omega_b, b_b))
    print((1+eclipse_phase(P_b*24*3600, a_b, R_star, R_b, i_b, e_b, omega_b, b_b)[1]-eclipse_phase(P_b*24*3600, a_b, R_star, R_b, i_b, e_b, omega_b, b_b)[0])*P_b*24,"h")
    print(eclipse(P_b*24*3600, a_b, R_star, R_b, i_b, phase_b, e_b, omega_b, b_b))
    
    k = 0
    while k < len(phase_b):
        if eclipse(P_b*24*3600, a_b, R_star, R_b, i_b, phase_b[k], e_b, omega_b, b_b) == True:
            print(phase_b[k])
        k += 1

    # print("TRAPPIST-1c")
    # phase_c = np.array([0, 0.25, 0.5, 0.75, 1])
    # b_c = eclipse_impact_parameter(a_c, i_c, e_c, R_star, omega_c)
    # print("Impact parameter:",b_c)
    # print("Transit depth:",transit_depth(R_c, R_star)*100, '%')
    # # print(np.sqrt((1+R_c/R_star)**2 - b_c**2)/(1-np.cos(i_c)**2))
    # # print("Total eclipse duration:",total_eclipse_duration(P_c*24, a_c, R_star, R_c, i_c, e_c, omega_c, b_c), 'h')
    # # print("Flat eclipse duration:",flat_eclipse_duration(P_c*24, a_c, R_star, R_c, i_c, e_c, omega_c, b_c), 'h')
    # print("Total transit duration:",total_transit_duration(P_c*24, a_c, R_star, R_c, i_c, e_c, omega_c, b_c), 'h')
    # print("Flat transit duration:",flat_transit_duration(P_c*24, a_c, R_star, R_c, i_c, e_c, omega_c, b_c), 'h')
    # print(eclipse_phase(P_c*24*3600, a_c, R_star, R_c, i_c, e_c, omega_c, b_c))
    # print((1+eclipse_phase(P_c*24*3600, a_c, R_star, R_c, i_c, e_c, omega_c, b_c)[1]-eclipse_phase(P_c*24*3600, a_c, R_star, R_c, i_c, e_c, omega_c, b_c)[0])*P_c*24,"h")
    # print(eclipse(P_c*24*3600, a_c, R_star, R_c, i_c, phase_c, e_c, omega_c, b_c))


if __name__ == "__main__":
    main()