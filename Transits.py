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


def total_transit_duration(P, a, R_star, R_planet, i):
    """
    Determines the total duration of an exoplanet transit (in s) for a circular orbit.

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

    :return: t_total
    :rtype: float
    """

    t_total = P/np.pi * np.arcsin(R_star/a * np.sqrt(((1+R_planet/R_star)**2 - (a/R_star * np.cos(i))**2) / (1-np.cos(i)**2)))
    return t_total


def flat_transit_duration(P, a, R_star, R_planet, i):
    """
    Determines the flat duration of an exoplanet transit (in s) for a circular orbit.

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

    :return: t_flat
    :rtype: float
    """

    t_flat = P/np.pi * np.arcsin(np.sin(total_transit_duration(P, a, R_star, R_planet, i)*np.pi/P) * np.sqrt(((1-R_planet/R_star)**2-(a/R_star*np.cos(i))**2)/((1+R_planet/R_star)**2-(a/R_star*np.cos(i))**2)))
    return t_flat


def eclipse_phase(P, a, R_star, R_planet, i):
    """
    Determines the phases of an exoplanet for which its secondary eclipse starts and ends for a circular orbit.

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

    :return: phase_eclipse_start, phase_eclipse_end
    :rtype: float
    """

    t_eclipse = total_transit_duration(P, a, R_star, R_planet, i)#+flat_transit_duration(P, a, R_star, R_planet, i))/2
    phase_eclipse_start = 1-t_eclipse/(2*P)
    phase_eclipse_end = t_eclipse/(2*P)

    return phase_eclipse_start, phase_eclipse_end


def eclipse(P, a, R_star, R_planet, i, phase):
    """
    Determines if an exoplanet is in eclipse or not.

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

    :return: in_eclipse
    :rtype: bool
    """

    phase_eclipse_start, phase_eclipse_end = eclipse_phase(P, a, R_star, R_planet, i)

    in_eclipse = (phase_eclipse_start < phase-np.trunc(phase)) + (phase-np.trunc(phase) < phase_eclipse_end)
    return in_eclipse


def main():
    phase_b = np.array([0, 0.25, 0.5, 0.75, 1])
    print(transit_depth(R_b, R_star)*100, '%')
    print("Total transit duration:",total_transit_duration(P_b*24, a_b, R_star, R_b, i_b), 'h')
    print("Flat transit duration:",flat_transit_duration(P_b*24, a_b, R_star, R_b, i_b), 'h')
    print(eclipse_phase(P_b*24*3600, a_b, R_star, R_b, i_b))
    print(eclipse(P_b*24*3600, a_b, R_star, R_b, i_b, phase_b))

if __name__ == "__main__":
    main()