#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# May 2025
# Phase curves with wavelength

import numpy as np
import matplotlib.pyplot as plt
from Orbital_motion import compute_true_anomaly
from Transits import eclipse
from Phase_curve_v1 import star_planet_separation, flux_star, flux_planet, luminosity_planet_dayside, phase_curve
from TRAPPIST1_parameters import *

def Planck_law(wavelength, T):
    """
    Determines the spectral radiance of a black body (in W/m^2 m^-1 sr^-1).

    :param wavelength: the wavelength (in m)
    :type wavelength: float

    :param T: the temperature (in K)
    :type T: float

    :return: B
    :rtype: float
    """

    h = 6.62607015e-34 # Planck constant in J.s
    c = 2.99792458e8 # speed of light in m/s
    k = 1.380649e-23 # Boltzmann constant in J/K

    B = (2*h*c**2)/(wavelength**5)*(1/(np.exp(h*c/(wavelength*k*T))-1))
    return B

# def Planck_law_integral_solid_angle(B):
#     """
#     Integrates the Planck law over all directions (in W/m^2 m^-1).

#     :param B: the spectral radiance (in W/m^2 m^-1 sr^-1)
#     :type B: float

#     :return: B_theta
#     :rtype: float
#     """

#     B_theta = np.pi*B
#     return B_theta

def flux_black_body_wave_unit(wavelength,T):
    """
    Determines the flux of a black body per unit wavelength (in W/m^2 m^-1).

    :param wavelength: the wavelength (in m)
    :type wavelength: float

    :param T: the temperature (in K)
    :type T: float

    :return: F
    :rtype: float
    """

    B = Planck_law(wavelength,T)
    # B_theta = Planck_law_integral_solid_angle(B)
    B_theta = B * np.pi
    return B_theta

def flux_black_body(wavelengths,T):
    """
    Determines the flux of a black body (in W/m^2).

    :param wavelengths: the wavelength array (in m)
    :type wavelengths: float

    :param T: the temperature (in K)
    :type T: float

    :return: F
    :rtype: float
    """

    F = np.trapz(flux_black_body_wave_unit(wavelengths,T),wavelengths)
    return F

def main():

    nb_points = 10000
    
    # Define MIRI 15 micron filter band

    lambda_min = 13e-6 # in m
    lambda_max = 17.3e-6 # in m

    wavelengths = np.linspace(lambda_min, lambda_max, nb_points)

    # Flux of star TRAPPIST-1 in the MIRI 15 micron filter band
    F_star = flux_black_body(wavelengths, T_eff_star)
    print("F_star = ", F_star)


if __name__ == "__main__":
    main()