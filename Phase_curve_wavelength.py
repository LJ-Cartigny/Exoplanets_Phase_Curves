#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# May 2025
# Phase curves with wavelength

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.constants import c, h, k, sigma
from tqdm import tqdm
import time

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

    # h = 6.62607015e-34 # Planck constant in J.s
    # c = 2.99792458e8 # speed of light in m/s
    # k = 1.380649e-23 # Boltzmann constant in J/K

    B = (2*h*c**2)/(wavelength**5)*(1/(np.exp(h*c/(wavelength*k*T))-1))

    return B


def flux_black_body(lambda_min, lambda_max,T):
    """
    Determines the flux of a black body (in W/m^2) over a range of wavelentgths.

    :param lambda_min: the minimum wavelength (in m)
    :type lambda_min: float

    :param lambda_max: the maximum wavelength (in m)
    :type lambda_max: float

    :param T: the temperature (in K)
    :type T: float

    :return: F
    :rtype: float
    """

    F = quad(lambda l: np.pi*Planck_law(l,T), lambda_min, lambda_max)[0]

    # nb_points = 10000
    # wavelengths = np.linspace(lambda_min, lambda_max, nb_points)
    # B = np.pi*Planck_law(wavelengths, T)
    # F = np.trapz(B, wavelengths)  # Integrate using the trapezoidal rule

    return F


def planet_equilibirium_temperature(T_star, R_star, d, albedo=0):
    """
    Determines the equilibrium temperature of the day side of a tidally locked planet (in K) without any redistribution.

    :param T_star: the effective temperature of the star (in K)
    :type T_star: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :param d: the distance between the star and the planet (in m)
    :type d: float

    :param albedo: the albedo of the planet (default: 0)
    :type albedo: float

    :return: T_eq
    :rtype: float
    """

    # L = 4 * np.pi * R_star**2 * sigma * T_eff_star**4
    # T_eq = ((1-albedo) * 2/3 * L / (4 * np.pi * d**2 * sigma))**0.25
    T_eq = T_star * (R_star / d)**0.5 * (2/3 * (1-albedo))**0.25

    return T_eq


def flux_ratio(R_planet, R_star, T_star, d, lambda_min, lambda_max):
    """
    Determines the flux ratio between the planet and the star (in ppm).

    :param F_planet: the flux of the planet (in W/m^2)
    :type F_planet: float

    :param F_star: the flux of the star (in W/m^2)
    :type F_star: float

    :param R_planet: the radius of the planet (in m)
    :type R_planet: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :return: F_ratio
    :rtype: float
    """

    T_planet = planet_equilibirium_temperature(T_star, R_star, d)
    F_ratio = (R_planet/R_star)**2 * flux_black_body(lambda_min, lambda_max, T_planet) / flux_black_body(lambda_min, lambda_max, T_star) * 1e6

    return F_ratio


def conversion_IS_to_mJy(F, wavelength, dist, R):
    """
    Converts the flux density (in W/m^2/m) to mJy.

    :param F: the flux density (in W/m^2/m)
    :type F: float

    :param wavelength: the wavelength (in m)
    :type wavelength: float

    :param dist: the distance of the object (in m)
    :type dist: float

    :param R: the radius of the object (in m)
    :type R: float

    :return: F_mJy
    :rtype: float
    """

    conv_factor = (dist/R)**-2 * wavelength**2/c * 1e26 * 1e3
    F_mJy = F * conv_factor

    return F_mJy


def flux_mJy(F, lambda_min, lambda_max, dist, R):
    """
    Compute the flux of an object in mJy over a range of wavelength.

    :param F: the flux density (in W/m^2/m)
    :type F: float

    :param lambda_min: the minimum wavelength (in m)
    :type lambda_min: float

    :param lambda_max: the maximum wavelength (in m)
    :type lambda_max: float

    :param dist: the distance of the object (in m)
    :type dist: float

    :param R: the radius of the object (in m)
    :type R: float

    :return: F_mJy
    :rtype: float
    """

    F_mJy = quad(lambda l: conversion_IS_to_mJy(F, l, dist, R), lambda_min, lambda_max)[0]

    return F_mJy

def conversion_mJy_to_IS(F_mJy, wavelength, dist, R):
    """
    Converts the flux density (in mJy) to W/m^2/m.

    :param F_mJy: the flux density (in mJy)
    :type F_mJy: float

    :param wavelength: the wavelength (in m)flux_T1_sphinx_cut *= QE
    :type wavelength: float

    :param dist: the distance of the object (in m)
    :type dist: float

    :param R: the radius of the object (in m)
    :type R: float

    :return: F
    :rtype: float
    """

    conv_factor = (dist/R)**2 * c/wavelength**2 * 1e-26 * 1e-3
    F = F_mJy * conv_factor

    return F

def flux_Wm2(F_mJy, lambda_min, lambda_max, dist, R):
    """
    Compute the flux of an object in W/m^2/m over a range of wavelength.

    :param F_mJy: the flux density (in mJy)
    :type F_mJy: float

    :param lambda_min: the minimum wavelength (in m)
    :type lambda_min: float

    :param lambda_max: the maximum wavelength (in m)
    :type lambda_max: float

    :param dist: the distance of the object (in m)
    :type dist: float

    :param R: the radius of the object (in m)
    :type R: float

    :return: F
    :rtype: float
    """

    F = quad(lambda l: conversion_mJy_to_IS(F_mJy, l, dist, R), lambda_min, lambda_max)[0]

    return F

def filter(filter_name):
    """
    Returns the filter band of the specified filter.

    :param filter_name: the name of the filter
    :type filter_name: str

    :return: filter_band
    :rtype: np.ndarray
    """

    if filter_name == "F560W":
        filter_band = np.loadtxt("miri_filter.csv", delimiter=',', unpack=True, skiprows=1, usecols=(0, 1))

    elif filter_name == "F770W":
        filter_band = np.loadtxt("miri_filter.csv", delimiter=',', unpack=True, skiprows=1, usecols=(0, 2))

    elif filter_name == "F1000W":
        filter_band = np.loadtxt("miri_filter.csv", delimiter=',', unpack=True, skiprows=1, usecols=(0, 3))

    elif filter_name == "F1130W":
        filter_band = np.loadtxt("miri_filter.csv", delimiter=',', unpack=True, skiprows=1, usecols=(0, 4))

    elif filter_name == "F1280W":
        filter_band = np.loadtxt("miri_filter.csv", delimiter=',', unpack=True, skiprows=1, usecols=(0, 5))

    elif filter_name == "F1500W":
        filter_band = np.loadtxt("miri_filter.csv", delimiter=',', unpack=True, skiprows=1, usecols=(0, 6))

    elif filter_name == "F1800W":
        filter_band = np.loadtxt("miri_filter.csv", delimiter=',', unpack=True, skiprows=1, usecols=(0, 7))

    elif filter_name == "F2100W":
        filter_band = np.loadtxt("miri_filter.csv", delimiter=',', unpack=True, skiprows=1, usecols=(0, 8))

    elif filter_name == "F2550W":
        filter_band = np.loadtxt("miri_filter.csv", delimiter=',', unpack=True, skiprows=1, usecols=(0, 9))


    return filter_band


def quantum_efficiency(filter_name, wavelength):
    """
    Returns the quantum efficiency of the specified filter at the given wavelength.

    :param filter_name: the name of the filter
    :type filter_name: str

    :param wavelength: the wavelength (in m)
    :type wavelength: float

    :return: QE
    :rtype: float
    """

    filter_band = filter(filter_name)
    wavelengths = filter_band[0,:]*1e-6 # Convert to m

    # lambda_min = np.min(wavelengths)
    # lambda_max = np.max(wavelengths)

    # i = 0
    # while wavelengths_T1_sphinx[i] < lambda_min:
    #     i += 1
    # i_min = i
    # while wavelengths_T1_sphinx[i] < lambda_max:
    #     i += 1
    # i_max = i

    # wavelengths_T1_sphinx_cut = wavelengths_T1_sphinx[i_min:i_max]

    # QE = np.interp(wavelength, wavelengths, filter_band[1,:]/np.max(filter_band[1,:]))

    QE = np.interp(wavelength, wavelengths, filter_band[1,:])

    # QE = np.interp(wavelengths_T1_sphinx_cut, wavelengths, filter_band[1,:])

    return QE


def flux_star_miri(filter_name):
    """
    Returns the flux of the star in the specified filter band.

    :param filter_name: the name of the filter
    :type filter_name: str

    :return: F_star
    :rtype: float
    """

    filter_band = filter(filter_name)
    wavelengths = filter_band[0,:]*1e-6 # Convert to m

    # print(np.min(wavelengths), np.max(wavelengths))
    # print(wavelengths.shape)

    # print(np.min(wavelengths_T1_sphinx), np.max(wavelengths_T1_sphinx))
    # print(wavelengths_T1_sphinx.shape)

    lambda_min = np.min(wavelengths)
    lambda_max = np.max(wavelengths)

    i = 0
    while wavelengths_T1_sphinx[i] < lambda_min:
        i += 1
    i_min = i
    while wavelengths_T1_sphinx[i] < lambda_max:
        i += 1
    i_max = i

    flux_T1_sphinx_cut = flux_T1_sphinx[i_min:i_max]
    wavelengths_T1_sphinx_cut = wavelengths_T1_sphinx[i_min:i_max]
    flux_interp = interp1d(wavelengths_T1_sphinx_cut, flux_T1_sphinx_cut, kind='linear', bounds_error=False, fill_value=0)

    # mid_point = (lambda_min + lambda_max) / 2
    # F_star_1 = quad(lambda l: quantum_efficiency(filter_name, l) * flux_interp(l), lambda_min, mid_point)[0]
    # F_star_2 = quad(lambda l: quantum_efficiency(filter_name, l) * flux_interp(l), mid_point, lambda_max)[0]
    # F_star = F_star_1 + F_star_2

    # print(flux_T1_sphinx.shape)

    # flux_T1_sphinx_interp = np.interp(wavelengths, wavelengths_T1_sphinx, flux_T1_sphinx)

    # print(wavelengths.shape)
    # print(flux_T1_sphinx_interp.shape)

    # F_star = quad(lambda l: quantum_efficiency(filter_name, l)*flux_interp(l), lambda_min, lambda_max, limit=100)[0]

    # F_star = quad(lambda l: Planck_law(l, T_eff_star)*quantum_efficiency(filter_name, l), wavelengths[0], wavelengths[-1])[0]
    # F_star = np.trapz(flux_T1_sphinx_interp*quantum_efficiency(filter_name,wavelengths), wavelengths)
    
    F_star = 0
    for l in tqdm(wavelengths_T1_sphinx_cut):
        time.sleep(0.1)
        QE = quantum_efficiency(filter_name, l)
        F_star += flux_interp(l) * QE

    F_star *= (wavelengths_T1_sphinx_cut[1]-wavelengths_T1_sphinx_cut[0])

    # F_star = np.sum(flux_T1_sphinx_cut*quantum_efficiency(filter_name))*(wavelengths_T1_sphinx[1]-wavelengths_T1_sphinx[0])

    return F_star

def flux_miri2(filter_name):
    """
    Returns the flux of the star in the specified filter band.

    :param filter_name: the name of the filter
    :type filter_name: str

    :return: F_star
    :rtype: float
    """

    filter_band = filter(filter_name)
    wavelengths_filter = filter_band[0,:]*1e-6 # Convert to m
    QE = filter_band[1,:]/np.max(filter_band[1,:])
    
    # The two arrays are not covering the same wevelength range so we need to deterùine the common range
    lambda_min_common = np.max([np.min(wavelengths_T1_sphinx), np.min(wavelengths_filter)])
    lambda_max_common = np.min([np.max(wavelengths_T1_sphinx), np.max(wavelengths_filter)])

    mask_common = (wavelengths_T1_sphinx >= lambda_min_common) & (wavelengths_T1_sphinx <= lambda_max_common)
    wavelengths_T1_sphinx_cut = wavelengths_T1_sphinx[mask_common]
    flux_T1_sphinx_cut = flux_T1_sphinx[mask_common]

    # Interpolate the QE values to match the wavelengths of the SPHINX spectrum
    interp_filter = interp1d(wavelengths_filter, QE, bounds_error=False, fill_value=0)
    QE_interp = interp_filter(wavelengths_T1_sphinx_cut)

    F_star = np.trapz(flux_T1_sphinx_cut * QE_interp, wavelengths_T1_sphinx_cut)

    return F_star


def main():

    nb_points = 10000

    # Comparison of the SPHINX and black body spectra of TRAPPIST-1

    # flux_Planck = Planck_law(wavelengths_T1_sphinx, T_eff_star)

    # plt.figure(figsize=(16,9))
    # plt.plot(wavelengths_T1_sphinx, flux_T1_sphinx, label="SPHINX")
    # plt.plot(wavelengths_T1_sphinx, flux_Planck, label="Black body")
    # plt.xlabel(r"Wavelength ($m$)")
    # plt.ylabel(r"Flux ($W/m^2/m$)")
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.xlim(np.min(wavelengths_T1_sphinx), np.max(wavelengths_T1_sphinx))
    # plt.title("Comparison of the SPHINX and black body spectra of TRAPPIST-1")
    # plt.legend()
    # plt.grid()
    # plt.show()
    
    # Define MIRI 15 micron filter band

    lambda_min_F1500 = 13e-6 # in m
    lambda_max_F1500 = 17.3e-6 # in m

    # QE = quantum_efficiency("F1500W", 15e-6)
    # print("QE = ", QE)

    #wavelengths = np.linspace(lambda_min, lambda_max, nb_points)

    # Flux of star TRAPPIST-1 in the MIRI 15 micron filter band
    F_star = flux_black_body(lambda_min_F1500, lambda_max_F1500, T_eff_star)
    print("F_star = ", F_star, "W/m^2")

    #F_star_mJy = conversion_IS_to_mJy(F_star, dist_system, R_star)
    # F_star_mJy = flux_mJy(F_star, lambda_min, lambda_max, dist_system, R_star)
    # print("F_star = ", F_star_mJy, "mJy (seen from Earth)")

    F_star_obs_mJy = 2.528
    #F_star_obs_Wm2m = conversion_mJy_to_IS(F_star_obs_mJy, dist_system, R_star)
    #print("F_star_obs = ", F_star_obs_Wm2m, "W/m^2/m (seen from Earth)")
    F_star_obs_Wm2 = flux_Wm2(F_star_obs_mJy, lambda_min_F1500, lambda_max_F1500, dist_system, R_star)
    print("F_star_obs = ", F_star_obs_Wm2, "W/m^2 (seen from Earth)")

    # F_star_sphinx = flux_star_miri("F1500W")
    # print("F_star_sphinx = ", F_star_sphinx, "W/m^2")

    F_star_sphinx = flux_miri2("F1500W")
    print("F_star_sphinx = ", F_star_sphinx, "W/m^2")


    # For TRAPPIST-1 b

    # T_eq_b = T_eff_star * (R_star / (2 * a_b))**0.5
    T_eq_b = planet_equilibirium_temperature(T_eff_star, R_star, a_b)
    print("T_eq_b = ", T_eq_b, "K")

    F_ratio_b = flux_ratio(R_b, R_star, T_eff_star, a_b, lambda_min_F1500, lambda_max_F1500)
    print("F_ratio_b = ", F_ratio_b, "ppm")
    print("12.8 µm filter:", flux_ratio(R_b, R_star, T_eff_star, a_b, 11.2e-6, 14.2e-6), "ppm")

    # F_b = flux_black_body(lambda_min, lambda_max, T_eq_b)
    # print("F_b/F_star_obs", F_b/F_star_obs_Wm2*1e6, "ppm")


    # For TRAPPIST-1 c

    # T_eq_c = T_eff_star * (R_star / (2 * a_c))**0.5
    # print("T_eq_c = ", T_eq_c)

    # F_c = 2/3* flux_black_body(lambda_min,lambda_max, T_eq_c)
    # print("F_c = ", F_c)

    # print("F_c/F_star = ", F_c/F_star* 1e6, "ppm")

    # T_eq_c = planet_equilibirium_temperature(T_eff_star, R_star, a_c)
    # print("T_eq_c = ", T_eq_c, "K")
    # F_ratio_c = flux_ratio(R_c, R_star, T_eff_star, a_c, lambda_min, lambda_max)
    # print("F_ratio_c = ", F_ratio_c, "ppm")


if __name__ == "__main__":
    main()