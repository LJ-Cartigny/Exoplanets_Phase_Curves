#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# May 2025
# Fluxes with wavelength

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.constants import c, h, k, sigma
# from tqdm import tqdm
# import time

from Orbital_motion import compute_true_anomaly
from Transits import eclipse
from Phase_curve_v1 import star_planet_separation, flux_star, flux_planet, luminosity_planet_dayside, phase_curve
from TRAPPIST1_parameters import *

def flux_model_interp(l,model='sphinx'):
    """
    Interpolates the flux of TRAPPIST-1 at a given wavelength from the SPHINX or PHOENIX model.

    :param l: the wavelength (in m)
    :type l: float

    :param model: the model to use: 'sphinx' or 'phoenix' (default: 'sphinx')
    :type model: str

    :return: F
    :rtype: float
    """

    # Interpolate the flux over wavelength
    if model == 'sphinx':
        flux_interp = interp1d(wavelengths_T1_sphinx, flux_T1_sphinx, bounds_error=False, fill_value=0)
    
    elif model == 'phoenix':
        flux_interp = interp1d(wavelengths_T1_phoenix, flux_T1_phoenix_mJy, bounds_error=False, fill_value=0)
    
    else:
        raise ValueError("Invalid model. Choose 'sphinx' or 'phoenix'.")

    # Compute the flux at the given wavelength
    F = flux_interp(l)

    return F

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

    return F


def planet_equilibirium_temperature(T_star, R_star, d, albedo=0., redistribution=0.):
    """
    Determines the equilibrium temperature of the day side of a tidally locked planet (in K).

    :param T_star: the effective temperature of the star (in K)
    :type T_star: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :param d: the distance between the star and the planet (in m)
    :type d: float

    :param albedo: the albedo of the planet (default: 0)
    :type albedo: float

    :param redistribution: the redistribution efficiency between the day side and night side (default: 0)
    :type redistribution: float

    :return: T_eq
    :rtype: float
    """

    T_eq = T_star * (R_star / d)**0.5 * ((2/3 - 5/12 * redistribution) * (1-albedo))**0.25

    return T_eq


def flux_ratio_black_body(R_planet, R_star, T_star, d, lambda_min, lambda_max):
    """
    Determines the flux ratio between the planet and the star as black bodies(in ppm).

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
    Compute the flux of an object in mJy over a range of wavelengths.

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


def flux_mJy_array(F_array, lambda_vals, lambda_min, lambda_max, dist, R):
    """
    Compute the integrated flux of an object in mJy over a given wavelength range.

    :param F_array: Array of flux densities (in W/m^2/m)
    :type F_array: array-like

    :param lambda_vals: Corresponding wavelengths for F_array (in m)
    :type lambda_vals: array-like

    :param lambda_min: Minimum wavelength for integration (in m)
    :type lambda_min: float

    :param lambda_max: Maximum wavelength for integration (in m)
    :type lambda_max: float

    :param dist: Distance to the object (in m)
    :type dist: float

    :param R: Radius of the object (in m)
    :type R: float

    :return: Integrated flux in mJy
    :rtype: float
    """
    print(f"Interpolating between {lambda_min} and {lambda_max}")
    print(f"Wavelength range: {np.min(lambda_vals)} – {np.max(lambda_vals)}")

    # Interpolate the flux over wavelength
    flux_interp = interp1d(lambda_vals, F_array, bounds_error=False, fill_value=0)

    # Define scalar integrand using the interpolated flux
    def integrand(l):
        return conversion_IS_to_mJy(flux_interp(l), l, dist, R)

    # Perform the integration
    F_mJy = quad(integrand, lambda_min, lambda_max)[0]

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
    Compute the flux of an object in W/m^2 over a range of wavelengths.

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

    else:
        raise ValueError("Invalid filter name. Choose from: F560W, F770W, F1000W, F1130W, F1280W, F1500W, F1800W, F2100W, F2550W.")

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

    QE = np.interp(wavelength, wavelengths, filter_band[1,:])/np.max(filter_band[1,:])

    return QE


def flux_star_miri(filter_name):
    """
    Returns the flux of the star TRAPPIST-1 in the specified MIRI filter band using the SPHINX model.

    :param filter_name: the name of the filter
    :type filter_name: str

    :return: F_star
    :rtype: float
    """

    filter_band = filter(filter_name)
    wavelengths_filter = filter_band[0,:]*1e-6 # Convert to m
    QE = filter_band[1,:]/np.max(filter_band[1,:])
    
    # The two arrays are not covering the same wavelength range so we need to determine the common range
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


def flux_planet_miri(filter_name, T_planet):
    """
    Returns the flux of the planet in the specified MIRI filter band.

    :param filter_name: the name of the filter
    :type filter_name: str

    :param T_planet: the temperature of the planet (in K)
    :type T_planet: float

    :return: F_planet_miri
    :rtype: float
    """

    filter_band = filter(filter_name)
    wavelengths_filter = filter_band[0,:]*1e-6 # Convert to m
    QE = filter_band[1,:]/np.max(filter_band[1,:])

    F_planet = np.pi*Planck_law(wavelengths_filter, T_planet) # in W/m^2/m

    F_planet_miri = np.trapz(F_planet*QE, wavelengths_filter)

    return F_planet_miri


def flux_ratio_miri(filter_name, R_planet, R_star, T_planet):
    """
    Returns the flux ratio between the planet and the star in the specified MIRI filter band (in ppm).

    :param filter_name: the name of the filter
    :type filter_name: str

    :param R_planet: the radius of the planet (in m)
    :type R_planet: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :param T_planet: the temperature of the planet (in K)
    :type T_planet: float

    :return: F_ratio_miri
    :rtype: float
    """

    F_star = flux_star_miri(filter_name)
    F_planet = flux_planet_miri(filter_name, T_planet)

    F_ratio_miri = (R_planet/R_star)**2 * F_planet / F_star * 1e6

    return F_ratio_miri


def integrate_flux_model_mJy(filter_name,model='sphinx'):
    """
    Integrates the flux (in mJy) of the SPHINX or PHOENIX model over the specified MIRI filter band.

    :param filter_name: the name of the filter
    :type filter_name: str

    :param model: the model to use: 'sphinx' or 'phoenix' (default: 'sphinx')
    :type model: str

    :return: F_miri
    :rtype: float
    """

    filter_band = filter(filter_name)
    
    spectrum_filter = flux_model_interp(filter_band[0,:]*1e-6,model=model)

    # Convert the flux density to mJy if necessary
    if model == 'sphinx':
        spectrum_filter_mJy = conversion_IS_to_mJy(spectrum_filter, filter_band[0,:]*1e-6, dist_system, R_star)
        spectrum_filter = spectrum_filter_mJy

    F_miri = 0
    norm_filter = 0

    for i in range(1, len(filter_band[0,:]), 1):
        norm_filter += filter_band[1,i] * (1./filter_band[0,i] - 1./filter_band[0,i-1])
        F_miri += spectrum_filter[i] * filter_band[1,i] * (1./filter_band[0,i] - 1./filter_band[0,i-1])

    # Normalize the flux
    F_miri /= norm_filter

    # Applying correction for the PHOENIX model
    if model == 'phoenix':
        if filter_name=='F1280W':
            F_miri /= 1.061
        if filter_name=='F1500W':
            F_miri /= 1.023
    
    return F_miri



def main():

    # Convert the SPHINX spectrum to mJy

    flux_T1_sphinx_mJy = conversion_IS_to_mJy(flux_T1_sphinx, wavelengths_T1_sphinx, dist_system, R_star)

    # Comparison of the SPHINX, PHOENIX and black body spectra of TRAPPIST-1

    flux_Planck = Planck_law(wavelengths_T1_sphinx, T_eff_star)
    flux_Planck_mJy = conversion_IS_to_mJy(flux_Planck, wavelengths_T1_sphinx, dist_system, R_star)

    # plt.figure(figsize=(16,9))
    # plt.plot(wavelengths_T1_sphinx*1e6, flux_T1_sphinx_mJy, label="SPHINX", alpha=0.4, linewidth=0.5)
    # plt.plot(wavelengths_T1_phoenix*1e6, flux_T1_phoenix_mJy, label='PHOENIX', alpha=0.4, linewidth=0.5)
    # plt.plot(wavelengths_T1_sphinx*1e6, flux_Planck_mJy, color='black', label="Black body")
    # plt.xlabel(r"Wavelength ($\mu m$)",size=15)
    # plt.ylabel(r"Flux ($mJy$)",size=15)
    # plt.xlim((10,20))
    # plt.ylim(0,7.5)
    # plt.title("Comparison of the SPHINX, PHOENIX and black body spectra of TRAPPIST-1",size=15)
    # plt.legend(fontsize=13)
    # plt.xticks(fontsize=13)
    # plt.yticks(fontsize=13)
    # plt.grid()
    # plt.savefig("Comparison_spectra.png", bbox_inches='tight')
    # plt.show()
    
    # Define MIRI 15 micron filter band

    # lambda_min_F1500 = 13e-6 # in m
    # lambda_max_F1500 = 17.3e-6 # in m

    l_eff_F1280 = 12.62
    l_eff_F1500 = 14.79
    flux_measured_12 = np.mean([3.425,3.428, 3.427, 3.428]) # From Ducrot et al. 2024
    flux_err_12 = 0.12
    flux_measured_15 = 2.589 # From Greene et al. 2023
    flux_err_15 = 0.078

    l = np.linspace(1, 20, 10000)*1e-6  # Wavelength range from 1 to 20 microns


    plt.figure(figsize=(16,9))
    plt.plot(wavelengths_T1_sphinx*1e6,flux_T1_sphinx_mJy, color='blue', alpha=0.4, label="SPHINX model")
    plt.plot(wavelengths_T1_phoenix*1e6, flux_T1_phoenix_mJy, color='red', alpha=0.4, label="PHOENIX model")
    # plt.scatter(l_eff_F1500, flux_measured_15, marker='*', color='green', label="Measured flux", zorder=5)
    plt.errorbar(l_eff_F1500, flux_measured_15, yerr=flux_err_15, fmt='.', color='green', markersize=5, elinewidth=2, capsize=5, label="F1500W measured flux", zorder=5)
    # plt.scatter(l_eff_F1280, flux_measured_12, marker='*', color='orange', label="Measured flux", zorder=5)
    plt.errorbar(l_eff_F1280, flux_measured_12, yerr=flux_err_12, fmt='.', color='orange', markersize=5, elinewidth=2, capsize=5, label="F1280W measured flux", zorder=5)
    plt.scatter(l_eff_F1500,integrate_flux_model_mJy('F1500W'),color='green',marker='s',label='F1500W simulation (SPHINX)', zorder=5)
    plt.scatter(l_eff_F1280,integrate_flux_model_mJy('F1280W'),color='orange',marker='s',label='F1280W simulation (SPHINX)', zorder=5)
    plt.scatter(l_eff_F1500,integrate_flux_model_mJy('F1500W',model='phoenix'),marker='x',color='green',label='F1500W simulation (PHOENIX)', zorder=5)
    plt.scatter(l_eff_F1280,integrate_flux_model_mJy('F1280W',model='phoenix'),marker='x',color='orange',label='F1280W simulation (PHOENIX)', zorder=5)
    plt.plot(l*1e6,quantum_efficiency('F1500W',l), color='green')
    plt.fill_between(l*1e6, quantum_efficiency('F1500W', l), color='green', alpha=0.2,label='F1500W filter')  # Fill under F1500W filter
    plt.plot(l*1e6,quantum_efficiency('F1280W',l), color ='orange')
    plt.fill_between(l*1e6, quantum_efficiency('F1280W', l), color='orange', alpha=0.2, label='F1280W filter')  # Fill under F1280W filter
    plt.xlabel(r"Wavelength ($\mu m$)", size=15)
    plt.ylabel(r"Flux ($mJy$)", size=15)
    plt.xlim(10,20)
    plt.ylim(0,4.5)
    plt.title("Comparisons of simulated and observed TRAPPIST-1 stellar fluxes", size=15)
    plt.legend(ncols=2, fontsize=13)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.grid()
    # plt.savefig("Flux_TRAPPIST1_mJy.png", bbox_inches='tight')
    plt.show()


    # Flux of star TRAPPIST-1 in the MIRI 15 micron filter band

    print("For star TRAPPIST-1:")

    # F_star_Planck = flux_black_body(lambda_min_F1500, lambda_max_F1500, T_eff_star)
    # print("F_star = ", F_star_Planck, "W/m^2 (as a black body)")

    #F_star_obs_mJy = 2.528
    print("F_star_obs_F1500_mJy = ", flux_measured_15, "mJy (observed value with MIRI F1500W filter)")

    F_star_sphinx_miri_F1500_mJy = integrate_flux_model_mJy("F1500W")
    print("F_star_sphinx_miri = ", F_star_sphinx_miri_F1500_mJy, "mJy (using the corrected SPHINX spectrum with MIRI F1500W filter)")
    print("Difference in % between SPHINX and observed flux: ", (F_star_sphinx_miri_F1500_mJy - flux_measured_15) / flux_measured_15 * 100, "%")

    F_star_phoenix_miri_F1500_mJy = integrate_flux_model_mJy("F1500W", model='phoenix')
    print("F_star_phoenix_miri = ", F_star_phoenix_miri_F1500_mJy, "mJy (using the PHOENIX spectrum with MIRI F1500W filter)")
    print("Difference in % between PHOENIX and observed flux: ", (F_star_phoenix_miri_F1500_mJy - flux_measured_15) / flux_measured_15 * 100, "%")
    

    print("\nF_star_obs_F1280_mJy = ", flux_measured_12, "mJy (observed value with MIRI F1280W filter)")
    # F_star_obs_Wm2 = flux_Wm2(F_star_obs_mJy, lambda_min_F1500, lambda_max_F1500, dist_system, R_star)
    # print("F_star_obs = ", F_star_obs_Wm2, "W/m^2 (seen from Earth)")
    # print("F_star_obs_mJy = ", flux_mJy(F_star_obs_Wm2, lambda_min_F1500, lambda_max_F1500, dist_system, R_star), "mJy (seen from Earth)")
    
    F_star_sphinx_miri_F1280_mJy = integrate_flux_model_mJy("F1280W")
    print("F_star_sphinx_miri = ", F_star_sphinx_miri_F1280_mJy, "mJy (using the corrected SPHINX spectrum with MIRI F1280W filter)")
    print("Difference in % between SPHINX and observed flux: ", (F_star_sphinx_miri_F1280_mJy - flux_measured_12) / flux_measured_12 * 100, "%")

    F_star_phoenix_miri_F1280_mJy = integrate_flux_model_mJy("F1280W", model='phoenix')
    print("F_star_phoenix_miri = ", F_star_phoenix_miri_F1280_mJy, "mJy (using the PHOENIX spectrum with MIRI F1280W filter)")
    print("Difference in % between PHOENIX and observed flux: ", (F_star_phoenix_miri_F1280_mJy - flux_measured_12) / flux_measured_12 * 100, "%")


    # # For TRAPPIST-1 b
    # print("\nFor planet TRAPPIST-1 b:")

    # T_eq_b = planet_equilibirium_temperature(T_eff_star, R_star, a_b)
    # print("T_eq_b = ", T_eq_b, "K")

    # F_b_miri = flux_planet_miri("F1500W", T_eq_b)
    # print("F_b = ", F_b_miri, "W/m^2 (F1500W MIRI filter)")

    # flux_ratio_miri_b = flux_ratio_miri("F1500W", R_b, R_star, T_eq_b)
    # print("F_ratio_b = ", flux_ratio_miri_b, "ppm (F1500W MIRI filter)")

    # F_b_miri_bis = flux_planet_miri("F1280W", T_eq_b)
    # print("F_b_bis = ", F_b_miri_bis, "W/m^2 (F1280W MIRI filter)")
    # flux_ratio_miri_b_bis = flux_ratio_miri("F1280W", R_b, R_star, T_eq_b)
    # print("F_ratio_b_bis = ", flux_ratio_miri_b_bis, "ppm (F1280W MIRI filter)")


    # # For TRAPPIST-1 c

    # print("\nFor planet TRAPPIST-1 c:")

    # T_eq_c = planet_equilibirium_temperature(T_eff_star, R_star, a_c)
    # print("T_eq_c = ", T_eq_c, "K")

    # F_c_miri = flux_planet_miri("F1500W", T_eq_c)
    # print("F_c = ", F_c_miri, "W/m^2 (F1500W MIRI filter)")

    # flux_ratio_miri_c = flux_ratio_miri("F1500W", R_c, R_star, T_eq_c)
    # print("F_ratio_c = ", flux_ratio_miri_c, "ppm (F1500W MIRI filter)")


    # # For TRAPPIST-1 d

    # print("\nFor planet TRAPPIST-1 d (without atmosphere):")

    # T_eq_d = planet_equilibirium_temperature(T_eff_star, R_star, a_d)
    # print("T_eq_d = ", T_eq_d, "K")
    # F_d_miri = flux_planet_miri("F1500W", T_eq_d)
    # print("F_d = ", F_d_miri, "W/m^2 (F1500W MIRI filter)")

    # flux_ratio_miri_d = flux_ratio_miri("F1500W", R_d, R_star, T_eq_d)
    # print("F_ratio_d = ", flux_ratio_miri_d, "ppm (F1500W MIRI filter)")

    # print("\nFor planet TRAPPIST-1 d (with atmosphere):")
    # T_eq_d_atm = planet_equilibirium_temperature(T_eff_star, R_star, a_d, redistribution = 1)
    # print("T_eq_d_atm = ", T_eq_d_atm, "K")
    # F_d_miri_atm = flux_planet_miri("F1500W", T_eq_d_atm)
    # print("F_d_atm = ", F_d_miri_atm, "W/m^2 (F1500W MIRI filter)")

    # flux_ratio_miri_d_atm = flux_ratio_miri("F1500W", R_d, R_star, T_eq_d_atm)
    # print("F_ratio_d_atm = ", flux_ratio_miri_d_atm, "ppm (F1500W MIRI filter)")


if __name__ == "__main__":
    main()