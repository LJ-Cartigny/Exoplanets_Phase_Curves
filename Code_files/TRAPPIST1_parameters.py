#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# February 2025
# Parameters for the TRAPPIST-1 system

"""
This module contains some parameters for the planets of the TRAPPIST-1 system.
"""

from Solar_System_constants import *
import numpy as np
from pathlib import Path
DATA_DIR = Path(__file__).resolve().parent


# Distance between the TRAPPIST-1 system and the Solar system

dist_system = 12.429888806540756 * 3.086e16  #: Distance between the TRAPPIST-1 system and the Solar system in meters (from NASA Exoplanet Archive)

# For star TRAPPIST-1 (using NASA Exoplanet Archive)

T_eff_star = 2566  #: Effective temperature of star TRAPPIST-1 in Kelvin (Agol et al. 2021)
R_star = 	0.1192 * R_Sun  #: Radius of star TRAPPIST-1 in meters (Agol et al. 2021)
M_star = 0.0898 * M_Sun  #: Mass of star TRAPPIST-1 in kilograms (Ducrot et al. 2020)
L_star = 10**(-3.26)*L_Sun  #: Luminosity of star TRAPPIST-1 in Watts (Ducrot et al. 2020)

# Load model spectra from files located next to this module. Fall back to empty arrays if unavailable (e.g. during Sphinx build).
try:
    wavelengths_T1_sphinx, flux_T1_sphinx = np.loadtxt(DATA_DIR / "sphinx_spectrum_T-1_aisha.txt", unpack=True, skiprows=1)  #: Wavelengths (m) and flux from SPHINX model spectrum of star TRAPPIST-1
    wavelengths_T1_sphinx = wavelengths_T1_sphinx * 1e-6  #: Wavelengths from the SPHINX model spectrum of star TRAPPIST-1 in meters
    flux_T1_sphinx = flux_T1_sphinx / 1.07  #: Stellar flux of TRAPPIST-1 in W/m^2/m from the SPHINX spectrum model (corrected by 7% to better match observations)
except Exception:
    wavelengths_T1_sphinx = np.array([])  #: Wavelengths (m) from SPHINX model spectrum of star TRAPPIST-1 (empty during docs build)
    flux_T1_sphinx = np.array([])  #: Flux from SPHINX model spectrum of star TRAPPIST-1 (empty during docs build)

try:
    wavelengths_T1_phoenix, flux_T1_phoenix_mJy = np.loadtxt(DATA_DIR / "TRAPPIST1_Phoenix_model.txt", unpack=True, skiprows=1)  #: Stellar flux of TRAPPIST-1 in mJy from the PHOENIX spectrum model
    wavelengths_T1_phoenix = wavelengths_T1_phoenix * 1e-6  #: Wavelengths from the PHOENIX model spectrum of star TRAPPIST-1 in meters
except Exception:
    wavelengths_T1_phoenix = np.array([])  #: Wavelengths (m) from Phoenix model spectrum of star TRAPPIST-1 (empty during docs build)
    flux_T1_phoenix_mJy = np.array([])  #: Flux (mJy) from Phoenix model spectrum of star TRAPPIST-1 (empty during docs build)


#For TRAPPIST-1 b (using NASA Exoplanet Archive)

a_b = 20.13 * R_star  #: Semi-major axis of TRAPPIST-1 b in meters (Ducrot et al. 2020)
P_b = 1.51088432  #: Orbital period of TRAPPIST-1 b in days (Ducrot et al. 2020)
i_b = 89.28 * np.pi/180  #: Orbital inclination of TRAPPIST-1 b in radians (Ducrot et al. 2020)
#i_b = 89.73 * np.pi/180 # in rad (from Agol et al. 2021)
omega_b = 336.86 * np.pi/180  #: Argument of periastron of TRAPPIST-1 b in radians (Grimm et al. 2018)
e_b = 0.00622  #: Orbital eccentricity of TRAPPIST-1 b (Grimm et al. 2018)
R_b = 1.116 * R_Earth  #: Radius of TRAPPIST-1 b in meters (Agol et al. 2021)


#For TRAPPIST-1 c (using NASA Exoplanet Archive)

a_c = 27.57 * R_star  #: Semi-major axis of TRAPPIST-1 c in meters (Ducrot et al. 2020)
P_c = 2.42179346  #: Orbital period of TRAPPIST-1 c in days (Ducrot et al. 2020)
i_c = np.radians(89.47)  #: Orbital inclination of TRAPPIST-1 c in radians (Ducrot et al. 2020)
omega_c = np.radians(282.45)  #: Argument of periastron of TRAPPIST-1 c in radians (Grimm et al. 2018)
e_c = 0.00654  #: Orbital eccentricity of TRAPPIST-1 c (Grimm et al. 2018)
R_c = 1.097 * R_Earth  #: Radius of TRAPPIST-1 c in meters (Agol et al. 2021)


#For TRAPPIST-1 d (using NASA Exoplanet Archive)

a_d = 38.85 * R_star  #: Semi-major axis of TRAPPIST-1 d in meters (Ducrot et al. 2020)
P_d = 4.04978035  #: Orbital period of TRAPPIST-1 d in days (Ducrot et al. 2020)
i_d = np.radians(89.65)  #: Orbital inclination of TRAPPIST-1 d in radians (Ducrot et al. 2020)
omega_d = np.radians(-8.73)  #: Argument of periastron of TRAPPIST-1 d in radians (Grimm et al. 2018)
e_d = 0.00837  #: Orbital eccentricity of TRAPPIST-1 d (Grimm et al. 2018)
R_d = 0.788 * R_Earth  #: Radius of TRAPPIST-1 d in meters (Agol et al. 2021)


#For TRAPPIST-1 e (using NASA Exoplanet Archive)

a_e = 51.0 * R_star  #: Semi-major axis of TRAPPIST-1 e in meters (Ducrot et al. 2020)
P_e = 6.09956479  #: Orbital period of TRAPPIST-1 e in days (Ducrot et al. 2020)
i_e = np.radians(89.663)  #: Orbital inclination of TRAPPIST-1 e in radians (Ducrot et al. 2020)
omega_e = np.radians(108.37)  #: Argument of periastron of TRAPPIST-1 e in radians (Grimm et al. 2018)
e_e = 0.000510  #: Orbital eccentricity of TRAPPIST-1 e (Grimm et al. 2018)
R_e = 0.920 * R_Earth  #: Radius of TRAPPIST-1 e in meters (Agol et al. 2021)


#For TRAPPIST-1 f (using NASA Exoplanet Archive)

a_f = 67.10 * R_star  #: Semi-major axis of TRAPPIST-1 f in meters (Ducrot et al. 2020)
P_f = 9.20659399  #: Orbital period of TRAPPIST-1 f in days (Ducrot et al. 2020)
i_f = np.radians(89.666)  #: Orbital inclination of TRAPPIST-1 f in radians (Ducrot et al. 2020)
omega_f = np.radians(368.81)  #: Argument of periastron of TRAPPIST-1 f in radians (Grimm et al. 2018)
e_f = 0.01007  #: Orbital eccentricity of TRAPPIST-1 f (Grimm et al. 2018)
R_f = 1.045 * R_Earth  #: Radius of TRAPPIST-1 f in meters (Agol et al. 2021)


#For TRAPPIST-1 g (using NASA Exoplanet Archive)

a_g = 81.7 * R_star  #: Semi-major axis of TRAPPIST-1 g in meters (Ducrot et al. 2020)
P_g = 12.35355570  #: Orbital period of TRAPPIST-1 g in days (Ducrot et al. 2020)
i_g = np.radians(89.698)  #: Orbital inclination of TRAPPIST-1 g in radians (Ducrot et al. 2020)
omega_g = np.radians(191.34)  #: Argument of periastron of TRAPPIST-1 g in radians (Grimm et al. 2018)
e_g = 0.00208  #: Orbital eccentricity of TRAPPIST-1 g (Grimm et al. 2018)
R_g = 1.129 * R_Earth  #: Radius of TRAPPIST-1 g in meters (Agol et al. 2021)


#For TRAPPIST-1 h (using NASA Exoplanet Archive)

a_h = 107.9 * R_star  #: Semi-major axis of TRAPPIST-1 h in meters (Ducrot et al. 2020)
P_h = 18.76727450  #: Orbital period of TRAPPIST-1 h in days (Ducrot et al. 2020)
i_h = np.radians(89.763)  #: Orbital inclination of TRAPPIST-1 h in radians (Ducrot et al. 2020)
omega_h = np.radians(338.92)  #: Argument of periastron of TRAPPIST-1 h in radians (Grimm et al. 2018)
e_h = 0.00567  #: Orbital eccentricity of TRAPPIST-1 h (Grimm et al. 2018)
R_h = 0.755 * R_Earth  #: Radius of TRAPPIST-1 h in meters (Agol et al. 2021)


__all__ = [
    "dist_system",
    "T_eff_star", "R_star", "M_star", "L_star",
    "wavelengths_T1_sphinx", "flux_T1_sphinx",
    "wavelengths_T1_phoenix", "flux_T1_phoenix_mJy",
    "a_b", "P_b", "i_b", "omega_b", "e_b", "R_b",
    "a_c", "P_c", "i_c", "omega_c", "e_c", "R_c",
    "a_d", "P_d", "i_d", "omega_d", "e_d", "R_d",
    "a_e", "P_e", "i_e", "omega_e", "e_e", "R_e",
    "a_f", "P_f", "i_f", "omega_f", "e_f", "R_f",
    "a_g", "P_g", "i_g", "omega_g", "e_g", "R_g",
    "a_h", "P_h", "i_h", "omega_h", "e_h", "R_h",
]