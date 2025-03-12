#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# February 2025
# Parameters for the TRAPPIST-1 system

from Solar_System_constants import *
import numpy as np

# For star TRAPPIST-1 (using NASA Exoplanet Archive)

T_eff_star = 2520 # in K (from Ducrot et al. 2020)
R_star = 0.1234 * R_Sun # in m (from Ducrot et al. 2020)
M_star = 0.0898 * M_Sun # in kg (from Ducrot et al. 2020)
L_star = 10**(-3.26)*L_Sun # in W (from Ducrot et al. 2020)


#For TRAPPIST-1 b (using NASA Exoplanet Archive)

a_b = 20.13 * R_star # in m (from Ducrot et al. 2020)
P_b = 1.51088432 # in days (from Ducrot et al. 2020)
i_b = 89.28 * np.pi/180 # in rad (from Ducrot et al. 2020)
#i_b = 89.73 * np.pi/180 # in rad (from Agol et al. 2021)
omega_b = 336.86 * np.pi/180 # in rad (from Grimm et al. 2018)
e_b = 0.00622 # (from Grimm et al. 2018)
R_b = 1.017 * R_Earth # in m (from Grimm et al. 2018)


#For TRAPPIST-1 c (using NASA Exoplanet Archive)

a_c = 27.57 * R_star # in m (from Ducrot et al. 2020)
P_c = 2.42179346 # in days (from Ducrot et al. 2020)
i_c = np.radians(89.47) # in rad (from Ducrot et al. 2020)
omega_c = np.radians(282.45) # in rad (from Grimm et al. 2018)
e_c = 0.00654 # (from Grimm et al. 2018)
R_c = 1.095 * R_Earth # in m (from Grimm et al. 2018)


#For TRAPPIST-1 d (using NASA Exoplanet Archive)

a_d = 38.85 * R_star # in m (from Ducrot et al. 2020)
P_d = 4.04978035 # in days (from Ducrot et al. 2020)
i_d = np.radians(89.65) # in rad (from Ducrot et al. 2020)
omega_d = np.radians(-8.73) # in rad (from Grimm et al. 2018)
e_d = 0.00837 # (from Grimm et al. 2018)
R_d = 0.784 * R_Earth # in m (from Grimm et al. 2018)


#For TRAPPIST-1 e (using NASA Exoplanet Archive)

a_e = 51.0 * R_star # in m (from Ducrot et al. 2020)
P_e = 6.09956479 # in days (from Ducrot et al. 2020)
i_e = np.radians(89.663) # in rad (from Ducrot et al. 2020)
omega_e = np.radians(108.37) # in rad (from Grimm et al. 2018)
e_e = 0.000510 # (from Grimm et al. 2018)
R_e = 0.910 * R_Earth # in m (from Grimm et al. 2018)


#For TRAPPIST-1 f (using NASA Exoplanet Archive)

a_f = 67.10 * R_star # in m (from Ducrot et al. 2020)
P_f = 9.20659399 # in days (from Ducrot et al. 2020)
i_f = np.radians(89.666) # in rad (from Ducrot et al. 2020)
omega_f = np.radians(368.81) # in rad (from Grimm et al. 2018)
e_f = 0.01007 # (from Grimm et al. 2018)
R_f = 1.046 * R_Earth # in m (from Grimm et al. 2018)


#For TRAPPIST-1 g (using NASA Exoplanet Archive)

a_g = 81.7 * R_star # in m (from Ducrot et al. 2020)
P_g = 12.35355570 # in days (from Ducrot et al. 2020)
i_g = np.radians(89.698) # in rad (from Ducrot et al. 2020)
omega_g = np.radians(191.34) # in rad (from Grimm et al. 2018)
e_g = 0.00208 # (from Grimm et al. 2018)
R_g = 1.148 * R_Earth # in m (from Grimm et al. 2018)


#For TRAPPIST-1 h (using NASA Exoplanet Archive)

a_h = 107.9 * R_star # in m (from Ducrot et al. 2020)
P_h = 18.76727450 # in days (from Ducrot et al. 2020)
i_h = np.radians(89.763) # in rad (from Ducrot et al. 2020)
omega_h = np.radians(338.92) # in rad (from Grimm et al. 2018)
e_h = 0.00567 # (from Grimm et al. 2018)
R_h = 0.773 * R_Earth # in m (from Grimm et al. 2018)