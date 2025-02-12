#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# February 2025
# Phase curves with an eyeball planet

import numpy as np
import matplotlib.pyplot as plt
from Orbital_motion import compute_true_anomaly
from Phase_curve_v1 import phase_angle, phase_function, star_planet_separation, surface_sphere, phase_curve, flux_star, flux_planet


def flux_star_zenith(L,d,theta):
    """
    Determines the flux received from a star (in W/m^2) at a distance d with a zenith angle theta.

    :param L: the star luminosity (in W)
    :type L: float

    :param d: the distance (in m)
    :type d: float

    :param theta: the zenith angle (in rad)
    :type theta: float

    :return: F
    :rtype: float
    """

    F = flux_star(L,d)*np.cos(theta)
    return F