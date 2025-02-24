#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# February 2025
# Exoplanet transit

import numpy as np
import matplotlib.pyplot as plt
from TRAPPIST1_parameters import *

def transit_depth(R_planet, R_star, a, i):
    """
    Determines the depth of an exoplanet transit.

    :param R_planet: the radius of the planet (in m)
    :type Rp: float

    :param R_star: the radius of the star (in m)
    :type R_star: float

    :param a: the semimajor axis (in m)
    :type a: float

    :param i: the inclination (in rad)
    :type i: float

    :return: d
    :rtype: float
    """

    d = (R_planet/R_star)**2 * (a/R_star) * np.cos(i)
    return d