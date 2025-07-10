#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# April 2025
# Simulations of the phase curves during JWST observations of TRAPPIST-1

import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from tqdm import tqdm
import time

from Phase_curve_TTV import phase_curve_simulation

def phase_curve_visit(planets, redistribution, filter, model, unit, nb_points=10000, Keplerian=True):
    """
    Simulates the phase curves of the TRAPPIST-1 planets during JWST visits.

    :param planets: the planets to simulate
    :type planets: str

    :param redistribution: the redistribution efficiency between the day side and night side (default: 0)
    :type redistribution: float

    :param filter: the MIRI filter to use
    :type filter: str

    :param model: the model to use for the stellar flux. If 'sphinx', the flux is computed using the SPHINX model. If 'phoenix', the flux is computed using the PHOENIX model.
    :type model: str

    :param unit: the unit of the phase curve. If 'ppm', the fluxes of the planets will be computed relatively to the stellar flux in ppm. If 'mJy', the planetary fluxes will be computed in absolute value in "mJy". Will be set automatically to 'mJy' if the model is 'phoenix'.
    :type unit: str

    :param nb_points: the number of points for the phase curves (default: 10000)
    :type nb_points: int

    :param Keplerian: whether to use the Keplerian periods or not (default: True)
    :type Keplerian: bool

    :rtype: None
    """

    program_ID, visit, t_start, t_end = np.loadtxt("JWST_Obs_times.txt", delimiter=',', skiprows=2, usecols=(0,1,2,3), unpack=True,dtype=str)

    t_start = Time(t_start, format='isot', scale='tdb')
    t_end = Time(t_end, format='isot', scale='tdb')

    t_start.format = 'jd'
    t_end.format = 'jd'

    nb_days = t_end - t_start

    t_start = t_start.jd
    t_end = t_end.jd
    nb_days = nb_days.jd

    t_start -= 2450000
    t_end -= 2450000

    print("Simulating the phase curves during the JWST visits...")
    for i in tqdm(range(len(t_start))):
        phase_curve_simulation(t_start[i], nb_days[i], nb_points=nb_points, planets=planets, redistribution=redistribution, filter=filter, model=model, unit=unit, Keplerian=Keplerian, plot=False,save_plot=True,save_txt=True)
        time.sleep(0.1)
