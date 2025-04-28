#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# April 2025
# Simulations of the phase curves during JWST observations of TRAPPIST-1b

import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

from Phase_curve_TTV import phase_curve_simulation

program_ID, visit, t_start, t_end = np.loadtxt("JWST_Obs_times.txt", delimiter=',', skiprows=2, unpack=True,dtype=str)

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

for i in range(len(t_start)):
    phase_curve_simulation(t_start[i], nb_days[i],save_plot=True,save_txt=True)
