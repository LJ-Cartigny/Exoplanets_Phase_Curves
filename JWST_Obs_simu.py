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

nb_points = 10000
planets='defgh'
redistribution = 0
Keplerian = True
filter = 'F1280W'
model = 'phoenix'
unit = 'mJy'

print("Simulating the phase curves during the JWST visits...")
for i in tqdm(range(len(t_start))):
    phase_curve_simulation(t_start[i], nb_days[i], nb_points=nb_points, planets=planets, redistribution=redistribution, filter=filter, model=model, unit=unit, Keplerian=Keplerian, plot=False,save_plot=True,save_txt=True)
    time.sleep(0.1)
