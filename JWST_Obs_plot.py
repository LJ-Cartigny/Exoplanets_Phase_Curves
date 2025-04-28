#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# April 2025
# Placement of the JWST observations on the phase curves

import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

from Phase_curve_TTV import phase_curve_simulation

program_ID, visit, t_start, t_end = np.loadtxt("JWST_Obs_times.txt", delimiter=',', skiprows=2, unpack=True,dtype=str)

t_start = Time(t_start, format='isot', scale='tdb')
t_end = Time(t_end, format='isot', scale='tdb')

t_start.format = 'jd'
t_end.format = 'jd'

t_start = t_start.jd
t_end = t_end.jd

t_start -= 2450000
t_end -= 2450000

nb_days = t_end[-1] - t_start[0] + 1
t0 = t_start[0]-0.5

# phase_curve_simulation(t0, nb_days,plot=False,save_plot=True,save_txt=True) # Uncomment if the simulation is not done yet


# Plot

plt.figure(figsize=(16, 9))

for p in "bcdefgh":  # Comment to not plot the individual planets
    t_simu, phase_simu = np.loadtxt("Phase_curve_TTV_output/phase_curve_"+p+"_TTV_"+str(t0)+".txt", delimiter=",", skiprows=1, unpack=True)
    plt.plot(t_simu, phase_simu, '--', label=p, linewidth=0.5)

t_total_simu, phase_curve_total_simu = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_TTV_"+str(t0)+".txt", delimiter=",",skiprows=1, unpack=True)

plt.plot(t_total_simu, phase_curve_total_simu, '--', color='grey', label="Total flux")


for i in range(len(t_start)):
    t0 = t_start[i]
    t_visit, phase_curve_total_visit = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_TTV_"+str(t0)+".txt", delimiter=',',skiprows=1, unpack=True)
    plt.plot(t_visit, phase_curve_total_visit, label=program_ID[i]+' visit '+visit[i] , linewidth=3)

plt.xlabel(r"Time ($BJD_{TBD} - 2450000$)")
plt.ylabel(r"$L_{planet}/L_{star}$ (ppm)")
plt.title("JWST Observations over the phase curves of TRAPPIST-1")
plt.legend()
plt.grid()
# plt.savefig("JWST_Obs_plots/JWST_Obs_phase_curves_GTO1177.png", dpi=300, bbox_inches='tight')
plt.show()