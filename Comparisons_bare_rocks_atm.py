#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# May 2025
# Comparison of TRAPPIST-1 phase curves with and without thick atmospheres


import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec

from Phase_curve_TTV import phase_curve_simulation

# Load the data

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

nb_days = np.max(t_end) - np.min(t_start) + 3
t0 = np.min(t_start)-1.5

print("t0 = ", t0)
print("t_end = ", t0+nb_days)
print("nb_days = ", nb_days)

# Simulation parameters

nb_points = 100000
Keplerian = True
planets = 'defgh'
filter = 'F1500W'
save_plots = False # Write True if you want to save the plots
do_simulation = False # Write True if the simulation hasn't been done yet
plot_individual_planets = True # Write True if you want to plot the individual planets as bare rocks to see their phases


if do_simulation:
    phase_curve_simulation(t0, nb_days, nb_points=nb_points, planets=planets, redistribution=1, filter=filter, Keplerian=Keplerian, plot=False,save_plot=True,save_txt=True)
    phase_curve_simulation(t0, nb_days, nb_points=nb_points, planets=planets, redistribution=0, filter=filter, Keplerian=Keplerian, plot=False,save_plot=True,save_txt=True)


# Plot

fig_comparison = plt.figure(figsize=(32, 18))

colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink']
j = 0

if plot_individual_planets:
    for p in planets:
        t_simu, phase_simu = np.loadtxt("Phase_curve_TTV_output/phase_curve_"+p+"_"+filter+"_"+str(t0)+".txt", delimiter=",", skiprows=1, unpack=True)
        # t_simu_atm, phase_simu_atm = np.loadtxt("Phase_curve_TTV_output/phase_curve_"+p+"_atm_"+filter+"_"+str(t0)+".txt", delimiter=",", skiprows=1, unpack=True)

        plt.plot(t_simu, phase_simu, label=p+" (bare rock)", linewidth=0.5)
        # plt.plot(t_simu_atm, phase_simu_atm, '--', color=colors[j], label=p+" (atmosphere)", linewidth=0.5)

        j += 1

t_total_simu, phase_curve_total_simu = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_"+filter+"_"+str(t0)+".txt", delimiter=",",skiprows=1, unpack=True)
t_total_simu_atm, phase_curve_total_simu_atm = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_atm_"+filter+"_"+str(t0)+".txt", delimiter=",",skiprows=1, unpack=True)

plt.plot(t_total_simu, phase_curve_total_simu, '--', color='grey', label="Total (bare rock)")
plt.plot(t_total_simu_atm, phase_curve_total_simu_atm, color='black', label="Total (atmospheres)")


j = -1

for i in range(len(t_start)):
    t0_visit= t_start[i]

    t_visit, phase_curve_total_visit = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_"+filter+"_"+str(t0_visit)+".txt", delimiter=',',skiprows=1, unpack=True)

    if i == 0 or program_ID[i] != program_ID[i-1]:
        j+=1
        plt.plot(t_visit, phase_curve_total_visit, color = colors[j], label=program_ID[i], linewidth=3)
    else:
        plt.plot(t_visit, phase_curve_total_visit, color = colors[j], linewidth=3)
    x_text = np.mean(t_visit)
    y_text = np.max(phase_curve_total_visit)
    plt.text(x_text, 1.05*y_text, "Visit "+visit[i], fontsize=12, ha='center', va='bottom', color = colors[j], bbox=dict(facecolor='white', alpha=0.6, edgecolor='white', boxstyle='square,pad=0.3'), zorder=10)

plt.xlim(t0, t0+nb_days)
plt.xlabel(r"Time ($BJD_{TBD} - 2450000$)")
plt.ylabel(r"$F_{planet}/F_{star}$ (ppm)")

plt.title("Comparison of TRAPPIST-1 phase curves with and without thick atmospheres")
plt.legend(ncols=2)
plt.grid()

if save_plots:
    plt.savefig("Comparisons_bare_rock_atm/phase_curve_comparison_"+planets+"_"+filter+"_"+str(t0)+".png", bbox_inches='tight')

plt.show()