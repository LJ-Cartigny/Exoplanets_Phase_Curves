#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# April 2025
# Placement of the JWST observations on the phase curves

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
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

nb_days = np.max(t_end) - np.min(t_start) + 3
t0 = np.min(t_start)-1.5

print("t0 = ", t0)
print("t_end = ", t0+nb_days)
print("nb_days = ", nb_days)

nb_points = 100000
Keplerian = True

phase_curve_simulation(t0, nb_days, nb_points=nb_points, Keplerian=Keplerian, plot=False,save_plot=True,save_txt=True) # Comment if the simulation is already done


# Overall plot

fig_overall = plt.figure(figsize=(32, 18))

for p in "bcdefgh":  # Comment to not plot the individual planets
    t_simu, phase_simu = np.loadtxt("Phase_curve_TTV_output/phase_curve_"+p+"_TTV_"+str(t0)+".txt", delimiter=",", skiprows=1, unpack=True)
    plt.plot(t_simu, phase_simu, '--', label=p, linewidth=0.5)

t_total_simu, phase_curve_total_simu = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_TTV_"+str(t0)+".txt", delimiter=",",skiprows=1, unpack=True)

colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink']
j = -1

line, = plt.plot(t_total_simu, phase_curve_total_simu, '--', color='grey', label="Total flux")


for i in range(len(t_start)):
    t0 = t_start[i]
    t_visit, phase_curve_total_visit = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_TTV_"+str(t0)+".txt", delimiter=',',skiprows=1, unpack=True)
    if i == 0 or program_ID[i] != program_ID[i-1]:
        j+=1
        plt.plot(t_visit, phase_curve_total_visit, color = colors[j], label=program_ID[i], linewidth=3)
    else:
        plt.plot(t_visit, phase_curve_total_visit, color = colors[j], linewidth=3)
    x_text = np.mean(t_visit)
    y_text = np.max(phase_curve_total_visit)
    plt.text(x_text, 1.05*y_text, "Visit "+visit[i], fontsize=12, ha='center', va='bottom', color = colors[j], bbox=dict(facecolor='white', alpha=0.6, edgecolor='white', boxstyle='square,pad=0.3'), zorder=10)

plt.xlabel(r"Time ($BJD_{TBD} - 2450000$)")
plt.ylabel(r"$L_{planet}/L_{star}$ (ppm)")
plt.title("JWST Observations over the phase curves of TRAPPIST-1 from Oct 2022 to Dec 2024")
plt.legend(loc='lower right', ncol = 2)
plt.grid()
# plt.savefig("JWST_Obs_plots/JWST_Obs_phase_curves_Oct2022-Dec2024.png", bbox_inches='tight') # Uncomment to save the figure
lines = plt.gca().get_lines()
texts = plt.gca().texts
ax_orig = plt.gca()
plt.show()



# Close-up on the observations

xlims = [(9879,9920),(10130,10150),(10270,10275),(10610,10651)] # Values found after looking at the first plot
widths = [xmax-xmin for (xmin, xmax) in xlims]
fig = plt.figure(figsize=(32, 18))
gs = gridspec.GridSpec(1, len(xlims), width_ratios=widths, wspace=0.05)

axes = []
for i in range(len(xlims)):
    if i == 0:
        ax = fig.add_subplot(gs[i])
    else:
        ax = fig.add_subplot(gs[i], sharey=axes[0])
    axes.append(ax)

for i, (ax, (xmin, xmax)) in enumerate(zip(axes, xlims)):
    for line in lines:
        x_data = line.get_xdata()
        y_data = line.get_ydata()
        mask = (x_data >= xmin) & (x_data <= xmax)
        ax.plot(x_data[mask], y_data[mask],color=line.get_color(), linewidth=line.get_linewidth(), linestyle=line.get_linestyle(),label=line.get_label())

    for txt in texts:
        x_txt, y_txt = txt.get_position()
        if xmin <= x_txt <= xmax:
            ax.text(x_txt, y_txt, txt.get_text(), fontsize=txt.get_fontsize(), fontstyle=txt.get_fontstyle(), ha='center', va='bottom', color=txt.get_color(), bbox=dict(facecolor='white', alpha=0.6, edgecolor='white', boxstyle='square,pad=0.3'), zorder=10)

    ax.set_xlim(xmin, xmax)
    ax.xaxis.set_major_locator(MultipleLocator(2.5))
    ax.ticklabel_format(style='plain', axis='x',useOffset=False)
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(labelleft=False, left=False, labelrotation=45)
    ax.grid(True)

axes[0].spines['left'].set_visible(True)
axes[0].tick_params(labelleft=True, left=True)

ysticks = axes[0].get_yticks()
for ax in axes[1:]:
    ax.set_yticks(ysticks)

d = .015
for i in range(len(axes)-1):
    kwargs = dict(transform=axes[i].transAxes, color='k', clip_on=False)
    axes[i].plot((1-d, 1+d), (-d, +d), **kwargs)
    axes[i].plot((1-d, 1+d), (1-d, 1+d), **kwargs)

    kwargs.update(transform=axes[i+1].transAxes)
    axes[i+1].plot((-d, +d), (-d, +d), **kwargs)
    axes[i+1].plot((-d, +d), (1-d, 1+d), **kwargs)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', ncol = 2)

fig.text(0.5, 0.04, r"Time ($BJD_{TBD} - 2450000$)", ha="center")
fig.text(0.04, 0.5, r"$L_{planet}/L_{star}$ (ppm)", va="center", rotation="vertical")
plt.subplots_adjust(wspace=0.05)
plt.suptitle("Close-up on JWST Observations over the phase curves of TRAPPIST-1 from Oct 2022 to Dec 2024", fontsize=20)

plt.tight_layout(rect=[0.05, 0.05, 1, 0.93])

# plt.savefig("JWST_Obs_plots/JWST_Obs_phase_curves_Oct2022-Dec2024_zoom.png", bbox_inches='tight') # Uncomment to save the figure
plt.show()