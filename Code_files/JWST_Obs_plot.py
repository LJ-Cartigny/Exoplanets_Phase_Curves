#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# April 2025
# Placement of the JWST observations on the phase curves

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
# from matplotlib import container
import matplotlib as mpl
from astropy.time import Time

from Phase_curve_TTV import phase_curve_simulation
from JWST_Obs_simu import phase_curve_visit

# Set the style for the plots

mpl.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 20,
    'axes.titlesize': 20,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 18,
    'figure.figsize': (16, 10),
    # 'lines.linewidth': 2,
    # 'grid.alpha': 0.5,
    # 'grid.linestyle': '--',
})

# Parameters

nb_points = 100000

Keplerian = True

planets = 'defgh'

redistribution = 0 # 0 for bare rocks, 1 for thick atmospheres (0 by default if comparison is True)

filter = 'F1280W'

unit = 'mJy' # 'ppm' or 'mJy' ('mJy' by default if plot_obs_points is True or model is 'phoenix')

model = 'phoenix' # 'phoenix' or 'sphinx'

save_plots = False # Write True if you want to save the plots

do_simulation = False # Write True if the simulation hasn't been done yet

plot_individual_planets = False # Write True if you want to plot the individual planets as bare rocks to see their phases

comparison = True # Write True if you want to compare the bare rock and thick atmosphere cases

plot_obs_points = True # Write True if you want to plot the observations points (in mJy) on the phase curves

points_offset = True # Write True if you want to add an offset to the observation points to place them closer to the phase curves (useful if the observations are too far from the phase curves)


# def plot_jwst_phase_curves(planets,
#                            filter,
#                            model,
#                            nb_points=10000,
#                            Keplerian=True,
#                            redistribution=0,
#                            unit='mJy',
#                            save_plots=False,
#                            do_simulation=False,
#                            plot_individual_planets=True,
#                            comparison=True,
#                            plot_obs_points=True,
#                            points_offset=False):
#     pass


# Simulations

# program_ID, visit, t_start, t_end, filter_obs, flux_obs, err_obs = np.loadtxt("JWST_Obs_times.txt", delimiter=',', skiprows=2, unpack=True,dtype=str)
program_ID, visit, t_start, t_end, filter_obs, flux_obs, err_obs = np.loadtxt("JWST_Obs_times.txt", usecols=(0,1,2,3,4,7,8), delimiter=',', skiprows=2, unpack=True,dtype=str) # With new stellar flux values

flux_obs = flux_obs.astype(float)
err_obs = err_obs.astype(float)

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

if comparison:
    redistribution = 0

if plot_obs_points or model=='phoenix':
    unit = 'mJy' # If we plot the observed fluxes or use the PHOENIX model, we use mJy as unit

if do_simulation:
    phase_curve_visit(planets,redistribution,filter,model,unit,Keplerian=Keplerian)
    print("Simulating the phase curves between t0 and t_end...")
    phase_curve_simulation(t0, nb_days, nb_points=nb_points, planets=planets, redistribution=redistribution, filter=filter, model=model, unit=unit, Keplerian=Keplerian, plot=False,save_plot=True,save_txt=True)
    
    if comparison:
        phase_curve_simulation(t0, nb_days, nb_points=nb_points, planets=planets, redistribution=1, filter=filter, model=model, unit=unit, Keplerian=Keplerian, plot=False,save_plot=True,save_txt=True)

if plot_obs_points:
    errorbar_data = []


# Overall plot

fig_overall = plt.figure()

if plot_individual_planets:
    for p in planets:
        if filter == None:
            t_simu, phase_simu = np.loadtxt("Phase_curve_TTV_output/phase_curve_"+p+"_bolometric_"+str(t0)+".txt", delimiter=",", skiprows=1, unpack=True)
        else:
            if redistribution == 0:
                t_simu, phase_simu = np.loadtxt("Phase_curve_TTV_output/phase_curve_"+p+"_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", delimiter=",", skiprows=1, unpack=True)
            else:
                t_simu, phase_simu = np.loadtxt("Phase_curve_TTV_output/phase_curve_"+p+"_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", delimiter=",", skiprows=1, unpack=True)
        plt.plot(t_simu, phase_simu, '--', label=p+" (bare rock)", linewidth=0.5)

if filter == None:
    t_total_simu, phase_curve_total_simu = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_bolometric_"+str(t0)+".txt", delimiter=",",skiprows=1, unpack=True)
else:
    if redistribution == 0:
        t_total_simu, phase_curve_total_simu = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", delimiter=",",skiprows=1, unpack=True)
    else:
        t_total_simu, phase_curve_total_simu = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", delimiter=",",skiprows=1, unpack=True)

if comparison:
    t_total_simu_atm, phase_curve_total_simu_atm = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", delimiter=",",skiprows=1, unpack=True)
    plt.plot(t_total_simu_atm, phase_curve_total_simu_atm, '--', color='grey', label="Total (atmospheres)")

colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink']
j = -1

if comparison:
    label_total = "Total (bare rocks)"
else:
    label_total = "Total flux"

line, = plt.plot(t_total_simu, phase_curve_total_simu, '-', color='black', label=label_total)


# Plot the visits

for i in range(len(t_start)):
    t0 = t_start[i]
    offset = 0
    if filter == None:
        t_visit, phase_curve_total_visit = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_bolometric_"+str(t0)+".txt", delimiter=',',skiprows=1, unpack=True)
    else:
        if redistribution == 0:
            t_visit, phase_curve_total_visit = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", delimiter=',',skiprows=1, unpack=True)
        else:
            t_visit, phase_curve_total_visit = np.loadtxt("Phase_curve_TTV_output/phase_curve_total_"+planets+"_atm_"+filter+"_"+model+"_"+unit+"_"+str(t0)+".txt", delimiter=',',skiprows=1, unpack=True)

    if i == 0 or program_ID[i] != program_ID[i-1]:
        j+=1
        # plt.plot(t_visit, phase_curve_total_visit, color = colors[j], label=program_ID[i], linewidth=3)
        if plot_obs_points and filter_obs[i]==filter and flux_obs[i] != np.nan:
            if points_offset and program_ID[i] == 'GTO_1279':
                if model == "phoenix":
                    offset = -0.16
                else:
                    offset = 0.04
            elif points_offset and program_ID[i] == 'GO_5191':
                if model == "phoenix":
                    offset = -0.19
                else:
                    offset = 0.01
            elif points_offset and (program_ID[i] == 'GTO_1177' or program_ID[i] == 'GO_2304'):
                if model == "phoenix":
                    offset = -0.215 # To be determined
                else:
                    offset = -0.13
            elif points_offset and program_ID[i] == 'GO_3077':
                if model == "phoenix":
                    offset = -0.15 # To be determined
                else:
                    offset = -0.06
            else:
                offset = 0
            plt.errorbar(np.mean(t_visit), flux_obs[i]+offset, yerr=err_obs[i], fmt='h', color=colors[j], markersize=5, elinewidth=2, capsize=5, label=program_ID[i]+" (observed)", zorder=10)
            errorbar_data.append((np.mean(t_visit), flux_obs[i]+offset, err_obs[i], dict(fmt='h',color=colors[j], markersize=5, elinewidth=2, capsize=5, label=program_ID[i]+" (observed)", zorder=10)))
    else:
        # plt.plot(t_visit, phase_curve_total_visit, color = colors[j], linewidth=3)
        if plot_obs_points and filter_obs[i]==filter and flux_obs[i] != np.nan:
            if points_offset and program_ID[i] == 'GTO_1279':
                if model == "phoenix":
                    offset = -0.16
                else:
                    offset = 0.04
            elif points_offset and program_ID[i] == 'GO_5191':
                if model == "phoenix":
                    offset = -0.19
                else:
                    offset = 0.01
            elif points_offset and (program_ID[i] == 'GTO_1177' or program_ID[i] == 'GO_2304'):
                if model == "phoenix":
                    offset = -0.215
                else:
                    offset = -0.13
            elif points_offset and program_ID[i] == 'GO_3077':
                if model == "phoenix":
                    offset = -0.15
                else:
                    offset = -0.06
            else:
                offset = 0
            plt.errorbar(np.mean(t_visit), flux_obs[i]+offset, yerr=err_obs[i], fmt='h', color=colors[j], markersize=5, elinewidth=2, capsize=5, zorder=10)
            errorbar_data.append((np.mean(t_visit), flux_obs[i]+offset, err_obs[i], dict(fmt='h',color=colors[j], markersize=5, elinewidth=2, capsize=5, zorder=10)))
    
    # Add visit label
    if plot_obs_points:
        x_text = np.mean(t_visit)
        y_text = flux_obs[i]+1.5*err_obs[i] + offset
        plt.text(x_text, y_text, "Visit "+visit[i], fontsize=18, ha='center', va='bottom', color = colors[j], bbox=dict(facecolor='white', alpha=0.6, edgecolor='white', boxstyle='square,pad=0.3'), zorder=10)
    else:
        x_text = np.mean(t_visit)
        y_text = np.max(phase_curve_total_visit)
        plt.text(x_text, y_text + 0.2 * np.ptp(phase_curve_total_visit), "Visit "+visit[i], fontsize=18, ha='center', va='bottom', color = colors[j], bbox=dict(facecolor='white', alpha=0.6, edgecolor='white', boxstyle='square,pad=0.3'), zorder=10)
    # plt.text(x_text, 1.05*y_text, "Visit "+visit[i], fontsize=12, ha='center', va='bottom', color = colors[j], bbox=dict(facecolor='white', alpha=0.6, edgecolor='white', boxstyle='square,pad=0.3'), zorder=10)

plt.xlabel(r"Time ($BJD_{TBD} - 2450000$)")
if unit == 'ppm':
    plt.ylabel(r"$F_{planet}/F_{star}$ (ppm)")
else:
    plt.ylabel(r"$F_{star}+F_{planet}$ (mJy)")

# if filter == None:  # No title for article figures
#     plt.title("JWST Observations over the phase curves of TRAPPIST-1 planets with bolometric fluxes from Oct 2022 to Dec 2024")
# else:
#     if redistribution == 1:
#         plt.title("JWST Observations over the phase curves of TRAPPIST-1 planets with atmospheres with MIRI "+filter+" filter using the "+model.upper()+" model from Oct 2022 to Dec 2024")
#     elif comparison:
#         plt.title("JWST Observations over the phase curves of TRAPPIST-1 planets with and without thick atmospheres with MIRI "+filter+" filter using the "+model.upper()+" model from Oct 2022 to Dec 2024")
#     else:
#         plt.title("JWST Observations over the phase curves of TRAPPIST-1 planets as bare rocks with MIRI "+filter+" filter using the "+model.upper()+" model from Oct 2022 to Dec 2024")

plt.legend(loc='lower right', ncol = 2)
plt.grid()

if save_plots:
    if filter == None:
        plt.savefig("JWST_Obs_plots/JWST_Obs_phase_curves_"+planets+"_bolometric_Oct2022-Dec2024.png", bbox_inches='tight')
    else:
        if plot_obs_points and comparison:
            # plt.savefig("Comparisons_obs_JWST/comparison_obs_"+planets+"_"+filter+"_"+model+"_Oct2022-Dec2024.png", bbox_inches='tight')
            plt.savefig("Article_figures/comparison_obs_"+planets+"_"+filter+"_"+model+"_Oct2022-Dec2024.png", bbox_inches='tight') # For article
        elif redistribution == 1:
            plt.savefig("JWST_Obs_plots/JWST_Obs_phase_curves_"+planets+"_atm_"+filter+"_"+model+"_"+unit+"_Oct2022-Dec2024.png", bbox_inches='tight')
        elif comparison:
            # plt.savefig("Comparisons_bare_rock_atm/comparison_"+planets+"_"+filter+"_"+model+"_"+unit+"_Oct2022-Dec2024.png", bbox_inches='tight')
            plt.savefig("Article_figures/comparison_"+planets+"_"+filter+"_"+model+"_"+unit+"_Oct2022-Dec2024.png", bbox_inches='tight') # For article
        else:
            plt.savefig("JWST_Obs_plots/JWST_Obs_phase_curves_"+planets+"_"+filter+"_"+model+"_"+unit+"_Oct2022-Dec2024.png", bbox_inches='tight')

lines = plt.gca().get_lines()
texts = plt.gca().texts
ax_orig = plt.gca()
plt.show()



# Close-up on the observations

if filter == 'F1500W' and plot_obs_points:
    xlims = [(9879,9920),(10270,10275)] # Values found after looking at the overall plot
elif filter == 'F1280W' and plot_obs_points:
    xlims = [(10130,10150),(10610,10651)]
else:
    xlims = [(9879,9920),(10130,10150),(10270,10275),(10610,10651)] # Values found after looking at the first plot
widths = [xmax-xmin for (xmin, xmax) in xlims]
fig = plt.figure()
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

    
    if plot_obs_points:
        for (x, y, err, fmt) in errorbar_data:
            if xmin <= x <= xmax:
                ax.errorbar(x, y, yerr=err, **fmt)

    for txt in texts:
        x_txt, y_txt = txt.get_position()
        if xmin <= x_txt <= xmax:
            ax.text(x_txt, y_txt, txt.get_text(), fontsize=txt.get_fontsize(), fontstyle=txt.get_fontstyle(), ha='center', va='bottom', color=txt.get_color(), bbox=dict(facecolor='white', alpha=0.6, edgecolor='white', boxstyle='square,pad=0.3'), zorder=10)

    ax.set_xlim(xmin, xmax)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.ticklabel_format(style='plain', axis='x',useOffset=False)
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%d'))
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

# Add dummy points for errorbar legend entries
for (_, _, _, fmt) in errorbar_data:
    label = fmt.get("label")
    if label:  # only once per label
        axes[0].plot([], [], fmt['fmt'], color=fmt['color'], label=label)
        fmt["label"] = None  # prevent duplicates


handles, labels = axes[0].get_legend_handles_labels()

if plot_obs_points == False:
    fig.legend(handles, labels, loc='upper right',  ncols = 2)
else:
    fig.legend(handles, labels, ncols = 2)

fig.text(0.5, 0.01, r"Time ($BJD_{TBD} - 2450000$)", ha="center")
if unit == 'ppm':
    fig.text(0.05, 0.5, r"$F_{p}/F_{s}$ (ppm)", va="center", rotation="vertical")
else:
    fig.text(0.05, 0.5, r"$F_{s}+F_{p}$ (mJy)", va="center", rotation="vertical")
plt.subplots_adjust(wspace=0.05)

# if filter == None: # No title for article figures
#     plt.suptitle("Close-up on JWST Observations over the phase curves of TRAPPIST-1 planets\n with bolometric fluxes from Oct 2022 to Dec 2024", fontsize=20)
# else:
#     if redistribution == 1:
#         plt.suptitle("Close-up on JWST Observations over the phase curves of TRAPPIST-1 planets\n with atmospheres with MIRI "+filter+" filter using the "+model.upper()+" model from Oct 2022 to Dec 2024", fontsize=20)
#     elif comparison:
#         plt.suptitle("Close-up on JWST Observations over the phase curves of TRAPPIST-1 planets\n with and without thick atmospheres with MIRI "+filter+" filter using the "+model.upper()+" model from Oct 2022 to Dec 2024")
#     else:
#         plt.suptitle("Close-up on JWST Observations over the phase curves of TRAPPIST-1 planets\n as bare rocks with MIRI "+filter+" filter using the "+model.upper()+" model from Oct 2022 to Dec 2024", fontsize=20)

# plt.tight_layout(rect=[0.05, 0.05, 1, 0.93])

if save_plots:
    if filter == None:
        plt.savefig("JWST_Obs_plots/JWST_Obs_phase_curves_"+planets+"_bolometric_Oct2022-Dec2024_zoom.png", bbox_inches='tight')
    else:
        if plot_obs_points and comparison:
            # plt.savefig("Comparisons_obs_JWST/comparison_obs_"+planets+"_"+filter+"_"+model+"_Oct2022-Dec2024_zoom.png", bbox_inches='tight')
            plt.savefig("Article_figures/comparison_obs_"+planets+"_"+filter+"_"+model+"_Oct2022-Dec2024_zoom.png", bbox_inches='tight') # For article
        elif redistribution == 1:
            plt.savefig("JWST_Obs_plots/JWST_Obs_phase_curves_"+planets+"_atm_"+filter+"_"+model+"_"+unit+"_Oct2022-Dec2024_zoom.png", bbox_inches='tight')
        elif comparison:
            # plt.savefig("Comparisons_bare_rock_atm/comparison_"+planets+"_"+filter+"_"+model+"_"+unit+"_Oct2022-Dec2024_zoom.png", bbox_inches='tight')
            plt.savefig("Article_figures/comparison_"+planets+"_"+filter+"_"+model+"_"+unit+"_Oct2022-Dec2024_zoom.png", bbox_inches='tight') # For article
        else:
            plt.savefig("JWST_Obs_plots/JWST_Obs_phase_curves_"+planets+"_"+filter+"_"+model+"_"+unit+"_Oct2022-Dec2024_zoom.png", bbox_inches='tight')

plt.show()