#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# March 2025
# Exoplanet transit timing variations

import numpy as np
import matplotlib.pyplot as plt
from TRAPPIST1_parameters import *

def period_TTV(P, transit_start, transit_end):
    """
    Computes the modified orbital periods (in days) of the planet due to the TTVs

    :param P: the initial period of the planet without TTVs (in days)
    :type P: float

    :param transit_start: the start of the transit (in days)
    :type transit_start: numpy.ndarray

    :param transit_end: the end of the transit (in days)
    :type transit_end: numpy.ndarray

    :return: P_TTV
    :rtype: numpy.ndarray
    """
    P_TTV = np.zeros(len(transit_start))
    for i in range(len(transit_start)-1):
        P_TTV[i] = ((transit_start[i+1]-transit_start[i]) + (transit_end[i+1]-transit_end[i]))/2
    P_TTV[-1] = P # the last period is set to the initial period

    return P_TTV

def transit_peak(transit_start, transit_end):
    """
    Computes the peak of the transit (in days)

    :param transit_start: the start of the transit (in days)
    :type transit_start: numpy.ndarray

    :param transit_end: the end of the transit (in days)
    :type transit_end: numpy.ndarray

    :return: transit_peak
    :rtype: numpy.ndarray
    """
    transit_peak = (transit_start + transit_end)/2

    return transit_peak



def main():

    TTV_file = "Files_TTV/T1_transit_timing_forecast_cycle5.csv"
    planet, epoch, transit_start, transit_end = np.loadtxt(TTV_file, delimiter=',', usecols=(0,1,2,4), unpack=True)

    #print(epoch.shape)

    i = 0
    while planet[i] == 1:
        i += 1
    j = i
    epoch_b = epoch[:i]
    transit_start_b = transit_start[:i]
    transit_end_b = transit_end[:i]

    while planet[i] == 2:
        i += 1
    epoch_c = epoch[j:i]
    transit_start_c = transit_start[j:i]
    transit_end_c = transit_end[j:i]

    j = i
    while planet[i] == 3:
        i += 1
    epoch_d = epoch[j:i]
    transit_start_d = transit_start[j:i]
    transit_end_d = transit_end[j:i]

    j = i
    while planet[i] == 4:
        i += 1
    epoch_e = epoch[j:i]
    transit_start_e = transit_start[j:i]
    transit_end_e = transit_end[j:i]

    j = i
    while planet[i] == 5:
        i += 1
    epoch_f = epoch[j:i]
    transit_start_f = transit_start[j:i]
    transit_end_f = transit_end[j:i]

    j = i
    while planet[i] == 6:
        i += 1
    epoch_g = epoch[j:i]
    transit_start_g = transit_start[j:i]
    transit_end_g = transit_end[j:i]

    epoch_h = epoch[i:]
    transit_start_h = transit_start[i:]
    transit_end_h = transit_end[i:]


    #print(len(epoch_b)+len(epoch_c)+len(epoch_d)+len(epoch_e)+len(epoch_f)+len(epoch_g)+len(epoch_h))
    
    P_b_TTV = period_TTV(P_b, transit_start_b, transit_end_b)
    transit_peak_b = transit_peak(transit_start_b, transit_end_b)
    # print(P_b)
    # print(P_b_TTV)
    # print(np.mean(P_b_TTV))
    np.savetxt("Files_TTV/TTV_b.txt", np.c_[epoch_b, P_b_TTV, transit_peak_b], delimiter=',', header="Epoch,Modified period,Transit peak", comments='')

    P_c_TTV = period_TTV(P_c, transit_start_c, transit_end_c)
    transit_peak_c = transit_peak(transit_start_c, transit_end_c)
    np.savetxt("Files_TTV/TTV_c.txt", np.c_[epoch_c, P_c_TTV, transit_peak_c], delimiter=',', header="Epoch,Modified period,Transit peak", comments='')

    P_d_TTV = period_TTV(P_d, transit_start_d, transit_end_d)
    transit_peak_d = transit_peak(transit_start_d, transit_end_d)
    np.savetxt("Files_TTV/TTV_d.txt", np.c_[epoch_d, P_d_TTV, transit_peak_d], delimiter=',', header="Epoch,Modified period,Transit peak", comments='')

    P_e_TTV = period_TTV(P_e, transit_start_e, transit_end_e)
    transit_peak_e = transit_peak(transit_start_e, transit_end_e)
    np.savetxt("Files_TTV/TTV_e.txt", np.c_[epoch_e, P_e_TTV, transit_peak_e], delimiter=',', header="Epoch,Modified period,Transit peak", comments='')

    P_f_TTV = period_TTV(P_f, transit_start_f, transit_end_f)
    transit_peak_f = transit_peak(transit_start_f, transit_end_f)
    np.savetxt("Files_TTV/TTV_f.txt", np.c_[epoch_f, P_f_TTV, transit_peak_f], delimiter=',', header="Epoch,Modified period,Transit peak", comments='')

    P_g_TTV = period_TTV(P_g, transit_start_g, transit_end_g)
    transit_peak_g = transit_peak(transit_start_g, transit_end_g)
    np.savetxt("Files_TTV/TTV_g.txt", np.c_[epoch_g, P_g_TTV, transit_peak_g], delimiter=',', header="Epoch,Modified period,Transit peak", comments='')

    P_h_TTV = period_TTV(P_h, transit_start_h, transit_end_h)
    transit_peak_h = transit_peak(transit_start_h, transit_end_h)
    np.savetxt("Files_TTV/TTV_h.txt", np.c_[epoch_h, P_h_TTV, transit_peak_h], delimiter=',', header="Epoch,Modified period,Transit peak", comments='')

if __name__ == "__main__":
    main()