#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# May 2025
# MIRI Filters

"""
This module is to plot MIRI filters' quantum efficiency.
"""

import numpy as np
import matplotlib.pyplot as plt



def plot_filters():
    """
    Load MIRI filter data from CSV and plot the quantum efficiency curves for F1280W and F1500W filters.
    The plot is saved as 'MIRI_filters.png' and also displayed.
    """
    # Load wavelength and quantum efficiency data for F1280W and F1500W from CSV
    wavelength, quantum_eff_1280, quantum_eff_1500 = np.loadtxt(
        "miri_filter.csv", delimiter=',', skiprows=2, usecols=(0,5,6), unpack=True
    )

    plt.figure(figsize=(16,9))
    plt.plot(wavelength, quantum_eff_1500, label="F1500W", color='blue')
    plt.fill_between(wavelength, quantum_eff_1500, color='blue', alpha=0.2)  # Fill under F1500W filter
    plt.plot(wavelength, quantum_eff_1280, label="F1280W", color='orange')
    plt.fill_between(wavelength, quantum_eff_1280, color='orange', alpha=0.2)  # Fill under F1280W filter
    plt.xlim(10,18)
    plt.xlabel("Wavelength (µm)", size=15)
    plt.ylabel("Quantum efficiency", size=15)
    plt.legend(fontsize=13)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.title("MIRI filters quantum efficiency", size=15)
    plt.grid()
    plt.savefig("MIRI_filters.png", bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    plot_filters()
