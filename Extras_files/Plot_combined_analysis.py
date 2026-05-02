#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# March 2026
# Plotting the combined analysis of the light curve

import numpy as np
import matplotlib.pyplot as plt

time, flux, flux_err = np.loadtxt('LC_combined_reduction_F1500W.txt', delimiter=',', unpack=True, skiprows=1)

time -= 2450000  # Convert to BJD - 2450000

plt.figure(figsize=(10, 6))
plt.errorbar(time, flux, yerr=flux_err, fmt='o', markersize=5, color='blue', ecolor='lightgray', elinewidth=1, capsize=2)
plt.xlabel(r'Time ($BJD_{TBD}$ - 2450000)', fontsize=14)
plt.ylabel('Flux (normalized)', fontsize=14)
plt.title('Combined Analysis', fontsize=16)
plt.grid(True)
plt.legend()
plt.tight_layout()
# plt.savefig('Combined_Analysis.png', bbox_inches='tight')
plt.show()