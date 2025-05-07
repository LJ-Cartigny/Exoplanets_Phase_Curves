#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# May 2025
# MIRI Filters

import numpy as np
import matplotlib.pyplot as plt


wavelength, quantum_eff = np.loadtxt("miri_filter.csv", delimiter=',', skiprows=2, usecols=(0,6), unpack=True)

plt.figure(figsize=(16,9))
plt.plot(wavelength, quantum_eff, label="F1500W", color='blue')
plt.xlim(12,18)
plt.xlabel("Wavelength (µm)")
plt.ylabel("Quantum efficiency")
plt.title("MIRI 15 micron filter quantum efficiency")
plt.grid()
# plt.savefig("MIRI_F1500W_filter.png")
plt.show()
