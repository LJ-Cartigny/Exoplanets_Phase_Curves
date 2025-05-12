#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# May 2025
# MIRI Filters

import numpy as np
import matplotlib.pyplot as plt


wavelength, quantum_eff_1280, quantum_eff_1500 = np.loadtxt("miri_filter.csv", delimiter=',', skiprows=2, usecols=(0,5,6), unpack=True)

plt.figure(figsize=(16,9))
plt.plot(wavelength, quantum_eff_1500, label="F1500W", color='blue')
plt.plot(wavelength, quantum_eff_1280, label="F1280W", color='orange')
plt.xlim(10,18)
plt.xlabel("Wavelength (µm)")
plt.ylabel("Quantum efficiency")
plt.title("MIRI filters quantum efficiency")
plt.grid()
# plt.savefig("MIRI_F1500W_filter.png")
plt.show()
