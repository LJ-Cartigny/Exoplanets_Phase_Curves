#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# May 2025
# Phoenix model

"""
Generates the PHOENIX spectrum of star TRAPPIST-1 using pysynphot.
"""

import pysynphot as S
import math as math

import numpy as np
import matplotlib.pyplot as plt

# TRAPPIST 1:
Teff = 2566 # Teff = 2566 +/- 26 K from Agol+ 2021 
#Fe_H = 0.04 # [Fe/H] = 0.04 +/- 0.08 from Agol+ 2021 
Fe_H = 0
Rs = 0.1192 * 6.957E10  # Stellar radius planet from Agol+ 2021  in cm
Ms = 0.0802 * 1.99E33 # Stellar mass in g
Jmag = 11.35
Kmag = 10.30
logg = math.log10(6.673E-8 * Ms / Rs**2 ) #cgs
dist = 12.467 * 3.086E18 # stellar distance in pc, converted to cm

# TRAPPIST-1 b planet
Rp = 1.116 * 6.378E8 # cm b planet from Agol+ 2021 
a = 0.01154 * 1.496E13 # cm b planet from Agol+ 2021 
e = 0.001   #Wang+ 2017 K2 (Note Luger+ 2017 K2 data give e = 0.001)
T14 = 0.6 * 0.9 * 3600 # transit duration (36 min x 0.9 eff)
P = 1.510826  # days orbital period b planet from Agol+ 2021 

star = S.Icat('phoenix', Teff, Fe_H, logg)  # Phoenix model
star_norm = star.renorm(Jmag, 'vegamag', S.ObsBandpass('johnson,j'))
star_mJy  = star.renorm(Jmag, 'vegamag', S.ObsBandpass('johnson,j'))

star_norm.convert('Micron')
star_norm.convert('flam')  # flam units: erg/s/cm^2/Ang
star_mJy.convert('Micron')
star_mJy.convert('mjy')  # convert flux to mJy for JWST ETC input



np.savetxt("TRAPPIST1_Phoenix_model.txt", np.column_stack((star_mJy.wave, star_mJy.flux)), header="Wavelength (microns), Flux (mJy)")

plt.figure(figsize=(32, 18))
plt.plot(star_mJy.wave, star_mJy.flux, label='Phoenix model')
plt.xlabel('Wavelength (microns)')
plt.ylabel('Flux (mJy)')
plt.title('Phoenix Model Spectrum for TRAPPIST-1')
plt.legend()
plt.grid()
plt.xlim(10,20)
plt.ylim(0,5)
plt.show()