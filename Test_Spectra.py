#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# May 2025
# Test on the PHOENIX and SPHINX spectra

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy import constants as const
import pysynphot as S
import math as math

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


# Comparaison measurement with aisha T-1 model
df_filters = pd.read_csv("miri_filter.csv", sep=',')
df_filters = df_filters.drop([0])
df_filters["#Wave"] = [float(x) for x in df_filters["#Wave"]]
df_filters["F1280W"] = [float(x) for x in df_filters["F1280W"]]
df_filters["#Wave"] = [float(x) for x in df_filters["#Wave"]]
df_filters["F1500W"] = [float(x) for x in df_filters["F1500W"]]
df_filters["#Wave"] = [float(x) for x in df_filters["#Wave"]]
df_filters["F1000W"] = [float(x) for x in df_filters["F1000W"]]
df_filters["#Wave"] = [float(x) for x in df_filters["#Wave"]]
df_filters["F1800W"] = [float(x) for x in df_filters["F1800W"]]

#WE LOAD JWST MIRI FILTER DATA AND RE-INTERPOLATE ON A REGULAR GRID
F1280W_raw=np.loadtxt('JWST_MIRI.F1280W.dat.txt')
F1500W_raw=np.loadtxt('JWST_MIRI.F1500W.dat.txt')
wv_1280W=np.arange(F1280W_raw[0,0],F1280W_raw[len(F1280W_raw[:,0])-1,0],50.) # in Angstroms
wv_1500W=np.arange(F1500W_raw[0,0],F1500W_raw[len(F1500W_raw[:,0])-1,0],50.) # in Angstroms
F1280W=np.zeros((len(wv_1280W),2),dtype='f')
F1500W=np.zeros((len(wv_1500W),2),dtype='f')
f_F1280W=interpolate.interp1d(F1280W_raw[:,0],F1280W_raw[:,1])
f_F1500W=interpolate.interp1d(F1500W_raw[:,0],F1500W_raw[:,1])
for i in range(0,len(wv_1280W),1):
    F1280W[i,0]=wv_1280W[i]
    F1280W[i,1]=f_F1280W(wv_1280W[i])
for i in range(0,len(wv_1500W),1):
    F1500W[i,0]=wv_1500W[i]
    F1500W[i,1]=f_F1500W(wv_1500W[i])

file_aisha = "sphinx_spectrum_T-1_aisha.txt"
spec_T1_sphinx = pd.read_csv(file_aisha,
                             skiprows=1, delim_whitespace=True, header=None)
spec_T1_sphinx.columns = ["#Wavelength(micrometer)", "Flux(W/m2/m)"]

R_star = 0.1192 * 6.957E10  # Stellar radius planet from Agol+ 2021  in cm
dist = 12.467 * 3.086E18 # stellar distance in pc, converted to cm
flux_sphinx_Jy = [1e26*spec_T1_sphinx["Flux(W/m2/m)"][i]*((spec_T1_sphinx["#Wavelength(micrometer)"][i]*1e-6**2)/const.c.value)*(R_star/dist)**2*spec_T1_sphinx["#Wavelength(micrometer)"][i]
                  for i in range(len(spec_T1_sphinx))]

# %matplotlib inline
plt.figure(figsize=(10,6))
ax = plt.gca()

l_eff_F1280 = 12.62 #12.83
l_eff_F1500 = 14.79 #15.09

#Observations
x_F1000W = df_filters["#Wave"]
y_F1000W = df_filters["F1000W"]
f_F1000W = interpolate.interp1d(x_F1000W,y_F1000W) 
x_F1280W = df_filters["#Wave"]
y_F1280W = df_filters["F1280W"]
f_F1280W = interpolate.interp1d(x_F1280W,y_F1280W)
x_F1500W = df_filters["#Wave"]
y_F1500W = df_filters["F1500W"]
f_F1500W = interpolate.interp1d(x_F1500W,y_F1500W) 
x_F1800W = df_filters["#Wave"]
y_F1800W = df_filters["F1800W"]
f_F1800W = interpolate.interp1d(x_F1500W,y_F1800W) 

flux_measured_12 = np.mean([3.425,3.428, 3.427, 3.428])
flux_measured_15 = 2.589 #2.54 #2.559 #2.844 #

plt.errorbar(l_eff_F1280, flux_measured_12 ,yerr=0.10284+0.024,
             xerr=[np.abs(l_eff_F1280-min(x_F1280W[np.where((y_F1280W>0.2))[0]]))],
             fmt='h', color='k', capsize=10, markeredgewidth=2, markersize=10, zorder=5, alpha=0.5)#ecolor ='k')
#plt.errorbar(12.62, 3.725,yerr=0.031, fmt='o', color='r', capsize=10, markeredgewidth=2, markersize=10, )#ecolor ='k')
plt.errorbar(l_eff_F1500, flux_measured_15,#np.mean([2.532,2.528, 2.529,2.520,2.531,2.546,2.547, 2.599])
                            yerr=0.079,xerr=[np.abs(l_eff_F1500-min(x_F1500W[np.where((y_F1500W>0.2))[0]]))], 
                                                                    fmt='o', color='k', alpha=0.5,capsize=10, markeredgewidth=2, markersize=10,zorder=5,label="measurements" )#ecolor ='k')

# models 
## sphinx
t = np.linspace(1,19,1000)
x_model_sphinx = spec_T1_sphinx["#Wavelength(micrometer)"]
y_model_sphinx = np.array(flux_sphinx_Jy)*1e3 /1.07
f_model_sphinx = interpolate.interp1d(x_model_sphinx,y_model_sphinx ) 
plt.plot(t, f_model_sphinx(t), color='blue', alpha=0.4, label= "Sphinx model")
spec1_F1500W=f_model_sphinx(F1500W[:,0]*1e-4)
cum1_F1500W=0.
norm_F1500W=0.
for i in range(1,len(F1500W[:,0]),1):
    norm_F1500W=norm_F1500W+F1500W[i,1]*(1./F1500W[i,0]-1./F1500W[i-1,0])
    cum1_F1500W+=spec1_F1500W[i]*F1500W[i,1]*(1./F1500W[i,0]-1./F1500W[i-1,0])
spec1_F1280W=f_model_sphinx(F1280W[:,0]*1e-4)
cum1_F1280W=0.
norm_F1280W=0.
for i in range(1,len(F1280W[:,0]),1):
    # norm_F1280W=norm_F1280W+F1500W[i,1]*(1./F1280W[i,0]-1./F1280W[i-1,0])
    norm_F1280W=norm_F1280W+F1280W[i,1]*(1./F1280W[i,0]-1./F1280W[i-1,0])
    cum1_F1280W+=spec1_F1280W[i]*F1280W[i,1]*(1./F1280W[i,0]-1./F1280W[i-1,0])
ax.scatter(l_eff_F1500,cum1_F1500W/norm_F1500W,#np.median(f_model_sphinx(t[np.where((f_F1500W(t)!=0))[0]])),
                marker='h',s=100,color='blue',edgecolors='k',linewidth=1.5 ,zorder=7,alpha=0.5)
ax.scatter(l_eff_F1280,cum1_F1280W/norm_F1280W,#np.median(f_model_sphinx(t[np.where((f_F1280W(t)!=0))[0]])),
                marker='h',s=100,color='blue',edgecolors='k',linewidth=1.5 ,zorder=7,alpha=0.5,label="Sphinx integrated MIRI band")
print("Fluw modeled at 12.8 corrected SPHINX", cum1_F1280W / norm_F1280W)
print("Fluw modeled at 15 corrected SPHINX", cum1_F1500W / norm_F1500W)
print("Difference in % between measurements and corrected Sphinx at 12.8 :", (cum1_F1280W/norm_F1280W/flux_measured_12)*100 - 100)
print("Difference in % between measurements and corrected Sphinx at 15 :", (cum1_F1500W/norm_F1500W/flux_measured_15)*100 - 100)

## phoenix
t = np.linspace(1,19,2000)
x_model_phoenix = star_mJy.wave
y_model_phoenix = star_mJy.flux
f_model_phoenix = interpolate.interp1d(x_model_phoenix,y_model_phoenix ) 
plt.plot(t, f_model_phoenix(t),color='orange', alpha=0.4,label="PHOENIX")

spec1_F1500W=f_model_phoenix(F1500W[:,0]*1e-4)
cum1_F1500W=0.
norm_F1500W=0.
for i in range(1,len(F1500W[:,0]),1):
    norm_F1500W=norm_F1500W+F1500W[i,1]*(1./F1500W[i,0]-1./F1500W[i-1,0])
    cum1_F1500W+=spec1_F1500W[i]*F1500W[i,1]*(1./F1500W[i,0]-1./F1500W[i-1,0])
spec1_F1280W=f_model_phoenix(F1280W[:,0]*1e-4)
cum1_F1280W=0.
norm_F1280W=0.
for i in range(1,len(F1280W[:,0]),1):
    # norm_F1280W=norm_F1280W+F1500W[i,1]*(1./F1280W[i,0]-1./F1280W[i-1,0])
    norm_F1280W=norm_F1280W+F1280W[i,1]*(1./F1280W[i,0]-1./F1280W[i-1,0])
    cum1_F1280W+=spec1_F1280W[i]*F1280W[i,1]*(1./F1280W[i,0]-1./F1280W[i-1,0])
    
ax.scatter(l_eff_F1500,cum1_F1500W/norm_F1500W,#np.median(f_model_phoenix(t[np.where((f_F1500W(t)!=0))[0]])),
                marker='h',s=100,color='orange',edgecolors='k',linewidth=1.5 ,zorder=7)
ax.scatter(l_eff_F1280,cum1_F1280W/norm_F1280W,#np.median(f_model_phoenix(t[np.where((f_F1280W(t)!=0))[0]])),
                marker='h',s=100,color='orange',edgecolors='k',linewidth=1.5 ,zorder=7,label="Phoenix integrated MIRI band")

print("Fluw modeled at 12.8 PHOENIX", cum1_F1280W / norm_F1280W)
print("Fluw modeled at 15 PHOENIX", cum1_F1500W / norm_F1500W)
print("Difference in % between measurements and PHOENIX at 12.8 :", (cum1_F1280W/norm_F1280W/flux_measured_12)*100 - 100)
print("Difference in % between measurements and PHOENIX at 15 :", (cum1_F1500W/norm_F1500W/flux_measured_15)*100 - 100)


ax.fill_between(t,f_F1280W(t)*4,color='pink',alpha=0.8,linewidth=2,)#label='F1280W')
ax.fill_between(t,f_F1500W(t)*4,color='brown',alpha=0.5,linewidth=2,)#label='F1500W')
plt.xlim(11,17.5)
plt.ylim(0,5.)
plt.ylabel("Flux (mJy)",fontsize=14)
plt.xlabel(r"Wavelength ($\mu m$)",fontsize=14)
ax.text(12.2, 0.35, "F1280W", fontsize=14, fontweight="bold",color='red')
ax.text(14.5, 0.35, "F1500W", fontsize=14, fontweight="bold",color='brown')
ax.tick_params(axis="y", labelsize=14, labelcolor="k")
ax.tick_params(axis='x', labelsize=14,labelcolor='k')
plt.legend(fontsize=12, loc=1)
# plt.savefig("compa_measures_to_models.png")
plt.show()