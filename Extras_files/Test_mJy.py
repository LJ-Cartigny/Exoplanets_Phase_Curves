#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# May 2025
# Test conversion mJy

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from TRAPPIST1_parameters import flux_T1_sphinx, wavelengths_T1_sphinx, dist_system, R_star
from Flux_wavelength import conversion_IS_to_mJy

#Observations
#WE LOAD JWST MIRI FILTER DATA AND RE-INTERPOLATE ON A REGULAR GRID
F1280W_raw=np.loadtxt('JWST_MIRI.F1280W.dat.txt')
F1500W_raw=np.loadtxt('JWST_MIRI.F1500W.dat.txt')
wv_1280W=np.arange(F1280W_raw[0,0],F1280W_raw[len(F1280W_raw[:,0])-1,0],50.) # in Angstroms
wv_1500W=np.arange(F1500W_raw[0,0],F1500W_raw[len(F1500W_raw[:,0])-1,0],50.) # in Angstroms
F1280W=np.zeros((len(wv_1280W),2),dtype='f')
F1500W=np.zeros((len(wv_1500W),2),dtype='f')
f_F1280W=interp1d(F1280W_raw[:,0],F1280W_raw[:,1])
f_F1500W=interp1d(F1500W_raw[:,0],F1500W_raw[:,1])
for i in range(0,len(wv_1280W),1):
    F1280W[i,0]=wv_1280W[i]
    F1280W[i,1]=f_F1280W(wv_1280W[i])
for i in range(0,len(wv_1500W),1):
    F1500W[i,0]=wv_1500W[i]
    F1500W[i,1]=f_F1500W(wv_1500W[i])

l_eff_F1280 = 12.62 
l_eff_F1500 = 14.79 
flux_measured_12 = np.mean([3.425,3.428, 3.427, 3.428])
flux_measured_15 = 2.589 

# plt.errorbar(l_eff_F1280, flux_measured_12 ,yerr=0.10284+0.024,
#              xerr=[np.abs(l_eff_F1280-min(x_F1280W[np.where((y_F1280W>0.2))[0]]))],
#              fmt='h', color='k', capsize=10, markeredgewidth=2, markersize=10, zorder=5, alpha=0.5)#ecolor ='k')
# plt.errorbar(l_eff_F1500, flux_measured_15,#np.mean([2.532,2.528, 2.529,2.520,2.531,2.546,2.547, 2.599])
#                             yerr=0.079,xerr=[np.abs(l_eff_F1500-min(x_F1500W[np.where((y_F1500W>0.2))[0]]))], 
#                                                                     fmt='o', color='k', alpha=0.5,capsize=10, markeredgewidth=2, markersize=10,zorder=5,label="measurements" )#ecolor ='k')

# models 
## sphinx
t = np.linspace(1,19,1000)
# x_model_sphinx = spec_T1_sphinx["#Wavelength(micrometer)"]
# y_model_sphinx = np.array(flux_sphinx_Jy)*1e3 #/1.07
x_model_sphinx = wavelengths_T1_sphinx * 1e6
y_model_sphinx = conversion_IS_to_mJy(flux_T1_sphinx, x_model_sphinx, dist_system, R_star) #/1.07
f_model_sphinx = interp1d(x_model_sphinx,y_model_sphinx ) 
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
    norm_F1280W=norm_F1280W+F1500W[i,1]*(1./F1280W[i,0]-1./F1280W[i-1,0])
    cum1_F1280W+=spec1_F1280W[i]*F1280W[i,1]*(1./F1280W[i,0]-1./F1280W[i-1,0])

plt.scatter(l_eff_F1500,cum1_F1500W/norm_F1500W,#np.median(f_model_sphinx(t[np.where((f_F1500W(t)!=0))[0]])),
                marker='h',s=100,color='blue',edgecolors='k',linewidth=1.5 ,zorder=7,alpha=0.5)
plt.scatter(l_eff_F1280,cum1_F1280W/norm_F1280W,#np.median(f_model_sphinx(t[np.where((f_F1280W(t)!=0))[0]])),
                marker='h',s=100,color='blue',edgecolors='k',linewidth=1.5 ,zorder=7,alpha=0.5,label="Sphinx integrated MIRI band")

print("Difference in % between measurements and Sphinx at 12.8 :", (cum1_F1280W/norm_F1280W/flux_measured_12)*100 - 100)
print("Difference in % between measurements and Sphinx at 15 :", (cum1_F1500W/norm_F1500W/flux_measured_15)*100 - 100)

plt.xlabel("Wavelength (µm)")
plt.ylabel("Flux (mJy)")
plt.legend()
plt.show()