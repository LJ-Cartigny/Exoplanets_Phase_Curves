#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# January 2025
# Phase curves

import numpy as np
import matplotlib.pyplot as plt
from Orbital_motion import compute_true_anomaly
from Transits import eclipse, eclipse_impact_parameter
from TRAPPIST1_parameters import *

def phase_angle(omega, nu, i):
    """
    Determines the phase angle of a planet from its orbital parameters (in rad).

    :param omega: the argument of pericentre (in rad)
    :type omega: float

    :param nu: the true anomaly (in rad)
    :type nu: float

    :param i: the inclination (in rad)
    :type i: float

    :return: alpha
    :rtype: float
    """

    alpha = np.arccos(np.sin(omega+nu)*np.sin(i))
    return alpha

def phase_function(alpha):
    """
    Determines the phase function of a Lambert sphere.

    :param alpha: the phase angle (in rad)
    :type alpha: float

    :return: g
    :rtype: float
    """

    g = (np.sin(alpha)+(np.pi-alpha)*np.cos(alpha))/np.pi
    return g

def star_planet_separation(a,e,nu):
    """
    Determines the distance between a planet and its star using its orbital parameters.

    :param a: the semimajor axis (in m)
    :type a: float

    :param e: the eccentricity
    :type e: float

    :param nu: the true anomaly (in rad)
    :type nu: float

    :return: r
    :rtype: float
    """

    r = (a*(1-e**2))/(1+e*np.cos(nu))
    return r

def flux_star(L,d):
    """
    Determines the flux received from a star (in W/m^2) at a distance d.

    :param L: the star luminosity (in W)
    :type L: float

    :param d: the distance (in m)
    :type d: float

    :return: F
    :rtype: float
    """

    F = L/(4*np.pi*d**2)
    return F

def flux_planet(F_star):
    """
    Determines the flux reemitted by a planet (in W/m^2) from the one it receives from its star considering the planet is a black body.

    :param F_star: the flux received by the planet from its star (in W/m^2)
    :type F_star: float

    :return: F_planet
    :rtype: float
    """

    F_planet = F_star*2/3
    return F_planet

def surface_sphere(R):
    """
    Determines the surface of a sphere of radius R.

    :param R: the radius (in m)
    :type R: float

    :return: S
    :rtype: float
    """

    S = 4*np.pi*R**2
    return S

def luminosity_planet_dayside(F_planet,R_planet):
    """
    Determines the luminosity of the dayside of a planet from the flux it reemits and its radius.

    :param F_planet: the flux reemitted by the planet's dayside (in W/m^2)
    :type F_planet: float

    :param R_planet: the planet radius (in m)
    :type R_planet: float

    :return: L_planet
    :rtype: float
    """

    L_planet = F_planet * surface_sphere(R_planet)/2
    return L_planet

def phase_curve(L_star, L_planet, R_star, R_planet, phase_planet, eclipse):
    """
    Determines the phase curve of a planet from its luminosity, its star's luminosity and its phase function expressed as the ratio between the planet and star's luminosities in ppm.

    :param L_star: the star luminosity (in W)
    :type L_star: float

    :param L_planet: the planet luminosity (in W)
    :type L_planet: float

    :param R_star: the star radius (in m)
    :type R_star: float

    :param R_planet: the planet radius (in m)
    :type R_planet: float

    :param phase_planet: the phase function of the planet
    :type phase_planet: float

    :param eclipse: True if the planet is in eclipse, False otherwise
    :type eclipse: bool

    :return: curve
    :rtype: float
    """

    curve = L_planet/L_star*phase_planet/(R_planet/R_star)**2*10**6 * (-1*eclipse+1)
    return curve


def main():
    
    t_end = 30 # simulation duration in days
    nb_points = 10000 # number of points in the time array

    t = np.linspace(0,t_end,nb_points) # time array in days


    # For TRAPPIST-1 b
    
    nu_b = compute_true_anomaly(0,e_b,P_b,t)

    alpha_b = phase_angle(omega_b,nu_b,i_b)
    phase_b = phase_function(alpha_b)
    b_b = eclipse_impact_parameter(a_b,i_b,e_b,R_star,omega_b)
    eclipse_b = eclipse(P_b,a_b,R_star,R_b,i_b,phase_b,e_b,omega_b,b_b)

    r_b = star_planet_separation(a_b,e_b,nu_b)

    flux_star_b = flux_star(L_star,r_b)
    flux_b = flux_planet(flux_star_b)
    L_b = luminosity_planet_dayside(flux_b,R_b)

    phase_curve_b = phase_curve(L_star,L_b,R_star,R_b,phase_b,eclipse_b)
    # np.savetxt("Phase_curve_v1_output/phase_curve_b_ter.txt",np.concatenate((t.reshape(nb_points,1),phase_curve_b.reshape(nb_points,1)),axis=1))


    #For TRAPPIST-1 c

    nu_c = compute_true_anomaly(0,e_c,P_c,t)
    alpha_c = phase_angle(omega_c,nu_c,i_c)
    phase_c = phase_function(alpha_c)
    b_c = eclipse_impact_parameter(a_c,i_c,e_c,R_star,omega_c)
    eclipse_c = eclipse(P_c,a_c,R_star,R_c,i_c,phase_c,e_c,omega_c,b_c)

    r_c = star_planet_separation(a_c,e_c,nu_c)

    flux_star_c = flux_star(L_star,r_c)
    flux_c = flux_planet(flux_star_c)
    L_c = luminosity_planet_dayside(flux_c,R_c)

    phase_curve_c = phase_curve(L_star,L_c,R_star,R_c,phase_c,eclipse_c)
    #np.savetxt("Phase_curve_v1_output/phase_curve_c_ter.txt",np.concatenate((t.reshape(nb_points,1),phase_curve_c.reshape(nb_points,1)),axis=1))

    
    #For TRAPPIST-1 d

    nu_d = compute_true_anomaly(0,e_d,P_d,t)
    alpha_d = phase_angle(omega_d,nu_d,i_d)
    phase_d = phase_function(alpha_d)
    b_d = eclipse_impact_parameter(a_d,i_d,e_d,R_star,omega_d)
    eclipse_d = eclipse(P_d,a_d,R_star,R_d,i_d,phase_d, e_d, omega_d, b_d)

    r_d = star_planet_separation(a_d,e_d,nu_d)

    flux_star_d = flux_star(L_star,r_d)
    flux_d = flux_planet(flux_star_d)
    L_d = luminosity_planet_dayside(flux_d,R_d)

    phase_curve_d = phase_curve(L_star,L_d,R_star,R_d,phase_d,eclipse_d)
    #np.savetxt("Phase_curve_v1_output/phase_curve_d.txt",np.concatenate((t.reshape(nb_points,1),phase_curve_d.reshape(nb_points,1)),axis=1))


    #For TRAPPIST-1 e

    nu_e = compute_true_anomaly(0,e_e,P_e,t)
    alpha_e = phase_angle(omega_e,nu_e,i_e)
    phase_e = phase_function(alpha_e)
    b_e = eclipse_impact_parameter(a_e,i_e,e_e,R_star,omega_e)
    eclipse_e = eclipse(P_e,a_e,R_star,R_e,i_e,phase_e, e_e, omega_e, b_e)

    r_e = star_planet_separation(a_e,e_e,nu_e)

    flux_star_e = flux_star(L_star,r_e)
    flux_e = flux_planet(flux_star_e)
    L_e = luminosity_planet_dayside(flux_e,R_e)

    phase_curve_e = phase_curve(L_star,L_e,R_star,R_e,phase_e,eclipse_e)
    #np.savetxt("Phase_curve_v1_output/phase_curve_e.txt",np.concatenate((t.reshape(nb_points,1),phase_curve_e.reshape(nb_points,1)),axis=1))


    #For TRAPPIST-1 f

    nu_f = compute_true_anomaly(0,e_f,P_f,t)
    alpha_f = phase_angle(omega_f,nu_f,i_f)
    phase_f = phase_function(alpha_f)
    b_f = eclipse_impact_parameter(a_f,i_f,e_f,R_star,omega_f)
    eclipse_f = eclipse(P_f,a_f,R_star,R_f,i_f,phase_f, e_f, omega_f, b_f)

    r_f = star_planet_separation(a_f,e_f,nu_f)

    flux_star_f = flux_star(L_star,r_f)
    flux_f = flux_planet(flux_star_f)
    L_f = luminosity_planet_dayside(flux_f,R_f)

    phase_curve_f = phase_curve(L_star,L_f,R_star,R_f,phase_f,eclipse_f)
    #np.savetxt("Phase_curve_v1_output/phase_curve_f.txt",np.concatenate((t.reshape(nb_points,1),phase_curve_f.reshape(nb_points,1)),axis=1))


    #For TRAPPIST-1 g

    nu_g = compute_true_anomaly(0,e_g,P_g,t)
    alpha_g = phase_angle(omega_g,nu_g,i_g)
    phase_g = phase_function(alpha_g)
    b_g = eclipse_impact_parameter(a_g,i_g,e_g,R_star,omega_g)
    eclipse_g = eclipse(P_g,a_g,R_star,R_g,i_g,phase_g, e_g, omega_g, b_g)

    r_g = star_planet_separation(a_g,e_g,nu_g)

    flux_star_g = flux_star(L_star,r_g)
    flux_g = flux_planet(flux_star_g)
    L_g = luminosity_planet_dayside(flux_g,R_g)

    phase_curve_g = phase_curve(L_star,L_g,R_star,R_g,phase_g,eclipse_g)
    #np.savetxt("Phase_curve_v1_output/phase_curve_g.txt",np.concatenate((t.reshape(nb_points,1),phase_curve_g.reshape(nb_points,1)),axis=1))


    #For TRAPPIST-1 h

    nu_h = compute_true_anomaly(0,e_h,P_h,t)
    alpha_h = phase_angle(omega_h,nu_h,i_h)
    phase_h = phase_function(alpha_h)
    b_h = eclipse_impact_parameter(a_h,i_h,e_h,R_star,omega_h)
    eclipse_h = eclipse(P_h,a_h,R_star,R_h,i_h,phase_h, e_h, omega_h, b_h)
    
    r_h = star_planet_separation(a_h,e_h,nu_h)

    flux_star_h = flux_star(L_star,r_h)
    flux_h = flux_planet(flux_star_h)
    L_h = luminosity_planet_dayside(flux_h,R_h)

    phase_curve_h = phase_curve(L_star,L_h,R_star,R_h,phase_h,eclipse_h)
    #np.savetxt("Phase_curve_v1_output/phase_curve_h.txt",np.concatenate((t.reshape(nb_points,1),phase_curve_h.reshape(nb_points,1)),axis=1))


    phase_curve_total = phase_curve_b + phase_curve_c + phase_curve_d + phase_curve_e + phase_curve_f + phase_curve_g + phase_curve_h
    #np.savetxt("Phase_curve_v1_output/phase_curve_total.txt",np.concatenate((t.reshape(nb_points,1),phase_curve_total.reshape(nb_points,1)),axis=1))



    # Plot

    plt.figure()
    plt.plot(t,phase_curve_b,label="b")
    plt.plot(t,phase_curve_c,label="c")
    # plt.plot(t,phase_curve_d,label="d")
    # plt.plot(t,phase_curve_e,label="e")
    # plt.plot(t,phase_curve_f,label="f")
    # plt.plot(t,phase_curve_g,label="g")
    # plt.plot(t,phase_curve_h,label="h")
    # plt.plot(t,phase_curve_total,label="Total")
    plt.xlabel("Time (days)")
    plt.ylabel("$L_{planet}/L_{star}$ (ppm)")
    plt.title("Phase curves of planets of TRAPPIST-1 as black bodies")
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()