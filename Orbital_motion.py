#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Louis-Julien Cartigny
# January 2025
# Orbital motion

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def kepler_equation(E,M,e):
    """
    Returns the Kepler equation: M = E - e sin(E)

    :param E: the eccentric anomaly
    :type E: float

    :param M: the mean anomaly
    :type M: float

    :param e: the eccentricity
    :type e: float

    :return: E - e*np.sin(E) - M
    :rtype: float
    """

    return E - e*np.sin(E) - M


def solve_kepler(M,e):
    """
    Solve the Kepler equation to find the eccentric anomaly E.

    :param E: the eccentric anomaly
    :type E: float

    :param M: the mean anomaly
    :type M: float

    :return: E
    :rtype: float
    """
    
    M = np.atleast_1d(M)
    E_ini = M # initial approximation
    E_solution = np.array([fsolve(kepler_equation, E0, args=(M_i,e))[0] for M_i, E0 in zip(M, E_ini)])

    return E_solution


def true_anomaly(E,e):
    """
    Computes the true anomaly from the eccentric anomaly and the eccentricity.

    :param E: the eccentric anomaly
    :type E: float

    :param e: the eccentricity
    :type e: float

    :return: nu
    :rtype: float
    """

    nu = 2 * np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    return nu


def compute_true_anomaly(nu_0, e, T, t, t_0=0):
    """
    Computes the true anomaly with respect to time.

    :param nu_0: the initial true anomaly (in rad)
    :type nu_0: float

    :param e: the eccentricity
    :type e: float

    :param T: the orbital period
    :type T: float

    :param t: the time passed
    :type t: float

    :param t_0: the initial time (default value: 0)
    :type t_0: float

    :return: nu_t
    :rtype: float
    """

    # Step 1: Computing the initial eccentric anomaly E_0
    E_0 = 2 * np.arctan2(np.sqrt(1-e)*np.sin(nu_0/2), np.sqrt(1+e)*np.cos(nu_0/2))

    # Step 2: Computing the initial mean anomaly
    M_0 = E_0 - e * np.sin(E_0)

    # Step 3: Computing the mean anomaly at time t
    n = 2*np.pi / T # mean motion
    M_t = M_0 + n*(t-t_0)

    # Step 4: Solving the Kepler equation for E_t
    E_t = solve_kepler(M_t,e)

    # Step 5: Computing the mean anomaly nu_t
    nu_t = true_anomaly(E_t,e)

    return nu_t



def main():
    # Example with Earth
    nu_0 = np.radians(30) # initial true anomaly of 30° in rad
    e = 0.07
    T = 365.25 # days
    t = np.linspace(0,500,500) # days

    nu_t = compute_true_anomaly(nu_0,e,T,t)
    #print(np.degrees(nu_t))

    print(nu_t.shape)
    print(t.shape)

    plt.figure()
    plt.plot(np.degrees(nu_t),t)
    plt.xlabel("$t$ (days)")
    plt.ylabel("True anomaly (°)")
    plt.show()




if __name__ == "__main__":
    main()