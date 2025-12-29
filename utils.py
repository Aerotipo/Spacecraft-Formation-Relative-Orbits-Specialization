"""
Docstring for Spacecraft-Formation-Relative-Orbits-Specialization.utils
This file contains the utility functions used in the following specialization:
    Title: Spacecraft Fformation Relative Orbits Specialization
    Institution: University of Colorado Boulder (on Coursera)
    Period: December 2025 - March 2026
All definitions, terminologies, symbols and conventions used here follow the courses conventions.

Created on December 2025
@author: marcos perez
"""

import numpy as np
from math import sin, cos, tan, asin, acos, atan, atan2, sqrt, pi

def trueanomaly_2_meananomaly(f,e):
    """Given the true anomaly f this function calculates the mean anomaly M using Kepler's equations.
    Inputs:
      f [float]: true anomaly in radians
      e [float]: eccentricity
    Output:
      M [float]: mean anomaly in radians
    """
    
    E = 2*atan(sqrt((1-e)/(1+e))*tan(f/2))

    return E - e*sin(E)


def meananomaly_2_trueanomaly(M,e):
    """Given the mean anomaly M, this function calculates the true anomaly f using Kepler's equations and
    a Newton-Raphson solver.
    Inputs:
      M [float]: mean anomaly in radians
      e [float]: eccentricity
    Output:
      f [float]: true anomaly in radians
    """
    #**********Newton-Raphson method
    eps = 1e-7    #residual tolerance
    E = M         #seed value (recommended)
    Residual = M - E +e*sin(E)  #residual to be minimized

    while abs(Residual) > eps:
        dResidual = e*cos(E) - 1
        E = E - Residual/dResidual
        Residual = M - E + e*sin(E)
    #**********
        
    return 2*atan(sqrt((1+e)/(1-e))*tan(E/2))

def f_dot(f,n,e,t_step):
    """Given the true anomaly f and the orbital parameters, this function calculates the rate f_dot
    numerically using a time step defined by t_step.
    Inputs:
      f [float]: the true anomaly at time t, in radians
      n [float]: the mean orbit rate n = sqrt(mu/a^3), in rad/s
      e [float]: eccentricity
      t_step [float]: time step for numerical derivative, in seconds
    Output:
      f_dot [float]: the true anomaly rate in rad/s
      """
    
    M_0 = trueanomaly_2_meananomaly(f,e)
    M_1 = M_0 + n*t_step
    f_1 = meananomaly_2_trueanomaly(M_1,e)
    f_dot = (f_1 - f)/t_step

    return f_dot


def Euler_313(angle1,angle2,angle3):
    """This program computes the matrix BN corresponding to a 3-1-3 Euler transformation.
    Inputs:
        angle1 [float]: the 1st rotation angle around axis 3, in degrees
        angle2 [float]: the 2nd rotation angle around axis 1, in degrees
        angle3 [float]: the 3rd rotation angle around axis 3, in degrees
    """
    a1 = angle1*pi/180
    a2 = angle2*pi/180
    a3 = angle3*pi/180
    
    C11 = cos(a3)*cos(a1)-sin(a3)*cos(a2)*sin(a1)
    C21 = -sin(a3)*cos(a1)-cos(a3)*cos(a2)*sin(a1)
    C31 = sin(a2)*sin(a1)
    C12 = cos(a3)*sin(a1)+sin(a3)*cos(a2)*cos(a1)
    C22 = -sin(a3)*sin(a1)+cos(a3)*cos(a2)*cos(a1)
    C32 = -sin(a2)*cos(a1)
    C13 = sin(a3)*sin(a2)
    C23 = cos(a3)*sin(a2)
    C33 = cos(a2)
    
    return np.array([[C11,C12,C13],
                     [C21,C22,C23],
                     [C31,C32,C33]])

def inverse_Euler_313(C):
    """
    Computes the Ohm, i, w values corresponding to an Euler 3-1-3 DCM given
    by the ndarray C (input).
    Input:
        C [ndarray]: 3x3 DCM for 3-1-3 Euler rotation
    Output: the rotation angles 
        Ohm (ascending node), in degrees 
        i (inclination), in degrees
        w (argument of perigree), in degrees
    """
    Ohm = atan2(C[2,0],-C[2,1])*180/pi
    i = acos(C[2,2])*180/pi
    w = atan2(C[0,2],C[1,2])*180/pi
    
    return Ohm, i, w

def orbitelement_2_cartesian(mu,a,e,i,ohm,w,f):
    """This function performs a coordinate transformation from the invariant orbit elements (elliptical orbit) and the true anomaly 
    to the equivalent position and velocity vector states in cartesian coordinates, expressed in ECI frame.
    Inputs:
      mu [float]: gravitational constant in SI units
      a [float]: semi-major axis of elliptical orbit in [m]
      e [float]: eccentricity
      i [float]: orbit inclination angle [rad]
      ohm [float]: ascending node [rad]
      w [float]: argument of periapses [rad]
      f [float]: initial true anomaly [rad]
    Outputs:
      r [ndarray]: the cartesian position components, expressed in ECI frame (inertial), in [m]
      v [ndarray]: the cartesian velocity components, expressed in ECI frame (inertial), in [m]
      """
    
    n = sqrt(mu/a**3)
    theta = w + f
    p = a*(1-e)*(1+e)
    r_norm = p/(1 + e*cos(f))
    fdot = f_dot(f,n,e,t_step=0.001)
    h = fdot*r_norm**2
    
    r_ix = cos(ohm)*cos(theta) - sin(ohm)*sin(theta)*cos(i)
    r_iy = sin(ohm)*cos(theta) + cos(ohm)*sin(theta)*cos(i)
    r_iz = sin(theta)*sin(i)
    r = r_norm*np.array([r_ix, r_iy, r_iz])

    v_ix = cos(ohm)*( sin(theta) + e*sin(w) ) + sin(ohm)*( cos(theta) + e*cos(w) )*cos(i)
    v_iy = sin(ohm)*( sin(theta) + e*sin(w) ) - cos(ohm)*( cos(theta) + e*cos(w) )*cos(i)
    v_iz = -( cos(theta) + e*cos(w) )*sin(i)
    v = -(mu/h)*np.array([v_ix, v_iy, v_iz])
    print(np.linalg.norm(np.cross(r,v)))

    return r, v
    

def cartesian_2_orbitelement(mu, rVEC, vVEC):
    """This function performs a coordinate transformation from cartesian ECI (intertial) to orbit elements.
    Inputs:
      mu [float]: is the gravitational constant expressed in SI units
      rVEC [ndarray]: is the position vector in inertial ECI frame coordinates
      vVEC [ndarray]: is the velocity vector in inertial ECI frame coordinates
    Outputs:
      sma [float]: semi-major axis of elliptical orbit in [m]
      ecc [float]: eccentricity
      inc [float]: orbit inclination angle [rad]
      AN  [float]: ascending node [rad]
      AP  [float]: argument of periapses [rad]
      f   [float]: true anomaly [rad]
    """

    r = np.linalg.norm(rVEC)
    vv = np.dot(vVEC,vVEC)
    hVEC = np.cross(rVEC,vVEC)
    h = np.linalg.norm(hVEC)
    sma = 1/(2/r - vv/mu)
    eVEC = np.cross(vVEC,hVEC)/mu - rVEC/r
    ecc = np.linalg.norm(eVEC)
    i_e = eVEC/ecc
    i_h = hVEC/h
    i_p = np.cross(i_h, i_e)
    i_r = rVEC/r
    PN = np.array([i_e,i_p,i_h])
    AN, inc, AP = inverse_Euler_313(PN)
    f = atan2(np.dot(np.cross(i_e,i_r),i_h),np.dot(i_e,i_r))
    
    return sma, ecc, inc*pi/180, AN*pi/180, AP*pi/180, f

