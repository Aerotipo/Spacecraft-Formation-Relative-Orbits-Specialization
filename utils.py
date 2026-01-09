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

def a_J2(r, mu, J2, r_eq):
    """This function generates the disturbance acceleration a_J2 produced by J2
    as a function of the position r in cartesian coordinates.
    Input: 
      r = np.array(x,y,z) --> satellite position
    Output: 
      a = np.array(a_x, a_y, a_z)"""
    x = r[0]
    y = r[1]
    z = r[2]
    rnorm = np.linalg.norm(r)
    ax = (1 - 5*(z/rnorm)**2)*(x/rnorm)
    ay = (1 - 5*(z/rnorm)**2)*(y/rnorm)
    az = (3 - 5*(z/rnorm)**2)*(z/rnorm)

    return (-3/2*J2*(mu/rnorm**2)*(r_eq/rnorm)**2)*np.array([ax,ay,az])


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
    eps = 1e-9    #residual tolerance
    E = M         #seed value (recommended)
    Residual = M - E +e*sin(E)  #residual to be minimized

    while abs(Residual) > eps:
        dResidual = e*cos(E) - 1
        E = E - Residual/dResidual
        Residual = M - E + e*sin(E)
    #**********
        
    return 2*atan(sqrt((1+e)/(1-e))*tan(E/2))

def compute_f_dot(f,n,e,t_step):
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
    p = a*(1-e**2)
    r_norm = p/(1 + e*cos(f))
    h = sqrt(mu*p)
        
    r_ix = cos(ohm)*cos(theta) - sin(ohm)*sin(theta)*cos(i)
    r_iy = sin(ohm)*cos(theta) + cos(ohm)*sin(theta)*cos(i)
    r_iz = sin(theta)*sin(i)
    r = r_norm*np.array([r_ix, r_iy, r_iz])

    v_ix = cos(ohm)*( sin(theta) + e*sin(w) ) + sin(ohm)*( cos(theta) + e*cos(w) )*cos(i)
    v_iy = sin(ohm)*( sin(theta) + e*sin(w) ) - cos(ohm)*( cos(theta) + e*cos(w) )*cos(i)
    v_iz = -( cos(theta) + e*cos(w) )*sin(i)
    v = -(mu/h)*np.array([v_ix, v_iy, v_iz])
    
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
    # f = atan(np.dot(np.cross(i_e,i_r),i_h)/np.dot(i_e,i_r))
    
    return sma, ecc, inc*pi/180, AN*pi/180, AP*pi/180, f


def Inertial_2_Hill_mapping(rc_N,vc_N,rd_N,vd_N):
    """This function maps the inertial (ECI) position and velocity of the chief and deputy satellites, 
    to the Orbital (LVLH) or relative Hills frame (centred in the chief).
    Inputs:
      rc_N [ndarray]: position vector of the chief satellite expressed in ECI cartesian coordinates, in [m];
      vc_N [ndarray]: velocity vector of the chief satellite expressed in ECI cartesian coordinates, [m/s];
      rd_N [ndarray]: position vector of the deputy satellite expressed in ECI cartesian coordinates, in [m];
      vd_N [ndarray]: velocity vector of the deputy satellite expressed in ECI cartesian coordinates, [m/s];
    Outputs:
      rho_H [ndarray]: relative position vector of deputy as seen by chief in Hills coordinate frame, in [m];
      rhoP_H [ndarray]: relative velocity vector (rho prime) of deputy as seen by chief in Hills coordinate frame, in [m/s]
    """
    h_N = np.cross(rc_N,vc_N)
    h_norm = np.linalg.norm(h_N)
    rc_norm = np.linalg.norm(rc_N)
    o_r = rc_N/rc_norm
    o_h = h_N/h_norm
    o_t = np.cross(o_h, o_r)
    ON = np.array([o_r, o_t, o_h])
    w_ON_H = (h_norm/rc_norm**2)*o_h

    rho_H = np.dot(ON, rd_N - rc_N)
    rhoP_H = np.dot(ON, vd_N - vc_N) - np.cross(w_ON_H, rho_H)

    return rho_H, rhoP_H



def Hill_2_Inertial_mapping(rc_N,vc_N,rho_H,rhoP_H):
    """This function maps the Orbital (LVLH) or relative Hills frame position and relative velocity of the deputy satellite
    with respect to chief, to the inertial (ECI) cartesian position and velocity vectors of the deputy.
    Inputs:
      rc_N [ndarray]: position vector of the chief satellite expressed in ECI cartesian coordinates, in [m];
      vc_N [ndarray]: velocity vector of the chief satellite expressed in ECI cartesian coordinates, [m/s];
      rho_H [ndarray]: relative position vector of deputy as seen by chief in Hills coordinate frame, in [m];
      rhoP_H [ndarray]: relative velocity vector (rho prime) of deputy as seen by chief in Hills coordinate frame, in [m/s]
    Outputs:
      rd_N [ndarray]: position vector of the deputy satellite expressed in ECI cartesian coordinates, in [m];
      vd_N [ndarray]: velocity vector of the deputy satellite expressed in ECI cartesian coordinates, [m/s];
    """
    h_N = np.cross(rc_N,vc_N)
    h_norm = np.linalg.norm(h_N)
    rc_norm = np.linalg.norm(rc_N)
    o_r = rc_N/rc_norm
    o_h = h_N/h_norm
    o_t = np.cross(o_h, o_r)
    NO = np.transpose(np.array([o_r, o_t, o_h]))
    w_ON_H = (h_norm/rc_norm**2)*o_h

    rd_N = np.dot(NO, rho_H) + rc_N
    vd_N = vc_N + np.dot(NO, rhoP_H + np.cross(w_ON_H,rho_H))

    return rd_N, vd_N 


def Rectilinear_2_Curvilinear_mapping(x,y,r,x_dot,y_dot, r_dot):
    """This function converts Hill-frame rectilinear coordinates into Hill-frame cylindrical coordinates
    for the linearized relative equations of motion (CWH equations)"""
    rd = sqrt((r + x)**2 + y**2)
    delta_r = rd - r
    s = r*atan(y/(r+x))
    rd_dot = 1/(sqrt((x + r)**2 + y**2))*((x + r)*(x_dot + r_dot) + y*y_dot )
    delta_r_dot = rd_dot - r_dot
    aux = y/(r+x)
    s_dot = r_dot*atan(aux) + ( r/(1 + (aux)**2) )*( y_dot/(r+x) - y*(r_dot + x_dot)/(r+x)**2 )
    
    return delta_r, delta_r_dot, s, s_dot


def Curvilinear_2_Rectilinear_mapping(delta_r, s, r, delta_r_dot, s_dot, r_dot):
    """This function converts Hill-frame curvilinear coordinates into Hill-frame rectilinear coordinates
    for the linearized relative equations of motion (CWH equations)"""
    
    delta_theta = s/r
    delta_theta_dot = s_dot/r - s*r_dot/r**2
    rd = r + delta_r
    rd_dot = r_dot + delta_r_dot
    x = rd*cos(delta_theta) - r
    y = rd*sin(delta_theta)
    x_dot = rd_dot*cos(delta_theta) - rd*sin(delta_theta)*delta_theta_dot - r_dot
    y_dot = rd_dot*sin(delta_theta) + rd*cos(delta_theta)*delta_theta_dot

    return x, x_dot, y, y_dot

