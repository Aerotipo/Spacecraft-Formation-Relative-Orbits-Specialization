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

# def compute_f_dot(f,n,e,t_step):
#     """Given the true anomaly f and the orbital parameters, this function calculates the rate f_dot
#     numerically using a time step defined by t_step.
#     Inputs:
#       f [float]: the true anomaly at time t, in radians
#       n [float]: the mean orbit rate n = sqrt(mu/a^3), in rad/s
#       e [float]: eccentricity
#       t_step [float]: time step for numerical derivative, in seconds
#     Output:
#       f_dot [float]: the true anomaly rate in rad/s
#       """
    
#     M_0 = trueanomaly_2_meananomaly(f,e)
#     M_1 = M_0 + n*t_step
#     f_1 = meananomaly_2_trueanomaly(M_1,e)
#     f_dot = (f_1 - f)/t_step

#     return f_dot


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


def Kepler_solver(mu, r_N_0, v_N_0, time_lapse, report = True):
    """ This function propagates the Keplerian orbit of a satellite during a time_lapse, from its initial position
    and velocity, to its final position and velocity, all in ECI cartesian coordinates.
    Inputs:
      mu [float]: is the gravitational constant expressed in SI units
      r_N_0 [ndarray]: is the initial position vector in inertial ECI frame coordinates, in [m]
      v_N_0 [ndarray]: is the initial velocity vector in inertial ECI frame coordinates, in [m]
      time_lapse [float]: is the elapsed time for the propagation, in seconds
      orbit_report (optional): if set to "true" it prints the values of the orbit elements
    Outputs:
      r_N_final [ndarray]: is the position vector in inertial ECI frame coordinates, in [m]
      v_N_final [ndarray]: is the velocity vector in inertial ECI frame coordinates, in [m]
    """
    sma, ecc, inc, AN, AP, f_initial = cartesian_2_orbitelement(mu, r_N_0, v_N_0)
    n = sqrt(mu/sma**3)
    M_initial = trueanomaly_2_meananomaly(f_initial, ecc)
    M_final = M_initial + sqrt(mu/sma**3)*time_lapse
    f_final = meananomaly_2_trueanomaly(M_final,ecc)
    r_N_final, v_N_final = orbitelement_2_cartesian(mu,sma,ecc,inc,AN,AP,f_final)

    if (report == True):
        print("Semi-major axis [km] = ", sma/1000)
        print("Eccentricity = ", ecc)
        print("Inclination [rad] = ", inc)
        print("Ascending Node [rad] = ", AN)
        print("Argument of Periapsis [rad] = ", AP)
        print("Initial true anomaly f [rad] = ", f_initial)
        print("chief orbital period: {} hours".format(2*pi/n/3600))
        print("chief final true anomaly f [rad] = ", f_final)
        print("\b")
    
    return r_N_final, v_N_final
    

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


def OE_2_LVLH_linear_mapping_matrix_A(mu, sma, theta, inc, q1, q2, AN):
    """ This function computes the [A_oe] matrix with the linear mapping from OE Relative Orbit Coordinates to
    Cartesian LVLH, that describe the deputy orbital elements relative to the chief orbital elements.
    Inputs:
      mu [float]: is the gravitational constant expressed in SI units
      sma [float]: the chief semi-major axis, in m;
      theta [float]: the chief total anomaly (theta = Argument of periapsis (AP) + true anomaly (f)) in rad;
      inc [float]: the chief orbit inclination, in rad;
      q1 [float]: q1 = ecc*cos(AP), where ecc is the eccentricity of the chief orbit, non-dimensional;
      q2 [float]: q2 = ecc*sin(AP), where ecc is the eccentricity of the chief orbit, non-dimensional;
      AN [float]: the chief Ascending Node (Ohmega)
    Outputs:
      A_oe [ndarray]: a 6x6 matrix mapping the two spaces as per the vectorial relation delta_X = [A_oe]*delta_oe
    """
    # Simplify parameters notation
    a = sma
    i = inc
    Ohmega = AN
    
    # Compute orbital parameters
    ecc = sqrt(q1**2 + q2**2)   #eccentricity
    if ecc == 0:
        print("Warning! Chief circular orbit can create singularities in w")
    AP = acos(q1/ecc)   # argument of periapsis (w)
    f = theta - AP    # true anomaly
    r = a*(1 - q1**2 - q2**2)/(1 + ecc*cos(f))
    p = a*(1 - q1**2 - q2**2)
    h = sqrt(mu*p)
    Vr = h/p*(q1*sin(theta) - q2*cos(theta))
    Vt = h/p*(1 + q1*cos(theta) + q2*sin(theta))

    # Compute Matrix elements
    A_oe = np.zeros([6,6])
    A_oe[0,0] = r/a
    A_oe[0,1] = Vr/Vt*r
    A_oe[0,3] = -r/p*(2*a*q1 + r*cos(theta))
    A_oe[0,4] = -r/p*(2*a*q2 + r*sin(theta))
    A_oe[1,1] = r
    A_oe[1,5] = r*cos(i)
    A_oe[2,2] = r*sin(theta)
    A_oe[2,5] = -r*cos(theta)*sin(i)
    A_oe[3,0] = -Vr/(2*a)
    A_oe[3,1] = (1/r - 1/p)*h
    A_oe[3,3] = (Vr*a*q1 + h*sin(theta))/p
    A_oe[3,4] = (Vr*a*q2 - h*cos(theta))/p
    A_oe[4,0] = -3*Vt/(2*a)
    A_oe[4,1] = -Vr
    A_oe[4,3] = (3*Vt*a*q1 + 2*h*cos(theta))/p
    A_oe[4,4] = (3*Vt*a*q2 + 2*h*sin(theta))/p
    A_oe[4,5] = Vr*cos(i)
    A_oe[5,2] = (Vt*cos(theta) + Vr*sin(theta))
    A_oe[5,5] = (Vt*sin(theta) - Vr*cos(theta))*sin(i)

    return A_oe


def LVLH_2_OE_linear_mapping_inversematrix_A(mu, sma, theta, inc, q1, q2, AN):
    """ This function computes the inverse inv[A_oe] matrix with the linear mapping from Cartesian LVLH to OE Relative 
    Orbit Coordinates, that describe the relative deputy orbital elements wrt to the chief orbital elements.
    Inputs:
      mu [float]: is the gravitational constant expressed in SI units
      sma [float]: the chief semi-major axis, in m;
      theta [float]: the chief total anomaly (theta = Argument of periapsis (AP) + true anomaly (f)) in rad;
      inc [float]: the chief orbit inclination, in rad;
      q1 [float]: q1 = ecc*cos(AP), where ecc is the eccentricity of the chief orbit, non-dimensional;
      q2 [float]: q2 = ecc*sin(AP), where ecc is the eccentricity of the chief orbit, non-dimensional;
      AN [float]: the chief Ascending Node (Ohmega)
    Outputs:
      inverse_A_oe [ndarray]: a 6x6 matrix mapping the two spaces as per the vectorial relation delta_oe = inv([A_oe])*delta_X
    """
    # Simplify parameters notation:
    a = sma
    i = inc
    Ohmega = AN
    
    # Compute orbital parameters:
    ecc = sqrt(q1**2 + q2**2)   #eccentricity
    AP = acos(q1/ecc)           # argument of periapsis (w)
    f = theta - AP              # true anomaly
    R = a*(1 - q1**2 - q2**2)/(1 + ecc*cos(f))
    p = a*(1 - q1**2 - q2**2)
    h = sqrt(mu*p)
    Vr = h/p*(q1*sin(theta) - q2*cos(theta))
    Vt = h/p*(1 + q1*cos(theta) + q2*sin(theta))

    # Compute non-dimensional matrix components:
    alpha = a/R
    vu = Vr/Vt
    rho = R/p
    k1 = alpha*(1/rho - 1)
    k2 = alpha*(1/rho)*vu**2

    # Compute Matrix elements:
    A_oe_inverse = np.zeros([6,6])
    A_oe_inverse[0,0] = 2*alpha*(2 + 3*k1 + 2*k2)
    A_oe_inverse[0,1] = -2*alpha*vu*(1 + 2*k1 + k2)
    A_oe_inverse[0,3] = 2*alpha**2*vu*p/Vt
    A_oe_inverse[0,4] = 2*a/Vt*(1 + 2*k1 + k2)
    A_oe_inverse[1,1] = 1/R
    A_oe_inverse[1,2] = 1/(tan(i)*R)*(cos(theta) + vu*sin(theta))
    A_oe_inverse[1,5] = -sin(theta)/(tan(i)*Vt)
    A_oe_inverse[2,2] = (sin(theta) - vu*cos(theta))/R
    A_oe_inverse[2,5] = cos(theta)/Vt
    A_oe_inverse[3,0] = 1/(rho*R)*(3*cos(theta) + 2*vu*sin(theta))
    A_oe_inverse[3,1] = -(1/R)*(vu**2/rho*sin(theta) + q1*sin(2*theta) - q2*cos(2*theta))
    A_oe_inverse[3,2] = -q2/(tan(i)*R)*(cos(theta) + vu*sin(theta))
    A_oe_inverse[3,3] = sin(theta)/(rho*Vt)
    A_oe_inverse[3,4] = (1/(rho*Vt))*(2*cos(theta) + vu*sin(theta))
    A_oe_inverse[3,5] = (q2*sin(theta))/(tan(i)*Vt)
    A_oe_inverse[4,0] = (1/(rho*R))*(3*sin(theta) - 2*vu*cos(theta))
    A_oe_inverse[4,1] = 1/R*(vu**2/rho*cos(theta) + q2*sin(2*theta) + q1*cos(2*theta))
    A_oe_inverse[4,2] = q1/(tan(i)*R)*(cos(theta) + vu*sin(theta))
    A_oe_inverse[4,3] = -cos(theta)/(rho*Vt)
    A_oe_inverse[4,4] = 1/(rho*Vt)*(2*sin(theta) - vu*cos(theta))
    A_oe_inverse[4,5] = -q1*sin(theta)/(tan(i)*Vt)
    A_oe_inverse[5,2] = -(cos(theta) + vu*sin(theta))/(R*sin(i))
    A_oe_inverse[5,5] = sin(theta)/(Vt*sin(i))

    return A_oe_inverse

