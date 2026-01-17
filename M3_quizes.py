import numpy as np 
import matplotlib.pyplot as plt
from math import sin, cos, tan, atan, pi, sqrt
from utils import *

# #===================================
# #quiz 3 - Relative Motion Mapping
# #===================================
# # exercise 1
# #===================================
# rc_N = np.array([-4893.268, 3864.478, 3169.646])*1000     #[m] chief satellite position in ECI
# vc_N = np.array([-3.91386, -6.257673, 1.59797])*1000      #[m/s] chief satellite velocity in ECI
# rd_N = np.array([-4892.98, 3863.073, 3170.619])*1000     #[m] chief satellite position in ECI
# # vd_N = np.array([-3.913302, -6.258661, 1.598199])*1000      #[m/s] chief satellite velocity in ECI
# vd_N = np.array([-3.912794, -6.258249, 1.598643])*1000      #[m/s] chief satellite velocity in ECI --> correct input

# rho_H, rhoP_H = Inertial_2_Hill_mapping(rc_N,vc_N,rd_N,vd_N)
# print("exercise 1 results")
# print("rho_H = ", rho_H/1000," km" )
# print("rhoP_H = ", rhoP_H/1000, " km/s")
# print("\b")

# # rd_N_check, vd_N_check = Hill_2_Inertial_mapping(rc_N,vc_N,rho_H,rhoP_H)

# # print("health check")
# # print("rd_N (check) = ", rd_N_check/1000," km" )
# # print("vd_N (check) = ", vd_N_check/1000, " km/s")
# # print("\b")

# #===================================
# # exercise 2
# #===================================
# rho_H = np.array([-0.537, 1.221, 1.106])*1000           #[m] deputy satellite relative pos in H-frame
# rhoP_H = np.array([0.000486, 0.001158, 0.0005590])*1000 #[m/s] deputy satellite relative vel in H-frame

# rd_N, vd_N = Hill_2_Inertial_mapping(rc_N,vc_N,rho_H,rhoP_H)

# print("exercise 2 results\b")
# print("rd_N = ", rd_N/1000," km" )
# print("vd_N = ", vd_N/1000, " km/s")
# print("\b")

# # rho_H_check, rhoP_H_check = Inertial_2_Hill_mapping(rc_N,vc_N,rd_N,vd_N)

# # print("health check exercise 2")
# # print("rho_H (check) = ", rho_H_check/1000," km" )
# # print("rhoP_H (check) = ", rhoP_H_check/1000, " km/s")

# #===================================
# # exercise 3
# #===================================
# x = 10          #km
# y = 500         #km
# r = 7000        #km
# x_dot = 0.1     #km/s
# y_dot = -0.1    #km/s
# r_dot = 0.05    #km/s

# delta_r, delta_r_dot, s, s_dot = Rectilinear_2_Curvilinear_mapping(x,y,r,x_dot,y_dot, r_dot)

# print("exercise 3 results")
# print("delta_r = ", delta_r)
# print("delta_r_dot = ", delta_r_dot)
# print("s = ", s)
# print("s_dot = ", s_dot)
# print("\b")

#===================================
# exercise 4
#===================================
# dr = 10         #km
# s =  500        #km
# r = 7000        #km
# dr_dot = 0.1    #km/s
# s_dot = -0.1    #km/s
# r_dot = 0.05    #km/s

# x, x_dot, y, y_dot = Curvilinear_2_Rectilinear_mapping(dr, s, r, dr_dot, s_dot, r_dot)

# print("exercise 4 results")
# print("x = ", x)
# print("x_dot = ", x_dot)
# print("y = ", y)
# print("y_dot = ", y_dot)
# print("\b")

# #==============================================
# #quiz 4 - General Relative Equations of Motion
# #==============================================
# # exercise 2 - Task 1
# #===================================
# # mu = 398600441.8e6   #Earth geocentric gravitational constant
# mu = 398600000.0e6   #Earth geocentric gravitational constant

# #Initial conditions of chief in ECI frame (N):
# rc_N_0 = np.array([-6685.20926, 601.51244, 3346.066634])*1000.   #m
# vc_N_0 = np.array([-1.74294, -6.70242, -2.27739])*1000.          #m/s
# #Initial conditions of deputy in Hill frame:
# rho_H_0 = np.array([-81.22301, 248.14201, 94.95904])*1000.     #[m] deputy satellite relative pos in H-frame
# rhoP_H_0 = np.array([0.47884, 0.14857, 0.13577])*1000.        #[m/s] deputy satellite relative vel in H-frame

# #Mapping deputy rd, vd from Hill frame to inertial:

# rd_N_0, vd_N_0 = Hill_2_Inertial_mapping(rc_N_0,vc_N_0,rho_H_0,rhoP_H_0)

# print("Results:")
# print("initial_rd_N = ", rd_N_0/1000.," km" )
# print("initial_vd_N = ", vd_N_0/1000., " km/s")
# print("\b")

# #Time domain
# time_lapse = 20   #seconds
# dt = 0.001          #time step
# t = np.arange(0, time_lapse + dt, dt)

# # #METHOD 1A --> Using cartesian ECI coordinates and a [first order] ODE integrator
# print("---------------------------------------------------------------------------------------")
# print("**** METHOD 1-A: Using cartesian ECI coordinates and a [first order] ODE integrator ****")
# print("\b")

# #history arrays creation
# rc_history = np.zeros([len(t)+1,3])  #chief satellite position history (trajectory) in inertial frame
# vc_history = np.zeros([len(t)+1,3])  #chief satellite velocity history in inertial frame
# rd_history = np.zeros([len(t)+1,3])  #deputy satellite position history (trajectory) in inertial frame
# vd_history = np.zeros([len(t)+1,3])  #deputy satellite velocity history in inertial frame

# #Main loop
# rc = rc_N_0.copy()
# vc = vc_N_0.copy()
# rd = rd_N_0.copy()
# vd = vd_N_0.copy()
# rc_history[0,:] = rc/1000.   #km
# vc_history[0,:] = vc/1000.   #km/s
# rd_history[0,:] = rd/1000.   #km
# vd_history[0,:] = vd/1000.   #km/s

# for i,ti in enumerate(t):
#     #chief integration of inertial equations of motion
#     rc_norm = np.linalg.norm(rc)
#     rc_dot_dot = -(mu/rc_norm**3)*rc
#     vc += rc_dot_dot*dt
#     rc += vc*dt
#     rc_history[i+1,:]=rc/1000.
#     vc_history[i+1,:]=vc/1000.
    
#     #deputy satellite integration of inertial equations of motion
#     rd_norm = np.linalg.norm(rd)
#     rd_dot_dot = -(mu/rd_norm**3)*rd
#     vd += rd_dot_dot*dt
#     rd += vd*dt
#     rd_history[i+1,:]=rd/1000.
#     vd_history[i+1,:]=vd/1000.
    
# print("chief position [in km] at {} seconds: rc_N = ".format(t[-1]),rc_history[-1])
# print("chief velocity [in km/s] at {} seconds: vc_N = ".format(t[-1]),vc_history[-1])
# print("deputy position [in km] at {} seconds: rd_N = ".format(t[-1]),rd_history[-1])
# print("deputy velocity [in km/s] at {} seconds: vd_N = ".format(t[-1]),vd_history[-1])
# print("\b")

# #Mapping back deputy position and velocity from inertial to Hill frame
# rho_H, rhoP_H = Inertial_2_Hill_mapping(rc,vc,rd,vd)
# print("Final position of deputy in Hill frame")
# print("final_rho_H = ", rho_H/1000.," km" )
# print("final_rhoP_H = ", rhoP_H/1000., " km/s")
# print("\b") 

# # #plot results
# # ax = plt.figure().add_subplot(projection='3d')
# # ax.plot(rc_history[:,0],rc_history[:,1],rc_history[:,2], 'b', label='chief trajectory (in km)')
# # ax.plot(rd_history[:,0],rd_history[:,1],rd_history[:,2], 'r', label='deputy trajectory (in km)')
# # ax.legend()
# # plt.show()

# #METHOD 1-B: Using cartesian ECI coordinates and Kepler equations for time of flight
# print("---------------------------------------------------------------------------------------")
# print("**** METHOD 1-B: Using cartesian ECI coordinates and Kepler equations ****")
# print("\b")

# # Chief satellite time of flight:
# print("Chief satellite time of flight:")
# print("\b")
# rc_N_f, vc_N_f = Kepler_solver(mu, rc_N_0, vc_N_0, time_lapse, report = False)

# # Deputy satellite time of flight:
# print("Deputy satellite time of flight:")
# print("\b")
# rd_N_f, vd_N_f = Kepler_solver(mu, rd_N_0, vd_N_0, time_lapse, report = False)

# print("final results (Kepler):")
# print("(Kepler) chief position [in km] at {0} seconds: rc_N = {1}".format(time_lapse,rc_N_f/1000.))
# print("(Kepler) chief velocity [in km/s] at {0} seconds: vc_N = {1}".format(time_lapse,vc_N_f/1000.))
# print("(Kepler) deputy position [in km] at {0} seconds: rd_N = {1}".format(time_lapse,rd_N_f/1000.))
# print("(Kepler) deputy velocity [in km/s] at {0} seconds: vd_N = {1}".format(time_lapse,vd_N_f/1000.))
# print("\b")

# #Mapping back deputy position and velocity from inertial to Hill frame
# rho_H_f, rhoP_H_f = Inertial_2_Hill_mapping(rc_N_f,vc_N_f,rd_N_f,vd_N_f)
# print("Final position of deputy in Hill frame:")
# print("Mapping: Inertial ==> Hill:")
# print("(Kepler) final_rho_H = ", rho_H_f/1000.," km" )
# print("(Kepler) final_rhoP_H = ", rhoP_H_f/1000., " km/s")
# print("\b") 


# # #METHOD 2 --> Using Non Linear Relative Orbit EOM (Hill frame) and a [first order] ODE integrator
# print("---------------------------------------------------------------------------------------")
# print("**** METHOD 2: Integrating the general relative equations of motion in the Hill Frame ****")
# print("\b")

# #Initial conditions (H-frame relative position & velocity)
# rc = rc_N_0.copy()
# vc = vc_N_0.copy()
# rho = rho_H_0.copy()
# rhoP = rhoP_H_0.copy()

# h = np.cross(rc,vc)
# o_hat_r = rc/np.linalg.norm(rc)   #radial direction (unit vector)
# x = rho[0]
# y = rho[1]
# z = rho[2]
# x_dot = rhoP[0]
# y_dot = rhoP[1]
# z_dot = rhoP[2]

# #history list creation
# x_history = np.zeros([len(t)+1,3])
# y_history = np.zeros([len(t)+1,3])
# z_history = np.zeros([len(t)+1,3])

# #Main loop
# x_history[0,:] = x
# y_history[0,:] = y
# z_history[0,:] = z 

# for i, ti in enumerate(t):
#     rc_norm = np.linalg.norm(rc)
#     rc_dot_norm = np.dot(vc,o_hat_r)
#     rd_norm = sqrt((rc_norm + x)**2 + y**2 + z**2)
#     fdot = np.linalg.norm(h)/rc_norm**2
#     #relative position propagation
#     x_dot_dot = (-mu/rd_norm**3)*(rc_norm + x) + 2*fdot*(y_dot - y*rc_dot_norm/rc_norm) + x*fdot**2 + mu/rc_norm**2
#     y_dot_dot = (-mu/rd_norm**3)*y - 2*fdot*(x_dot - x*rc_dot_norm/rc_norm) + y*fdot**2
#     z_dot_dot = (-mu/rd_norm**3)*z
#     x_dot += x_dot_dot*dt
#     x += x_dot*dt
#     y_dot += y_dot_dot*dt
#     y += y_dot*dt
#     z_dot += z_dot_dot*dt
#     z += z_dot*dt
#     #chief satellite orbit propagation
#     rc_dot_dot = -(mu/rc_norm**3)*rc
#     vc += rc_dot_dot*dt
#     rc += vc*dt
#     h = np.cross(rc,vc)
#     o_hat_r = rc/np.linalg.norm(rc)
#     #history update
#     x_history[i+1,:] = x/1000
#     y_history[i+1,:] = y/1000
#     z_history[i+1,:] = z/1000

# rho_H_final = np.array([x,y,z])
# rhoP_H_final = np.array([x_dot,y_dot,z_dot])
# rc_N_final = rc.copy()
# vc_N_final = vc.copy()

# print("final results (Relative ODE):")
# print("Final position of deputy in Hill frame")
# print("(Relative ODE) final_rho_H = ", rho_H_final/1000," km" )
# print("(Relative ODE) final_rhoP_H = ", rhoP_H_final/1000, " km/s")
# print("\b")
# print("(Relative ODE) chief position [in km] at {0} seconds: rc_N = {1}".format(time_lapse,rc_N_final/1000))
# print("(Relative ODE) chief velocity [in km/s] at {0} seconds: vc_N = {1}".format(time_lapse,vc_N_final/1000))

# #Map back from Hill frame to inertial
# rd_N_final, vd_N_final = Hill_2_Inertial_mapping(rc_N_final,vc_N_final,rho_H_final,rhoP_H_final)
# print("Mapping: Hill ==> Inertial:")
# print("(Relative ODE) deputy position [in km] at {0} seconds: rd_N = {1}".format(time_lapse,rd_N_final/1000))
# print("(Relative ODE) deputy velocity [in km/s] at {0} seconds: vd_N = {1}".format(time_lapse,vd_N_final/1000))
# print("\b") 

# #==============================================
# #quiz 5 - Linearized General Relative EOM
# #==============================================
# # exercise 2
# #===================================
# # mu = 398600441.8e6   #Earth geocentric gravitational constant
# mu = 398600000.0e6   #Earth geocentric gravitational constant

# #Initial conditions of chief in ECI frame (N):
# rc_N_0 = np.array([-6685.20926, 601.51244, 3346.066634])*1000   #m
# vc_N_0 = np.array([-1.74294, -6.70242, -2.27739])*1000          #m/s
# #Initial conditions of deputy in Hill frame
# rho_H_0 = np.array([-81.22301, 248.14201, 94.95904])*1000     #[m] deputy satellite relative pos in H-frame
# rhoP_H_0 = np.array([0.47884, 0.14857, 0.13577])*1000        #[m/s] deputy satellite relative vel in H-frame

# #Time domain
# time_lapse = 2000     #seconds
# dt = 0.001          #time step
# t = np.arange(0, time_lapse + dt, dt)

# #METHOD 3 --> Using linearized Relative Orbit EOM and a [first order] ODE integrator
# print("---------------------------------------------------------------------------------------")
# print("**** METHOD 3: Integrating the linearized relative equations of motion in the Hill Frame ****")
# print("\b")

# #Initial conditions (H-frame relative position & velocity)
# rc = rc_N_0.copy()
# vc = vc_N_0.copy()
# rho = rho_H_0.copy()
# rhoP = rhoP_H_0.copy()

# h = np.cross(rc,vc)
# o_hat_r = rc/np.linalg.norm(rc)   #radial direction (unit vector)
# x = rho[0]
# y = rho[1]
# z = rho[2]
# x_dot = rhoP[0]
# y_dot = rhoP[1]
# z_dot = rhoP[2]

# #history list creation
# x_history = np.zeros([len(t)+1,3])
# y_history = np.zeros([len(t)+1,3])
# z_history = np.zeros([len(t)+1,3])

# #Main loop
# x_history[0,:] = x
# y_history[0,:] = y
# z_history[0,:] = z 

# for i, ti in enumerate(t):
#     rc_norm = np.linalg.norm(rc)
#     rc_dot_norm = np.dot(vc,o_hat_r)
#     fdot = np.linalg.norm(h)/rc_norm**2
#     p = np.dot(h,h)/mu
#     #relative position propagation
#     x_dot_dot = x*fdot**2*(1 + 2*rc_norm/p) + 2*fdot*(y_dot - y*rc_dot_norm/rc_norm)
#     y_dot_dot = y*fdot**2*(1 - rc_norm/p) - 2*fdot*(x_dot - x*rc_dot_norm/rc_norm)
#     z_dot_dot = -(rc_norm/p*fdot**2)*z
#     x_dot += x_dot_dot*dt
#     x += x_dot*dt
#     y_dot += y_dot_dot*dt
#     y += y_dot*dt
#     z_dot += z_dot_dot*dt
#     z += z_dot*dt
#     #chief satellite orbit propagation
#     rc_dot_dot = -(mu/rc_norm**3)*rc
#     vc += rc_dot_dot*dt
#     rc += vc*dt
#     h = np.cross(rc,vc)
#     o_hat_r = rc/np.linalg.norm(rc)
#     #history update
#     x_history[i+1,:] = x/1000
#     y_history[i+1,:] = y/1000
#     z_history[i+1,:] = z/1000

# rho_H_final = np.array([x,y,z])
# rhoP_H_final = np.array([x_dot,y_dot,z_dot])

# print("Final position of deputy in Hill frame (Linealized Relative ODE)")
# print("final_rho_H (relative ODE) = ", rho_H_final/1000," km" )
# print("final_rhoP_H (relative ODE) = ", rhoP_H_final/1000, " km/s")
# print("\b") 

#==============================================
#quiz 6 - CWH Equations
#==============================================
# exercise 1
#===================================
# mu = 398600000.0e6   #Earth geocentric gravitational constant
# sma = 6800*1000      # [m], chief orbit semi-major axis

# #Chief orbit rate:
# n = sqrt(mu/sma**3)

# #Initial relative conditions:
# rho_H_0 = np.array([1.299038, -1.0000, 0.3213938])*1000     #[m] deputy satellite relative pos in H-frame
# rhoP_H_0 = np.array([-0.000844437, -0.00292521, -0.000431250])*1000        #[m/s] deputy satellite relative vel in H-frame

# #Time domain
# time_lapse = 1300     #seconds
# dt = 0.001         #time step
# t = np.arange(0, time_lapse + dt, dt)

# print("---------------------------------------------------------------------------------------")
# print("**** Numerical integration of CWH equations of motion in the Hill Frame ****")
# print("\b")

# #Initial conditions (H-frame relative position & velocity)
# rho = rho_H_0.copy()
# rhoP = rhoP_H_0.copy()

# x = rho[0]
# y = rho[1]
# z = rho[2]
# x_dot = rhoP[0]
# y_dot = rhoP[1]
# z_dot = rhoP[2]

# #history list creation
# rho_H_history = np.zeros([len(t)+1,3])

# #Main loop
# rho_H_history[0,:] = rho/1000


# for i, ti in enumerate(t):
#     #relative position propagation
#     x_dot_dot = 2*n*y_dot + 3*n**2*x
#     y_dot_dot = -2*n*x_dot
#     z_dot_dot = -n**2*z
#     x_dot += x_dot_dot*dt
#     x += x_dot*dt
#     y_dot += y_dot_dot*dt
#     y += y_dot*dt
#     z_dot += z_dot_dot*dt
#     z += z_dot*dt
#     #history update
#     rho_H_history[i+1,:] = np.array([x,y,z])/1000   #km
    
# rho_H_final = np.array([x,y,z])
# rhoP_H_final = np.array([x_dot,y_dot,z_dot])

# print("Final position of deputy in Hill frame (CWH equations) at {} seconds".format(t[-1]))
# print("final_rho_H (CWH equations) = ", rho_H_final/1000," km" )
# print("final_rhoP_H (CWH equations) = ", rhoP_H_final/1000, " km/s")
# print("\b")

# #plot results
# ax = plt.figure().add_subplot(projection='3d')
# ax.plot(rho_H_history[:,0],rho_H_history[:,1],rho_H_history[:,2], 'r', label='deputy relative trajectory (in km)')
# ax.legend()
# plt.show()


# #===================================
# # exercise 2
# #===================================
# mu = 398600000.0e6   #Earth geocentric gravitational constant
# time_lapse = 2000    #seconds

# # chief circular motion in inertial frame:
# rc = 7500*1000         # [m], chief orbit semi-major axis
# n = sqrt(mu/rc**3)    # rad/s

# #deputy relative motion in Hill-frame:
# y_off = 200*1000    #[m], Hill-frame deputy relative position

# #Initial conditions of chief (inertial frame)
# fc_initial = 0                   #initial true anomaly (arbitrarily defined)
# rc_N_0 = rc*np.array([cos(fc_initial), sin(fc_initial), 0])
# vc_N_0 = n*rc*np.array([cos(fc_initial+pi/2), sin(fc_initial+pi/2), 0])
# # #Initial conditions (H-frame relative position & velocity)
# rho_H_0 = np.array([0, y_off, 0])
# rhoP_H_0 = np.array([0, 0, 0])

# #Mapping deputy rd, vd from Hill frame to inertial:

# rd_N_0, vd_N_0 = Hill_2_Inertial_mapping(rc_N_0,vc_N_0,rho_H_0,rhoP_H_0)

# print("Results:")
# print("initial_rd_N = ", rd_N_0/1000.," km" )
# print("initial_vd_N = ", vd_N_0/1000., " km/s")
# print("\b")

# #Time domain
# time_lapse = 2000     #seconds
# dt = 0.001          #time step
# t = np.arange(0, time_lapse + dt, dt)

# # #METHOD 1A --> Using cartesian ECI coordinates and a [first order] ODE integrator
# print("---------------------------------------------------------------------------------------")
# print("**** METHOD 1-A: Using cartesian ECI coordinates and a [first order] ODE integrator ****")
# print("\b")

# #history arrays creation
# rc_history = np.zeros([len(t)+1,3])  #chief satellite position history (trajectory) in inertial frame
# vc_history = np.zeros([len(t)+1,3])  #chief satellite velocity history in inertial frame
# rd_history = np.zeros([len(t)+1,3])  #deputy satellite position history (trajectory) in inertial frame
# vd_history = np.zeros([len(t)+1,3])  #deputy satellite velocity history in inertial frame

# #Main loop
# rc = rc_N_0.copy()
# vc = vc_N_0.copy()
# rd = rd_N_0.copy()
# vd = vd_N_0.copy()
# rc_history[0,:] = rc/1000.   #km
# vc_history[0,:] = vc/1000.   #km/s
# rd_history[0,:] = rd/1000.   #km
# vd_history[0,:] = vd/1000.   #km/s

# for i,ti in enumerate(t):
#     #chief integration of inertial equations of motion
#     rc_norm = np.linalg.norm(rc)
#     rc_dot_dot = -(mu/rc_norm**3)*rc
#     vc += rc_dot_dot*dt
#     rc += vc*dt
#     rc_history[i+1,:]=rc/1000.
#     vc_history[i+1,:]=vc/1000.
    
#     #deputy satellite integration of inertial equations of motion
#     rd_norm = np.linalg.norm(rd)
#     rd_dot_dot = -(mu/rd_norm**3)*rd
#     vd += rd_dot_dot*dt
#     rd += vd*dt
#     rd_history[i+1,:]=rd/1000.
#     vd_history[i+1,:]=vd/1000.
    
# print("chief position [in km] at {} seconds: rc_N = ".format(t[-1]),rc_history[-1])
# print("chief velocity [in km/s] at {} seconds: vc_N = ".format(t[-1]),vc_history[-1])
# print("deputy position [in km] at {} seconds: rd_N = ".format(t[-1]),rd_history[-1])
# print("deputy velocity [in km/s] at {} seconds: vd_N = ".format(t[-1]),vd_history[-1])
# print("\b")

# #Mapping back deputy position and velocity from inertial to Hill frame
# rho_H, rhoP_H = Inertial_2_Hill_mapping(rc,vc,rd,vd)
# print("Final position of deputy in Hill frame")
# print("final_rho_H = ", rho_H/1000.," km" )
# print("final_rhoP_H = ", rhoP_H/1000., " km/s")
# print("\b") 


# #=====================================================================
# #quiz 12 - Linear Mapping between OE Differences and LVLH Coordinates
# #=====================================================================
# # exercise 3
# #=============
# # mu = 398600000.0e6  # Earth geocentric gravitational constant in SI units
# mu = 398600         # [km^3/s^2] Earth geocentric gravitational constant
# sma = 7500          # [m] semi-major axis
# theta = 13*pi/180   # [rad] full anomaly
# inc = 22*pi/180     # [rad] inclination
# q1 = 0.00707107
# q2 = 0.00707107
# AN = 70*pi/180      # [rad] ascending node

# A_matrix = OE_2_LVLH_linear_mapping_matrix_A(mu, sma, theta, inc, q1, q2, AN)

# print("[A_eo] = ", A_matrix)

# #=============
# # exercise 4
# #=============

# A_inverse = LVLH_2_OE_linear_mapping_inversematrix_A(mu, sma, theta, inc, q1, q2, AN)

# print("inverse([A_eo]) = ", A_inverse)

# # Health check:
# print("health check: [A]*inv([A]) = ", np.dot(A_matrix, A_inverse))

# #=====================================================================
# #quiz 14 - Formation Barycenter
# #=====================================================================
# # exercise 1
# #=============
# mu = 398600000.0e6  # Earth geocentric gravitational constant in SI units
# a = 7000*1000       # [m] semi-major axis
# e = 0.01            # eccentricity
# i = 78*pi/180       # [rad] inclination
# ohm = 120*pi/180    # [rad] Ascending Node
# w = 33*pi/180       # [rad] Argument of Periapsis
# f = 45*pi/180       # [rad] True Anomaly
# df = 1*pi/180       # [rad] deputy phase difference

# r1_N, v1_N = orbitelement_2_cartesian(mu,a,e,i,ohm,w,f)
# r2_N, v2_N = orbitelement_2_cartesian(mu,a,e,i,ohm,w,f+df)
# rCCM_N = (r1_N + r2_N)/2

# print("rCCM_N = ", rCCM_N/1000, "[km]")

# #=============
# # exercise 2
# #=============
# rOECM_N, _ = orbitelement_2_cartesian(mu,a,e,i,ohm,w,f+df/2)

# print("rEOCM_N = ", rOECM_N/1000, "[km]")

# #=============
# # exercise 3
# #=============
# # Step 1: inertial position and velocity vectors of the Cartesian center of mass location
# rCCM_N = (r1_N + r2_N)/2
# vCCM_N = (v1_N + v2_N)/2

# # Step 2: map these into equivalent classical orbit elements to determine the semi-major axis.
# sma, _, _, _, _, _ = cartesian_2_orbitelement(mu, rCCM_N, vCCM_N)    #sma, ecc, inc, AN, AP, f

# # Step 3: Compute the corresponding orbit period
# n = sqrt(mu/sma**3)
# period = 2*pi/n     # in seconds

# # Step 4: Compare this period to the true period of this formation Cartesian barycenter location.
# n_true = sqrt(mu/a**3)
# period_true = 2*pi/n_true  #in seconds

# print("absolute difference in periods in seconds: ", abs(period - period_true))


# #=====================================================================
# #quiz 15 - Linearized Differential Mean Anomaly
# #=====================================================================
# # exercise 1
# #=============
# mu = 398600         # [km^3/s^2] Earth geocentric gravitational constant
# e = 0.01
# a = 7500            # [km]
# da = 50             # [km]
# dM_0 = 15*pi/180    # [radians]
# delta_t = 14400     # [seconds]


# dM_t_exact = dM_0 + (sqrt(mu/(a + da)**3) - sqrt(mu/a**3))*delta_t

# print("Exact Solution:")
# print("d_M(t)_exact = ", dM_t_exact*180/pi, "degrees")

# dM_t_linear = dM_0 - 3/2*( sqrt(mu/a**3)*delta_t )*da/a

# print("Linearized Solution:")
# print("d_M(t)_linearized = ", dM_t_linear*180/pi, "degrees")


#=====================================================================
#quiz 16 - OE Differential Solutions
#=====================================================================
# exercise 1
# #=============
mu = 398600000.0e6      # Earth geocentric gravitational constant in SI units
# Chief orbit
a = 10000*1000          # [m] chief semi-major axis
e = 0.2                 # chief eccentricity
i = 37*pi/180           # [rad] chief inclination
ohm = 40*pi/180         # [rad] chief ascending node
w = 65*pi/180           # [rad] chief argument of periapsis
f_0 = 10*pi/180         # [rad] chief true anomaly

# Deputy relative motion (orbital elements)
delta_a = 0              # [km]
delta_e = 0.0001
delta_i = 0.001
delta_ohm = 0.001
delta_w = 0.001
delta_M_0 = -0.001       # [rad]

# METHOD 1: Non-linear mapping between orbit elements and Hill state vector X
print("\b")
print("---------- METHOD 1: Non-linear mapping between orbit elements and Hill state vector X ----------")
print("\b")
# Chief initial conditions:
rc_N_0, vc_N_0 = orbitelement_2_cartesian(mu,a,e,i,ohm,w,f_0)
Mc_0 = trueanomaly_2_meananomaly(f_0,e)
# Deputy initial conditions
Md_0 = Mc_0 + delta_M_0
fd_0 = meananomaly_2_trueanomaly(Md_0,e+delta_e)
rd_N_0, vd_N_0 = orbitelement_2_cartesian(mu,a+delta_a,e+delta_e,i+delta_i,ohm+delta_ohm,w+delta_w,fd_0)
# Map to Hill Frame initial conditions
rho_H_0, rhoP_H_0 = Inertial_2_Hill_mapping(rc_N_0, vc_N_0, rd_N_0, vd_N_0)
print("Method 1:   rho_0_H = ", rho_H_0/1000, "km")

# Chief final conditions
f_t = f_0 + 60*pi/180
rc_N_t, vc_N_t = orbitelement_2_cartesian(mu,a,e,i,ohm,w,f_t)
Mc_t = trueanomaly_2_meananomaly(f_t,e)
# Deputy final conditions
delta_M = delta_M_0 - 3/2*( trueanomaly_2_meananomaly(f_t,e) - trueanomaly_2_meananomaly(f_0,e) )*delta_a/a
Md_t = Mc_t + delta_M
fd_t = meananomaly_2_trueanomaly(Md_t,e+delta_e)
rd_N_t, vd_N_t = orbitelement_2_cartesian(mu,a+delta_a,e+delta_e,i+delta_i,ohm+delta_ohm,w+delta_w,fd_t)
# Map to Hill Frame final conditions
rho_H_t, rhoP_H_t = Inertial_2_Hill_mapping(rc_N_t, vc_N_t, rd_N_t, vd_N_t)
print("Method 1:   rho_t_H = ", rho_H_t/1000, "km")

# METHOD 2: Using the EO difference solution for general elliptic chief orbits
print("\b")
print("---------- METHOD 2: Using the EO difference solution for general elliptic chief orbits ----------")
print("\b")
# Initial conditions
f = f_0
theta = w + f
delta_M = delta_M_0
r = a*(1 - e**2)/(1 + e*cos(f))
eta = sqrt(1 - e**2)
f_u = atan2(e*delta_M,(-eta*delta_e))
# f_v = atan2(eta*delta_e,(e*delta_M))
theta_w = atan2(delta_i,(-sin(i)*delta_ohm))

delta_u = sqrt((e**2*delta_M**2)/eta**2 + delta_e**2)
delta_w = sqrt(delta_i**2 + (sin(i))**2*delta_ohm**2)

u = delta_a/a - e*delta_e/(2*eta**2) + (delta_u/eta**2)*(cos(f - f_u) + e/2*cos(2*f - f_u))
v = ((1 + 0.5*e**2)*delta_M/eta**3 + delta_w + cos(i)*delta_ohm) - (delta_u/eta**2)*(2*sin(f - f_u) + e/2*sin(2*f - f_u))
w = delta_w*cos(theta - theta_w)

rho_0_H = r*np.array([u,v,w])
print("Method 2:   rho_0_H = ", rho_0_H/1000, "[km]")

# Final conditions ( f2 = f_0 + 60 degrees)
f = f_0 + 60*pi/180
theta = w + f
delta_M = delta_M_0 - 3/2*( trueanomaly_2_meananomaly(f,e) - trueanomaly_2_meananomaly(f_0,e) )*delta_a/a
r = a*(1 - e**2)/(1 + e*cos(f))
eta = sqrt(1 - e**2)
f_u = atan2(e*delta_M,(-eta*delta_e))
f_v = atan2(eta*delta_e,(e*delta_M))
theta_w = atan2(delta_i,(-sin(i)*delta_ohm))

delta_u = sqrt(e**2*delta_M**2/eta**2 + delta_e**2)
delta_w = sqrt(delta_i**2 + (sin(i))**2*delta_ohm**2)

u = delta_a/a - e*delta_e/(2*eta**2) + (delta_u/eta**2)*(cos(f - f_u) + e/2*cos(2*f - f_u))
v = ((1 + 0.5*e**2)*delta_M/eta**3 + delta_w + cos(i)*delta_ohm) - (delta_u/eta**2)*(2*sin(f - f_u) + e/2*sin(2*f - f_u))
w = delta_w*cos(theta - theta_w)

rho_t_H = r*np.array([u,v,w])
print("Method 2:   rho_t_H = ", rho_t_H/1000, "[km]")

# METHOD 3: Linearized mapping about chief orbit elements using matrix [A(oe)]
print("\b")
print("---------- METHOD 3: Linearized mapping about chief orbit elements using matrix [A(oe)] ----------")
print("\b")

# Initial conditions
theta_0 = w + f_0
q1 = e*cos(w)
q2 = e*sin(w)
A_oe_0 = OE_2_LVLH_linear_mapping_matrix_A(mu,a,theta_0,i,q1,q2,ohm)
# differential orbit elements array
delta_theta_0 = delta_w + (fd_0 - f_0)
delta_q1 = delta_e*cos(w) - e*sin(w)*delta_w
delta_q2 = delta_e*sin(w) + e*cos(w)*delta_w
delta_oe_0 = np.array([delta_a,delta_theta_0,delta_i,delta_q1,delta_q2,delta_ohm])
X_0 = np.dot(A_oe_0,delta_oe_0)
Xrho_0_H = X_0[:3]
print("Method 3:   rho_0_H = ", Xrho_0_H/1000, "[km]")

# Final conditions
theta_t = w + f_0 + 60*pi/180
A_oe_t = OE_2_LVLH_linear_mapping_matrix_A(mu,a,theta_t,i,q1,q2,ohm)
# differential orbit elements array
delta_theta_t = delta_w + (fd_t - f_t)
delta_oe_t = np.array([delta_a,delta_theta_t,delta_i,delta_q1,delta_q2,delta_ohm])
X_t = np.dot(A_oe_t,delta_oe_t)
Xrho_t_H = X_t[:3]
print("Method 3:   rho_0_H = ", Xrho_t_H/1000, "[km]")
