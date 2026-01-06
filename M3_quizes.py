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

# # #===================================
# # # exercise 3
# # #===================================
# x = 10          #km
# y = 500         #km
# r = 7000        #km
# x_dot = 0.1     #km/s
# y_dot = -0.1    #km/s
# r_dot = 0.05    #km/s


# rd = sqrt((r + x)**2 + y**2)
# delta_r = rd - r
# s = r*atan(y/(r+x))
# rd_dot = 1/(sqrt((x + r)**2 + y**2))*((x + r)*(x_dot + r_dot) + y*y_dot )
# delta_r_dot = rd_dot - r_dot
# aux = y/(r+x)
# s_dot = r_dot*atan(aux) + ( r/(1 + (aux)**2) )*( y_dot/(r+x) - y*(r_dot + x_dot)/(r+x)**2 )

# print("exercise 3 results")
# print("delta_r = ", delta_r)
# print("delta_r_dot = ", delta_r_dot)
# print("s = ", s)
# print("s_dot = ", s_dot)
# print("\b")

# # #===================================
# # # exercise 4
# # #===================================
# dr = 10         #km
# s =  500        #km
# r = 7000        #km
# dr_dot = 0.1    #km/s
# s_dot = -0.1    #km/s
# r_dot = 0.05    #km/s

# delta_theta = s/r
# delta_theta_dot = s_dot/r - s*r_dot/r**2
# rd = r + dr
# rd_dot = r_dot + dr_dot

# x = rd*cos(delta_theta) - r
# y = rd*sin(delta_theta)
# x_dot = rd_dot*cos(delta_theta) - rd*sin(delta_theta)*delta_theta_dot - r_dot
# y_dot = rd_dot*sin(delta_theta) + rd*cos(delta_theta)*delta_theta_dot

# print("exercise 4 results")
# print("x = ", x)
# print("x_dot = ", x_dot)
# print("y = ", y)
# print("y_dot = ", y_dot)
# print("\b")

#==============================================
#quiz 4 - General Relative Equations of Motion
#==============================================
# exercise 2 - Task 1
#===================================
# mu = 398600441.8e6   #Earth geocentric gravitational constant
mu = 398600000.0e6   #Earth geocentric gravitational constant

#Initial conditions (chief & deputy satellites) in ECI frame (N):
rc_N_0 = np.array([-6685.20926, 601.51244, 3346.066634])*1000   #m
vc_N_0 = np.array([-1.74294, -6.70242, -2.27739])*1000          #m/s
rho_H_0 = np.array([-81.22301, 248.14201, 94.95904])*1000     #[m] deputy satellite relative pos in H-frame
rhoP_H_0 = np.array([0.47884, 0.14857, 0.13577])*1000        #[m/s] deputy satellite relative vel in H-frame

#Mapping deputy rd, vd from Hill frame to inertial:

rd_N_0, vd_N_0 = Hill_2_Inertial_mapping(rc_N_0,vc_N_0,rho_H_0,rhoP_H_0)

print("Results:")
print("initial_rd_N = ", rd_N_0/1000," km" )
print("initial_vd_N = ", vd_N_0/1000, " km/s")
print("\b")

#Time domain
time_lapse = 2000     #seconds
dt = 0.001          #time step
t = np.arange(0, time_lapse + dt, dt)

# #METHOD 1A --> Using cartesian ECI coordinates and a [first order] ODE integrator
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
# rc_history[0,:] = rc/1000   #km
# vc_history[0,:] = vc/1000   #km/s
# rd_history[0,:] = rd/1000   #km
# vd_history[0,:] = vd/1000   #km/s

# for i,ti in enumerate(t):
#     #chief integration of inertial equations of motion
#     rc_norm = np.linalg.norm(rc)
#     rc_dot_dot = -(mu/rc_norm**3)*rc
#     vc += rc_dot_dot*dt
#     rc += vc*dt
#     rc_history[i+1,:]=rc/1000
#     vc_history[i+1,:]=vc/1000
    
#     #deputy satellite integration of inertial equations of motion
#     rd_norm = np.linalg.norm(rd)
#     rd_dot_dot = -(mu/rd_norm**3)*rd
#     vd += rd_dot_dot*dt
#     rd += vd*dt
#     rd_history[i+1,:]=rd/1000
#     vd_history[i+1,:]=vd/1000
    
# print("chief position [in km] at {} seconds: rc_N = ".format(t[-1]),rc_history[-1])
# print("chief velocity [in km/s] at {} seconds: vc_N = ".format(t[-1]),vc_history[-1])
# print("deputy position [in km] at {} seconds: rd_N = ".format(t[-1]),rd_history[-1])
# print("deputy velocity [in km/s] at {} seconds: vd_N = ".format(t[-1]),vd_history[-1])
# print("\b")

# #Mapping back deputy position and velocity from inertial to Hill frame
# rho_H, rhoP_H = Inertial_2_Hill_mapping(rc,vc,rd,vd)
# print("Final position of deputy in Hill frame")
# print("final_rho_H = ", rho_H/1000," km" )
# print("final_rhoP_H = ", rhoP_H/1000, " km/s")
# print("\b") 

#plot results
# ax = plt.figure().add_subplot(projection='3d')
# ax.plot(rc_history[:,0],rc_history[:,1],rc_history[:,2], 'b', label='chief trajectory (in km)')
# ax.plot(rd_history[:,0],rd_history[:,1],rd_history[:,2], 'r', label='deputy trajectory (in km)')
# ax.legend()
# plt.show()

#METHOD 1-B: Using cartesian ECI coordinates and Kepler equations
print("---------------------------------------------------------------------------------------")
print("**** METHOD 1-B: Using cartesian ECI coordinates and Kepler equations ****")
print("\b")

# Chief satellite time of flight:
print("Chief satellite time of flight:")
print("\b")
sma, ecc, inc, AN, AP, f_initial = cartesian_2_orbitelement(mu, rc_N_0, vc_N_0)
n = sqrt(mu/sma**3)
print("chief semi-major axis [km] = ", sma/1000)
print("chief eccentricity = ", ecc)
print("chief inclination [rad] = ", inc)
print("chief Ascending Node [rad] = ", AN)
print("chief Argument of Periapsis [rad] = ", AP)
print("chief initial true anomaly f [rad] = ", f_initial)
print("chief orbital period: {} hours".format(2*pi/n/3600))
M_initial = trueanomaly_2_meananomaly(f_initial, ecc)
M_final = M_initial + sqrt(mu/sma**3)*time_lapse
f_final = meananomaly_2_trueanomaly(M_final,ecc)
# if f_final < 0:
    # f_final += 2*pi
print("chief final true anomaly f [rad] = ", f_final)
print("\b")
rc_N_f, vc_N_f = orbitelement_2_cartesian(mu,sma,ecc,inc,AN,AP,f_final)

# Deputy satellite time of flight:
print("Deputy satellite time of flight:")
print("\b")
sma, ecc, inc, AN, AP, f_initial = cartesian_2_orbitelement(mu, rd_N_0, vd_N_0)
n = sqrt(mu/sma**3)
print("deputy semi-major axis [km] = ", sma/1000)
print("deputy eccentricity = ", ecc)
print("deputy inclination [rad] = ", inc)
print("deputy Ascending Node [rad] = ", AN)
print("deputy Argument of Periapsis [rad] = ", AP)
print("deputy initial true anomaly f [rad] = ", f_initial)
print("deputy orbital period: {} hours".format(2*pi/n/3600))
M_initial = trueanomaly_2_meananomaly(f_initial, ecc)
M_final = M_initial + sqrt(mu/sma**3)*time_lapse
f_final = meananomaly_2_trueanomaly(M_final,ecc)
# if f_final < 0:
    # f_final += 2*pi
print("deputy final true anomaly f [rad] = ", f_final)
print("\b")
rd_N_f, vd_N_f = orbitelement_2_cartesian(mu,sma,ecc,inc,AN,AP,f_final)

print("chief position [in km] at {0} seconds: rc_N = {1}".format(time_lapse,rc_N_f/1000))
print("chief velocity [in km/s] at {0} seconds: vc_N = {1}".format(time_lapse,vc_N_f/1000))
print("deputy position [in km] at {0} seconds: rd_N = {1}".format(time_lapse,rd_N_f/1000))
print("deputy velocity [in km/s] at {0} seconds: vd_N = {1}".format(time_lapse,vd_N_f/1000))
print("\b")

#Mapping back deputy position and velocity from inertial to Hill frame
rho_H_final, rhoP_H_final = Inertial_2_Hill_mapping(rc_N_f,vc_N_f,rd_N_f,vd_N_f)
print("Final position of deputy in Hill frame (Kepler)")
print("final_rho_H (Kepler) = ", rho_H_final/1000," km" )
print("final_rhoP_H (Kepler) = ", rhoP_H_final/1000, " km/s")
print("\b") 


#METHOD 2 --> Using Hill frame differential equations and a [first order] ODE integrator
print("---------------------------------------------------------------------------------------")
print("**** METHOD 2: Integrating the general relative equations of motion in the Hill Frame ****")
print("\b")

#Initial conditions (H-frame relative position & velocity)
rc = rc_N_0.copy()
vc = vc_N_0.copy()
rho = rho_H_0.copy()
rhoP = rhoP_H_0.copy()

h = np.cross(rc,vc)
o_hat_r = rc/np.linalg.norm(rc)   #radial direction (unit vector)
x = rho[0]
y = rho[1]
z = rho[2]
x_dot = rhoP[0]
y_dot = rhoP[1]
z_dot = rhoP[2]

#history list creation
x_history = np.zeros([len(t)+1,3])
y_history = np.zeros([len(t)+1,3])
z_history = np.zeros([len(t)+1,3])

#Main loop
x_history[0,:] = x
y_history[0,:] = y
z_history[0,:] = z 

for i, ti in enumerate(t):
    rc_norm = np.linalg.norm(rc)
    rc_dot_norm = np.dot(vc,o_hat_r)
    rd_norm = sqrt((rc_norm + x)**2 + y**2 + z**2)
    fdot = np.linalg.norm(h)/rc_norm**2
    #relative position propagation
    x_dot_dot = (-mu/rd_norm**3)*(rc_norm + x) + 2*fdot*(y_dot - y*rc_dot_norm/rc_norm) + x*fdot**2 + mu/rc_norm**2
    y_dot_dot = (-mu/rd_norm**3)*y - 2*fdot*(x_dot - x*rc_dot_norm/rc_norm) + y*fdot**2
    z_dot_dot = (-mu/rd_norm**3)*z
    x_dot += x_dot_dot*dt
    x += x_dot*dt
    y_dot += y_dot_dot*dt
    y += y_dot*dt
    z_dot += z_dot_dot*dt
    z += z_dot*dt
    #chief satellite orbit propagation
    rc_dot_dot = -(mu/rc_norm**3)*rc
    vc += rc_dot_dot*dt
    rc += vc*dt
    h = np.cross(rc,vc)
    o_hat_r = rc/np.linalg.norm(rc)
    #history update
    x_history[i+1,:] = x/1000
    y_history[i+1,:] = y/1000
    z_history[i+1,:] = z/1000

rho_H_final = np.array([x,y,z])
rhoP_H_final = np.array([x_dot,y_dot,z_dot])

print("Final position of deputy in Hill frame (Relative ODE)")
print("final_rho_H (relative ODE) = ", rho_H_final/1000," km" )
print("final_rhoP_H (relative ODE) = ", rhoP_H_final/1000, " km/s")
print("\b") 
