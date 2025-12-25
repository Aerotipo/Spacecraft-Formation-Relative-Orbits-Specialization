import numpy as np 
import matplotlib.pyplot as plt
from math import sin, cos, tan, atan, pi, sqrt

#===================================
#quiz 4 - Conic Section Properties
#===================================
# a = 8000  #km
# e = 0.23  #eccentricity (ellipse)
# f = 295*pi/180 # true anomaly angle
# rp = a*(1-e)
# ra = a*(1+e)
# p = rp*(1+e)
# b = a*sqrt(1-e**2)
# E = 2*atan(sqrt((1-e)/(1+e))*tan(f/2))*180/pi

# print(a,e,f,rp,ra,p,b,E+360)

#================================================
#quiz 5 - Orbit Differential Equations of Motion
#================================================
# Exercise 3
#================================================
# mu = 398600441.8e6   #Earth geocentric gravitational constant

# #Initial conditions
# r = np.array([2466.69, 5941.54, 3282.71])*1000   #m
# v = np.array([-6.80822, 1.04998, 3.61939])*1000  #m/s

# #Time domain
# dt = 0.001
# t = np.arange(0,60*60+dt,dt)

# #history list creation
# r_history = np.zeros([len(t)+1,3])  #history of position vector (trajectory)
# v_history = np.zeros([len(t)+1,3])  #history of velocity vector
# h_history = np.zeros([len(t)+1,3])  #history of angular momentum (per unit mass) --> constant for unperturbed orbits

# #Main loop
# r_history[0,:] = r
# v_history[0,:] = v
# h_history[0,:] = np.cross(r,v)

# for i,ti in enumerate(t):
#     r_norm = np.linalg.norm(r)
#     r_dot_dot = -(mu/r_norm**3)*r
#     v += r_dot_dot*dt
#     r += v*dt
#     r_history[i+1,:]=r
#     v_history[i+1,:]=v
#     h_history[i+1,:]=np.cross(r,v)

# print('initial h',h_history[0,:])
# print('final h',h_history[-1,:])
# print('rVec_N (km) = ',r/1000)
# print('vVec_N (km/s) = ',v/1000)

# #plot results
# ax = plt.figure().add_subplot(projection='3d')
# ax.plot(r_history[:,0],r_history[:,1],r_history[:,2], label='orbit trajectory')
# ax.legend()
# plt.show()

#================================================
# Exercise 4
#================================================
mu = 398600441.8e6   #Earth geocentric gravitational constant
J2 = 1082.6310e-6    #-
r_eq = 6378.14e3     #m

#Function definition
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

#Initial conditions (satellitte 1 & 2)
r1 = np.array([-6685.20926, 601.51244, 3346.066634])*1000  #m
v1 = np.array([-1.74294, -6.70242, -2.27739])*1000  #m/s
r2 = np.array([-6685.21657, 592.52839, 3345.6716])*1000  #m
v2 = np.array([-1.74283, -6.70475, -2.27334])*1000  #m/s

#Time domain
dt = 0.001
t = np.arange(0, 4848+dt, dt)

#history list creation
r1_history = np.zeros([len(t)+1,3])  #satellite 1 history of position vector (trajectory)
v1_history = np.zeros([len(t)+1,3])  #satellite 2 history of velocity vector
r2_history = np.zeros([len(t)+1,3])  #satellite 1 history of position vector (trajectory)
v2_history = np.zeros([len(t)+1,3])  #satellite 2 history of velocity vector
h_history = np.zeros([len(t)+1,3])   #history of angular momentum (per unit mass) --> constant for unperturbed orbits

#Main loop
r1_history[0,:] = r1
v1_history[0,:] = v1
r2_history[0,:] = r2
v2_history[0,:] = v2
h_history[0,:] = np.cross(r1,v1)

for i,ti in enumerate(t):
    a_J2_1 = a_J2(r1, mu, J2, r_eq)
    r1_norm = np.linalg.norm(r1)
    r1_dot_dot = -(mu/r1_norm**3)*r1 + a_J2_1
    v1 += r1_dot_dot*dt
    r1 += v1*dt
    r1_history[i+1,:]=r1
    v1_history[i+1,:]=v1
    h_history[i+1,:]=np.cross(r1,v1)

    a_J2_2 = a_J2(r2, mu, J2, r_eq)
    r2_norm = np.linalg.norm(r2)
    r2_dot_dot = -(mu/r2_norm**3)*r2 + a_J2_2
    v2 += r2_dot_dot*dt
    r2 += v2*dt
    r2_history[i+1,:]=r2
    v2_history[i+1,:]=v2
    
print("satellite 1 position [in km] at {} seconds: r1 = ".format(t[-1]),r1/1000)
print("satellite 1 velocity [in km/s] at {} seconds: v1 = ".format(t[-1]),v1/1000)
print("satellite 2 position [in km] at {} seconds: r2 = ".format(t[-1]),r2/1000)
print("satellite 2 velocity [in km/s] at {} seconds: v2 = ".format(t[-1]),v2/1000)

#plot results
ax = plt.figure().add_subplot(projection='3d')
ax.plot(r1_history[:,0],r1_history[:,1],r1_history[:,2], 'b', label='satellite 1 trajectory')
ax.plot(r2_history[:,0],r2_history[:,1],r2_history[:,2], 'r', label='satellite 2 trajectory')
ax.legend()
plt.show()