import numpy as np 
# import matplotlib.pyplot as plt
from math import sin, cos, tan, atan, pi, sqrt

#===================================
#quiz 1 - Lagrange Matrix
#===================================
# exercise 4
#===================================
# k1 = 3
# k3 = 1
# w = sqrt(k1)

# #Time domain
# dt = 0.0001
# t = np.arange(0,10+dt,dt)

# #Initial conditions
# e1 = 1
# e2 = 0

# for ti in t:
#     x = e1*cos(w*ti) + e2/w*sin(w*ti)
#     ad = -k3*x**3
#     e1_dot = -(1/w*sin(w*ti))*ad
#     e2_dot = cos(w*ti)*ad
#     e1 += e1_dot*dt
#     e2 += e2_dot*dt

# print("e1 at", t[-1], "s = ", e1)
# print("e2 at", t[-1], "s = ", e2)

# #===================================
# #quiz 2 - Lagrangian Brackets
# #===================================
# # exercise 2
# #===================================
# k = 3       #spring constant
# g = 9.81    #m/s, gravity acceleration
# w = sqrt(k) #natural frequency of the "imagined" spring-mass system

# #Time domain
# dt = 0.0001
# t = np.arange(0,10+dt,dt)

# #Initial conditions
# e1 = 1
# e2 = 0

# for ti in t:
#     x = e1*cos(w*ti) + e2/w*sin(w*ti)
#     ad = -g
#     e1_dot = -(1/w*sin(w*ti))*ad
#     e2_dot = cos(w*ti)*ad
#     e1 += e1_dot*dt
#     e2 += e2_dot*dt

# print("e1 at", t[-1], "s = ", e1)
# print("e2 at", t[-1], "s = ", e2)

#===================================
#quiz 4 - Poisson Brackets
#===================================
# exercise 2
#===================================
m = 1       #mass
k = 3       #spring constant
g = 9.81    #m/s, gravity acceleration
w = sqrt(k/m) #natural frequency of the "imagined" spring-mass system

#Time domain
dt = 0.0001
t = np.arange(0,20+dt,dt)

#Initial conditions
e1 = 1
e2 = 0

for ti in t:
    x = e1*cos(w*ti) + e2/w*sin(w*ti)
    ad = -g
    e1_dot = -(1/w*sin(w*ti))*ad
    e2_dot = cos(w*ti)*ad
    e1 += e1_dot*dt
    e2 += e2_dot*dt

print("e1 at", t[-1], "s = ", e1)
print("e2 at", t[-1], "s = ", e2)
