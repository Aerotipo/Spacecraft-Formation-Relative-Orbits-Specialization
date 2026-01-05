import numpy as np 
from math import sin, cos, tan, atan, pi, sqrt
from utils import *

#===================================
#quiz 3 - Relative Motion Mapping
#===================================
# exercise 1
#===================================
rc_N = np.array([-4893.268, 3864.478, 3169.646])*1000     #[m] chief satellite position in ECI
vc_N = np.array([-3.91386, -6.257673, 1.59797])*1000      #[m/s] chief satellite velocity in ECI
rd_N = np.array([-4892.98, 3863.073, 3170.619])*1000     #[m] chief satellite position in ECI
# vd_N = np.array([-3.913302, -6.258661, 1.598199])*1000      #[m/s] chief satellite velocity in ECI
vd_N = np.array([-3.912794, -6.258249, 1.598643])*1000      #[m/s] chief satellite velocity in ECI --> correct input

rho_H, rhoP_H = Inertial_2_Hill_mapping(rc_N,vc_N,rd_N,vd_N)  #---> not working
print("exercise 1 results")
print("rho_H = ", rho_H/1000," km" )
print("rhoP_H = ", rhoP_H/1000, " km/s")
print("\b")

# rd_N_check, vd_N_check = Hill_2_Inertial_mapping(rc_N,vc_N,rho_H,rhoP_H)  #-- working fine

# print("health check")
# print("rd_N (check) = ", rd_N_check/1000," km" )
# print("vd_N (check) = ", vd_N_check/1000, " km/s")
# print("\b")

#===================================
# exercise 2
#===================================
rho_H = np.array([-0.537, 1.221, 1.106])*1000           #[m] deputy satellite relative pos in H-frame
rhoP_H = np.array([0.000486, 0.001158, 0.0005590])*1000 #[m/s] deputy satellite relative vel in H-frame

rd_N, vd_N = Hill_2_Inertial_mapping(rc_N,vc_N,rho_H,rhoP_H)  #--> working well

print("exercise 2 results\b")
print("rd_N = ", rd_N/1000," km" )
print("vd_N = ", vd_N/1000, " km/s")
print("\b")

# rho_H_check, rhoP_H_check = Inertial_2_Hill_mapping(rc_N,vc_N,rd_N,vd_N)

# print("health check exercise 2")
# print("rho_H (check) = ", rho_H_check/1000," km" )
# print("rhoP_H (check) = ", rhoP_H_check/1000, " km/s")

# #===================================
# # exercise 3
# #===================================
x = 10          #km
y = 500         #km
r = 7000        #km
x_dot = 0.1     #km/s
y_dot = -0.1    #km/s
r_dot = 0.05    #km/s


rd = sqrt((r + x)**2 + y**2)
delta_r = rd - r
s = r*atan(y/(r+x))
rd_dot = 1/(sqrt((x + r)**2 + y**2))*((x + r)*(x_dot + r_dot) + y*y_dot )
delta_r_dot = rd_dot - r_dot
aux = y/(r+x)
s_dot = r_dot*atan(aux) + ( r/(1 + (aux)**2) )*( y_dot/(r+x) - y*(r_dot + x_dot)/(r+x)**2 )

print("exercise 3 results")
print("delta_r = ", delta_r)
print("delta_r_dot = ", delta_r_dot)
print("s = ", s)
print("s_dot = ", s_dot)
print("\b")

# #===================================
# # exercise 4
# #===================================
dr = 10         #km
s =  500        #km
r = 7000        #km
dr_dot = 0.1    #km/s
s_dot = -0.1    #km/s
r_dot = 0.05    #km/s

delta_theta = s/r
delta_theta_dot = s_dot/r - s*r_dot/r**2
rd = r + dr
rd_dot = r_dot + dr_dot

x = rd*cos(delta_theta) - r
y = rd*sin(delta_theta)
x_dot = rd_dot*cos(delta_theta) - rd*sin(delta_theta)*delta_theta_dot - r_dot
y_dot = rd_dot*sin(delta_theta) + rd*cos(delta_theta)*delta_theta_dot

print("exercise 4 results")
print("x = ", x)
print("x_dot = ", x_dot)
print("y = ", y)
print("y_dot = ", y_dot)
print("\b")