# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 11:01:22 2019

@author: Stijn
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 14:49:46 2019

@author: Stijn
"""

import detailed_design.parameters as p
import numpy as np
from matplotlib import pyplot as plt

# Runway Input Parameters
g = p.g
mu = 0.3

# Atmosphere Input Parameters
temperature0 = 288.15
temperature_gradient = -0.0065
gamma_at = 1.4
rho0 = 1.225
R = 287
p0 = 101325


# Aircraft Input Parameters
S = p.S
CL = p.C_L_max_land
CD0 = 0.055
CD = CD0 + CL**2 / (np.pi * p.A * p.e)
T = p.T_TO
W = p.MTOW
gamma_d = np.arcsin(CD/CL)

altitude = np.arange(0, 1, 1)

LFlength = []


for i in range(len(altitude)):
    altitude_actual = altitude[i]
    temperature, pressure, rho, speed_of_sound = p.atmosphere_calc(altitude_actual, temperature0, temperature_gradient, g, R, gamma_at)
    rho = rho*rho0
    
    
    # Landing Input Parameters
    gamma = 3 * (np.pi/180)
    dn = 0.15
    h_screen = 50 * 0.3048
    V_stall = np.sqrt(W*2/(CL * rho *S))
    landing_factor = 10/6
    
    
    # Landing distance computation from Jenkinson
    V_approach = 1.3*V_stall    #ref CS25
    V_TD = 1.15*V_stall                 #ref Jenkinson
    V_flare = (V_approach + V_TD)/2
    R_flare = V_flare**2/(g*dn)
    h_flare = R_flare*(1-np.cos(gamma_d))
    x_approach = (h_screen - h_flare)/np.tan(gamma_d)
    x_flare = R_flare*gamma_d
    x_airborne = (x_approach + x_flare)
    
    t_roll = 2
    x_roll = t_roll*V_TD
     
    V_B = 0.95*V_TD
    CLg = 1
    CDg = (mu*(W-CLg*0.5*rho*S*V_B**2))/(0.5*rho*S*V_B**2)
    Dgmax = 0.4*W
    Z = (CDg - mu*CLg)/(mu*(W/(0.5*rho*S*V_B**2)))
    x_b1 = V_TD**2/(2*g*mu*Z)*np.log((1+Z)/(1+Z*(V_B**2/V_TD**2)))
    
    x_b2 = W/(CDg*rho*g*S)*np.log((CDg*0.5*rho*V_B**2*S)/(Dgmax)+1)
    
    x_groundrun = x_b1+x_b2
    
    x_total = x_airborne + x_groundrun
#    print(x_airborne, x_groundrun, x_total)
#    print(x_total*landing_factor)
    LFlength.append(x_total*landing_factor)

#plt.plot(altitude, LFlength)

# ATR 72 verification
S = 61
CL = 2.1
CD0 = 0.055
CD = CD0 + CL**2 / (np.pi * 12 * p.e)
T = p.T_TO
W = 21350*g
gamma_d = np.arcsin(CD/CL)

# Landing Input Parameters
gamma = 3 * (np.pi/180)
dn = 0.15
h_screen = 50 * 0.3048
V_stall = np.sqrt(W*2/(CL * rho *S))
landing_factor = 10/6


# Landing distance computation from Jenkinson
V_approach = 1.3*V_stall    #ref CS25
V_TD = 1.15*V_stall                 #ref Jenkinson
V_flare = (V_approach + V_TD)/2
R_flare = V_flare**2/(g*dn)
h_flare = R_flare*(1-np.cos(gamma_d))
x_approach = (h_screen - h_flare)/np.tan(gamma_d)
x_flare = R_flare*gamma_d
x_airborne_atr = (x_approach + x_flare)

t_roll = 2
x_roll = t_roll*V_TD
 
V_B = 0.95*V_TD
CLg = 1
CDg = (mu*(W-CLg*0.5*rho*S*V_B**2))/(0.5*rho*S*V_B**2)
Dgmax = 0.4*W
Z = (CDg - mu*CLg)/(mu*(W/(0.5*rho*S*V_B**2)))
x_b1 = V_TD**2/(2*g*mu*Z)*np.log((1+Z)/(1+Z*(V_B**2/V_TD**2)))

x_b2 = W/(CDg*rho*g*S)*np.log((CDg*0.5*rho*V_B**2*S)/(Dgmax)+1)

x_groundrun_atr = x_b1+x_b2

x_total_atr = x_airborne + x_groundrun
x_atr_LFL = 1210
print(1 - x_total_atr*landing_factor/x_atr_LFL)