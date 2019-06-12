# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 14:49:46 2019

@author: Stijn
"""

import parameters as p
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

altitude = np.arange(0, 1501, 1)

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
    landing_factor = 1.67
    
    
    # Landing distance computation from Jenkinson
    V_approach = 1.3*V_stall    #ref CS25
    V_TD = 1.15*V_stall                 #ref Jenkinson
    V_flare = (V_approach + V_TD)/2
    R_flare = V_flare**2/(g*dn)
    h_flare = R_flare*gamma**2/2
    x_approach = (h_screen - h_flare)/np.tan(gamma)
    x_flare = R_flare*gamma
    x_tot_airborne = (x_approach + x_flare)*landing_factor
    
    t_roll = 2
    x_roll = t_roll*V_TD
    K_T = -mu
    K_A = -(rho*CD)/(2*W/S)
    x_brake = -(1/(2*g*K_A))*np.log((K_T + K_A*V_TD**2)/K_T)
    x_tot_ground = (x_roll + x_brake)*landing_factor
    
    x_total = x_tot_airborne + x_tot_ground
#    print(x_tot_airborne, x_tot_ground, x_total)
    
    
    # Landing distance computation from Howe
    L_L = ((0.38/mu) + (5.59/(mu/0.38)) + 20.6*np.tan(gamma))
    LandingLength = 25.55/np.tan(gamma) + 4.5*V_approach + 0.0255*L_L*V_approach**2
#    print(LandingLength)
    LFlength.append([x_total, LandingLength])
    
plt.plot(LFlength)