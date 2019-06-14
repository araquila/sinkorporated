# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:30:57 2019

@author: Stijn
"""

import parameters as p
import numpy as np
import matplotlib.pyplot as plt

from atmosphere import atmosphere_calc

# Atmosphere Input Parameters
t0 = p.temperature0
t_gradient = p.temperature_gradient
atR = p.R
atgamma = p.gamma
g = p.g

# Aircraft Input Parameters
S = p.S
#S = 10
A = p.A
e = p.e
CD0 = p.Cd0
#CD0 = 0.0012
CL = p.C_L_cruise
CD = CD0 + CL**2/(np.pi*A*e)
T = p.tot_thrust
W = p.MTOW


range_cruise_tbp = p.range_cruise
eff_cruise_tbp = p.eff_cruise
cp_cruise_tbp = p.cp_cruise
#LD_cruise_tbp = p.LD_ratio
LD_cruise_tbp = 22

# Fuel fractions from Roskam for regional tbp
f1_tbp = 0.990      # W_1 / W_TO (Engine start, warm-up)
f2_tbp = 0.995      # W_2 / W_1 (Taxi)
f3_tbp = 0.995      # W_3 / W_2 (Take-off)
f4_tbp = 0.985      # W_4 / W_3 (Climb)
f5_tbp = None       # W_5 / W_4 (Cruise)

# Calculation of weight during descent
f5_tbp = 1/np.exp(range_cruise_tbp/((eff_cruise_tbp/(g*cp_cruise_tbp))*LD_cruise_tbp))
W = W*f1_tbp*f2_tbp*f3_tbp*f4_tbp*f5_tbp
#W = 4000
gammalist = []
gammalist2 = []
RDlist = []
Vh_dlist = []
V_dlist = []
altitudelist = []
gammalist3 = []
RDlist3 = []
Vh_dlist3 = []
V_dlist3 = []
altitudelist3 = []
altitude = [8000, 7500, 7000, 6500,  6000, 5500, 5000, 4500,  4000, 3500,  3000, 2500, 2000, 1500, 1000, 500, 0]
#altitude = [2000]
for j in range(len(altitude)):
    altitude1 = altitude[j]
    temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude1, t0, t_gradient, g, atR, atgamma)
    rho = rho*p.rho0
    CLlist = np.arange(0.01, 1.5, 0.01)
    gamma_d = 0.004
    T_idle = 0.0*T
    RD = []
    Vh_d = []
    V_d = []
    gammalist = []
    H1 = [0]
    #Rate of descent
    for i in range(len(CLlist)):
        for iteration1 in range(1):
            CL_d = CLlist[i]
            CD_d = CD0 + CL_d**2/(np.pi*A*e)
            gamma_d = np.arctan((CD_d/CL_d))
            for iteration2 in range(10):
                gamma_d = np.arctan(CD_d/CL_d - T_idle/W/np.cos(gamma_d))
#            print(gamma_d)
        gammalist.append(gamma_d*180/np.pi)
        RD_gamma = np.sqrt((2*W)/(S*rho*CL_d)*np.cos(gamma_d))*(CD_d/CL_d - T_idle/(W*np.cos(gamma_d)))*(np.cos(gamma_d))
        RD.append(RD_gamma)
        Vh_d.append(RD_gamma/np.tan(gamma_d))
        V_d.append(np.sqrt((W*2)/(S*rho*CL_d)*np.cos(gamma_d)))
        if RD_gamma >= 4.95 and RD_gamma < 5.05:
            altitudelist.append(altitude[j])
            gammalist2.append(gamma_d*180/np.pi)
            RDlist.append(RD_gamma)
            V_dlist.append(V_d[-1])
        elif RD_gamma >= 7.4 and RD_gamma < 7.6:
            altitudelist3.append(altitude[j])
            gammalist3.append(gamma_d*180/np.pi)
            RDlist3.append(RD_gamma)
            V_dlist3.append(V_d[-1])
        
#    gammalist = gammalist[::-1]
#    RD = RD[::-1]
#    Vh_d = Vh_d[::-1]
#    V_d = V_d[::-1]
#    plt.figure()
#    plt.plot(gammalist, V_d)
#    plt.xlim(-1, 10)
#    plt.ylim(0, 200)
#    plt.plot(gammalist, RD)
##    plt.xlim(-1, 10)
##    plt.ylim(0, 10)
#    plt.show()
#    
for f1 in range(5):
    for f in range(1, 17):
        if altitudelist[f] == altitudelist[f-1]:
            altitudelist.remove(altitudelist[f])
            gammalist2.remove(gammalist2[f])
            RDlist.remove(RDlist[f])
            V_dlist.remove(V_dlist[f])
for f in range(1, 18):
        if altitudelist[f] == altitudelist[f-1]:
            altitudelist.remove(altitudelist[f])
            gammalist2.remove(gammalist2[f])
            RDlist.remove(RDlist[f])
            V_dlist.remove(V_dlist[f])

x = len(altitudelist3)
for f1 in range(5):
    for f in range(1, 17):
        if altitudelist3[f] == altitudelist3[f-1]:
            altitudelist3.remove(altitudelist3[f])
            gammalist3.remove(gammalist3[f])
            RDlist3.remove(RDlist3[f])
            V_dlist3.remove(V_dlist3[f])
for f in range(1, 17):
        if altitudelist3[f] == altitudelist3[f-1]:
            altitudelist3.remove(altitudelist3[f])
            gammalist3.remove(gammalist3[f])
            RDlist3.remove(RDlist3[f])
            V_dlist3.remove(V_dlist3[f])
