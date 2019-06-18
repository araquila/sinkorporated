# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:30:57 2019

@author: Stijn
"""

import parameters as p
import numpy as np
import scipy as sc
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
gammalist1 = []
RDlist = []
Vh_dlist = []
V_dlist = []
altitudelist = []
gammalist2 = []
RDlist2 = []
Vh_dlist2 = []
V_dlist2 = []
altitudelist2 = []
gammalist3 = []
RDlist3 = []
Vh_dlist3 = []
V_dlist3 = []
altitudelist3 = []
gammalist4 = []
RDlist4 = []
Vh_dlist4 = []
V_dlist4 = []
altitudelist4 = []
#altitude = [8000, 7500, 7000, 6500,  6000, 5500, 5000, 4500,  4000, 3500,  3000, 2500, 2000, 1500, 1000, 500, 0]
altitude = np.arange(9000, -1, -100)
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
            gammalist1.append(gamma_d*180/np.pi)
            RDlist.append(RD_gamma)
            V_dlist.append(V_d[-1])
        elif RD_gamma >= 5.94 and RD_gamma < 6.08:
            altitudelist2.append(altitude[j])
            gammalist2.append(gamma_d*180/np.pi)
            RDlist2.append(RD_gamma)
            V_dlist2.append(V_d[-1])
        elif RD_gamma >= 6.89 and RD_gamma < 7.08:
            altitudelist3.append(altitude[j])
            gammalist3.append(gamma_d*180/np.pi)
            RDlist3.append(RD_gamma)
            V_dlist3.append(V_d[-1])
        elif RD_gamma >= 7.8 and RD_gamma < 8.05:
            altitudelist4.append(altitude[j])
            gammalist4.append(gamma_d*180/np.pi)
            RDlist4.append(RD_gamma)
            V_dlist4.append(V_d[-1])
        
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
z = 0
for f1 in range(10):
    x = len(altitudelist)
    for f in range(1, x):
        if altitudelist[f-z] == altitudelist[f-1-z]:
            altitudelist.remove(altitudelist[f-z])
            z = z+1
            gammalist1.remove(gammalist1[f-z])
            RDlist.remove(RDlist[f-z])
            V_dlist.remove(V_dlist[f-z])


z= 0
for f1 in range(10):
    x = len(altitudelist)
    for f in range(1, x):
        if altitudelist3[f-z] == altitudelist3[f-1-z]:
            altitudelist3.remove(altitudelist3[f-z])
            z= z+1
            gammalist3.remove(gammalist3[f-z])
            RDlist3.remove(RDlist3[f-z])
            V_dlist3.remove(V_dlist3[f-z])
           
z = 0
for f1 in range(5):
    x = len(altitudelist)
    for f in range(1, x):
        if altitudelist2[f-z] == altitudelist2[f-1-z]:
            altitudelist2.remove(altitudelist2[f-z])
            z = z+1
            gammalist2.remove(gammalist2[f-z])
            RDlist2.remove(RDlist2[f-z])
            V_dlist2.remove(V_dlist2[f-z])

z = 0
for f1 in range(5):
    x = len(altitudelist)
    for f in range(1, x):
        if altitudelist4[f-z] == altitudelist4[f-1-z]:
            altitudelist4.remove(altitudelist4[f-z])
            gammalist4.remove(gammalist4[f-z])
            RDlist4.remove(RDlist4[f-z])
            V_dlist4.remove(V_dlist4[f-z])


for i in range(1, len(gammalist1)):
    gammalist1[i] = gammalist1[i] - (gammalist1[i] - gammalist1[i-1])/2
    gammalist2[i] = gammalist2[i] - (gammalist2[i] - gammalist2[i-1])/2
    gammalist3[i] = gammalist3[i] - (gammalist3[i] - gammalist3[i-1])/2
    gammalist4[i] = gammalist4[i] - (gammalist4[i] - gammalist4[i-1])/2
    V_dlist[i] = V_dlist[i] - (V_dlist[i] - V_dlist[i-1])/2
    V_dlist2[i] = V_dlist2[i] - (V_dlist2[i] - V_dlist2[i-1])/2
    V_dlist3[i] = V_dlist3[i] - (V_dlist3[i] - V_dlist3[i-1])/2
    V_dlist4[i] = V_dlist4[i] - (V_dlist4[i] - V_dlist4[i-1])/2

for i in range(1, len(gammalist1)):
    gammalist1[i] = gammalist1[i] - (gammalist1[i] - gammalist1[i-1])/2
    gammalist2[i] = gammalist2[i] - (gammalist2[i] - gammalist2[i-1])/2
    gammalist3[i] = gammalist3[i] - (gammalist3[i] - gammalist3[i-1])/2
    gammalist4[i] = gammalist4[i] - (gammalist4[i] - gammalist4[i-1])/2
    V_dlist[i] = V_dlist[i] - (V_dlist[i] - V_dlist[i-1])/2
    V_dlist2[i] = V_dlist2[i] - (V_dlist2[i] - V_dlist2[i-1])/2
    V_dlist3[i] = V_dlist3[i] - (V_dlist3[i] - V_dlist3[i-1])/2
    V_dlist4[i] = V_dlist4[i] - (V_dlist4[i] - V_dlist4[i-1])/2
    
for i in range(1, len(gammalist1)):
    gammalist1[i] = gammalist1[i] - (gammalist1[i] - gammalist1[i-1])/2
    gammalist2[i] = gammalist2[i] - (gammalist2[i] - gammalist2[i-1])/2
    gammalist3[i] = gammalist3[i] - (gammalist3[i] - gammalist3[i-1])/2
    gammalist4[i] = gammalist4[i] - (gammalist4[i] - gammalist4[i-1])/2
    V_dlist[i] = V_dlist[i] - (V_dlist[i] - V_dlist[i-1])/2
    V_dlist2[i] = V_dlist2[i] - (V_dlist2[i] - V_dlist2[i-1])/2
    V_dlist3[i] = V_dlist3[i] - (V_dlist3[i] - V_dlist3[i-1])/2
    V_dlist4[i] = V_dlist4[i] - (V_dlist4[i] - V_dlist4[i-1])/2
    
# Plot lines
plt.figure()
plt.plot(gammalist1, altitudelist, label="ROD = 5 m/s")
plt.plot(gammalist2, altitudelist2, label="ROD = 6 m/s")
plt.plot(gammalist3, altitudelist3, label="ROD = 7 m/s")
plt.plot(gammalist4, altitudelist4, label="ROD = 8 m/s")

# Label of the axes
plt.xlabel("Path angle [deg]", size="large")
plt.ylabel("Altitude [m]", size="large")

# Limits of the axes
plt.xlim([2,4.5])
plt.ylim([0,9000])

# Create Legend
plt.legend(loc="best", fontsize="large")

# Show plot
plt.show()

# Plot lines
plt.figure()
plt.plot(V_dlist, altitudelist, label="ROD = 5 m/s")
plt.plot(V_dlist2, altitudelist2, label="ROD = 6 m/s")
plt.plot(V_dlist3, altitudelist3, label="ROD = 7 m/s")
plt.plot(V_dlist4, altitudelist4, label="ROD = 8 m/s")

# Label of the axes
plt.xlabel("Airspeed [m/s]", size="large")
plt.ylabel("Altitude [m]", size="large")

# Limits of the axes
plt.xlim([85, 160])
plt.ylim([0,9000])

# Create Legend
plt.legend(loc="best", fontsize="large")

# Show plot
plt.show()


#Calculate distance and time
time1 = 9000/5/60
time2 = 9000/6/60
time3 = 9000/7/60
time4 = 9000/8/60

Vhd1 = []
Vhd2 =[]
Vhd3 = []
Vhd4 = []
for i in range(len(altitude)):
    Vhd1.append(np.sqrt(V_dlist[i]**2 - 5**2))
    Vhd2.append(np.sqrt(V_dlist2[i]**2 - 6**2))
    Vhd3.append(np.sqrt(V_dlist3[i]**2 - 7**2))
    Vhd4.append(np.sqrt(V_dlist4[i]**2 - 8**2))


dist1 = sum(Vhd1)/len(Vhd1)*time1*60/1000
dist2 = sum(Vhd2)/len(Vhd2)*time2*60/1000
dist3 = sum(Vhd3)/len(Vhd3)*time3*60/1000
dist4 = sum(Vhd4)/len(Vhd4)*time4*60/1000

data = [ ["ROD = 5", time1  ,dist1],["ROD = 6", time2  , dist2], ["ROD = 7" , time3 , dist3], ["ROD = 8" , time4 , dist4]]