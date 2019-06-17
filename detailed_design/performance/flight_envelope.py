# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 14:31:11 2019

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
rho0 = p.rho0
altitude1 = 8000
temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude1, t0, t_gradient, g, atR, atgamma)
rho = rho*rho0

# Aircraft Input Parameters
S = p.S
W = p.MTOW
A = p.A
e = p.e
CD0 = p.Cd0
CLmaxprofile = 1.5
CLmaxto = (p.C_L_max_TO + CLmaxprofile)/2*0.85
CLmaxld = (p.C_L_max_land + CLmaxprofile)/2*0.85
CLmaxclean = CLmaxprofile*0.85
#CDcruise = CD0 + CLcruise**2/(np.pi*A*e)

V_S = np.sqrt((2*W)/(rho*CLmaxclean*S))
V_S0 = np.sqrt((2*W)/(rho*CLmaxto*S))
V_S1 = np.sqrt((2*W)/(rho*CLmaxld*S))

n11 = 2.1 + 21000/(W/g*2.20462 + 10000)
n22 = -1
n1 = 1.5*n11
n2 = 1.5*n22
#n1 = n11
#n2 = n22

V_A = V_S * np.sqrt(n1/1.5)
V_F = max(1.4*V_S, 2*V_S0)
V_F = 1.6*V_S1
V_FE = V_F
VCmin = 4.77*np.sqrt(W/S)*0.514444
VC = p.V_cruise
VD = 1.4*VCmin

Vdiagram = []
ndiagram = []
V1range = np.arange(0, V_A, 0.1)
V2range = []
n1range = []
n2range = []
nflap = []
Vflap = []
Vtotal = []
ntotal = []
for i in range(len(V1range)):
    Lift = 0.5*rho0*S*CLmaxclean*0.85*V1range[i]**2
    if Lift/W <= n1:
        n1range.append(Lift/W)
        ndiagram.append(Lift/W)
        Vdiagram.append(V1range[i])
        Vtotal.append(V1range[i])
        ntotal.append(Lift/W)
    if -Lift/W >= n2:
        n2range.append(-Lift/W)
        V2range.append(V1range[i])
    
for i in range(len(V1range)):
    Liftflap = 0.5*rho0*S*CLmaxld*0.85*V1range[i]**2
    if Liftflap/W <= 2:
        nflap.append(Liftflap/W)
        Vflap.append(V1range[i])
    
Vflap.append(V_F)
nflap.append(2)
Vflap.append(V_F)
nflap.append(0)
Vflap.append(0)
nflap.append(0)


Vdiagram.append(VD)
ndiagram.append(n1range[-1])
Vdiagram.append(VD)
ndiagram.append(0)
Vdiagram.append(VC)
ndiagram.append(n2range[-1])

for i in range(len(n2range)):
    Vdiagram.append(V2range[len(n2range)-1-i])
    ndiagram.append(n2range[len(n2range)-1-i])


#plt.plot(Vdiagram, ndiagram)
Ude_1 = 20
Ude_C = 15.24*1.3
Ude_D = 7.62*1.3
a0 = 0.11 
b = p.b

a = a0*A/(2+np.sqrt(4+A**2))*180/np.pi
mu_g = 2*(W/g/S)/(rho0*S/b*a)
K_g = 0.88*mu_g/(5.3 + mu_g)

n_1 = 1 + 0.5*rho0*V_S*K_g*Ude_1/(W/S)*a -0.1
n_C = 1 + 0.5*rho0*VC*K_g*Ude_C/(W/S)*a -0.12
n_D = 1 + 0.5*rho0*VD*K_g*Ude_D/(W/S)*a 
n_1down = 1 - 0.5*rho0*V_S*K_g*Ude_1/(W/S)*a
n_Cdown = 1 - 0.5*rho0*VC*K_g*Ude_C/(W/S)*a +0.12
n_Ddown = 1 - 0.5*rho0*VD*K_g*Ude_D/(W/S)*a

ngustup = [1, n_C, n_D, 1]
Vgust = [0, VC, VD, 1]
ngustdown = [1, n_Cdown, n_Ddown, 1]



plt.plot(Vdiagram, ndiagram, label = "Wing flaps up")
plt.plot(Vflap, nflap, label = "Wing flaps down")
plt.plot(Vgust, ngustup, label = "Positive gust")
plt.plot(Vgust, ngustdown, label = "Negative gust")

# Label of the axes
plt.xlabel("Airspeed [m/s]", size="large")
plt.ylabel("Load factor [-]", size="large")

# Limits of the axes
plt.xlim([0, 220])
plt.ylim([-2,4.5])

# Create Legend
plt.legend(loc="best", fontsize="large")

# Show plot
plt.show()