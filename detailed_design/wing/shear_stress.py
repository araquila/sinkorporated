# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:01:19 2019

@author: Stijn
"""
from matplotlib import pyplot as plt
import parameters as p
from wing.section_properties import I_zz_wingbox, I_yy_wingbox, width_wingbox, height_wingbox
from wing.wing_deflection2 import read_aero_data, CallForces


V_cruise = p.V_cruise
rho_cruise = p.rho
lengthdata = 50
Lift, Chord, Yle, Drag = read_aero_data("wing/aquiladata1.txt", lengthdata, V_cruise, rho_cruise)
Frx, Fry, Fs, Mr, momentyi, momentzi, shearyi, shearzi, vii = CallForces(Lift, Chord, Yle, Drag, np.ones(len(Lift)+1)*10**(-4), 70*10**9, 0.2, 0.4, 0.4)

lengthdata = len(Lift)
b = p.b

di = b/2/lengthdata
xi = np.zeros(lengthdata + 1)
for i in range(len(xi)):
    xi[i] = i*di


hi = []
bi = []
t = 0.002
Izz = []
Iyy = []
for i in range(len(xi)):
    hi.append(height_wingbox(xi[i]))
    bi.append(width_wingbox(xi[i]))
    Izz.append(I_zz_wingbox(xi[i]))
    Iyy.append(I_yy_wingbox(xi[i]))
    
#t = 0.002
#xi = [0, 1]
#hi = [0.2, 0.2]
#bi = [2, 2]
#Izz = [10**-5, 10**-5]
#Iyy = [10**-5, 10**-5]

sheary_root = 87096.88196835331
sheary_1 = -11901.700560542904
shearz_root = 0
shearz_1 = 0
sheary = [sheary_root, sheary_1]
shearz = [shearz_root, shearz_1]
T = [0, 0]
taumax = []

for i in range(len(xi)):
    step = 0.001
    s1 = np.arange(0, bi[i]/2 + step , step)
    s2 = np.arange(0, hi[i] + step, step)
    s3 = np.arange(0, bi[i]/2 + step, step)
    dz1 = (1/8*t*bi[i]**3 + 1/4*t*bi[i]*hi[i]**2 + 3/8*t*bi[i]**2*hi[i])
    dzeta = 1/Iyy[i]*(1/4*t*bi[i]*hi[i]**2 + 1/8*t*bi[i]**2*hi[i] - hi[i]/(hi[i]+bi[i])*dz1 + 1/4*bi[i]**3*t - 1/8*t*bi[i]**3 + 1/2*t*bi[i]**2*hi[i] + 1/8*t*bi[i]**3 - bi[i]/(bi[i]+hi[i])*dz1)
    eta = 0
    ShC = [eta, dzeta]
    
    qs0 = -sheary[i]/Izz[i]*(1/8*bi[i]*t + 1/4*hi[i]*t + 1/12*hi[i]**2*t/bi[i]) - shearz[i]/Iyy[i]*(1/8*bi[i]**2/hi[i]*t + 3/8*bi[i]*t + 1/4*hi[i]*t) + shearz[i]*dzeta/(2*hi[i]*bi[i])
    qT = T[i]/(2*(hi[i]*bi[i]))
    
    qb12 = sheary[i]/Izz[i]*(1/2*hi[i]*t*s1) + shearz[i]/Iyy[i]*(1/2*t*(s1)**2)
    qb23 = -sheary[i]/Izz[i]*(1/2*t*(s2)**2 - 1/2*t*hi[i]*s2 - 1/4*hi[i]*bi[i]*t) - shearz[i]/Iyy[i]*(-1/2*t*bi[i]*s2 - 1/8*t*bi[i]**2)
    qb34 = -sheary[i]/Izz[i]*(hi[i]/2*t*s3 - 1/4*hi[i]*bi[i]*t) - shearz[i]/Iyy[i]*(-bi[i]/2*t*s3 + 1/2*t*(s3)**2 - 1/2*t*bi[i]*hi[i] - 1/8*t*bi[i]**2)
    
    q12 = qb12 + qs0 + qT
    q23 = qb23 + qs0 + qT
    q34 = qb34 + qs0 + qT
    
    tau12 = q12/t
    tau23 = q23/t
    tau34 = q34/t
    
    maxtau12 = max(tau12), min(tau12)
    maxtau23 = max(tau23), min(tau23)
    maxtau34 = max(tau34), min(tau34)
    
    maxtau = max(maxtau12[0], maxtau23[0], maxtau34[0]), min(maxtau12[1], maxtau23[1], maxtau34[1])
    taumax.append(max(abs(maxtau[0]), abs(maxtau[1])))

#plt.figure()
#plt.plot(s1, q12)
#plt.plot(s2, q23)
#plt.plot(s3, q34)
#plt.show()