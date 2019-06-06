# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:01:19 2019

@author: Stijn
"""
from matplotlib import pyplot as plt
import numpy as np
import parameters as p
from wing.section_properties import I_zz_wingbox, I_yy_wingbox, width_wingbox, height_wingbox
from wing.wing_deflection2 import read_aero_data, CallForces


tot_thrust = 20000
V_cruise = p.V_cruise
rho_cruise = p.rho
perc_engine = 0.15
perc_strut = 0.5
perc_pod = 0.5

di = p.b/2/lengthdata
xi = np.zeros(lengthdata + 1)
for i in range(len(xi)):
    xi[i] = i*di
    
hi = []
bi = []
t = p.t_sheet
Izz = []
Iyy = []
for i in range(len(xi)):
    hi.append(height_wingbox(xi[i]))
    bi.append(width_wingbox(xi[i]))
    Izz.append(I_zz_wingbox(xi[i]))
    Iyy.append(I_yy_wingbox(xi[i]))


lengthdata = 50
Lift, Chord, Yle, Drag = read_aero_data("wing/aquiladata1.txt", lengthdata, V_cruise, rho_cruise)
Frx, Fry, Fs, Mrz, Frz, Fsz, Mry, momentyi, momentzi, shearyi, shearzi, vyi, vny, vzi, vnz, xi, theta = CallForces(Lift, Yle, Drag, tot_thrust, Iyy, Izz , p.E_al2014, perc_engine, perc_strut, perc_pod)
lengthdata = len(Lift)
b = p.b



    
#t = 0.002
#xi = [0, 1]
#hi = [0.2, 0.2]
#bi = [2, 2]
#Izz = [10**-5, 10**-5]
#Iyy = [10**-5, 10**-5]

#sheary_root = 87096.88196835331
#sheary_1 = -11901.700560542904
#shearz_root = 0
#shearz_1 = 0
#sheary = [sheary_root, sheary_1]
#shearz = [shearz_root, shearz_1]

sheary = shearyi
shearz = shearzi

T = np.zeros(51)
taumax = []
Drag.append(0)

q12list = []
q23list = []
q34list = []
q45list = []
q56list = []
q61list = []
qs0list = []
qTlist = []

s1list = []
s2list = []
s3list = []

Shear_center = []
for i in range(len(xi)):
    step = 0.001
    s1 = np.arange(0, bi[i]/2 + step , step)
    s2 = np.arange(0, hi[i] + step, step)
    s3 = np.arange(0, bi[i]/2 + step, step)
    s1list.append(s1)
    s2list.append(s2)
    s3list.append(s3)
    
    dz1 = (1/8*t*bi[i]**3 + 1/4*t*bi[i]*hi[i]**2 + 3/8*t*bi[i]**2*hi[i])
    dzeta = 1/Iyy[i]*(1/4*t*bi[i]*hi[i]**2 + 1/8*t*bi[i]**2*hi[i] - hi[i]/(hi[i]+bi[i])*dz1 + 1/4*bi[i]**3*t - 1/8*t*bi[i]**3 + 1/2*t*bi[i]**2*hi[i] + 1/8*t*bi[i]**3 - bi[i]/(bi[i]+hi[i])*dz1)
    eta = 0
    ShC = [eta, dzeta]
    Shear_center.append(ShC)
    Torsion = Drag[i]*dzeta
    
    
    qs0 = -sheary[i]/Izz[i]*(1/8*bi[i]*t + 1/4*hi[i]*t + 1/12*hi[i]**2*t/bi[i]) - shearz[i]/Iyy[i]*(1/8*bi[i]**2/hi[i]*t + 3/8*bi[i]*t + 1/4*hi[i]*t) + shearz[i]*dzeta/(2*hi[i]*bi[i])
    qT = Torsion/(2*(hi[i]*bi[i]))
    
    qb12 = sheary[i]/Izz[i]*(1/2*hi[i]*t*s1) + shearz[i]/Iyy[i]*(1/2*t*(s1)**2)
    qb23 = -sheary[i]/Izz[i]*(1/2*t*(s2)**2 - 1/2*t*hi[i]*s2 - 1/4*hi[i]*bi[i]*t) - shearz[i]/Iyy[i]*(-1/2*t*bi[i]*s2 - 1/8*t*bi[i]**2)
    qb34 = -sheary[i]/Izz[i]*(hi[i]/2*t*s3 - 1/4*hi[i]*bi[i]*t) - shearz[i]/Iyy[i]*(-bi[i]/2*t*s3 + 1/2*t*(s3)**2 - 1/2*t*bi[i]*hi[i] - 1/8*t*bi[i]**2)
    qb45 = qb34[::-1]
    qb56 = qb23[::-1]
    qb61 = qb12[::-1]
    
    
    q12 = qb12 + qs0 + qT
    q23 = qb23 + qs0 + qT
    q34 = qb34 + qs0 + qT
    q45 = qb45 + qs0 + qT
    q56 = qb56 + qs0 + qT
    q61 = qb61 + qs0 + qT
    
    q12list.append(q12)
    q23list.append(q23)
    q34list.append(q34)
    q45list.append(q45)
    q56list.append(q56)
    q61list.append(q61)
    qs0list.append(qs0)
    qTlist.append(qT)
    
    
    
    tau12 = q12/t
    tau23 = q23/t
    tau34 = q34/t
    
    maxtau12 = max(tau12), min(tau12)
    maxtau23 = max(tau23), min(tau23)
    maxtau34 = max(tau34), min(tau34)
    
    maxtau = max(maxtau12[0], maxtau23[0], maxtau34[0]), min(maxtau12[1], maxtau23[1], maxtau34[1])
    taumax.append(max(abs(maxtau[0]), abs(maxtau[1]))/10**6)


#Plot on a single line
#sectionplot = 0
#plt.figure()
#plt.plot(s1list[sectionplot], q12list[sectionplot])
#plt.plot(s2list[sectionplot] + s1list[sectionplot][-1], q23list[sectionplot])
#plt.plot(s3list[sectionplot] + s1list[sectionplot][-1] + s2list[sectionplot][-1], q34list[sectionplot])
#plt.plot(s3list[sectionplot] + s1list[sectionplot][-1] + s2list[sectionplot][-1] + s3list[sectionplot][-1], q45list[sectionplot])
#plt.plot(s2list[sectionplot] + s1list[sectionplot][-1] + s2list[sectionplot][-1] + 2*s3list[sectionplot][-1], q56list[sectionplot])
#plt.plot(s1list[sectionplot] + s1list[sectionplot][-1] + 2*s2list[sectionplot][-1] + 2*s3list[sectionplot][-1], q61list[sectionplot])
#plt.show()


#Plot as a "box"
sectionplot = 0
plt.figure()
plt.plot(-s1list[sectionplot], q12list[sectionplot])
plt.plot(-s2list[sectionplot] - s1list[sectionplot][-1], q23list[sectionplot])
plt.plot(s3list[sectionplot] - s1list[sectionplot][-1], q34list[sectionplot])
plt.plot(s3list[sectionplot], q45list[sectionplot])
plt.plot(s2list[sectionplot] + s3list[sectionplot][-1], q56list[sectionplot])
plt.plot(-s1list[sectionplot] + s3list[sectionplot][-1], q61list[sectionplot])
plt.show()