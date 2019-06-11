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

def max_shear_stress(Lift, Drag, AeroMoment, Chord, shearyi, shearzi, height_along_wingbox, width_along_wingbox, Izz_along_wingbox, Iyy_along_wingbox):
        
    hi = height_along_wingbox
    bi = width_along_wingbox
    Izz = Izz_along_wingbox
    Iyy = Iyy_along_wingbox
    
    t = p.t_sheet
    
    sheary = shearyi
    shearz = shearzi
#    sheary = np.zeros(len(shearzi))
    
    taumax = []
    Drag.append(0)
    Lift.append(0)
    AeroMoment.append(0)
    Chord.append(p.tip_chord)
    
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
    cp_perc = 0.25
    for i in range(len(Iyy)):
        step = 0.001
        s1 = np.arange(0, bi[i]/2 + step , step)
        s2 = np.arange(0, hi[i] + step, step)
        s3 = np.arange(0, bi[i]/2 + step, step)
        s4 = np.arange(bi[i]/2 , bi[i]+step, step)
        s1list.append(s1)
        s2list.append(s2)
        s3list.append(s3)
        
        dz1 = (1/8*t*bi[i]**3 + 1/4*t*bi[i]*hi[i]**2 + 3/8*t*bi[i]**2*hi[i])
        dzeta = 1/Iyy[i]*(1/4*t*bi[i]*hi[i]**2 + 1/8*t*bi[i]**2*hi[i] - hi[i]/(hi[i]+bi[i])*dz1 + 1/4*bi[i]**3*t - 1/8*t*bi[i]**3 + 1/2*t*bi[i]**2*hi[i] + 1/8*t*bi[i]**3 - bi[i]/(bi[i]+hi[i])*dz1)
        eta = 0
    
        ShC = [eta, dzeta]
        Shear_center.append(ShC)
        Torsion = Drag[i]*dzeta - Lift[i]*(bi[i]/2 - (cp_perc-0.15)*Chord[i]) - AeroMoment[i]
        
        
        qs0 = -sheary[i]/Izz[i]*(1/8*bi[i]*t + 1/4*hi[i]*t + 1/12*hi[i]**2*t/bi[i]) - shearz[i]/Iyy[i]*(1/8*bi[i]**2/hi[i]*t + 3/8*bi[i]*t + 1/4*hi[i]*t) + shearz[i]*dzeta/(2*hi[i]*bi[i])
        qT = Torsion/(2*(hi[i]*bi[i]))
#        print(Torsion)
        qb12y = sheary[i]/Izz[i]*(1/2*hi[i]*t*s1) 
        qb23y = -sheary[i]/Izz[i]*(1/2*t*(s2)**2 - 1/2*t*hi[i]*s2 - 1/4*hi[i]*bi[i]*t) 
        qb34y = -sheary[i]/Izz[i]*(hi[i]/2*t*s3 - 1/4*hi[i]*bi[i]*t) 
        qb12z = shearz[i]/Iyy[i]*(1/2*t*(s1)**2)
        qb23z = - shearz[i]/Iyy[i]*(-1/2*t*bi[i]*s2 - 1/8*t*bi[i]**2)
        qb34z =  - shearz[i]/Iyy[i]*(-bi[i]/2*t*s3 + 1/2*t*(s3)**2 - 1/2*t*bi[i]*hi[i] - 1/8*t*bi[i]**2)
        
        #qb45 = -sheary[i]/Izz[i]*(hi[i]/2*t*s4 + 1/4*hi[i]*bi[i]*t) - shearz[i]/Iyy[i]*(-bi[i]/2*t*s4 + 1/2*t*(s4)**2 + 1/2*t*bi[i]*hi[i] + 1/8*t*bi[i]**2)
             
        qb45y = -qb34y[::-1]
        qb56y = -qb23y[::-1]
        qb61y = -qb12y[::-1]
        qb45z = qb34z[::-1]
        qb56z = qb23z[::-1]
        qb61z = qb12z[::-1]
        
        qset = qb34z[-1]/4
        for n in range(len(qb12z)):
            qb12z[n] = qb12z[n] - qset
            qb34z[n] = qb34z[n] - qset
            qb45z[n] = qb45z[n] - qset
            qb61z[n] = qb61z[n] - qset
        for n in range(len(qb23z)):
            qb23z[n] = qb23z[n] - qset
            qb56z[n] = qb56z[n] - qset
        
        q12 = qb12y + qb12z + qs0  + qT
        q23 = qb23y + qb23z + qs0 + qT
        q34 = qb34y + qb34z + qs0 + qT
        q45 = qb45y + qb45z + qs0 + qT
        q56 = qb56y + qb56z + qs0 + qT
        q61 = qb61y + qb61z + qs0 + qT
        
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
#    sectionplot = 0
#    plt.figure()
#    plt.plot(-s1list[sectionplot], q12list[sectionplot])
#    plt.plot(-s2list[sectionplot] - s1list[sectionplot][-1], q23list[sectionplot])
#    plt.plot(s3list[sectionplot] - s1list[sectionplot][-1], q34list[sectionplot])
#    plt.plot(s3list[sectionplot], q45list[sectionplot])
#    plt.plot(s2list[sectionplot] + s3list[sectionplot][-1], q56list[sectionplot])
#    plt.plot(-s1list[sectionplot] + s3list[sectionplot][-1], q61list[sectionplot])
#    plt.show()
    maximum_stress_along_wingbox = taumax 
    #returns shear stress in MPa
    return maximum_stress_along_wingbox



#TESTING
lengthdata = 100
tot_thrust = 20000
V_cruise = p.V_cruise
rho_cruise = p.rho
perc_engine = p.engine_pos_perc
perc_strut = p.strut_pos_perc
perc_pod = p.pod_pos_perc

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



Lift, Chord, Yle, Drag, AeroMoment = read_aero_data("wing/datastrut4.txt", lengthdata, V_cruise, rho_cruise)
Frx, Fry, Fs, Mrz, Frz, Fsz, Mry, momentyi, momentzi, shearyi, shearzi, vyi, vny, vzi, vnz, xi, theta = CallForces(Lift, Yle, Drag, tot_thrust, Iyy, Izz , p.E_al2014, perc_engine, perc_strut, perc_pod)
lengthdata = len(Lift)
b = p.b

taumax = max_shear_stress(Lift, Drag, AeroMoment, Chord, shearyi, shearzi, hi, bi, Izz, Iyy)
#TESTING
