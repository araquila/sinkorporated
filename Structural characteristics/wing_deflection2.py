# -*- coding: utf-8 -*-
"""
Created on Wed May 29 12:31:33 2019

@author: Stijn
"""

from math import *
import numpy as np
import parameters as p
#from wing_deflection1 import DefForces

def DefForces(Iavg, Li, offsetmax):
    g = 9.80665
    
    d_fuselage_outside = 2.8
    b = 29.76
    S = 49.2
    c_root = 2.36
    c_tip = 0.4*c_root
    x_strut = 0.4*b/2
    x_pod = x_strut
    x_engine = 0.2*b/2
    theta = atan(d_fuselage_outside/x_strut)
    
    E = 70*10**9
    EI = E*I
    vset = offsetmax
    
    #Weights
    W_empty = 9301 * g                  
    W_engine = (2.02 + 9.68)/100*W_empty
    W_wing = 13.85/100*W_empty
    W_pod = 5.09/100*W_empty + 1576 * g
    
    
    di = b/2/len(Li)
    xi = np.zeros(len(Li)+1)
    ci = np.zeros(len(Li)+1)
    Si = np.zeros(len(Li)+1)
    Wi = np.zeros(len(Li)+1)
    FyW = np.zeros(len(Li))
    FyL = np.zeros(len(Li))
    MoW = np.zeros(len(Li))
    MoL = np.zeros(len(Li))
    defW = np.zeros(len(Li))
    defL = np.zeros(len(Li))
    xset = b/2
    for i in range(len(Li)+1):
        xi[i] = i*di
        ci[i] = c_root - (c_root-c_tip)/(b/2)*xi[i]
    
    for i in range(len(Li)):
        Si[i] = (ci[i]+ci[i+1])/2*di
        Wi[i] = W_wing/S*Si[i]
        FyW[i] = Wi[i]*(xi[i+1]-xi[i])
        FyL[i] = Li[i]*(xi[i+1]-xi[i])
        MoW[i] = Wi[i]*(xi[i+1]-xi[i])*(xi[i]+di/2)
        MoL[i] = Li[i]*(xi[i+1]-xi[i])*(xi[i]+di/2)
        defW[i] = Wi[i]/24*(xset-xi[i])**4 - Wi[i]/24*(xset-xi[i+1])**4
        defL[i] = -Li[i]/24*(xset-xi[i])**4 + Li[i]/24*(xset-xi[i+1])**4
    
    A = [[1, 0, -cos(theta), 0], [0, 1, -sin(theta), 0], [0,0, -sin(theta)*x_strut, 1], [0, -1/6*(xset)**3, sin(theta)/6*(xset-x_strut)**3, 1/2*(xset)**2]]
    B = [[0], [W_engine + W_pod + sum(FyW) - sum(FyL)], [W_engine*x_engine + W_pod*x_pod + sum(MoW) - sum(MoL)], [-W_engine/6*(xset-x_engine)**3 - W_pod/6*(xset-x_pod)**3 + sum(defL) + sum(defW) - EI*vset]]
    C = np.matmul(np.linalg.inv(A),B)
    Frx = C[0][0]
    Fry = C[1][0]
    Fs = C[2][0]
    Mr = C[3][0]
    return Frx, Fry, Fs, Mr


g = 9.80665

#wing characteristics
d_fuselage_outside = 2.8
b = 29.76
S = 49.2
c_root = 2.36
c_tip = 0.4*c_root
x_strut = 0.4*b/2
x_pod = x_strut
x_engine = 0.2*b/2
theta = atan(d_fuselage_outside/x_strut)

#wing box characteristics
E = 70*10**9
Iavg = 10**(-5)
vset = 0.1*b/2

#Weights
W_empty = 9301 * g                  
W_engine = (2.02 + 9.68)/100*W_empty
W_wing = 13.85/100*W_empty
W_pod = 5.09/100*W_empty + 1576 * g


Li = [1000, 1000, 990, 990, 900, 800, 600, 300]

Frx, Fry, Fs, Mr = DefForces(Iavg, Li, vset)

di = b/2/len(Li)
xi = np.zeros(len(Li)+1)
ci = np.zeros(len(Li)+1)
vi = np.zeros(len(Li)+1)
Si = np.zeros(len(Li)+1)
Wi = np.zeros(len(Li)+1)

for i in range(len(Li)+1):
    xi[i] = i*di
    ci[i] = c_root - (c_root-c_tip)/(b/2)*xi[i]
    vi[i] = vset/(b/2)**2*xi[i]**2
    
for i in range(len(Li)):
    Si[i] = (ci[i]+ci[i+1])/2*di
    Wi[i] = W_wing/S*Si[i]
    FyW[i] = Wi[i]*(xi[i+1]-xi[i])
    FyL[i] = Li[i]*(xi[i+1]-xi[i])
    MoW[i] = Wi[i]*(xi[i+1]-xi[i])*(xi[i]+di/2)
    MoL[i] = Li[i]*(xi[i+1]-xi[i])*(xi[i]+di/2)

Ii = np.zeros(len(Li))
for i in range(len(Li)):
    xset = xi[i+1]
    vset = vi[i+1]-vi[i]

    defW = []
    defL = []
    
    for j in range(len(Li)):
        if xset < xi[j]:
            break
        elif xset >= xi[j] and xset >= xi[j+1]:
            defW.append(Wi[j]/24*(xset-xi[j])**4 - Wi[j]/24*(xset-xi[j+1])**4)
            defL.append(-Li[j]/24*(xset-xi[j])**4 + Li[j]/24*(xset-xi[j+1])**4)
        elif xset >= xi[j] and xset < xi[j+1]:
            defW.append(Wi[j]/24*(xset-xi[j])**4)
            defL.append(-Li[j]/24*(xset-xi[j])**4)
    #print(defL)
    if xset >= x_engine and xset >= x_strut and xset >= x_pod:
        Ii[i] = -(Mr/2*xset**2 - Fry/6*xset**3 + W_engine/6*(xset - x_engine)**3 + Fs*cos(theta)/6*(xset - x_strut)**3 + W_pod/6*(xset - x_pod)**3 + sum(defL) + sum(defW))/(E*vset)
    elif xset >= x_engine and xset >= x_strut and xset < x_pod:
        Ii[i] = -(Mr/2*xset**2 - Fry/6*xset**3 + W_engine/6*(xset - x_engine)**3 + Fs*cos(theta)/6*(xset - x_strut)**3 + sum(defL) + sum(defW))/(E*vset)
    elif xset >= x_engine and xset < x_strut and xset < x_pod:
        Ii[i] = -(Mr/2*xset**2 - Fry/6*xset**3 + W_engine/6*(xset - x_engine)**3 + sum(defL) + sum(defW))/(E*vset)
    elif xset < x_engine and xset < x_strut and xset < x_pod:
        Ii[i] = -(Mr/2*xset**2 - Fry/6*xset**3 + sum(defL) + sum(defW))/(E*vset)