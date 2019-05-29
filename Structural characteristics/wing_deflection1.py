# -*- coding: utf-8 -*-
"""
Created on Wed May 29 10:05:29 2019

@author: Stijn
"""
from math import *
import numpy as np

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
EI = 10**6

#Weights
W_empty = 9301 * g                  
W_engine = (2.02 + 9.68)/100*W_empty
W_wing = 13.85/100*W_empty
W_pod = 5.09/100*W_empty + 1576 * g


Li = [1000, 1000, 990, 990, 900, 800, 600, 300]
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
B = [[0], [W_engine + W_pod + sum(FyW) - sum(FyL)], [W_engine*x_engine + W_pod*x_pod + sum(MoW) - sum(MoL)], [-W_engine/6*(xset-x_engine)**3 - W_pod/6*(xset-x_pod)**3 + sum(defL) + sum(defW)]]
C = np.matmul(np.linalg.inv(A),B)