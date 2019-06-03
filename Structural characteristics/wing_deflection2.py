# -*- coding: utf-8 -*-
"""
Created on Wed May 29 12:31:33 2019

@author: Stijn
"""

from math import *
import numpy as np
import parameters as p
from matplotlib import pyplot as plt
#from wing_deflection1 import DefForces

def DefForces(Iavg, Li, offsetmax):
    g = 9.80665
    
    d_fuselage_outside = 2.84
    b = 29.76
    S = 49.2
    c_root = 2.36
    c_tip = 0.4*c_root
    x_strut = 0.4*b/2
    x_pod = x_strut
    x_engine = 0.2*b/2
    theta = atan(d_fuselage_outside/x_strut)
    
    E = 70*10**9
    EI = E*Iavg
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
    print(W_engine*x_engine + W_pod*x_pod + sum(MoW) - sum(MoL) + Fs*sin(theta)*x_strut - Mr)
    return Frx, Fry, Fs, Mr


#set initial parameters
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
Iavg = 10**(-6)
offsetmax = 0.1*b/2

#Weights
W_empty = 9301 * g                  
W_engine = (2.02 + 9.68)/100*W_empty
W_wing = 13.85/100*W_empty
W_pod = 5.09/100*W_empty + 1576 * g


#Li = [1000, 1000, 990, 990, 900, 800, 600, 300]
Li = np.ones(50)*1000

#set forces
Frx, Fry, Fs, Mr = DefForces(Iavg, Li, offsetmax)
Forces = [Frx, Fry, Fs, Mr]

#for iter in range(3):
#   
#    di = b/2/len(Li)
#    xi = np.zeros(len(Li)+1)
#    ci = np.zeros(len(Li)+1)
#    vi = np.zeros(len(Li)+1)
#    Si = np.zeros(len(Li)+1)
#    Wi = np.zeros(len(Li)+1)
#    
#    for i in range(len(Li)+1):
#        xi[i] = i*di
#        ci[i] = c_root - (c_root-c_tip)/(b/2)*xi[i]
#        vi[i] = offsetmax/(b/2)**2*xi[i]**2
#        
#    for i in range(len(Li)):
#        Si[i] = (ci[i]+ci[i+1])/2*di
#        Wi[i] = W_wing/S*Si[i]
#    
#    Ii = np.zeros(len(Li))
#    for i in range(len(Li)):
#        xset = xi[i+1]
#        vset = vi[i+1]
#    
#        defW = []
#        defL = []
#        
#        for j in range(len(Li)):
#            if xset < xi[j]:
#                break
#            elif xset >= xi[j] and xset >= xi[j+1]:
#                defW.append(Wi[j]/24*(xset-xi[j])**4 - Wi[j]/24*(xset-xi[j+1])**4)
#                defL.append(-Li[j]/24*(xset-xi[j])**4 + Li[j]/24*(xset-xi[j+1])**4)
#            elif xset >= xi[j] and xset < xi[j+1]:
#                defW.append(Wi[j]/24*(xset-xi[j])**4)
#                defL.append(-Li[j]/24*(xset-xi[j])**4)
#        #print(defL)
#        if xset >= x_engine and xset >= x_strut and xset >= x_pod:
#            Ii[i] = -(Mr/2*xset**2 - Fry/6*xset**3 + W_engine/6*(xset - x_engine)**3 + Fs*cos(theta)/6*(xset - x_strut)**3 + W_pod/6*(xset - x_pod)**3 + sum(defL) + sum(defW))/(E*vset)
#        elif xset >= x_engine and xset >= x_strut and xset < x_pod:
#            Ii[i] = -(Mr/2*xset**2 - Fry/6*xset**3 + W_engine/6*(xset - x_engine)**3 + Fs*cos(theta)/6*(xset - x_strut)**3 + sum(defL) + sum(defW))/(E*vset)
#        elif xset >= x_engine and xset < x_strut and xset < x_pod:
#            Ii[i] = -(Mr/2*xset**2 - Fry/6*xset**3 + W_engine/6*(xset - x_engine)**3 + sum(defL) + sum(defW))/(E*vset)
#        elif xset < x_engine and xset < x_strut and xset < x_pod:
#            Ii[i] = -(Mr/2*xset**2 - Fry/6*xset**3 + sum(defL) + sum(defW))/(E*vset)
#    
#    #iterate
#    xend = xi[-1]
#    Iavg = sum(Ii)/len(Ii)
#    offsetmax = -(Mr/2*xset**2 - Fry/6*xset**3 + W_engine/6*(xset - x_engine)**3 + Fs*cos(theta)/6*(xset - x_strut)**3 + W_pod/6*(xset - x_pod)**3 + sum(defL) + sum(defW))/(E*Iavg)
#    Frx, Fry, Fs, Mr = DefForces(Iavg, Li, offsetmax)
#    Forces = [Frx, Fry, Fs, Mr]
#    
#    
#print(Ii)
#print(Iavg)
#print(offsetmax)
#print(Forces)



def SumzM(list1, xi, elementset):
    list2 = []
    for i in range(len(list1)):
        if xi[i] < elementset:
            outcome = -list1[i]/2*(elementset - xi[i])**2 + list1[i]/2*(elementset - xi[i+1])**2
        else:
            outcome = 0
        list2.append(outcome)
    return sum(list2)

def SumzV(list1, xi, elementset):
    list2 = []
    for i in range(len(list1)):
        if xi[i] < elementset:
            outcome = list1[i]*xi[-1]/len(Li)
        else:
            outcome = 0
        list2.append(outcome)
    return sum(list2)

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
Iavg = 10**(-6)
offsetmax = 0.1*b/2

#Weights
W_empty = 9301 * g                  
W_engine = (2.02 + 9.68)/100*W_empty
W_wing = 13.85/100*W_empty
W_pod = 5.09/100*W_empty + 1576 * g


#Li = [1000, 1000, 990, 990, 900, 800, 600, 300]
Li = np.ones(50)*10
Ii = np.ones(51)*10**(-5)

#set forces
Frx, Fry, Fs, Mr = DefForces(Iavg, Li, offsetmax)
Forces = [Frx, Fry, Fs, Mr]

di = b/2/len(Li)
xi = np.zeros(len(Li)+1)
ci = np.zeros(len(Li)+1)
vi = np.zeros(len(Li)+1)
Si = np.zeros(len(Li)+1)
Wi = np.zeros(len(Li)+1)

for i in range(len(Li)+1):
    xi[i] = i*di
    ci[i] = c_root - (c_root-c_tip)/(b/2)*xi[i]
    vi[i] = offsetmax/(b/2)**2*xi[i]**2
    
for i in range(len(Li)):
    Si[i] = (ci[i]+ci[i+1])/2*di
    Wi[i] = W_wing/S*Si[i]

    
momentLi = []
momentWi = []
for i in range(len(Li)+1):
    xset = xi[i]
    momentL = SumzM(Li, xi, xset)
    momentW = SumzM(Wi, xi, xset)
    momentLi.append(momentL)
    momentWi.append(momentW)

funcMr = []
funcFy = []
funcFs = []
funcWe = []
funcWp = []
funcmL = []
funcmW = []
varMr1 = []
varRy1 = []
varFs1 = np.zeros(len(Li)+1)
varWe1 = np.zeros(len(Li)+1)
varWp1 = np.zeros(len(Li)+1)
varmL1 = []
varmW1 = []
varMr2 = np.zeros(len(Li)+1)
varRy2 = np.zeros(len(Li)+1)
varFs2 = np.zeros(len(Li)+1)
varWe2 = np.zeros(len(Li)+1)
varWp2 = np.zeros(len(Li)+1)
varmL2 = np.zeros(len(Li)+1)
varmW2 = np.zeros(len(Li)+1)
for i in range(len(Li)+1):
    varMr1.append(-di/(2*E)*(sum(1/Ii[:i])+sum(1/Ii[:(i-1)])))
    varRy1.append(di/(2*E)*(sum(xi[:i]/Ii[:i])+sum(xi[:(i-1)]/Ii[:(i-1)])))
    varmL1.append(-di/(2*E)*(sum(momentLi[:i]/Ii[:i])+sum(momentLi[:(i-1)]/Ii[:(i-1)])))
    varmW1.append(-di/(2*E)*(sum(momentWi[:i]/Ii[:i])+sum(momentWi[:(i-1)]/Ii[:(i-1)])))
    if xi[i] >= x_engine and xi[i] < x_strut and xi[i] < x_pod:
        varWe1[i] = (W_engine*(-di/(2*E)*(sum((xi[:i]-x_engine)/Ii[:i])+sum((xi[:(i-1)]-x_engine)/Ii[:(i-1)]))))
    elif xi[i] >= x_strut and xi[i] >= x_engine and xi[i] < x_pod:
        varFs1[i] = (-sin(theta)*di/(2*E)*(sum((xi[:i]-x_strut)/Ii[:i])+sum((xi[:(i-1)]-x_strut)/Ii[:(i-1)])))
        varWe1[i] = (W_engine*(-di/(2*E)*(sum((xi[:i]-x_engine)/Ii[:i])+sum((xi[:(i-1)]-x_engine)/Ii[:(i-1)]))))
    elif xi[i] >= x_pod and xi[i] >= x_engine and xi[i] >= x_strut:
        varWp1[i] = (W_pod*(-di/(2*E)*(sum((xi[:i]-x_pod)/Ii[:i])+sum((xi[:(i-1)]-x_pod)/Ii[:(i-1)]))))
        varFs1[i] = (-sin(theta)*di/(2*E)*(sum((xi[:i]-x_strut)/Ii[:i])+sum((xi[:(i-1)]-x_strut)/Ii[:(i-1)])))
        varWe1[i] = (W_engine*(-di/(2*E)*(sum((xi[:i]-x_engine)/Ii[:i])+sum((xi[:(i-1)]-x_engine)/Ii[:(i-1)]))))
        


for i in range(1,len(Li)+1):
    varMr2[i] = di/2*(sum(varMr1[:i]) + sum(varMr1[:(i-1)]))
    varRy2[i] = di/2*(sum(varRy1[:i]) + sum(varRy1[:(i-1)]))
    varmL2[i] = di/2*(sum(varmL1[:i]) + sum(varmL1[:(i-1)]))
    varmW2[i] = di/2*(sum(varmW1[:i]) + sum(varmW1[:(i-1)]))
    varWe2[i] = di/2*(sum(varWe1[:i]) + sum(varWe1[:(i-1)]))
    varWp2[i] = di/2*(sum(varWp1[:i]) + sum(varWp1[:(i-1)]))
    varFs2[i] = di/2*(sum(varFs1[:i]) + sum(varFs1[:(i-1)]))
    
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
vset = b/2*0.002
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
#    defW[i] = Wi[i]/24*(xset-xi[i])**4 - Wi[i]/24*(xset-xi[i+1])**4
#    defL[i] = -Li[i]/24*(xset-xi[i])**4 + Li[i]/24*(xset-xi[i+1])**4

Amatrix = A = [[1, 0, -cos(theta), 0], [0, 1, -sin(theta), 0], [0,0, -sin(theta)*x_strut, 1], [0, varRy2[-1], varFs2[-1] , varMr2[-1]]]
B = [[0], [W_engine + W_pod + sum(FyW) - sum(FyL)], [W_engine*x_engine + W_pod*x_pod + sum(MoW) - sum(MoL)], [ vset - varWe2[-1] - varWp2[-1] - varmL2[-1] - varmW2[-1]]]
C = np.matmul(np.linalg.inv(A),B)

Frx = C[0][0]
Fry = C[1][0]
Fs = C[2][0]
Mr = C[3][0]


#moment diagram
momenti = []
for i in range(len(Li)+1):
    xset = xi[i]
    momentL = SumzM(Li, xi, xset)
    momentW = SumzM(Wi, xi, xset)
    if xset >= x_engine and xset >= x_strut and xset >= x_pod:
        momi = -(Mr - Fry*xset + W_engine*(xset - x_engine) + Fs*cos(theta)*(xset - x_strut) + W_pod*(xset - x_pod) + momentL + momentW)
    elif xset >= x_engine and xset >= x_strut and xset < x_pod:
        momi = -(Mr - Fry*xset + W_engine*(xset - x_engine) + Fs*cos(theta)*(xset - x_strut) + momentL + momentW)
    elif xset >= x_engine and xset < x_strut and xset < x_pod:
        momi = -(Mr - Fry*xset + W_engine*(xset - x_engine) + momentL + momentW)
    elif xset < x_engine and xset < x_strut and xset < x_pod:
        momi = -(Mr - Fry*xset + momentL + momentW)
    momenti.append(momi)

#shear diagram
sheari = []
for i in range(len(Li)+1):
    xset = xi[i]
    shearL = SumzV(Li, xi, xset)
    shearW = SumzV(Wi, xi, xset)
    if xset >= x_engine and xset >= x_strut and xset >= x_pod:
        shri = Fry - W_engine - W_pod - Fs*sin(theta) + shearL + shearW
    elif xset >= x_engine and xset >= x_strut and xset < x_pod:
        shri = Fry - W_engine - Fs*sin(theta) + shearL + shearW
    elif xset >= x_engine and xset < x_strut and xset < x_pod:
        shri = Fry - W_engine + shearL + shearW
    elif xset < x_engine and xset < x_strut and xset < x_pod:
        shri = Fry + shearL + shearW
    sheari.append(shri)
    
thetai = [0]
vi1 = [0]
for i in range(1,(len(Li)+1)):
    theti = thetai[i-1] + 1/2*(momenti[i]/(E*Ii[i]) + momenti[i-1]/(E*Ii[i-1]))*(di)
    thetai.append(theti)
    vi2 = vi1[i-1] + 1/2*(thetai[i] + thetai[i-1])*(di)
    vi1.append(vi2)

#plt.plot(xi, vi1)    

plt.plot(xi, vi1)
    