# -*- coding: utf-8 -*-
"""
Created on Wed May 29 12:31:33 2019

@author: Stijn
"""

from math import *
import numpy as np
#import parameters as p
from matplotlib import pyplot as plt
#from wing_deflection1 import DefForces

def feet_to_meter(length_in_feet):
    length_in_meter = length_in_feet / 3.28084
    return length_in_meter

def feetsquared_to_metersquared(area_in_feetsquared):
    area_in_metersquared = area_in_feetsquared / 10.7639
    return area_in_metersquared


def read_aero_data(datafile, lengthdata, V_cruise, rho_cruise):
#    V_cruise = 184.84
#    rho_cruise = 0.5258
    
    filename = datafile
    
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    lengthlines = len(lines)
    
#    lengthdata = 50
    skipheader = 20
    table = np.genfromtxt(filename,delimiter=None , skip_header=skipheader, skip_footer=(lengthlines-lengthdata-skipheader-5))
    
    ji = []
    Yle = []
    Chord = []
    Area = []
    c_cl = []
    ai = []
    cl_norm = []
    cl = []
    cd = []
    cdv = []
    cm_c4 = []
    cm_LE = []
    CP_xc = []
    for i in range(lengthdata):
        ji.append(table[i][0])
        Yle.append(feet_to_meter(table[i][1]))
        Chord.append(feet_to_meter(table[i][2]))
        Area.append(feetsquared_to_metersquared(table[i][3]))
        c_cl.append(table[i][4])
        ai.append(table[i][5])
        cl_norm.append(table[i][6])
        cl.append(table[i][7])
        cd.append(table[i][8])
        cdv.append(table[i][9])
        cm_c4.append(table[i][10])
        cm_LE.append(table[i][11])
        CP_xc.append(table[i][12])
    
    Lift = []
    for i in range(len(ji)):
        Lift.append(1/2*rho_cruise*Area[i]*cl[i]*V_cruise**2)

    return Lift, Chord, Yle



#def DefForces(Iavg, Li, offsetmax):
#    g = 9.80665
#    
#    d_fuselage_outside = 2.84
#    b = 29.76
#    S = 49.2
#    c_root = 2.36
#    c_tip = 0.4*c_root
#    x_strut = 0.4*b/2
#    x_pod = x_strut
#    x_engine = 0.2*b/2
#    theta = atan(d_fuselage_outside/x_strut)
#    
#    E = 70*10**9
#    EI = E*Iavg
#    vset = offsetmax
#    
#    #Weights
#    W_empty = 9301 * g                  
#    W_engine = (2.02 + 9.68)/100*W_empty
#    W_wing = 13.85/100*W_empty
#    W_pod = 5.09/100*W_empty + 1576 * g
#    
#    
#    di = b/2/len(Li)
#    xi = np.zeros(len(Li)+1)
#    ci = np.zeros(len(Li)+1)
#    Si = np.zeros(len(Li)+1)
#    Wi = np.zeros(len(Li)+1)
#    FyW = np.zeros(len(Li))
#    FyL = np.zeros(len(Li))
#    MoW = np.zeros(len(Li))
#    MoL = np.zeros(len(Li))
#    defW = np.zeros(len(Li))
#    defL = np.zeros(len(Li))
#    xset = b/2
#    for i in range(len(Li)+1):
#        xi[i] = i*di
#        ci[i] = c_root - (c_root-c_tip)/(b/2)*xi[i]
#    
#    for i in range(len(Li)):
#        Si[i] = (ci[i]+ci[i+1])/2*di
#        Wi[i] = W_wing/S*Si[i]
#        FyW[i] = Wi[i]*(xi[i+1]-xi[i])
#        FyL[i] = Li[i]*(xi[i+1]-xi[i])
#        MoW[i] = Wi[i]*(xi[i+1]-xi[i])*(xi[i]+di/2)
#        MoL[i] = Li[i]*(xi[i+1]-xi[i])*(xi[i]+di/2)
#        defW[i] = Wi[i]/24*(xset-xi[i])**4 - Wi[i]/24*(xset-xi[i+1])**4
#        defL[i] = -Li[i]/24*(xset-xi[i])**4 + Li[i]/24*(xset-xi[i+1])**4
#    
#    A = [[1, 0, -cos(theta), 0], [0, 1, -sin(theta), 0], [0,0, -sin(theta)*x_strut, 1], [0, -1/6*(xset)**3, sin(theta)/6*(xset-x_strut)**3, 1/2*(xset)**2]]
#    B = [[0], [W_engine + W_pod + sum(FyW) - sum(FyL)], [W_engine*x_engine + W_pod*x_pod + sum(MoW) - sum(MoL)], [-W_engine/6*(xset-x_engine)**3 - W_pod/6*(xset-x_pod)**3 + sum(defL) + sum(defW) - EI*vset]]
#    C = np.matmul(np.linalg.inv(A),B)
#    
#    Frx = C[0][0]
#    Fry = C[1][0]
#    Fs = C[2][0]
#    Mr = C[3][0]
#    print(W_engine*x_engine + W_pod*x_pod + sum(MoW) - sum(MoL) + Fs*sin(theta)*x_strut - Mr)
#    return Frx, Fry, Fs, Mr


##set initial parameters
#g = 9.80665
#
##wing characteristics
#d_fuselage_outside = 2.8
#b = 29.76
#S = 49.2
#c_root = 2.36
#c_tip = 0.4*c_root
#x_strut = 0.4*b/2
#x_pod = x_strut
#x_engine = 0.2*b/2
#theta = atan(d_fuselage_outside/x_strut)
#
##wing box characteristics
#E = 70*10**9
#Iavg = 10**(-6)
#offsetmax = 0.1*b/2
#
##Weights
#W_empty = 9301 * g                  
#W_engine = (2.02 + 9.68)/100*W_empty
#W_wing = 13.85/100*W_empty
#W_pod = 5.09/100*W_empty + 1576 * g
#
#
#
##Li = [1000, 1000, 990, 990, 900, 800, 600, 300]
#Li = np.ones(50)*20000*g/di/50
#
##set forces
#Frx, Fry, Fs, Mr = DefForces(Iavg, Li, offsetmax)
#Forces = [Frx, Fry, Fs, Mr]

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



def SumzM(list1, xi, elementset, typeLW):
    list2 = []
    for i in range(len(list1)):
        if xi[i] < elementset and typeLW == 1:
            outcome = list1[i]/2*(elementset - xi[i])**2 - list1[i]/2*(elementset - xi[i+1])**2
        elif xi[i] < elementset and typeLW == 2:
            outcome = -list1[i]/2*(elementset - xi[i])**2 + list1[i]/2*(elementset - xi[i+1])**2
        else:
            outcome = 0
        list2.append(outcome)
    return sum(list2)

def SumzV(list1, xi, elementset, typeLW):
    list2 = []
    for i in range(len(list1)):
        if xi[i] < elementset and typeLW == 1:
            outcome = list1[i]*xi[-1]/len(Li)
        if xi[i] < elementset and typeLW == 2:
            outcome = -list1[i]*xi[-1]/len(Li)
        else:
            outcome = 0
        list2.append(outcome)
    return sum(list2)

V_cruise = 184.84
rho_cruise = 0.5258
Lift, Chord, Yle = read_aero_data("aquiladata1.txt", 50, V_cruise, rho_cruise)

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
#W_empty = 0                 
#W_engine = 0
#W_pod = 0

Li = np.ones(50)*20000*g/di/50
#Li = np.ones(50)*1000
Ii = np.ones(len(Li)+1)*10**(-4)

di = b/2/len(Li)
xi = np.zeros(len(Li)+1)
dii = np.zeros(len(Li)+1)
ci = np.zeros(len(Li)+1)
vi = np.zeros(len(Li)+1)
Si = np.zeros(len(Li))
Wi = np.zeros(len(Li))

for i in range(len(Li)+1):
    xi[i] = i*di
    ci[i] = c_root - (c_root-c_tip)/(b/2)*xi[i]

#    vi[i] = offsetmax/(b/2)**2*xi[i]**2
    
for i in range(len(Li)):
    dii[i+1] = xi[i+1] - xi[i]
    Si[i] = (ci[i]+ci[i+1])/2*di
    Wi[i] = W_wing/S*Si[i]/di

    
momentLi = []
momentWi = []
for i in range(len(Li)+1):
    xset = xi[i]
    momentL = SumzM(Li, xi, xset, 1)
    momentW = SumzM(Wi, xi, xset, 2)
    momentLi.append(momentL)
    momentWi.append(momentW)

xmr = [0]
xry = [0]
xfs = [0]
xwe = [0]
xwp = [0]
xL = [0]
xW = [0]
vmr = [0]
vry = [0]
vfs = [0]
vwe = [0]
vwp = [0]
vL = [0]
vW = [0]
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
for i in range(1,len(Li)+1):
    varMr1.append(di/(2*E)*(sum(1/Ii[:i])+sum(1/Ii[:(i-1)])))
    varRy1.append(di/(2*E)*(sum(xi[:i]/Ii[:i])+sum(xi[:(i-1)]/Ii[:(i-1)])))
    varmL1.append(di/(2*E)*(sum(momentLi[:i]/Ii[:i]) + sum(momentLi[:(i-1)]/Ii[:(i-1)])))
    varmW1.append(-di/(2*E)*(sum(momentWi[:i]/Ii[:i]) + sum(momentWi[:(i-1)]/Ii[:(i-1)])))
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
    
for i in range(1,len(Li)+1):
    xmr.append(-(dii[i]+dii[i-1])/(2*E*Ii[i]) + xmr[i-1])
    xry.append((dii[i]+dii[i-1])/(2*E*Ii[i])*xi[i] + xry[i-1])
    xL.append(momentLi[i]*(dii[i]+dii[i-1])/(2*E*Ii[i]) + xL[i-1])
    xW.append(momentWi[i]*(dii[i]+dii[i-1])/(2*E*Ii[i]) + xW[i-1])
    if xi[i] >= x_engine and xi[i] < x_strut and xi[i] < x_pod:
        xwe.append(-(dii[i]+dii[i-1])/(2*E*Ii[i])*(xi[i]-x_engine) + xwe[i-1])
        xwp.append(0)
        xfs.append(0)
    elif xi[i] >= x_strut and xi[i] >= x_engine and xi[i] < x_pod:
        xwe.append(-(dii[i]+dii[i-1])/(2*E*Ii[i])*(xi[i]-x_engine) + xwe[i-1])
        xfs.append(-sin(theta)*(dii[i]+dii[i-1])/(2*E*Ii[i])*(xi[i]-x_strut) + xfs[i-1])
        xwp.append(0)
    elif xi[i] >= x_strut and xi[i] >= x_engine and xi[i] >= x_pod:
        xwe.append(-(dii[i]+dii[i-1])/(2*E*Ii[i])*(xi[i]-x_engine) + xwe[i-1])
        xfs.append(-sin(theta)*(dii[i]+dii[i-1])/(2*E*Ii[i])*(xi[i]-x_strut) + xfs[i-1])
        xwp.append(-(dii[i]+dii[i-1])/(2*E*Ii[i])*(xi[i]-x_pod) + xwp[i-1])
    elif xi[i] < x_strut and xi[i] < x_engine and xi[i] < x_pod:
        xfs.append(0)
        xwe.append(0)
        xwp.append(0)

for i in range(1,len(Li)+1):
    vmr.append((xmr[i]+xmr[i-1])/2*dii[i] + vmr[i-1])
    vry.append((xry[i]+xry[i-1])/2*dii[i] + vry[i-1])
    vL.append((xL[i]+xL[i-1])/2*dii[i] + vL[i-1])
    vW.append((xW[i]+xW[i-1])/2*dii[i] + vW[i-1])
    vwe.append((xwe[i]+xwe[i-1])/2*dii[i] + vwe[i-1])
    vwp.append((xwp[i]+xwp[i-1])/2*dii[i] + vwp[i-1])
    vfs.append((xfs[i]+xfs[i-1])/2*dii[i] + vfs[i-1])

    
di = b/2/len(Li)
xi = np.zeros(len(Li)+1)
ci = np.zeros(len(Li)+1)
FyW = np.zeros(len(Li))
FyL = np.zeros(len(Li))
MoW = np.zeros(len(Li))
MoL = np.zeros(len(Li))
defW = np.zeros(len(Li))
defL = np.zeros(len(Li))
xset = b/2
vset = b/2*0.1

for i in range(len(Li)+1):
    xi[i] = i*di
    ci[i] = c_root - (c_root-c_tip)/(b/2)*xi[i]

for i in range(len(Li)):
    FyW[i] = Wi[i]*(xi[i+1]-xi[i])
    FyL[i] = Li[i]*(xi[i+1]-xi[i])
    MoW[i] = Wi[i]*(xi[i+1]-xi[i])*(xi[i]+(xi[i+1]-xi[i])/2)
    MoL[i] = Li[i]*(xi[i+1]-xi[i])*(xi[i]+(xi[i+1]-xi[i])/2)

A = [[1, 0, -cos(theta), 0], [0, 1, -sin(theta), 0], [0, -(xi[-1]), sin(theta)*(xi[-1] - x_strut), 1], [0, -vry[-1], -vfs[-1] , -vmr[-1]]]
B = [[0], [W_engine + W_pod + sum(FyW) - sum(FyL)], [-W_engine*(xi[-1] - x_engine) - W_pod*(xi[-1] - x_pod) - sum(MoW) + sum(MoL)], [ -vset + vwe[-1]*W_engine + vwp[-1]*W_pod + vL[-1] + vW[-1]]]
C = np.matmul(np.linalg.inv(A),B)

Frx = C[0][0]
Fry = C[1][0]
Fs = C[2][0]
Mr = C[3][0]



#moment diagram
momenti = []
for i in range(len(Li)+1):
    xset = xi[i]
    if xset >= x_engine and xset >= x_strut and xset >= x_pod:
        momi = -Mr + Fry*xset - W_engine*(xset - x_engine) - Fs*sin(theta)*(xset - x_strut) - W_pod*(xset - x_pod) + momentLi[i] + momentWi[i]
        momenti.append(momi)
    elif xset >= x_engine and xset >= x_strut and xset < x_pod:
        momi = -Mr + Fry*xset - W_engine*(xset - x_engine) - Fs*sin(theta)*(xset - x_strut) + momentLi[i] + momentWi[i]
        momenti.append(momi)
    elif xset >= x_engine and xset < x_strut and xset < x_pod:
        momi = -Mr + Fry*xset - W_engine*(xset - x_engine) + momentLi[i] + momentWi[i]
        momenti.append(momi)
    elif xset < x_engine and xset < x_strut and xset < x_pod:
        momi = -Mr + Fry*xset + momentLi[i] + momentWi[i]
        momenti.append(momi)
    

#shear diagram
sheari = []
for i in range(len(Li)+1):
    xset = xi[i]
    shearL = SumzV(Li, xi, xset, 1)
    shearW = SumzV(Wi, xi, xset, 2)
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
vii = [0]
vn = [0]
for i in range(1,(len(Li)+1)):
    thetai.append(thetai[i-1] + 1/2*(momenti[i]/(E*Ii[i]) + momenti[i-1]/(E*Ii[i-1]))*(dii[i]))
    
for i in range(1,(len(Li)+1)):
    vii.append(vii[i-1] + 1/2*(thetai[i] + thetai[i-1])*dii[i])

for n in range(1, len(Li)+1):
    vn.append(vmr[n]*Mr + vry[n]*Fry + vwe[n]*W_engine + vfs[n]*Fs + vwp[n]*W_pod + vW[n]+ vL[n] )
#plt.plot(xi, vi1)    

plt.plot(xi, vii)
    