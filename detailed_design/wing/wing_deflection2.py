# -*- coding: utf-8 -*-
"""
Created on Wed May 29 12:31:33 2019

@author: Stijn
"""

from math import *
import numpy as np
#import parameters as p
from matplotlib import pyplot as plt

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

def SumzML(list1, xi, elementset):
    list2 = []
    for i in range(len(list1)):
        if xi[i] < elementset:
            outcome = list1[i]*(elementset - xi[i])
        else:
            outcome = 0
        list2.append(outcome)
    return sum(list2)

def SumzV(list1, xi, elementset, typeLW):
    list2 = []
    for i in range(len(list1)):
        if xi[i] < elementset and typeLW == 1:
            outcome = list1[i]
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

#UPDATE THESE VALUES
d_fuselage_outside = 2.8
S = 49.209
b = 31.372
taper = 0.4
c_root = 2.241
c_tip = taper*c_root
x_strut = 0.4*b/2
x_pod = x_strut
x_engine = 0.2*b/2
theta = atan(d_fuselage_outside/x_strut)

#wing box characteristics
E = 70*10**9
#Iavg = 10**(-6)
#offsetmax = 0.1*b/2

#Weights
W_empty = 9301 * g                  
W_engine = (2.02 + 9.68)/100*W_empty
W_wing = 13.85/100*W_empty
W_pod = 5.09/100*W_empty + 1576 * g
#W_empty = 0                 
#W_engine = 0
#W_pod = 0
Li = Lift
di = b/2/len(Li)
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
    momentL = SumzML(Li, Yle, xset)
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

    
for i in range(1,len(Li)+1):
    if i == 50:
        xmr.append(-(dii[i])/(2*E*Ii[i]) + xmr[i-1])
        xry.append((dii[i])/(2*E*Ii[i])*xi[i] + xry[i-1])
        xL.append(momentLi[i]*(dii[i])/(2*E*Ii[i]) + xL[i-1])
        xW.append(momentWi[i]*(dii[i])/(2*E*Ii[i]) + xW[i-1])
        if xi[i] >= x_engine and xi[i] < x_strut and xi[i] < x_pod:
            xwe.append(-(dii[i])/(2*E*Ii[i])*(xi[i]-x_engine) + xwe[i-1])
            xwp.append(0)
            xfs.append(0)
        elif xi[i] >= x_strut and xi[i] >= x_engine and xi[i] < x_pod:
            xwe.append(-(dii[i])/(2*E*Ii[i])*(xi[i]-x_engine) + xwe[i-1])
            xfs.append(-sin(theta)*(dii[i])/(2*E*Ii[i])*(xi[i]-x_strut) + xfs[i-1])
            xwp.append(0)
        elif xi[i] >= x_strut and xi[i] >= x_engine and xi[i] >= x_pod:
            xwe.append(-(dii[i])/(2*E*Ii[i])*(xi[i]-x_engine) + xwe[i-1])
            xfs.append(-sin(theta)*(dii[i])/(2*E*Ii[i])*(xi[i]-x_strut) + xfs[i-1])
            xwp.append(-(dii[i])/(2*E*Ii[i])*(xi[i]-x_pod) + xwp[i-1])
        elif xi[i] < x_strut and xi[i] < x_engine and xi[i] < x_pod:
            xfs.append(0)
            xwe.append(0)
            xwp.append(0)
    else:
        xmr.append(-(dii[i]+dii[i+1])/(2*E*Ii[i]) + xmr[i-1])
        xry.append((dii[i]+dii[i+1])/(2*E*Ii[i])*xi[i] + xry[i-1])
        xL.append(momentLi[i]*(dii[i]+dii[i+1])/(2*E*Ii[i]) + xL[i-1])
        xW.append(momentWi[i]*(dii[i]+dii[i+1])/(2*E*Ii[i]) + xW[i-1])
        if xi[i] >= x_engine and xi[i] < x_strut and xi[i] < x_pod:
            xwe.append(-(dii[i]+dii[i+1])/(2*E*Ii[i])*(xi[i]-x_engine) + xwe[i-1])
            xwp.append(0)
            xfs.append(0)
        elif xi[i] >= x_strut and xi[i] >= x_engine and xi[i] < x_pod:
            xwe.append(-(dii[i]+dii[i+1])/(2*E*Ii[i])*(xi[i]-x_engine) + xwe[i-1])
            xfs.append(-sin(theta)*(dii[i]+dii[i+1])/(2*E*Ii[i])*(xi[i]-x_strut) + xfs[i-1])
            xwp.append(0)
        elif xi[i] >= x_strut and xi[i] >= x_engine and xi[i] >= x_pod:
            xwe.append(-(dii[i]+dii[i+1])/(2*E*Ii[i])*(xi[i]-x_engine) + xwe[i-1])
            xfs.append(-sin(theta)*(dii[i]+dii[i+1])/(2*E*Ii[i])*(xi[i]-x_strut) + xfs[i-1])
            xwp.append(-(dii[i]+dii[i+1])/(2*E*Ii[i])*(xi[i]-x_pod) + xwp[i-1])
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
    FyL[i] = Li[i]
    MoW[i] = Wi[i]*(xi[i+1]-xi[i])*(xi[i]+(xi[i+1]-xi[i])/2)
    MoL[i] = Li[i]*(xi[-1] - Yle[i])

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
    shearL = sum(Li[:(i-1)])
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

#plt.plot(xi, vii)
    