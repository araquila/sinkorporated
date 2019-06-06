# -*- coding: utf-8 -*-
"""
Created on Wed May 29 12:31:33 2019

@author: Stijn
"""

from math import *
import numpy as np
import parameters as p
from matplotlib import pyplot as plt


#    V_cruise = 184.84
#    rho_cruise = 0.5258
#   skipheader = 50
def read_aero_data(datafile, lengthdata, V_cruise, rho_cruise):
#    V_cruise = 184.84
#    rho_cruise = 0.5258
    def feet_to_meter(length_in_feet):
        length_in_meter = length_in_feet / 3.28084
        return length_in_meter

    def feetsquared_to_metersquared(area_in_feetsquared):
        area_in_metersquared = area_in_feetsquared / 10.7639
        return area_in_metersquared
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
    Drag = []
    for i in range(len(ji)):
        Lift.append(1/2*rho_cruise*Area[i]*cl[i]*V_cruise**2)
        Drag.append(1/2*rho_cruise*Area[i]*cd[i]*V_cruise**2)

    return Lift, Chord, Yle, Drag











def CallForces(Lift, Yle, Drag, tot_thrust, Iyy, Izz, E, perc_engine, perc_strut, perc_pod):
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
    
#    V_cruise = 184.84
#    rho_cruise = 0.5258
#    Lift, Chord, Yle, Drag = read_aero_data("aquiladata1.txt", 50, V_cruise, rho_cruise)
    
    g = 9.80665
    
    #wing characteristics
    
    #UPDATE THESE VALUES
    d_fuselage_outside = p.d_fuselage_outside
    S = p.S
    b = p.b
    taper = p.taper
    c_root = p.root_chord
    c_tip = taper*c_root
    x_strut = perc_strut*b/2
    x_pod = perc_pod*b/2
    x_engine = perc_engine*b/2
    theta = atan(d_fuselage_outside/x_strut)
    
    thrust_per_engine = tot_thrust/2
    #wing box characteristics
    #E = 70*10**9
    #Iavg = 10**(-6)
    #offsetmax = 0.1*b/2
    
    #Weights
    W_empty = 9301 * g
    W_engine = (2.02 + 9.68)/100*W_empty/2
    W_wing = 13.85/100*W_empty/2
    W_pod = (5.09/100*W_empty + 1576 * g)/2
    #W_empty = 0                 
    #W_engine = 0
    #W_pod = 0
    Li = Lift
    di = b/2/len(Li)
    #Ii = np.ones(len(Li)+1)*10**(-4)
    Ii = Izz
    
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

        
    for i in range(len(Li)):
        dii[i+1] = xi[i+1] - xi[i]
        Si[i] = (ci[i]+ci[i+1])/2*di
        Wi[i] = W_wing/S*Si[i]/di
    
        
    momentLi = []
    momentWi = []
    momentDi = []
    for i in range(len(Li)+1):
        xset = xi[i]
        momentL = SumzML(Li, Yle, xset)
        momentD = SumzML(Drag, Yle, xset)
        momentW = SumzM(Wi, xi, xset, 2)
        momentLi.append(momentL)
        momentWi.append(momentW)
        momentDi.append(momentD)
    
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
    
    
    ymr = [0]
    yrz = [0]
    yfsz = [0]
    yT = [0]
    yD = [0]
    vymr = [0]
    vyrz = [0]
    vyfsz = [0]
    vyT = [0]
    vD = [0]
    
    
    for i in range(1,len(Drag)+1):
        if i == 50:
            ymr.append((dii[i])/(2*E*Iyy[i]) + ymr[i-1])
            yrz.append((dii[i])/(2*E*Iyy[i])*xi[i] + yrz[i-1])
            yD.append(-momentDi[i]*(dii[i])/(2*E*Iyy[i]) + yD[i-1])
            yT.append((dii[i])/(2*E*Iyy[i])*(xi[i]-x_engine) + yT[i-1])
            yfsz.append((dii[i])/(2*E*Iyy[i])*(xi[i]-x_strut) + yfsz[i-1])

        else:
            ymr.append((dii[i]+dii[i+1])/(2*E*Iyy[i]) + ymr[i-1])
            yrz.append((dii[i]+dii[i+1])/(2*E*Iyy[i])*xi[i] + yrz[i-1])
            yD.append(-momentDi[i]*(dii[i]+dii[i+1])/(2*E*Iyy[i]) + yD[i-1])
            if xi[i] >= x_engine and xi[i] < x_strut:
                yT.append((dii[i]+dii[i+1])/(2*E*Iyy[i])*(xi[i]-x_engine) + yT[i-1])
                yfsz.append(0)
            elif xi[i] >= x_strut and xi[i] >= x_engine:
                yT.append((dii[i]+dii[i+1])/(2*E*Iyy[i])*(xi[i]-x_engine) + yT[i-1])
                yfsz.append((dii[i]+dii[i+1])/(2*E*Iyy[i])*(xi[i]-x_strut) + yfsz[i-1])
            elif xi[i] < x_strut and xi[i] < x_engine:
                yT.append(0)
                yfsz.append(0)
    
    
    for i in range(1,len(Li)+1):
        vymr.append((ymr[i]+ymr[i-1])/2*dii[i] + vymr[i-1])
        vyrz.append((yrz[i]+yrz[i-1])/2*dii[i] + vyrz[i-1])
        vD.append((yD[i]+yD[i-1])/2*dii[i] + vD[i-1])
        vyfsz.append((yfsz[i]+yfsz[i-1])/2*dii[i] + vyfsz[i-1])
        vyT.append((yT[i]+yT[i-1])/2*dii[i] + vyT[i-1])
    
    
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
    
    vset_end = b/2*0.1
    vset_strut = 0
    
    for i in range(len(Li)+1):
        xi[i] = i*di
        ci[i] = c_root - (c_root-c_tip)/(b/2)*xi[i]
    
    for i in range(len(Li)):
        FyW[i] = Wi[i]*(xi[i+1]-xi[i])
        FyL[i] = Li[i]
        MoW[i] = Wi[i]*(xi[i+1]-xi[i])*(xi[i]+(xi[i+1]-xi[i])/2)
        MoL[i] = Li[i]*(xi[-1] - Yle[i])
    
    indexend = len(Li)
    for i in range(len(Li)):
        if xi[i] >= x_strut:
            indexstrut = i
            break

    
    A = [[1, 0, -cos(theta), 0], [0, 1, -sin(theta), 0], [0, -(xi[-1]), sin(theta)*(xi[-1] - x_strut), 1], [0, -vry[indexstrut], -vfs[indexstrut] , -vmr[indexstrut]]]
    B = [[0], [W_engine + W_pod + sum(FyW) - sum(FyL)], [-W_engine*(xi[-1] - x_engine) - W_pod*(xi[-1] - x_pod) - sum(MoW) + sum(MoL)], [ -vset_strut + vwe[indexstrut]*W_engine + vwp[indexstrut]*W_pod + vL[indexstrut] + vW[indexstrut]]]
    C = np.matmul(np.linalg.inv(A),B)
    
    Frx = C[0][0]
    Fry = C[1][0]
    Fs = C[2][0]
    Mr = C[3][0]
    
    
    
    MoD = []
    MtD = []
    for i in range(len(Drag)):
        MoD.append(Drag[i] * Yle[i])
        MtD.append(Drag[i] * (xi[-1] - Yle[i]))
    
    
    A1 = [[1, 1, 0], [0, (x_strut), 1], [vyrz[indexstrut], vyfsz[indexstrut], vymr[indexstrut]]]
    B1 = [[-thrust_per_engine + sum(Drag)],[-thrust_per_engine*(xi[-1] - x_engine) + sum(MoD)], [vset_strut - vD[indexstrut] - vyT[indexstrut]*thrust_per_engine]]
    C1 = np.matmul(np.linalg.inv(A1), B1)
    Frz = C1[0][0]
    Fsz = C1[1][0]
    Mry = C1[2][0]    
    
    
    #moment diagram
    momentzi = []
    for i in range(len(Li)+1):
        xset = xi[i]
        if xset >= x_engine and xset >= x_strut and xset >= x_pod:
            momi = -Mr + Fry*xset - W_engine*(xset - x_engine) - Fs*sin(theta)*(xset - x_strut) - W_pod*(xset - x_pod) + momentLi[i] + momentWi[i]
            momentzi.append(momi)
        elif xset >= x_engine and xset >= x_strut and xset < x_pod:
            momi = -Mr + Fry*xset - W_engine*(xset - x_engine) - Fs*sin(theta)*(xset - x_strut) + momentLi[i] + momentWi[i]
            momentzi.append(momi)
        elif xset >= x_engine and xset < x_strut and xset < x_pod:
            momi = -Mr + Fry*xset - W_engine*(xset - x_engine) + momentLi[i] + momentWi[i]
            momentzi.append(momi)
        elif xset < x_engine and xset < x_strut and xset < x_pod:
            momi = -Mr + Fry*xset + momentLi[i] + momentWi[i]
            momentzi.append(momi)
        
    
    #shear diagram
    shearyi = []
    for i in range(len(Li)+1):
        xset = xi[i]
        shearL = sum(Li[:(i)])
        shearW = SumzV(Wi, xi, xset, 2)
        if xset >= x_engine and xset >= x_strut and xset >= x_pod:
            shri = Fry - W_engine - W_pod - Fs*sin(theta) + shearL + shearW
        elif xset >= x_engine and xset >= x_strut and xset < x_pod:
            shri = Fry - W_engine - Fs*sin(theta) + shearL + shearW
        elif xset >= x_engine and xset < x_strut and xset < x_pod:
            shri = Fry - W_engine + shearL + shearW
        elif xset < x_engine and xset < x_strut and xset < x_pod:
            shri = Fry + shearL + shearW
        shearyi.append(shri)
        
    shearzi = []
    for i in range(len(Drag)+1):
        xset = xi[i]
        shearD = sum(Drag[:i])
        if xset >= x_engine and xset >= x_strut:
            shrzi = Frz + thrust_per_engine + Fsz - shearD
            shearzi.append(-shrzi)
        elif xset >= x_engine and xset < x_strut:
            shrzi = Frz + thrust_per_engine - shearD 
            shearzi.append(-shrzi)
        elif xset < x_engine and xset < x_strut:
            shrzi = Frz - shearD
            shearzi.append(-shrzi)
        
    momentyi = []
    for i in range(len(Drag)+1):
        xset = xi[i]
        momentD = SumzML(Drag, Yle, xset)
#        print(momentD)
        if xset >= x_engine and xset >= x_strut:
            momi = -Mry - Frz*xset - thrust_per_engine*(xset - x_engine) - Fsz*(xset - x_strut) + momentD
            momentyi.append(momi)
        elif xset >= x_engine and xset < x_strut:
            momi = -Mry - Frz*xset - thrust_per_engine*(xset - x_engine) + momentD
            momentyi.append(momi)
        elif xset < x_engine and xset < x_strut:
            momi = -Mry - Frz*xset + momentD
            momentyi.append(momi)
#    plt.plot(momentyi)
    
    thetayi = [0]
    vzi = [0]
    thetazi = [0]
    vyi = [0]
    vny = [0]
    vnz = [0]
    for i in range(1,(len(Li)+1)):
        thetayi.append(thetayi[i-1] + 1/2*(momentzi[i]/(E*Ii[i]) + momentzi[i-1]/(E*Ii[i-1]))*(dii[i]))
        
    for i in range(1,(len(Li)+1)):
        vyi.append(vyi[i-1] + 1/2*(thetayi[i] + thetayi[i-1])*dii[i])
    
    for i in range(1,(len(Drag)+1)):
        thetazi.append(thetazi[i-1] + 1/2*(momentyi[i]/(E*Iyy[i]) + momentyi[i-1]/(E*Iyy[i-1]))*(dii[i]))
        
    for i in range(1,(len(Drag)+1)):
        vzi.append(vzi[i-1] + 1/2*(thetazi[i] + thetazi[i-1])*dii[i])
    
    
    for n in range(1, len(Li)+1):
        vny.append(vmr[n]*Mr + vry[n]*Fry + vwe[n]*W_engine + vfs[n]*Fs + vwp[n]*W_pod + vW[n]+ vL[n] )
    
    for n in range(1, len(Drag)+1):
        vnz.append(-(vymr[n]*Mry + vyrz[n]*Frz + vyT[n]*thrust_per_engine + vyfsz[n]*Fsz + vD[n] ))
#    plt.plot(xi, vnz)
#    plt.plot(xi, vzi)    
    
    #plt.plot(xi, vii)
    Mrz = Mr
    return Frx, Fry, Fs, Mrz, Frz, Fsz, Mry, momentyi, momentzi, shearyi, shearzi, vyi, vny, vzi, vnz, xi, theta



#TESTING
tot_thrust = 20000
V_cruise = p.V_cruise
rho_cruise = p.rho
perc_engine = 0.15
perc_strut = 0.5
perc_pod = 0.5

lengthdata = 50
Lift, Chord, Yle, Drag = read_aero_data("wing/aquiladata1.txt", lengthdata, V_cruise, rho_cruise)
Frx, Fry, Fs, Mrz, Frz, Fsz, Mry, momentyi, momentzi, shearyi, shearzi, vyi, vny, vzi, vnz, xi, theta = CallForces(Lift, Yle, Drag, tot_thrust, np.ones(len(Lift)+1)*10**(-4), np.ones(len(Lift)+1)*10**(-4) , 70*10**9, perc_engine, perc_strut, perc_pod)
Forces1 = [['Frx = ', Frx], ['Fry = ', Fry], ['Fs = ', Fs], ['Mr = ', Mrz]]
Forces2 = [['Frz = ', Frz], ['Fsz = ', Fsz], ['Mry = ', Mry]]
#TESTING