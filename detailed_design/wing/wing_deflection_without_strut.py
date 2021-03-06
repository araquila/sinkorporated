# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 13:38:27 2019

@author: Stijn
"""

from math import *
import numpy as np
import parameters as p
from matplotlib import pyplot as plt

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
    table = np.genfromtxt(filename,delimiter=None , skip_header=skipheader, skip_footer=(lengthlines-lengthdata-skipheader-12))
    
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
    AeroMoment = []
    for i in range(len(ji)):
#        Lift.append(1/2*rho_cruise*Area[i]*cl[i]*V_cruise**2)
#        Drag.append(1/2*rho_cruise*Area[i]*cd[i]*V_cruise**2)
#        AeroMoment.append(1/2*rho_cruise*Area[i]*cm_c4[i]*V_cruise**2*Chord[i])
##        
        Lift.append(1.5*2.5*1/2*rho_cruise*Area[i]*cl[i]*V_cruise**2)
        Drag.append(1.5*2.5*1/2*rho_cruise*Area[i]*cd[i]*V_cruise**2)
        AeroMoment.append(1.5*2.5*1/2*rho_cruise*Area[i]*cm_c4[i]*V_cruise**2*Chord[i])

    return Lift, Chord, Yle, Drag, AeroMoment


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
    x_strut = p.strut_pos_perc*p.b/2
    x_pod = p.pod_pos_perc*p.b/2
    x_engine = p.engine_pos_perc*p.b/2
    theta = atan(d_fuselage_outside/x_strut)
    
    thrust_per_engine = tot_thrust/2
    #wing box characteristics
    #E = 70*10**9
    #Iavg = 10**(-6)
    #offsetmax = 0.1*b/2
    
    #Weights
    W_empty = p.W_empty
    W_engine = p.W_engine+0.5*p.W_fuel_system
    W_wing = p.W_wing
    W_pod = p.W_fuel+p.W_pod+0.5*p.W_fuel_system
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
        if i == 100:
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
        if i == 100:
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
    vset_strut = 0.0
    
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
    
    
    A1 = [[1, 1, 0], [xi[-1], (xi[-1] - x_strut), 1], [vyrz[indexstrut], vyfsz[indexstrut],vymr[indexstrut]]]
    B1 = [[-thrust_per_engine + sum(Drag)],[-thrust_per_engine*(xi[-1] - x_engine) + sum(MtD)], [vset_strut - vD[indexstrut] - vyT[indexstrut]*thrust_per_engine]]
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
        thetazi.append((thetazi[i-1] + 1/2*(momentyi[i]/(E*Iyy[i]) + momentyi[i-1]/(E*Iyy[i-1]))*(dii[i])))
        
    for i in range(1,(len(Drag)+1)):
        vzi.append((vzi[i-1] + 1/2*(thetazi[i] + thetazi[i-1])*dii[i]))
    
    for i in range(1,(len(Drag)+1)):
        vzi[i] = -vzi[i]
    
    
    for n in range(1, len(Li)+1):
        vny.append(vmr[n]*Mr + vry[n]*Fry + vwe[n]*W_engine + vfs[n]*Fs + vwp[n]*W_pod + vW[n]+ vL[n] )
    
    for n in range(1, len(Drag)+1):
        vnz.append((vymr[n]*Mry + vyrz[n]*Frz + vyT[n]*thrust_per_engine + vyfsz[n]*Fsz + vD[n] ))
#    plt.plot(xi, vnz)
#    plt.plot(xi, vzi)    
    
    #plt.plot(xi, vii)
    Mrz = Mr
    return Frx, Fry, Fs, Mrz, Frz, Fsz, Mry, momentyi, momentzi, shearyi, shearzi, vyi, vny, vzi, vnz, xi, theta


def CallForcesWithoutStrut(Lift, Yle, Drag, tot_thrust, Iyy, Izz, E, perc_engine, perc_strut, perc_pod):
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
    x_strut = p.strut_pos_perc*p.b/2
    x_pod = p.pod_pos_perc*p.b/2
    x_engine = p.engine_pos_perc*p.b/2
    theta = atan(d_fuselage_outside/x_strut)
    
    thrust_per_engine = tot_thrust/2
    #wing box characteristics
    #E = 70*10**9
    #Iavg = 10**(-6)
    #offsetmax = 0.1*b/2
    
    #Weights
    W_empty = p.W_empty
    W_engine = p.W_engine+0.5*p.W_fuel_system
    W_wing = p.W_wing
    W_pod = p.W_fuel+p.W_pod+0.5*p.W_fuel_system
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
        if i == 100:
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
        if i == 100:
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
    vset_strut = 0.0
    
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

    
    Awithout = [[1, 0,  0], [0, 1, 0], [0, -(xi[-1]),  1]]
    Bwithout = [[0], [W_engine + W_pod + sum(FyW) - sum(FyL)], [-W_engine*(xi[-1] - x_engine) - W_pod*(xi[-1] - x_pod) - sum(MoW) + sum(MoL)]]
    Cwithout = np.matmul(np.linalg.inv(Awithout),Bwithout)
    
    Frx = Cwithout[0][0]
    Fry = Cwithout[1][0]
    Mr = Cwithout[2][0]
    
    MoD = []
    MtD = []
    for i in range(len(Drag)):
        MoD.append(Drag[i] * Yle[i])
        MtD.append(Drag[i] * (xi[-1] - Yle[i]))
    
    
    
    
    A1without = [[1, 0], [xi[-1],  1]]
    B1without = [[-thrust_per_engine + sum(Drag)],[-thrust_per_engine*(xi[-1] - x_engine) + sum(MtD)]]
    C1without = np.matmul(np.linalg.inv(A1without), B1without)
    Frz = C1without[0][0]
    Mry = C1without[1][0]    

    Fs = 0
    Fsz = 0
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
        thetazi.append((thetazi[i-1] + 1/2*(momentyi[i]/(E*Iyy[i]) + momentyi[i-1]/(E*Iyy[i-1]))*(dii[i])))
        
    for i in range(1,(len(Drag)+1)):
        vzi.append((vzi[i-1] + 1/2*(thetazi[i] + thetazi[i-1])*dii[i]))
    
    
    
    for n in range(1, len(Li)+1):
        vny.append(vmr[n]*Mr + vry[n]*Fry + vwe[n]*W_engine + vfs[n]*Fs + vwp[n]*W_pod + vW[n]+ vL[n] )
    
    for n in range(1, len(Drag)+1):
        vnz.append((vymr[n]*Mry + vyrz[n]*Frz + vyT[n]*thrust_per_engine + vD[n] ))
#    plt.plot(xi, vnz)
#    plt.plot(xi, vzi)    
    
    #plt.plot(xi, vii)
    Mrz = Mr
    
    return Frx, Fry, Fs, Mrz, Frz, Fsz, Mry, momentyi, momentzi, shearyi, shearzi, vyi, vny, vzi, vnz, xi, theta


#TESTING
tot_thrust = p.tot_thrust
V_cruise = p.V_cruise
rho_cruise = p.rho
perc_engine = p.engine_pos_perc
strut_pos_perc = 0
perc_pod = p.pod_pos_perc

import wing.section_properties as sp
import wing.wing_deflection2 as wd2 
from wing.shear_stress import max_shear_stress

### DISCRETIZATION OF THE WINGBOX ###
n = 101
x_pos = np.linspace(0,p.b/2,n)


### OBTAIN CROSSECTIONAL PROPERTIES ###
Izz_list = []
Iyy_list = []
first_moment_of_area_list = []
area_list = []
y_max_list = []
hi = []
bi = []


for x in x_pos:
    Izz_list.append(sp.I_zz_wingbox(x))
    Iyy_list.append(sp.I_yy_wingbox(x))
    first_moment_of_area_list.append(sp.first_moment_of_area_y(x))
    area_list.append(sp.cross_sectional_area(x))
    y_max_list.append(sp.y_max(x))
    hi.append(sp.height_wingbox(x))
    bi.append(sp.width_wingbox(x))



lengthdata = 100
Lift, Chord, Yle, Drag, AeroMoment = read_aero_data("wing/datastrut5.txt", lengthdata, V_cruise, rho_cruise)
Frx, Fry, Fs, Mrz, Frz, Fsz, Mry, momentyi, momentzi, shearyi, shearzi, vyi, vny, vzi, vnz, xi, theta = CallForces(Lift, Yle, Drag, tot_thrust, Iyy_list, Izz_list , 70*10**9, perc_engine, perc_strut, perc_pod)
Forces1 = [['Frx = ', Frx], ['Fry = ', Fry], ['Fs = ', Fs], ['Mr = ', Mrz]]
Forces2 = [['Frz = ', Frz], ['Fsz = ', Fsz], ['Mry = ', Mry]]
#TESTING


### NORMAL STRESS CALCULATOR ###
def normal_stress(x,y,moment_z,moment_y,normal_force,I_zz,I_yy,area):
    """Returns bending moment and normal force stress""" 
    
    z = sp.width_wingbox(x)/2
    
    #-my/I according to the formula, makes sense because for a positive Mz the top skin will be in compression
    moment_z_upperskin = -moment_z*(sp.height_wingbox(x)-abs(y))/I_zz
    moment_z_lowerskin = -moment_z*y/I_zz
    
    moment_y_rightflange = -moment_y*z/I_yy
    moment_y_leftflange = -moment_y*-z/I_yy
    
    area = sp.cross_sectional_area(x)
    
    normal_force_stress = normal_force/area
    
    normal_ru = moment_z_upperskin + moment_y_rightflange +     normal_force_stress 
    normal_lu = moment_z_upperskin + moment_y_leftflange +     normal_force_stress 
    normal_rl = moment_z_lowerskin + moment_y_rightflange +     normal_force_stress 
    normal_ll = moment_z_lowerskin + moment_y_leftflange +  normal_force_stress 
    
#    if max(normal_ru,normal_lu,normal_rl,normal_ll) == normal_ru:
#        print("Max tension at right upper corner")
#    elif max(normal_ru,normal_lu,normal_rl,normal_ll) == normal_lu:
#        print("Max tension at left upper corner")
#    elif max(normal_ru,normal_lu,normal_rl,normal_ll) == normal_rl:
#        print("Max tension at right lower corner")
#    elif max(normal_ru,normal_lu,normal_rl,normal_ll) == normal_ll:
#        print("Max tension at left lower corner")
#        
#    
#    if min(normal_ru,normal_lu,normal_rl,normal_ll) == normal_ru:
#        print("Max compression at right upper corner")
#    elif min(normal_ru,normal_lu,normal_rl,normal_ll) == normal_lu:
#        print("Max compression at left upper corner")
#    elif min(normal_ru,normal_lu,normal_rl,normal_ll) == normal_rl:
#        print("Max compression at right lower corner")
#    elif min(normal_ru,normal_lu,normal_rl,normal_ll) == normal_ll:
#        print("Max compression at left lower corner")
#        
#        
#    print("x: ",x)
#    print("Neutral axis: y =",sp.centroid_y(x),"(",sp.centroid_y(x)/sp.height_wingbox(x)*100,"%)",)
#    print("Maximum tension: ",max(normal_ru,normal_lu,normal_rl,normal_ll)/10**6,"MPa")
#    print("Maximum compression: ", min(normal_ru,normal_lu,normal_rl,normal_ll)/10**6,"MPa")
#    print("Normal force stress: ",normal_force_stress,"MPa")
#    print("")
    
    return normal_ru,normal_lu,normal_rl,normal_ll

 
def von_mises(x,y,moment_z,moment_y,normal_force,I_zz,I_yy,area, tau_max):
    """Returns maximum Von Mises stress"""
    
    z = sp.width_wingbox(x)/2
    
    #-my/I according to the formula, makes sense because for a positive Mz the top skin will be in compression
    moment_z_upperskin = -moment_z*(sp.height_wingbox(x)-abs(y))/I_zz
    moment_z_lowerskin = -moment_z*y/I_zz
    
    moment_y_rightflange = -moment_y*z/I_yy
    moment_y_leftflange = -moment_y*-z/I_yy
    
    area = sp.cross_sectional_area(x)
    
    normal_force_stress = normal_force/area
    
    stress_x_lower = moment_z_lowerskin + normal_force_stress
    stress_x_upper = moment_z_upperskin + normal_force_stress
    
    stress_y_right = moment_y_rightflange
    stress_y_left = moment_y_leftflange
    
    
    vm_ll_1 = (stress_x_lower + stress_y_left)/2 + np.sqrt(((stress_x_lower-stress_y_left)/2)**2 + tau_max**2)
    vm_ll_2 = (stress_x_lower + stress_y_left)/2 - np.sqrt(((stress_x_lower-stress_y_left)/2)**2 + tau_max**2)
    vm_ll = np.sqrt(vm_ll_1**2 + vm_ll_2**2 -vm_ll_1*vm_ll_2 + 3*tau_max**2)
    
    vm_lr_1 = (stress_x_lower + stress_y_right)/2 + np.sqrt(((stress_x_lower-stress_y_right)/2)**2 + tau_max**2)
    vm_lr_2 = (stress_x_lower + stress_y_right)/2 - np.sqrt(((stress_x_lower-stress_y_right)/2)**2 + tau_max**2)
    vm_lr = np.sqrt(vm_lr_1**2 + vm_lr_2**2 -vm_lr_1*vm_lr_2 + 3*tau_max**2)
    
    vm_ur_1 = (stress_x_upper + stress_y_right)/2 + np.sqrt(((stress_x_upper-stress_y_right)/2)**2 + tau_max**2)
    vm_ur_2 = (stress_x_upper + stress_y_right)/2 - np.sqrt(((stress_x_upper-stress_y_right)/2)**2 + tau_max**2)
    vm_ur = np.sqrt(vm_ur_1**2 + vm_ur_2**2 -vm_ur_1*vm_ur_2 + 3*tau_max**2)
    
    vm_ul_1 = (stress_x_upper + stress_y_left)/2 + np.sqrt(((stress_x_upper-stress_y_left)/2)**2 + tau_max**2)
    vm_ul_2 = (stress_x_upper + stress_y_left)/2 - np.sqrt(((stress_x_upper-stress_y_left)/2)**2 + tau_max**2)
    vm_ul = np.sqrt(vm_ul_1**2 + vm_ul_2**2 -vm_ul_1*vm_ul_2 + 3*tau_max**2)
  
    
    return vm_ll,vm_lr,vm_ur,vm_ul
 


def skin_buckling_stress(x):
    """Returns compressive stress at which buckling will occur"""
    
    k = 4 # to be determined based on the "final" stiffener and rib spacing
    b = sp.top_spacing
    t = p.t_sheet
    
    return k*np.pi**2*p.E_sheet/(12*(1-p.poisson_ratio_al2014**2)) * (t/b)**2


def shear_skin_buckling_stress(x):
    """Returns compressive stress at which buckling will occur"""
    
    k = 5.35 # to be determined based on the "final" stiffener and rib spacing
    b = sp.top_spacing
    t = p.t_sheet
    
    return k*np.pi**2*p.E_sheet/(12*(1-p.poisson_ratio_al2014**2)) * (t/b)**2



#print("Skin buckling limit: ", skin_buckling_stress(0))

def critical_column_buckling(x):
    """Critical column buckling stress on the panel""" 
    
    l_eff = p.rib_spacing
    
    return (np.pi**2*p.E_sheet*sp.I_yy_wingbox(x))/(l_eff**2*sp.cross_sectional_area(x))/10**6


def column_buckling_stiffener(x):
    """Critical column buckling stress on the stiffener"""
    
    l_eff = p.rib_spacing
    
    return (np.pi**2*p.E_compressive_2099*sp.I_yy_hat)/(l_eff**2*sp.A_hat)/10**6


def critical_crippling_stiffener(x):
    """Critical crippling stress for aluminium top stiffener in MPa"""
    
    alpha = 0.8
    n = 0.6
    yield_stress = p.ultimate_yield_strength_2099
    
    
    def stress_cc(K,b):
        return K*(np.pi**2*p.E_compressive_2099*(sp.t_hat/b)**2/(12*(1-p.poisson_ratio_al2014**2)))
    
    t_sheet = sp.t_hat
    
    #areas
    area_a = sp.a*t_sheet
    area_b = (sp.b-t_sheet)*t_sheet
    area_c = (sp.c-t_sheet)*t_sheet
    
    #critical crippling stress
    cc_a = stress_cc(4,sp.a)
    cc_b = stress_cc(0.425,sp.b)
    cc_c = stress_cc(4,sp.c)
    
    ratio_a = alpha*(cc_a/yield_stress)**(1-n)
    ratio_b = alpha*(cc_b/yield_stress)**(1-n)
    ratio_c = alpha*(cc_c/yield_stress)**(1-n)

    total_crippling = yield_stress*(2*(ratio_a*area_a+ratio_b*area_b)+ratio_c*area_c) / (2*area_a + 2*area_b + area_c)
    
    return total_crippling


def critical_panel_buckling(x):
    """Returns critical panel buckling stress"""
    C = 6.98
    v = p.poisson_ratio_al2014
    
    crippling = critical_crippling_stiffener(x)
    
    we = p.t_sheet * np.sqrt(C*np.pi**2/(12*(1-v**2))) * np.sqrt(p.E_sheet/crippling)
    
    
    return we


def crack_length_sheet(stress):
    """Fast fracture crack length for given stress"""
    
    Kic  = p.fracture_toughness_2195
    a = Kic/((stress*10**6)**2 * np.pi)
    
    return a
    
### CALCULATE MAXIMUM SHEAR STRESS PER SECTION ###
tau_max = max_shear_stress(Lift, Drag, AeroMoment, Chord, shearyi, shearzi, hi, bi, Izz_list, Iyy_list)
 
### CALCULATE SHEAR AND NORMAL STRESS ###
normal_ru_list= []
normal_lu_list= []
normal_rl_list= []
normal_ll_list= []

vm_ll_list = []
vm_lr_list = []
vm_ur_list = []
vm_ul_list = []


for i in range(len(x_pos)):
    if x_pos[i] < p.x_strut:
        F_normal = -Frx
    else:
        F_normal = 0
    
    normal_ru, normal_lu, normal_rl, normal_ll = normal_stress(x_pos[i],y_max_list[i],momentzi[i],momentyi[i],F_normal,Izz_list[i],Iyy_list[i],area_list[i])
    normal_ru_list.append(normal_ru/10**6)
    normal_lu_list.append(normal_lu/10**6)
    normal_rl_list.append(normal_rl/10**6)
    normal_ll_list.append(normal_ll/10**6)
    
    vm_ll,vm_lr,vm_ur,vm_ul = von_mises(x_pos[i],y_max_list[i],momentzi[i],momentyi[i],F_normal,Izz_list[i],Iyy_list[i],area_list[i],tau_max[i])
    vm_ll_list.append(vm_ll/10**6)
    vm_lr_list.append(vm_lr/10**6)
    vm_ul_list.append(vm_ul/10**6)
    vm_ur_list.append(vm_ur/10**6)
    
    
    
    

# Error correction
error_lower = normal_ll_list[-1]
error_upper = normal_lu_list[-1]
    
for i in range(len(x_pos)):
    normal_ru_list[i] -= error_upper
    normal_lu_list[i] -= error_upper
    normal_rl_list[i] -= error_lower
    normal_ll_list[i] -= error_lower
    
 
#plt.plot(x_pos,vm_ll_list,marker='x')
#plt.plot(x_pos,vm_lr_list)
#plt.plot(x_pos,vm_ul_list)
#plt.plot(x_pos,vm_ur_list)
#plt.show()


plt.figure()
plt.xlim([0,p.b/2])

#### MOMENT AND SHEAR DIAGRAM ###
#plt.figure(2,figsize = (8,6))
#plt.xlabel('Location along the length of the wingbox [m]',size='large')
#plt.ylabel('Moment in z-axis [kNm]',size='large')
#momentzi = np.array(momentzi)
#plt.plot(x_pos, momentzi/1000, 'b') 

##plt.figure(3,figsize = (8,6))
#plt.xlabel('Location along the length of the wingbox [m]',size='large')
#shearyi = np.array(shearyi)
#plt.ylabel('Shear force in y-direction [kN]',size='large')
#plt.plot(x_pos, shearyi/1000,'r')   

#plt.figure(4,figsize = (8,6))
#plt.xlabel('Location along the length of the wingbox [m]',size='large')
#momentyi = np.array(momentyi)
#plt.ylabel('Moment in around y-axis [kNm]',size='large')
#plt.plot(x_pos, momentyi/1000,'b')

#plt.figure(5,figsize = (8,6))
#plt.xlabel('Location along the length of the wingbox [m]',size='large')
#shearzi = np.array(shearzi)
#plt.ylabel('Shear force in z-direction [kN]',size='large')
#plt.plot(x_pos, shearzi/1000,'r') 

#### PLOT SHEAR STRESS ###
#plt.figure(6,figsize = (8,6))
plt.xlabel('Location along the length of the wingbox [m]',size='large')
plt.ylabel('Shear stress [MPa]',size='large')
plt.plot(x_pos, tau_max,'r')

## PLOT NORMAL STRESS AT THE FOUR CORNERS
plt.figure(7,figsize = (8,6))
plt.xlabel('Location along the length of the wingbox [m]',size='large')
plt.ylabel('Normal stress [MPa]',size='large')
plt.plot(x_pos, normal_ru_list, 'r', label='Right-up')
plt.plot(x_pos, normal_lu_list, 'g', label='Left-up')
plt.plot(x_pos, normal_rl_list, 'b', label='Right-bottom')
plt.plot(x_pos, normal_ll_list, 'y', label='Left-bottom')
plt.legend(loc="best", fontsize="large")
plt.xlim([0,p.b/2])


plt.show()

max_compressive_stress = min(normal_lu_list) #*p.safety_factor_compression
max_tensile_stress = max(normal_rl_list) #*p.safety_factor_tension


print("Top stiffeners: ",p.n_upper_skin_wingbox)
print("Bottom stiffeners: ",p.n_lower_skin_wingbox)
print("Skin thickness: ",p.t_sheet*1000, "mm")
print("Total weight: ",sp.total_weight,"kg")
print("")


print("Max shear: ",max(tau_max),"MPa")
print("Max compressive: ",max_compressive_stress,"MPa")
print("Max tensile: ",max_tensile_stress,"MPa")
print("Max Von Mises: ", max(max(vm_ll_list),max(vm_lr_list),max(vm_ul_list),max(vm_ur_list)))
print("")

print("Skin buckling limit: ",skin_buckling_stress(strut_pos_perc*p.b/2)/10**6,"MPa")

print("Panel column buckling limit: ",critical_column_buckling(strut_pos_perc*p.b/2),"MPa")

print("Stiffener column buckling limit: ",column_buckling_stiffener(strut_pos_perc*p.b/2),"MPa")

print("Critical crippling stress of the hat stiffener: ",critical_crippling_stiffener(strut_pos_perc*p.b/2)/10**6,"MPa")

print("Shear stress buckling limit: ",shear_skin_buckling_stress(strut_pos_perc*p.b/2)/10**6,"MPa (still to be determined where the maximum shear stress occurs spanwise)")

print("")
print("Tests: ")
#print("Crack length: ",crack_length_sheet(max(max_tensile_stress,max_compressive_stress)),"m")

if abs(max_compressive_stress) < abs(skin_buckling_stress(strut_pos_perc*p.b/2)/10**6):
    print("Skin buckling passed")
else:
    print("Failure on skin buckling")
    
if abs(max_compressive_stress) < abs(critical_column_buckling(strut_pos_perc*p.b/2)):
    print("Panel column buckling passed")
else:
    print("Failure on panel column buckling")
    
if abs(max_compressive_stress) < abs(column_buckling_stiffener(strut_pos_perc*p.b/2)):
    print("Stiffener column buckling passed")
else:
    print("Failure on stiffener column buckling")
    
if abs(max_compressive_stress)*10**6 < abs(critical_crippling_stiffener(strut_pos_perc*p.b/2)):
    print("Cripple limit of the stiffener passed")
else:
    print("Failure on stiffener crippling")
    
if shear_skin_buckling_stress(strut_pos_perc*p.b/2)/4 > max(tau_max):
    print("Shear buckling passed")
else:
    print("Failure on shear buckling")

print("")

if (max_tensile_stress > p.ult_yield_strength_2195/10**6) :
    print("Failure on tensile yielding at the strut",max_tensile_stress/(p.ult_yield_strength_2195/10**6)*100)
else:
    print("Max tensile stress",max_tensile_stress/(p.ult_yield_strength_2195/10**6)*100,"% of the yield strength")

    
if abs(max_compressive_stress) > (p.ult_compressive_strength_2195/10**6):
    print("Failure on compressive yielding at the strut",-max_compressive_stress/(p.ult_compressive_strength_2195/10**6)*100)
else:
    print("Max compressive stress",-max_compressive_stress/(p.ult_compressive_strength_2195/10**6)*100,"% of the yield strength")

if p.ultimate_shear_stress_al2024 > max(tau_max)*10**6:
    print("Shear limit of the material passed: ",max(tau_max)*10**6/p.ult_shear_strength_2195 *100,"% of the max")
else: 
    print("Shear failure of the material")


plt.figure()
# Plot lines
plt.plot(xi, vyi)
plt.xlabel("Location along the length of the wing", size="large")
plt.ylabel("Wing deflection in y [m]", size="large")

# Limits of the axes
#plt.xlim([0,10])
#plt.ylim([0,20])

# Create Legend
#plt.legend(loc="best", fontsize="large")

# Show plot
plt.show()

plt.figure()
# Plot lines
plt.plot(xi, vzi)
plt.xlabel("Location along the length of the wing", size="large")
plt.ylabel("Wing deflection in z [m]", size="large")

# Limits of the axes
#plt.xlim([0,10])
#plt.ylim([0,20])

# Create Legend
#plt.legend(loc="best", fontsize="large")

# Show plot
plt.show()
