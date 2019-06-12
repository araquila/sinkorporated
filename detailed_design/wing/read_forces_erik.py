# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:07:56 2019

@author: Stijn
"""
import numpy as np
from matplotlib import pyplot as plt

def feet_to_meter(length_in_feet):
    length_in_meter = length_in_feet / 3.28084
    return length_in_meter

def feetsquared_to_metersquared(area_in_feetsquared):
    area_in_metersquared = area_in_feetsquared / 10.7639
    return area_in_metersquared


#    V_cruise = 184.84
#    rho_cruise = 0.5258
    



V_cruise = 184.84
rho_cruise = 0.5258

z = [0,1,2,3]

for a in z:

    if a == 0:
        skipheader = 20
        lengthdata = 100
        name = 'wing 1'
        random = 12
    elif a == 1:
        skipheader = 135
        lengthdata = 100
        name = 'wing 2'
        random = 9
    elif a == 2:
        skipheader = 250
        lengthdata = 20
        name = 'strut 1'
        random = 4
    elif a == 3:
        skipheader = 285
        lengthdata = 20
        name = 'strut 2'
        random = 0
    elif a == 4:
        skipheader = 20
        lengthdata = 100
        name = 'wing 1'
        random = 5
    filename = 'datastrut4.txt'

    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    lengthlines = len(lines)
    
    #    lengthdata = 50
    
    table = np.genfromtxt(filename,delimiter=None , skip_header=skipheader, skip_footer=(lengthlines-lengthdata-skipheader-random))
    
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
    Cd = []
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
        cd012 = 0.00617
        cd015 = 0.00657
        e = 0.95
        A = 20.
#        if a == 0:
#            AR = (15.686/lengthdata)**2/feetsquared_to_metersquared(table[i][3])
#            print(AR)
        cd0 = cd015-(cd015-cd012)/(i+1)
        Cd.append(cd0+table[i][7]**2/np.pi/e/A)
    
    if a == 0:
        plt.figure()
        plt.plot(Yle,Cd)
        plt.show()
        
    Lift = []
    Drag = []
    for j in range(len(ji)):
        Lift.append(1/2*rho_cruise*Area[j]*cl[j]*V_cruise**2)
    for l in range(len(ji)):
        Drag.append(1/2*rho_cruise*Area[l]*Cd[l]*V_cruise**2)
    lift = sum(Lift)
    drag = sum(Drag)
    
    print(lift)
    print(drag)
    print(lift/drag)


