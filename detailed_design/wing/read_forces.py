# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:07:56 2019

@author: Stijn
"""
import numpy as np

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

V_cruise = 184.84
rho_cruise = 0.5258

Lift, Chord, Yle = read_aero_data("aquiladata1.txt", 50, V_cruise, rho_cruise)
