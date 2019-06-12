import matplotlib.pyplot as plt
# =============================================================================
## ## Constants ##
#import sys
#import os
#
#sys.path.append(os.getcwd())
#from detailed_design.parameters import *
# =============================================================================
#import parameters as p
import numpy as np
import parameters as p
#constants
n_passenger = p.n_passenger
n_seats_abreast = p.n_seats_abreast
n_rows = n_passenger / n_seats_abreast
seat_pitch = p.seat_pitch
#Dimensions
l_fuselage = p.l_fuselage

#MASSES
W_fuel = 1576
M_OEW = 9785
M_fusgroup = 5710
M_winggroup = M_OEW - M_fusgroup
M_total_cargo = p.M_total_cargo
weight_passenger = p.M_passenger
cargo_passenger = p.M_cargo

#CG
CG_fusgroup = 0.42 * l_fuselage
CG_cargo = p.l_nose + p.l_cabin - 1
cg_margin = 0.02
CG_winggroup = np.linspace(8,14,100)
xlemac = []
minCG = []
maxCG = []
CGmacmin = []
CGmacmax = []
for j in range(len(CG_winggroup)):
    ## CARGO ##

    CG_fuel = CG_winggroup[j]
    CG_OEW = (CG_fusgroup * M_fusgroup + CG_winggroup[j] * M_winggroup) / (M_OEW)
    CG=[CG_OEW]
    weight = [M_OEW]


    CG_cargo_loaded = (CG_OEW * M_OEW + CG_cargo * M_total_cargo) / (M_OEW + M_total_cargo)
    CG.append(CG_cargo_loaded)
    weight.append(M_OEW + M_total_cargo)

    ## PASSENGERS ##
    #AISLE LOADING#
    x_rowaft = [p.l_nose + p.l_lavatory + (n_rows) * seat_pitch]
    x_rowfront = [p.l_nose + p.l_lavatory +seat_pitch]

    #start loading aft#
    for i in range(int(n_rows-1)):
        x_rowaft.append(x_rowaft[-1] - seat_pitch)
        x_rowfront.append(x_rowfront[-1] + seat_pitch)

    for i in range(int(n_rows)):
        weight.append(weight[-1] + (n_seats_abreast / 2) * weight_passenger)
        CG.append((weight[-2]*CG[-1] + (n_seats_abreast / 2) * weight_passenger * x_rowaft[i])/weight[-1])

    for i in range(int(n_rows)):
        weight.append(weight[-1] + (n_seats_abreast / 2) * weight_passenger)
        CG.append((weight[-2]*CG[-1]+(n_seats_abreast / 2) * weight_passenger * x_rowaft[i])/weight[-1])

    for i in range(2):
        weight.append(weight[-1] + W_fuel / 2)
        CG.append((weight[-2] * CG[-1] + (W_fuel / 2) * CG_fuel) / (weight[-1]))

    weight.append(weight[1])
    CG.append(CG_cargo_loaded)
    for i in range(int(n_rows)):
        weight.append(weight[-1] + (n_seats_abreast / 2) * weight_passenger)
        CG.append((weight[-2]*CG[-1]+(n_seats_abreast / 2) * weight_passenger * x_rowfront[i])/weight[-1])

    for i in range(int(n_rows)):
        weight.append(weight[-1] + (n_seats_abreast / 2) * weight_passenger)
        CG.append((weight[-2]*CG[-1]+(n_seats_abreast / 2) * weight_passenger * x_rowfront[i])/weight[-1])

    CG_mostfor = (1 - cg_margin) * min(CG)
    CG_mostaft = (1 + cg_margin) * max(CG)
    CG = np.array(CG)

# =============================================================================

# =============================================================================
    minCG.append(CG_mostfor)
    maxCG.append(CG_mostaft)
    x_lemac = CG_winggroup[j] - 0.35 * p.MAC
    x_lemac_l_f = x_lemac / l_fuselage
    CG_MACmin = (CG_mostfor - x_lemac) / p.MAC
    CG_MACmax = (CG_mostaft - x_lemac) / p.MAC
    xlemac.append(x_lemac_l_f)
    CGmacmin.append(CG_MACmin)
    CGmacmax.append(CG_MACmax)
    CG = (CG - x_lemac) / p.MAC
    
    print("Most forward CG =", CG_mostfor)
    print("Most afterward CG =", CG_mostaft)
    plt.scatter(CG,weight)
    plt.show()


    #if CG_winggroup[j] > 10 and CG_winggroup[j] < 11:
        #plt.scatter(CG, weight)
        #plt.show()
    #print("Most forward CG =", CG_mostfor)
    #print("Most afterward CG =", CG_mostaft)
    
   
# =============================================================================
plt.xlim([-0.5,1])

# =============================================================================
plt.plot(CGmacmin, xlemac)
plt.plot(CGmacmax, xlemac)

plt.show()
