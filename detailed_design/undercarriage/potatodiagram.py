import matplotlib.pyplot as plt
## Constants ##
import sys
import os

sys.path.append(os.getcwd())
from detailed_design.parameters import *

#seat_pitch = 0.7874
weight_passenger = 80
cargo_passenger = 23
CG_OEW = 15.13
#OEW = 12162.4
W_fuel = 1576
OEW = 9785
CG = [CG_OEW]
weight = [OEW]
CG_wing = 0.4 * l_fuselage
n_rows = n_passenger / n_seats_abreast
CG_fuel = 14.5
## CARGO ##

CG_cargo = 18.24

CG_cargo_loaded = (CG_OEW * OEW + CG_cargo * M_total_cargo) / (OEW + M_total_cargo)
CG.append(CG_cargo_loaded)
weight.append(OEW + M_total_cargo)

## PASSENGERS ##
#AISLE LOADING#
x_rowaft = [16.74]
x_rowfront = [16.74 - (n_rows - 1) * seat_pitch]

#start loading aft#
for i in range(int(n_rows-1)):
    x_rowaft.append(x_rowaft[-1] - seat_pitch)
    x_rowfront.append(x_rowfront[-1] + seat_pitch)

for i in range(int(n_rows)):
    weight.append(weight[-1] + (n_seats_abreast / 2) * weight_passenger)
    CG.append((weight[-2]*CG[-1]+(n_seats_abreast / 2) * weight_passenger * x_rowaft[i])/weight[-1])

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

print("Most forward CG =", min(CG))
print("Most afterward CG =", max(CG))

plt.scatter(CG,weight)
plt.show()
