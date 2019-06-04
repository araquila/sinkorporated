import parameters as p
import numpy as np

# Insert Scrap Angle Here
scrap_angle = 15

# Nose Landing Gear Parameters
load_nl = 0.08
l_nl = 2        # m

# Centre of Gravity Calculation
x_cg = 0.71 * p.MAC + 0.4283 * p.l_fuselage


#Calculate Main Landing Gear Parameters
load_ml = 1- load_nl
l_ml = ((load_nl * (x_cg-l_nl)) / load_ml) + x_cg
distance_to_tail = p.l_fuselage - p.l_tail - l_ml

np.radians(scrap_angle)
