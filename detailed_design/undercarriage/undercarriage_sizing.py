import parameters as p
import numpy as np

# Load on Nose Landing Gear
load_nl = 0.08
load_ml = 1- load_nl

x_cg = 0.71 * p.MAC + 0.4283 * p.l_fuselage

# Position of the Nose Landing Gear
l_nl = 2        #m

l_ml = ((load_nl * (x_cg-l_nl)) / load_ml) + x_cg
