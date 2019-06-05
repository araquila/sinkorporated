from scipy.integrate import quad
import numpy as np

# Required turn angle in 1.5 seconds
angle = 30/180*np.pi

# Maximum aileron deflection
aileron_max = 20

MTOW = 173185.74
rho = 0.5258
S = 50
# Speed (cruise or landing? The latter is crucial)
V = 49    # np.sqrt(MTOW/(0.5*2.1*rho * S ))
print(V)
# Geometry
b = 31.3
root_chord = 2.24
taper = 0.4
tip_chord = root_chord * taper

# Aerodynamics
Cd0 = 0.008
Cl_alpha = 2 * np.pi
tau = 0.57

# Functions according to ADSEE
def chord(y):
    return root_chord - (root_chord-tip_chord)/(b/2) * y
def fun1(y):
    return chord(y) * y
def fun2(y):
    return y * y * chord(y)

c1, err = quad(fun1,0.8*b/2,0.95*b/2)
Cl_delta_aileron = 2 * Cl_alpha * tau / (S * b) * c1

c2, err = quad(fun2,0,b/2)
Clp = - 4 * (Cl_alpha + Cd0) / (S * b) * c2

P = - Cl_delta_aileron / Clp * aileron_max * (2 * V) / b

# Required time to turn; this should be smaller than 1.5
t = angle / P

print(t)
