from scipy.integrate import quad
import numpy as np
import parameters as p

# Required turn angle in 1.5 seconds
angle = np.deg2rad(30)

# Maximum aileron deflection
aileron_max = 20

rho = 0.5258
# Speed (cruise or landing? The latter is crucial)
#V = 49    # np.sqrt(p.MTOW/(0.5*2.1*rho * S ))
#print(V)
# Geometry
#b = 31.3
#root_chord = 2.24
#taper = 0.4
#tip_chord = root_chord * taper

# Aerodynamics
#Cd0 = 0.008
#Cl_alpha = 2 * np.pi
#tau = 0.57

# Functions according to ADSEE
def chord(y):
    return p.root_chord - (p.root_chord-p.tip_chord)/(p.b/2) * y
def fun1(y):
    return chord(y) * y
def fun2(y):
    return y * y * chord(y)

c1, err = quad(fun1,0.8*p.b/2,0.95*p.b/2)
Cl_delta_aileron = 2 * p.Cl_alpha * p.tau / (p.S * p.b) * c1

c2, err = quad(fun2,0,p.b/2)
Clp = - 4 * (p.Cl_alpha + p.Cd0) / (p.S * p.b) * c2

P = - Cl_delta_aileron / Clp * aileron_max * (2 * p.V_cruise) / p.b

# Required time to turn; this should be smaller than 1.5
t = angle / P

print(t)
