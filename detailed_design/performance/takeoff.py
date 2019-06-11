import parameters as p
import numpy as np

# Runway Input Parameters
g = p.g
rho = p.rho0
mu = 0.02

# Aircraft Input Parameters
S = p.S
CL = p.C_L_max_TO
CD0 = p.Cd0
CD = CD0 + CL**2 / (np.pi * p.A * p.e)
T = p.T_TO
W = p.MTOW

# Take-Off Input Parameters
V_stall = p.V_stall
gamma = 4 * (np.pi/180)
n_rotation = 1.15
h_screen = 15.2

V_min = V_stall
V_LOF = 1.05 * V_min
V_bar = V_LOF / (np.sqrt(2))

# Calculate Forces
D = CD * 0.5 * rho * V_bar**2 * S
L = CL * 0.5 * rho * V_bar**2 * S

# Calculate Acceleration and Rotation
a = (g/W) * (T - D - mu * (W - L))
R = (V_LOF**2) / ((n_rotation-1)*g)

# Calculate Distances
s_ground = (V_LOF**2) / (2*a)
s_rotation = R * np.sin(gamma)
s_climb = (h_screen - R*(1-np.cos(gamma))) / np.tan(gamma)
s_TO = (s_ground + s_rotation + s_climb) * 1.15

print("Take-Off Length Required:", np.round(s_TO, decimals=2), "m")