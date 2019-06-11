import parameters as p
import numpy as np


g = p.g
rho = p.rho0
S = p.S
mu = 0.02

gamma = 4 * (np.pi/180)

CL = p.C_L_max_TO
CD = p.Cd0 + CL**2 / (np.pi * p.A * p.e)

V_stall = p.V_stall
V_LOF = 1.05 * V_stall
V_bar = V_LOF / (np.sqrt(2))

D = CD * 0.5 * rho * V_bar**2 * S
L = CL * 0.5 * rho * V_bar**2 * S
T = 40000
W = p.MTOW

a = (g/W) * (T - D - mu * (W - L))

s_ground = (V_LOF**2) / (2*a)
s_rotation = (V_LOF**2) / (0.15*g) * np.sin(gamma)
s_climb = (15.2 - (1 - np.cos(gamma)) * (V_LOF**2/(0.15*g))) / np.tan(gamma)

s_TO = (s_ground + s_rotation + s_climb) * 1.15
