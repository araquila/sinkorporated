import parameters as p
import numpy as np
import matplotlib.pyplot as plt

from atmosphere import atmosphere_calc

# Atmosphere Input Parameters
t0 = p.temperature0
t_gradient = p.temperature_gradient
atR = p.R
atgamma = p.gamma
g = p.g

# Aircraft Input Parameters
S = p.S
CD0 = p.Cd0
T = p.T_TO
P = p.P_TO
W = p.MTOW
CL_max = 1.8

h_range = np.linspace(0, 15000, 3)

for altitude in h_range:
    temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, t0, t_gradient, g, atR, atgamma)
    rho = 1.225*rho
    V_stall = np.sqrt((W/S)*(2/rho)*(1/CL_max))
    V = np.linspace(V_stall, 300, 1000)
    CL = W / (0.5 * rho * V**2 * S)
    CD_L = CL**2 / (np.pi * p.A * p.e)
    CD = CD0 + CD_L
    D = CD * 0.5 * rho * V**2 * S
    Pr = D * V
    Pa = np.ones(len(V)) * P
    plt.plot(V, Pa)
    plt.plot(V, Pr)
  
plt.show()