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

h_range = np.linspace(0, 10000, 6)

for altitude in h_range:
    temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, t0, t_gradient, g, atR, atgamma)
    V = np.linspace(20, 300, 1000)
    CL = W / (0.5 * rho * V**2 * S)
    CD = CD0 + (CL**2 / (np.pi * p.A * p.e))
    D = CD * 0.5 * rho * V**2 * S
    Pr = D * V
    Pa = np.ones(len(V)) * P
    plt.plot(V, Pa)
    plt.plot(V, Pr)
  
plt.show()