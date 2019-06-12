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
A = p.A
e = p.e
CD0 = p.Cd0
T = p.T_TO
P = p.P_TO
W = p.MTOW

alt = [0]
for altitude in alt:
    temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, t0, t_gradient, g, atR, atgamma)
    rho = 1.225 * rho
    pressure = 101325 * pressure

    CLmax = 2

    V_stall = np.sqrt((W/S)*(2/rho)*(1/CLmax))

    V = np.linspace(V_stall, V_stall+150, 1001)
    k1 = (1 / (np.pi * A * e))
    CL = W / (0.5*rho*V**2*S)
    CD = CD0 + k1 * CL**2
    D = CD * 0.5 * rho * V**2 * S
    Pr = D*V
    Pr = W * np.sqrt((W/S)*(2/rho)*(CD**2/CL**3))
    Pa = np.ones(len(V))*P

    D0 = CD0 * 0.5 * rho * V**2 * S
    DL = k1 * CL**2 * 0.5 * rho * V**2 * S

    #plt.plot(V, D0)
    #plt.plot(V, DL)
    plt.plot(V, D*V)
    plt.plot(V, Pa)
plt.show()



ROC = (Pa-Pr)/W
ROC_max = np.max(ROC)

#plt.plot(V, Pr/1e6, label="Power Required")
#plt.plot(V, Pa/1e6, label="Power Available")
#plt.xlabel("Velocity [m/s]")
#plt.ylabel("Power [MW]")
#plt.legend()
#plt.show()

print("Maximum ROC:", np.round(ROC_max, decimals=2), "m/s")
