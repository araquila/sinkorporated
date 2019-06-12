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

alt = [8000, 9000, 13000]

Vmin = []
Vmax = []
CLlist = []
CDlist = []
ROClist = []
for altitude in alt:
    # Calculate density at certain height
    temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, t0, t_gradient, g, atR, atgamma)
    rho = 1.225 * rho
    pressure = 101325 * pressure

    # Calculate stall speed
    CLmax = 2
    V_stall = np.sqrt((W/S)*(2/rho)*(1/CLmax))

    # Determine Power Required and Power Available
    V = np.linspace(V_stall, V_stall+150, 1001)
    k1 = (1 / (np.pi * A * e))
    CL = W / (0.5*rho*V**2*S)
    CD = CD0 + k1 * CL**2
    D = CD * 0.5 * rho * V**2 * S
    Pr = D*V
    Pr = W * np.sqrt((W/S)*(2/rho)*(CD**2/CL**3))
    Pa = np.ones(len(V))*P

    for i in range(len(Pr)):
        if Pr[i] > Pa[i]:
            index_Vmax = i
            break

    for i in range(len(V)):
        if V[i] > 0.6*speed_of_sound:
            test = i
            break
    print(CL[test])
    print(D[test])
    
    Vmin.append(V_stall)
    Vmax.append(V[index_Vmax])    
    CLlist.append(CL)
    CDlist.append(CD)    
    
    ROC = (Pa-Pr)/W
    ROC_max = np.max(ROC)
    ROClist.append(ROC_max)

#plt.plot(ROClist, alt)
#plt.show()
#plt.plot(Vmin, alt)
#plt.plot(Vmax, alt)
#plt.show()
    
#plt.scatter(Vmin, alt)
#plt.scatter(Vmax, alt)



