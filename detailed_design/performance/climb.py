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

alt = np.linspace(0, 12825, 100)

Vmin = []
Vmax = []
machmin = []
machmax = []
ROClist = []
V_ROC_max = []
M_ROC_max = []
real_ROC = []

for altitude in alt:
# Aircraft Input Parameters
    S = p.S
    A = p.A
    e = p.e
    CD0 = p.Cd0
    P = p.P_TO
    W = p.MTOW

    # Calculate density at certain height
    temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, t0, t_gradient, g, atR, atgamma)
    rho = 1.225 * rho
    pressure = 101325 * pressure

    P = 1.25*P*(rho/1.225)
      
    # Calculate stall speed
    CLmax = 1.43
    V_stall = np.sqrt((W/S)*(2/rho)*(1/CLmax))

    # Determine Power Required and Power Available
    V = np.linspace(10, 0.6*speed_of_sound, 1001)
    
    for i in range(len(V)):
        if V[i] > 0.3*speed_of_sound:
            test = i
            break
        
    P = np.hstack([np.linspace(0, P, test),np.ones(1001-test)*P])
    Mach = np.linspace(10/speed_of_sound, 0.6, 1001)
    k1 = (1 / (np.pi * A * e))
    CL = W / (0.5*rho*V**2*S)
    CD = CD0 + k1 * CL**2
    Pr = W * np.sqrt((W/S)*(2/rho)*(CD**2/CL**3))
    Pa = np.ones(len(V))*P    
    
    climbpossible = []
    
    for i in range(len(Pr)):
        if Pr[i] < Pa[i]:
            climbpossible.append(V[i])
    
    if climbpossible[0] > V_stall:
        Vmin.append(climbpossible[0])
        machmin.append(climbpossible[0]/speed_of_sound)
    else:
        Vmin.append(V_stall)
        machmin.append(V_stall/speed_of_sound)
    Vmax.append(climbpossible[-1])    
    machmax.append(climbpossible[-1]/speed_of_sound)
    
    ROC = (Pa-Pr)/W
    ROC_max = np.max(ROC)
    idxVROCmax = list(ROC).index(ROC_max)
    ROClist.append(ROC_max)
    V_ROC_max.append(V[idxVROCmax])
    M_ROC_max.append(V[idxVROCmax]/speed_of_sound)
    
    real_ROC.append(ROClist[-1]/(1+0.567*(M_ROC_max[-1]**2)))
    
#plt.plot(Vmin, alt, label="Minimum Velocity")
#plt.plot(Vmax, alt, label="Maximum Velocity")
#plt.xlabel("Velocity [m/s]")
#plt.ylabel("Altitude [m]")
#plt.legend()
#plt.show()

print("ROC at sea level:", np.round(ROClist[0], decimals=3), "m/s")

for i in range(len(real_ROC)):
    if real_ROC[i] > 9:
        real_ROC[i] = 9

plt.plot(alt, real_ROC)
rrr = 0
zzz = 0
for i in range(len(alt)):
    if alt[i] > 8000:
        lim = i
        break
    
for i in range(lim):
    rrr = rrr+(alt[i+1]-alt[i])/real_ROC[i]
    zzz = zzz + V_ROC_max[i]
    
print("Time needed to climb:", np.round(rrr/60, decimals=2), "min")
print("Horizontal distance covered:", np.round((rrr*(zzz/(len(alt)-lim)))/1000, decimals=2), "km")