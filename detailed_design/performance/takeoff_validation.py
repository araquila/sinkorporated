import parameters as p
import numpy as np
import matplotlib.pyplot as plt

from atmosphere import atmosphere_calc

# Runway Input Parameters
g = p.g
rho0 = p.rho0
mu = 0.02

altrange = np.linspace(0, 1500, 100)
dist_TO = []

for altitude in altrange:
    temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, p.temperature0, p.temperature_gradient, p.g, p.R, p.gamma)
    rho = 1.225 * rho
    pressure = 101325 * pressure
    
    
    # Aircraft Input Parameters
    S = 61
    CL = 1.6
    CD0 = 0.02
    CD = CD0 + CL**2 / (np.pi * p.A * p.e)
    P = 3.689e6
    W = 223668
    
    # Take-Off Input Parameters
    V_stall = np.sqrt((W/S)*(2/rho)*(1/CL))
    gamma = 4 * (np.pi/180)
    n_rotation = 1.15
    h_screen = 15.2
    
    V_min = V_stall
    V_LOF = 1.05 * V_min
    V_bar = V_LOF / (np.sqrt(2))
    
    T = (0.9 * P) / V_bar
    
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
    
    dist_TO.append(s_TO)

    #print("Take-Off Length Required:", np.round(s_TO, decimals=2), "m")
    
plt.plot(altrange, dist_TO)
plt.xlim([altrange[0],altrange[-1]])
plt.ylim([dist_TO[0],dist_TO[-1]])
plt.xlabel("Runway height above sea level [m]")
plt.ylabel("Take-Off distance required [m]")

print((dist_TO[0] / 1175)*100)