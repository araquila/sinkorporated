import numpy as np
from matplotlib import pyplot as plt
import parameters as p

# Safety factor
safety = 2.5

# Weight and lift
design_load = p.MTOW * safety
L_wing = design_load / 2                # [N]

# Calculate ap.b/2solute distance
strut_pos = p.b/2 * p.strut_pos_perc          # [m]
engine_pos = p.b/2 * p.engine_pos_perc        # [m]

# Calculate angle p.b/2etween strut and wing
theta = np.arctan(p.d_fuselage_outside/strut_pos) # [radians]

# Calculate forces
L_wing_dist_start = L_wing / p.b/2 * (0.54/0.4) 
L_wing_dist_end = L_wing / p.b/2 * (0.33/0.4)

L_gradient = (L_wing_dist_end - L_wing_dist_start)/p.b/2

L_1 = L_wing_dist_end * p.b/2
L_2 = (L_wing_dist_start - L_wing_dist_end) * p.b/2 * 0.5

# Wing weight distrip.b/2ution
p.W_wing_dist = p.W_wing / p.b/2

# Set reaction moment
M = 100000

# Calc strut force and vertical reaction force
F_strut = -(-0.5 * p.b/2 * L_1 - 1/3 * p.b/2 * L_2 + p.W_engine * engine_pos + (p.W_pod + p.W_fuel) * strut_pos + p.W_wing * 0.5 * p.b/2 + M) / (strut_pos * np.sin(theta))

A_strut = (F_strut/p.ult_stress_carbon)*10000
print('Radius carbon strut: ' + str(np.sqrt(A_strut/np.pi)) + '[cm]')

F_y = p.W_wing + p.W_engine + p.W_fuel + p.W_pod + np.sin(theta) * F_strut  - L_1 - L_2

### SHEAR DIAGRAM #### 

x_shear = 0
V_list = [0]  
x_shear_list = [p.b/2]                             
while (p.b/2-x_shear) > 0:
    if (p.b/2 - x_shear) <= engine_pos:
        L_wing_dist_x_shear = L_wing_dist_end - L_gradient * x_shear
        L_x_shear1 = L_wing_dist_end * x_shear
        L_x_shear2 = (L_wing_dist_x_shear - L_wing_dist_end) * x_shear * 0.5
        V_x = p.W_wing_dist * x_shear - L_x_shear1 - L_x_shear2 + p.W_engine + p.W_pod + p.W_fuel + F_strut * np.sin(theta)
        V_list.append(V_x)
    elif (p.b/2 - x_shear) > engine_pos and (p.b/2 - x_shear) <= strut_pos:
        L_wing_dist_x_shear = L_wing_dist_end - L_gradient * x_shear
        L_x_shear1 = L_wing_dist_end * x_shear
        L_x_shear2 = (L_wing_dist_x_shear - L_wing_dist_end) * x_shear * 0.5
        V_x = p.W_wing_dist * x_shear - L_x_shear1 - L_x_shear2 + F_strut * np.sin(theta) + p.W_pod + p.W_fuel
        V_list.append(V_x)
    else:
        L_wing_dist_x_shear = L_wing_dist_end - L_gradient * x_shear
        L_x_shear1 = L_wing_dist_end * x_shear
        L_x_shear2 = (L_wing_dist_x_shear - L_wing_dist_end) * x_shear * 0.5
        V_x = p.W_wing_dist * x_shear - L_x_shear1 - L_x_shear2
        V_list.append(V_x)
    x_shear_list.append(p.b/2 - x_shear)
    x_shear += 0.1

x_shear_list.append(0)
V_list.append(0)

x_moment = 0
M_list = []  
x_moment_list = []                              
while (p.b/2-x_moment) > 0:
    if (p.b/2 - x_moment) < engine_pos:
        L_wing_dist_x_moment = L_wing_dist_end - L_gradient * x_moment
        L_x_moment1 = L_wing_dist_end * x_moment
        L_x_moment2 = (L_wing_dist_x_moment - L_wing_dist_end) * x_moment * 0.5
        M_x = p.W_wing_dist * x_moment * (x_moment/2) - L_x_moment1 * (x_moment/2) - L_x_moment2 * (x_moment*1/3) + p.W_engine * (x_moment - (p.b/2 - engine_pos)) +(p.W_pod + p.W_fuel + F_strut * np.sin(theta)) * (x_moment - (p.b/2 - strut_pos))
        M_list.append(-M_x)
    elif (p.b/2 - x_moment) > engine_pos and (p.b/2 - x_moment) <= strut_pos:
        L_wing_dist_x_moment = L_wing_dist_end - L_gradient * x_moment
        L_x_moment1 = L_wing_dist_end * x_moment
        L_x_moment2 = (L_wing_dist_x_moment - L_wing_dist_end) * x_moment * 0.5
        M_x = p.W_wing_dist * x_moment * (x_moment/2) - L_x_moment1 * (x_moment/2) - L_x_moment2 * (x_moment*1/3) + (F_strut * np.sin(theta) + p.W_pod  + p.W_fuel) * (x_moment - (p.b/2 - strut_pos))
        M_list.append(-M_x)
    else:
        L_wing_dist_x_moment = L_wing_dist_end - L_gradient * x_moment
        L_x_moment1 = L_wing_dist_end * x_moment
        L_x_moment2 = (L_wing_dist_x_moment - L_wing_dist_end) * x_moment * 0.5
        M_x = p.W_wing_dist * x_moment * (x_moment/2) - L_x_moment1 * (x_moment/2) - L_x_moment2 * (x_moment*1/3)
        M_list.append(-M_x)
    x_moment_list.append(p.b/2 - x_moment)
    x_moment += 0.01
    
M_list.append(0)
x_moment_list.append(0)

plt.figure(1,figsize = (8,6))
plt.xlabel('Location along span [m]',fontsize=13)
plt.ylabel('Shear force [N]',fontsize=13)
plt.plot(x_shear_list, V_list, 'r')

plt.figure(2,figsize = (8,6))
plt.xlabel('Location along span [m]',fontsize=13)
plt.ylabel('Moment [Nm]',fontsize=13)
plt.plot(x_moment_list, M_list, 'b')