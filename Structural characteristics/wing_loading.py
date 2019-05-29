import numpy as np
from matplotlib import pyplot as plt

# Constants
g = 9.80665
ult_stress_carbon = 600e6
safety = 2.5

strutwing = True
conventional = False

if strutwing:
    # Weight and lift
    W_wing = 1288 * g
    W_fuel = 7736.30
    MTOW = 173185.74 * safety
    W_cruise_start = MTOW                   # W_cruise_start calculated using fuel fraction
    L_wing = W_cruise_start / 2             # [N]
    W_pod = 157.55 * g                      # [N]
    W_engine = 450 * g                      # [N]
    
    # Aircraft parameters
    b = 29.76/2                             # [m]            
    strut_pos_perc = 0.5                    # % of span
    engine_pos_perc = 0.27                  # % of span
    d_fuselage =  2.84                      # [m]
    
    # Calculate absolute distance
    strut_pos = b * strut_pos_perc          # [m]
    engine_pos = b * engine_pos_perc        # [m]
    
    # Calculate angle between strut and wing
    theta = np.arctan(d_fuselage/strut_pos) # [radians]
    
    # Calculate forces
    L_wing_dist_start = L_wing / b * (0.54/0.4) 
    L_wing_dist_end = L_wing / b * (0.33/0.4)
    
    L_gradient = (L_wing_dist_end - L_wing_dist_start)/b
    
    L_1 = L_wing_dist_end * b
    L_2 = (L_wing_dist_start - L_wing_dist_end) * b * 0.5
    
    # Wing weight distribution
    W_wing_dist = W_wing / b
    
    # Set reaction moment
    M = 100000
    
    # Calc strut force and vertical reaction force
    F_strut = -(-0.5 * b * L_1 - 1/3 * b * L_2 + W_engine * engine_pos + (W_pod + W_fuel) * strut_pos + W_wing * 0.5 * b + M) / (strut_pos * np.sin(theta))
    
    A_strut = (F_strut/ult_stress_carbon)*10000
    print('Radius carbon strut: ' + str(np.sqrt(A_strut/np.pi)) + '[cm]')
    
    F_y = W_wing + W_engine + W_fuel + W_pod + np.sin(theta) * F_strut  - L_1 - L_2
    
    ### SHEAR DIAGRAM #### 
    
    x_shear = 0
    V_list = [0]  
    x_shear_list = [b]                             
    while (b-x_shear) > 0:
        if (b - x_shear) <= engine_pos:
            L_wing_dist_x_shear = L_wing_dist_end - L_gradient * x_shear
            L_x_shear1 = L_wing_dist_end * x_shear
            L_x_shear2 = (L_wing_dist_x_shear - L_wing_dist_end) * x_shear * 0.5
            V_x = W_wing_dist * x_shear - L_x_shear1 - L_x_shear2 + W_engine + W_pod + W_fuel + F_strut * np.sin(theta)
            V_list.append(V_x)
        elif (b - x_shear) > engine_pos and (b - x_shear) <= strut_pos:
            L_wing_dist_x_shear = L_wing_dist_end - L_gradient * x_shear
            L_x_shear1 = L_wing_dist_end * x_shear
            L_x_shear2 = (L_wing_dist_x_shear - L_wing_dist_end) * x_shear * 0.5
            V_x = W_wing_dist * x_shear - L_x_shear1 - L_x_shear2 + F_strut * np.sin(theta) + W_pod + W_fuel
            V_list.append(V_x)
        else:
            L_wing_dist_x_shear = L_wing_dist_end - L_gradient * x_shear
            L_x_shear1 = L_wing_dist_end * x_shear
            L_x_shear2 = (L_wing_dist_x_shear - L_wing_dist_end) * x_shear * 0.5
            V_x = W_wing_dist * x_shear - L_x_shear1 - L_x_shear2
            V_list.append(V_x)
        x_shear_list.append(b - x_shear)
        x_shear += 0.1
    
    x_shear_list.append(0)
    V_list.append(0)
    
    x_moment = 0
    M_list = []  
    x_moment_list = []                              
    while (b-x_moment) > 0:
        if (b - x_moment) < engine_pos:
            L_wing_dist_x_moment = L_wing_dist_end - L_gradient * x_moment
            L_x_moment1 = L_wing_dist_end * x_moment
            L_x_moment2 = (L_wing_dist_x_moment - L_wing_dist_end) * x_moment * 0.5
            M_x = W_wing_dist * x_moment * (x_moment/2) - L_x_moment1 * (x_moment/2) - L_x_moment2 * (x_moment*1/3) + W_engine * (x_moment - (b - engine_pos)) +(W_pod + W_fuel + F_strut * np.sin(theta)) * (x_moment - (b - strut_pos))
            M_list.append(-M_x)
        elif (b - x_moment) > engine_pos and (b - x_moment) <= strut_pos:
            L_wing_dist_x_moment = L_wing_dist_end - L_gradient * x_moment
            L_x_moment1 = L_wing_dist_end * x_moment
            L_x_moment2 = (L_wing_dist_x_moment - L_wing_dist_end) * x_moment * 0.5
            M_x = W_wing_dist * x_moment * (x_moment/2) - L_x_moment1 * (x_moment/2) - L_x_moment2 * (x_moment*1/3) + (F_strut * np.sin(theta) + W_pod  + W_fuel) * (x_moment - (b - strut_pos))
            M_list.append(-M_x)
        else:
            L_wing_dist_x_moment = L_wing_dist_end - L_gradient * x_moment
            L_x_moment1 = L_wing_dist_end * x_moment
            L_x_moment2 = (L_wing_dist_x_moment - L_wing_dist_end) * x_moment * 0.5
            M_x = W_wing_dist * x_moment * (x_moment/2) - L_x_moment1 * (x_moment/2) - L_x_moment2 * (x_moment*1/3)
            M_list.append(-M_x)
        x_moment_list.append(b - x_moment)
        x_moment += 0.01
        
    M_list.append(0)
    x_moment_list.append(0)
        
if conventional:
        # Weight and lift
    W_wing = 1708.14 * g
    W_fuel = 19670
    MTOW = 182209 * safety
    W_cruise_start = MTOW                   # W_cruise_start calculated using fuel fraction
    L_wing = W_cruise_start / 2             # [N]
    W_engine = 544 * g                      # [N]
    
    # Aircraft parameters
    b = 24/2                                # [m]            
    engine_pos_perc = 0.33                  # % of span
    
    # Calculate absolute distance
    engine_pos = b * engine_pos_perc        # [m]
    
    # Calculate forces
    L_wing_dist_start = L_wing / b * (0.54/0.4) 
    L_wing_dist_end = L_wing / b * (0.33/0.4)
    
    L_gradient = (L_wing_dist_end - L_wing_dist_start)/b
    
    L_1 = L_wing_dist_end * b
    L_2 = (L_wing_dist_start - L_wing_dist_end) * b * 0.5
    
    # Wing weight distribution
    W_wing_dist = (W_wing + W_fuel) / b
    
    # Calc reaction moment and vertical reaction force
    M_reaction = -0.5 * b * L_1 - 1/3 * b * L_2 + W_engine * engine_pos + (W_wing + W_fuel) * 0.5 * b
    
    F_y = W_wing + W_engine + W_fuel - L_1 - L_2
    
    ### SHEAR DIAGRAM #### 
    
    x_shear = 0
    V_list = [0]  
    x_shear_list = [b]                             
    while (b-x_shear) > 0:
        if (b - x_shear) <= engine_pos:
            L_wing_dist_x_shear = L_wing_dist_end - L_gradient * x_shear
            L_x_shear1 = L_wing_dist_end * x_shear
            L_x_shear2 = (L_wing_dist_x_shear - L_wing_dist_end) * x_shear * 0.5
            V_x = W_wing_dist * x_shear - L_x_shear1 - L_x_shear2 + W_engine
            V_list.append(V_x)
        else:
            L_wing_dist_x_shear = L_wing_dist_end - L_gradient * x_shear
            L_x_shear1 = L_wing_dist_end * x_shear
            L_x_shear2 = (L_wing_dist_x_shear - L_wing_dist_end) * x_shear * 0.5
            V_x = W_wing_dist * x_shear - L_x_shear1 - L_x_shear2
            V_list.append(V_x)
        x_shear_list.append(b - x_shear)
        x_shear += 0.1
    
    x_shear_list.append(0)
    V_list.append(0)
    
    x_moment = 0
    M_list = []  
    x_moment_list = []                              
    while (b-x_moment) > 0:
        if (b - x_moment) < engine_pos:
            L_wing_dist_x_moment = L_wing_dist_end - L_gradient * x_moment
            L_x_moment1 = L_wing_dist_end * x_moment
            L_x_moment2 = (L_wing_dist_x_moment - L_wing_dist_end) * x_moment * 0.5
            M_x = W_wing_dist * x_moment * (x_moment/2) - L_x_moment1 * (x_moment/2) - L_x_moment2 * (x_moment*1/3) + W_engine * (x_moment - (b - engine_pos))
            M_list.append(-M_x)
        else:
            L_wing_dist_x_moment = L_wing_dist_end - L_gradient * x_moment
            L_x_moment1 = L_wing_dist_end * x_moment
            L_x_moment2 = (L_wing_dist_x_moment - L_wing_dist_end) * x_moment * 0.5
            M_x = W_wing_dist * x_moment * (x_moment/2) - L_x_moment1 * (x_moment/2) - L_x_moment2 * (x_moment*1/3)
            M_list.append(-M_x)
        x_moment_list.append(b - x_moment)
        x_moment += 0.1
    
    x_moment_list.append(0)
    M_list.append(0)
    
plt.figure(1,figsize = (8,6))
plt.xlabel('Location along span [m]',fontsize=13)
plt.ylabel('Shear force [N]',fontsize=13)
plt.plot(x_shear_list, V_list, 'r')

plt.figure(2,figsize = (8,6))
plt.xlabel('Location along span [m]',fontsize=13)
plt.ylabel('Moment [Nm]',fontsize=13)
plt.plot(x_moment_list, M_list, 'b')