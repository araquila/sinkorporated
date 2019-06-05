import parameters as p
import numpy as np
from matplotlib import pyplot as plt

radius = p.d_fuselage_outside / 2 

# Choose which material
sigma_yield = p.yield_stress_al2024
density = p.density_aluminum

# Import strut and reaction forces from the wing
alpha = np.arctan((p.strut_pos_perc * p.b/2)/p.d_fuselage_outside)          # Angle of the strut with fuselage


# Constant weight distribution of fuselage weight
W_fuselage_dist = p.W_fuselage / p.l_fuselage

# Required thicknesses due to pressurisation
def t_circ(sigma_yield,p_in,p_out,R):
    t = p.safetyfactor_fuselage*(p_in-p_out)*R/sigma_yield
    return t

def t_long(sigma_yield,p_in,p_out,R):
    t = p.safetyfactor_fuselage*(p_in-p_out)*R/(2*sigma_yield)
    return t

t_circ = t_circ(sigma_yield,p.pressure_inside,p.pressure_outside,radius)
t_long = t_long(sigma_yield,p.pressure_inside,p.pressure_outside,radius)

# Print minimum required thickness
print('Required thickness [mm]: ' + str(max(t_circ * 1000,t_long * 1000)))

# Shear 
x_shear = 0
V_list = []
x_shear_list = []

while x_shear < p.l_fuselage:
    if x_shear < p.xLEMAC:
        V = x_shear * W_fuselage_dist
    else:
        V = x_shear * W_fuselage_dist - np.cos(alpha) * p.F_strut - p.R_y
    V_list.append(V)
    x_shear_list.append(x_shear)
    x_shear += 0.1
    
# Moment
x_moment = 0
M_list = []
x_moment_list = []

while x_moment < p.l_fuselage:
    if x_moment < p.xLEMAC:
        M = x_moment**2 * W_fuselage_dist * 0.5
    else:
        M = x_moment**2 * W_fuselage_dist * 0.5 - (np.cos(alpha) * p.F_strut + p.R_y) * (x_moment - p.xLEMAC)
    M_list.append(M)
    x_moment_list.append(x_moment)
    x_moment += 0.1    

plt.figure(1,figsize = (8,6))
plt.xlabel('Location along span [m]',fontsize=13)
plt.ylabel('Shear force [N]',fontsize=13)
plt.plot(x_shear_list, V_list, 'r')

plt.figure(2,figsize = (8,6))
plt.xlabel('Location along span [m]',fontsize=13)
plt.ylabel('Moment [Nm]',fontsize=13)
plt.plot(x_moment_list, M_list, 'b')
