import parameters as p
import numpy as np
from matplotlib import pyplot as plt
from wing.main_strutbox import Fry, Fs, Mrz

radius = p.d_fuselage_outside / 2 

# Choose which material
sigma_yield = p.yield_al_2060_T8E30
#density = p.density_aluminum

## Import strut and reaction forces from the wing
#lengthdata = 100
#Lift, Chord, Yle, Drag, AeroMoment = wd2.read_aero_data("wing/datastrut4.txt", lengthdata, p.V_cruise, p.rho)
#Frx, Fry, Fs, Mrz, Frz, Fsz, Mry, momentyi, momentzi, shearyi, shearzi, vyi, vny, vzi, vnz, xi, theta = wd2.CallForces(Lift, Yle, Drag, p.tot_thrust, Iyy_list, Izz_list,70*10**9, p.engine_pos_perc, p.strut_pos_perc, p.pod_pos_perc)
#
alpha = np.arctan((p.strut_pos_perc * p.b/2)/p.d_fuselage_outside)        # Angle of the strut with fuselage
F_strut_x = Fs * np.sin(alpha)

# Constant weight distribution of fuselage weight
W_fuselage_dist = (p.MTOW - 2*p.W_wing - 2*p.W_pod - 2*p.W_engine - 10000) / p.l_fuselage

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
        V = x_shear * W_fuselage_dist - np.cos(alpha) * Fs - Fry
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
        M = x_moment**2 * W_fuselage_dist * 0.5 - (np.cos(alpha) * Fs + Fry) * (x_moment - p.xLEMAC)
    M_list.append(M)
    x_moment_list.append(x_moment)
    x_moment += 0.1    

plt.figure(10,figsize = (8,6))
plt.xlabel('Location along the length of the fuselage [m]',fontsize=13)
plt.ylabel('Shear force in the y-direction [N]',fontsize=13)
plt.plot(x_shear_list, V_list, 'r')

plt.figure(11,figsize = (8,6))
plt.xlabel('Location along the length of the fuselage [m]',fontsize=13)
plt.ylabel('Moment around the x-axis [Nm]',fontsize=13)
plt.plot(x_moment_list, M_list, 'b')
