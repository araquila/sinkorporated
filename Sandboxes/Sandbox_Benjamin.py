import math
import matplotlib.pyplot as plt
import numpy as np
from skaero.atmosphere import coesa
tbp=True
jet=False
h, T, p, rho= coesa.table(0)
if tbp==True:
    v=40
if jet==True:
    v=200


a=math.sqrt(1.4*287*T)
M_cruise=v/a
delta_mach=0.03
C_L_cruise=0.3
half_chord_sweep=3*math.pi/180
t_c_ratio = min(0.18, (np.cos(half_chord_sweep)**3 * (1 - (M_cruise + delta_mach) * np.cos(half_chord_sweep)) - 0.115 * C_L_cruise**1.5) / np.cos(half_chord_sweep)**2)
print(t_c_ratio)
AR=11
S=55
b = np.sqrt(AR * S)
sweep=0
taper = 0.2 * (2 - sweep)
root_chord = (2 * S)/((1 + taper) * b)


MAC=root_chord * (2 / 3) * ((1 + taper + taper**2) / (1 + taper))
reynolds=(rho*v)/0.0000181206
print(reynolds)
print(MAC)
A=11
M_infinity=0.8
sweep_chord_0_5=10*math.pi/180
e_airfoil=0.95
beta_aero=math.sqrt(1-M_infinity**2)
C_L_slope = (2*math.pi*A)/(2+math.sqrt(4+(A*beta_aero/e_airfoil)**2*(1+(math.tan(sweep_chord_0_5)**2/beta_aero**2))))
