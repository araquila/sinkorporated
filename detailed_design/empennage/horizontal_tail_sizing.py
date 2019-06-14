import numpy as np
import matplotlib.pyplot as plt

import parameters as p
import undercarriage.potatodiagram as pt
import undercarriage.CGlocation as cg

M_cruise = p.M_cruise               # Cruise Mach number [-]
D_fuselage = p.d_fuselage_outside   # Diameter of the fuselage [m]
b = p.b                             # Span of the aircraft [m]
S = p.S                             # Wing surface area [m^2]
c_root = p.root_chord               # Wing root chord [m]

# Horizontal Tail Parameters
A_h = p.A_h                         # Aspect ratio of the horizontal tail [-]
sweep_h = p.sweep_h * (np.pi/180)   # Sweep of the horizontal tail [deg]

# Main Wing Parameters
A_w = p.A                           # Aspect ratio of the main wing [-]
sweep_w = - (1-p.taper)/(2*b)       # Sweep of the main wing [deg]

# General Aircraft Parameters
x_ac_w = p.x_ac_w                   # Position of the aerodynamic center [m]
l_h = p.l_h                         # Tail arm of the horizontal tail [-]
c = p.MAC                           # MAC length of the main wing [m]
Vh_V2 = 1                           # Correction factor of velocity over the tail, should be 1 for T-tail [-]
xLEMAC = p.xLEMAC                   # Position of the leading edge of the MAC [m]


C_L_h = -0.35*(A_h)**(1/3)          # Lift coefficient of the horizontal tail, for a fixed tail [-]
C_L_Aminh = 1.4 - C_L_h             # Lift coefficient of the entire aircraft without horizontal tail [-]
C_m_ac = -0.5                       # Moment coefficient around the aerodyamic center [-]
downwash = 4/(A_w+2)                # Downwash effect on the horizontal tail [-]

SM = 0.05                           # Stability margin as a percentage of MAC [-]

# DATCOM Estimation Methods
beta = np.sqrt(1-M_cruise**2)
C_L_alpha_h = (2*np.pi*A_h) / (2 + np.sqrt(4 + ((A_h*beta)/0.95)**2 * (1 + (np.tan(sweep_h)/beta**2)))) 
C_L_alpha_w = (2*np.pi*A_w) / (2 + np.sqrt(4 + ((A_w*beta)/0.95)**2 * (1 + (np.tan(sweep_w)/beta**2))))
C_L_alpha_Aminh = C_L_alpha_w * (1 + 2.15 * (D_fuselage / b)) * ((S - c_root * D_fuselage) / S) + (np.pi / 2) * (D_fuselage**2 / S)
x_ac = x_ac_w - (1.8/C_L_alpha_Aminh) * ((D_fuselage**2 * xLEMAC) / (S*c))

# Make array of tail surfaces
Sh_S = np.linspace(0, 0.6, 100)

# Stability Equations
x_cg_neutral_stability = x_ac + (C_L_alpha_h / C_L_alpha_Aminh) * (1 - downwash) * (l_h / c) * Sh_S * Vh_V2
x_cg_neutral_stability = (x_cg_neutral_stability - xLEMAC) / c
x_cg_stability = x_cg_neutral_stability - SM

# Controllability Equation
x_cg_controllability = x_ac - (C_m_ac / C_L_Aminh) + (C_L_h / C_L_Aminh) * (l_h / c) * Sh_S * Vh_V2 
x_cg_controllability = (x_cg_controllability - xLEMAC) / c

# Generate Scissor Plot
#plt.plot(x_cg_neutral_stability, Sh_S)
plt.plot(x_cg_stability, Sh_S)
plt.plot(x_cg_controllability, Sh_S)
plt.xlim([-0.5,1])
plt.ylim([0,0.6])
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(x_cg_stability, Sh_S, color='b')
ax1.plot(x_cg_controllability, Sh_S, color='b')
ax1.set_xlabel('x_cg/MAC')
ax1.set_xlim([-0.5, 1.5])
ax1.set_ylim([0, 0.6])
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('S_h/S', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(cg.CGmacmin, cg.xlemac, color='r')
ax2.plot(cg.CGmacmax, cg.xlemac, color='r')
ax2.set_ylabel('x_lemac/l_fus', color='r')
ax2.tick_params('y', colors='r')
ax2.set_xlim([-0.5, 1.5])
ax2.set_ylim([0.35, 0.54])

fig.tight_layout()
plt.show()
