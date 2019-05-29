import numpy as np
import matplotlib.pyplot as plt

M_cruise = 0.6                  # Cruise Mach number [-]
D_fuselage = 2.8                # Diameter of the fuselage [m]
b = 31                          # Span of the aircraft [m]
S = 50                          # Wing surface area [m^2]
c_root = 1.8                    # Wing root chord [m]

# Horizontal Tail Parameters
A_h = 4                         # Aspect ratio of the horizontal tail [-]
sweep_h = 15                    # Sweep of the horizontal tail [deg]

# Main Wing Parameters
A_w = 20                        # Aspect ratio of the main wing [-]
sweep_w = 0                     # Sweep of the main wing [deg]

# General Aircraft Parameters
x_ac = 10.2                     # Position of the aerodynamic center [m]
l_h = 11                        # Tail arm of the horizontal tail [-]
c = 2                           # Chord length of the main wing [m]
Vh_V2 = 1                       # Correction factor of velocity over the tail, should be 1 for T-tail [-]

xLEMAC = 10                     # Position of the leading edge of the MAC [m]
MAC = 1.8                       # Length of the MAC [m]

C_L_h = -0.35*(A_h)**(1/3)      # Lift coefficient of the horizontal tail, for a fixed tail [-]
C_L_Aminh = 1.93                # Lift coefficient of the entire aircraft without horizontal tail [-]
C_m_ac = -0.5                   # Moment coefficient around the aerodyamic center [-]
downwash = 4/(A_w+2)            # Downwash effect on the horizontal tail [-]

SM = 0.05                       # Stability margin as a percentage of MAC [-]


# DATCOM Estimation Methods
beta = np.sqrt(1-M_cruise**2)
C_L_alpha_h = (2*np.pi*A_h) / (2 + np.sqrt(4 + ((A_h*beta)/0.95)**2 * (1 + (np.tan(np.radians(sweep_h))/beta**2)))) 
C_L_alpha_w = (2*np.pi*A_w) / (2 + np.sqrt(4 + ((A_w*beta)/0.95)**2 * (1 + (np.tan(np.radians(sweep_w))/beta**2))))
C_L_alpha_Aminh = C_L_alpha_w * (1 + 2.15 * (D_fuselage / b)) * ((S - c_root * D_fuselage) / S) + (np.pi / 2) * (D_fuselage**2 / S)

# Make array of tail surfaces
Sh_S = np.linspace(0, 0.6, 100)

# Stability Equations
x_cg_neutral_stability = x_ac + (C_L_alpha_h / C_L_alpha_Aminh) * (1 - downwash) * (l_h / c) * Sh_S * Vh_V2
x_cg_neutral_stability = (x_cg_neutral_stability - xLEMAC) / MAC
x_cg_stability = x_cg_neutral_stability - SM

# Controllability Equation
x_cg_controllability = x_ac - (C_m_ac / C_L_Aminh) + (C_L_h / C_L_Aminh) * (l_h / c) * Sh_S * Vh_V2 
x_cg_controllability = (x_cg_controllability - xLEMAC) / MAC

# Generate Scissor Plot
plt.plot(x_cg_neutral_stability, Sh_S)
plt.plot(x_cg_stability, Sh_S)
plt.plot(x_cg_controllability, Sh_S)
plt.xlim([-0.5,1])
plt.ylim([0,0.5])
plt.show()
