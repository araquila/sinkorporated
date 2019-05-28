import numpy as np
import matplotlib.pyplot as plt

# General Aircraft Parameters
x_ac = 10                       # Position of the aerodynamic center [m]
l_h = 11                        # Tail arm of the horizontal tail [-]
c = 2                           # Chord length of the main wing [m]
Vh_V2 = 1                       # Correction factor of velocity over the tail, should be 1 for T-tail [-]

xLEMAC = 10                     # Position of the leading edge of the MAC [m]
MAC = 1.8                       # Length of the MAC [m]

C_L_alpha_h = 4.07              # Lift curve of the horizontal tail [rad^-1]
C_L_alpha_Aminh = 5.35          # Lift curve of the entire aircraft without horizontal tail [rad^-1]
C_L_h = -0.56                   # Lift coefficient of the horizontal tail [-]
C_L_Aminh = 1.93                # Lift coefficient of the entire aircraft without horizontal tail [-]
C_m_ac = -0.5                   # Moment coefficient around the aerodyamic center [-]
downwash = 0.022                # Downwash effect on the horizontal tail [-]

SM = 0.05                       # Stability margin as a percentage of MAC [-]

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
plt.show()
