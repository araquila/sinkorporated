import numpy as np
import parameters as p

# Aircraft Parameters
l_h = p.l_h            # [m]
l_v = p.l_v            # [m]
S = p.S              # [m2]
b = p.b              # [m]
c = p.MAC             # [m]
T = 31
y_engine = 4

# Horizontal Tail Parameters
A_h = p.A_h             # [-]
taper_h = p.taper_h       # [-]
sweep_h = p.sweep_h        # [deg]
V_h = p.V_h             # [-]

# Vertical Tail Parameters
A_v = p.A_v           # [-]
taper_v = p.taper_v       # [-]
sweep_v = p.sweep_v        # [deg]
V_v = p.V_v    
V_v2 = 0.379      # [-]

# Determination of the Horizontal Tail Surface
#S_h = (V_h * S * c) / l_h
S_h = 0.2745 * S

# Determination of the Vertical Tail Surface
#S_v = (V_v * S * b) / l_v
S_v = 0.21 * S
#S_v2 = (V_v2 * y_engine * T) / l_v


# Planform of the Vertical Tail
b_v = np.sqrt(S_v * A_v)
c_root_v = (2 * (b_v/A_v)) / (1 + taper_v)
c_tip_v = c_root_v * taper_v
MAC_v = (2/3) * c_root_v * ((1 + taper_v + taper_v**2)/(1 + taper_v))

# Planform of the Horizontal Tail
b_h = np.sqrt(S_h * A_h)
c_root_h = (2 * (b_h/A_h)) / (1 + taper_h)
c_tip_h = c_root_h * taper_h
MAC_h = (2/3) * c_root_h * ((1 + taper_h + taper_h**2)/(1 + taper_h))

print("HORIZONTAL TAIL PARAMETERS")
print("Root Chord")