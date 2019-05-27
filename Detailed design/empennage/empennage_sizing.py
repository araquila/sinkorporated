import numpy as np

# Aircraft Parameters
l_h = 11            # [m]
l_v = 11            # [m]
S = 50              # [m2]
b = 30              # [m]
c = 1.7             # [m]

# Horizontal Tail Parameters
A_h = 5             # [-]
taper_h = 0.4       # [-]
sweep_h = 15        # [deg]
V_h = 1             # [-]

# Vertical Tail Parameters
A_v = 1.5           # [-]
taper_v = 0.5       # [-]
sweep_v = 30        # [deg]
V_v = 0.05          # [-]

# Determination of the Horizontal Tail Surface
S_h = (V_h * S * c) / l_h

# Determination of the Vertical Tail Surface
S_v = (V_v * S * b) / l_v

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