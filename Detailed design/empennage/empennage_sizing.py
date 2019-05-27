# Aircraft Parameters
l_h = 11            # [m]
l_v = 11            # [m]
S = 50              # [m2]
b = 30              # [m]
c = 1.7             # [m]

# Horizontal Tail Parameters
A_h = 4             # [-]
taper_h = 0.5       # [-]
sweep_h = 15        # [deg]

# Vertical Tail Parameters
A_v = 1.5           # [-]
taper_v = 0.6       # [-]
sweep_v = 30        # [deg]

# Volume should be in the range of 0.04-0.08 according to statistics
V_v = 0.05          # [-]

# Determination of the Horizontal Tail Surface
# Here should be the stability function that calculates the required Sh/S
S_h = 11.61         # [m2]

# Determination of the Vertical Tail Surface
# This is purely based on statistics
S_v = (V_v * S * b) / l_v

