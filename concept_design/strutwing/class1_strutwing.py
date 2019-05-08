import numpy as np
import fuel_fraction as ff

# INPUT VARIABLES
W_payload = 60000. # [N]
A = 19.5 # [-]
e = 0.8 # [-]
C_fe = 0.0030 # [-]
S_wet = 5.5 # [m^2]
S = 1. # [m^2]
C_D_0 = (S_wet/S) * C_fe
LD = np.sqrt((np.pi*A*e)/(4*C_D_0))
print(LD)

# BASED ON STATISTICS: CALCULATION OF THE OEW AND MTOW_jet
payload_ratio = 0.22 # = PL/MTOW
emptyweight_ratio = 0.636 # = OEW/MTOW
MTOW_tbp = W_payload / payload_ratio
OEW_tbp = MTOW_tbp * emptyweight_ratio

# CALCULATE FUEL WEIGHT
ff_tbp = ff.fuel_fraction(LD_cruise_tbp = LD, tbp = True)
W_fuel = (1-ff_tbp)*MTOW_tbp

print(MTOW_tbp)
print(W_payload + W_fuel + OEW_tbp)
