import numpy as np
import fuel_fraction as ff

# INPUT VARIABLES
W_payload = 6120. # [kg]
A = 19.5 # [-]
e = 0.8 # [-]
C_fe = 0.0030 # [-]
S_wet = 5.5 # [m^2]
S = 1. # [m^2]
C_D_0 = (S_wet/S) * C_fe
LD_cruise_tbp = np.sqrt((np.pi*A*e)/(4*C_D_0))
LD_cruise_jet = (3/4) * np.sqrt((np.pi*A*e)/(3*C_D_0))


# BASED ON STATISTICS: CALCULATION OF THE OEW AND MTOW_jet
payload_ratio_tbp = 0.22 # = PL/MTOW
emptyweight_ratio_tbp = 0.615 # = OEW/MTOW
payload_ratio_jet = 0.22 # = PL/MTOW
emptyweight_ratio_jet = 0.570 # = OEW/MTOW
MTOW_tbp = W_payload / payload_ratio_tbp
OEW_tbp = MTOW_tbp * emptyweight_ratio_tbp
MTOW_jet = W_payload / payload_ratio_jet
OEW_jet = W_payload * emptyweight_ratio_jet


# CALCULATE FUEL WEIGHT
ff_tbp = ff.fuel_fraction(LD_cruise_tbp = LD_cruise_tbp, tbp = True)
ff_jet = ff.fuel_fraction(LD_cruise_jet = LD_cruise_jet, jet = True)
W_fuel_tbp = (1-ff_tbp)*MTOW_tbp
W_fuel_jet = (1-ff_jet)*MTOW_jet

print(MTOW_tbp)
print(W_payload + W_fuel_tbp + OEW_tbp)
