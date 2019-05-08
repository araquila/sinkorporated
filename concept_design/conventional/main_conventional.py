# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from fuel_fraction import fuel_fraction
from class1_conventional import Weights

# Gravitional constant
g = 9.8065

# Initial mass and fractions
M_payload = 10
M_crew = 10
f_trapped_fuel = 0.003      # Range 0.001-0.005
M_empty_tbp = 10000
M_empty_jet = 15000

# Convert to weights
W_payload = M_payload * g
W_crew = M_crew * g
W_empty_tbp = M_empty_tbp
W_empty_jet = M_empty_jet

## Initial jet and tbp aircraft parameters
C_fe = 1
S = 1
S_wet = 5 * S

# Jet
A_jet = 1
e_jet = 1
cj_loiter_jet = 0.5         # (0.4-0.6) [lbs/lbs/hr]
cj_cruise_jet = 0.6         # (0.5-0.9) [lbs/lbs/hr]

# Tbp
A_tbp = 1
e_tbp = 1
eff_cruise_tbp = 0.85       # [-]
eff_loiter_tbp = 0.77       # [-]
cp_cruise_tbp = 0.5         # (0.4-0.6) [lbs/hp/hr]
cp_loiter_tbp = 0.6         # (0.5-0.7) [lbs/hp/hr]

for iter in range(1):
    MTOW_jet, OEW_jet, W_fuel_jet =  Weights(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, jet = True)
    MTOW_tbp, OEW_tbp, W_fuel_tbp = Weights(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, tbp = True)
