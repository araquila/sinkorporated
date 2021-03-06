# CLASS 1 ESTIMATION OF THE STRUTTED AIRCRAFT SIZING PROGRAM

# Import modules
from class1_conventional import Weights

# Gravitional constant
g = 9.8065

# Passengers and crew
n_passenger = 60
M_passenger = 102 #(including luggage)
n_crew = 4
M_crew_member = 100

# Initial mass and fractions
M_payload = n_passenger * M_passenger
M_crew = n_crew * M_crew_member
f_trapped_fuel = 0.003      # Range 0.001-0.005
M_empty_tbp = 13000
M_empty_jet = 15000

# Convert to weights
W_payload = M_payload * g
W_crew = M_crew * g
W_empty_tbp = M_empty_tbp * g
W_empty_jet = M_empty_jet * g

## Initial jet and tbp aircraft parameters
C_fe = 0.003
S = 1
S_wet = 5 * S

# Jet
A_jet = 19.5
e_jet = 0.85
cj_loiter_jet = 0.5         # (0.4-0.6) [lbs/lbs/hr]
cj_cruise_jet = 0.6         # (0.5-0.9) [lbs/lbs/hr]

# Tbp
A_tbp = 19.5
e_tbp = 0.85
eff_cruise_tbp = 0.85       # [-]
eff_loiter_tbp = 0.77       # [-]
cp_cruise_tbp = 0.5         # (0.4-0.6) [lbs/hp/hr]
cp_loiter_tbp = 0.6         # (0.5-0.7) [lbs/hp/hr]

for iter in range(50):
    MTOW_jet, OEW_jet, W_fuel_jet =  Weights(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, jet = True)
    MTOW_tbp, OEW_tbp, W_fuel_tbp = Weights(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, tbp = True)
    W_empty_jet = 0.54*MTOW_jet
    W_empty_tbp = 0.52*MTOW_tbp

print(MTOW_tbp/g, OEW_tbp/g, W_fuel_tbp/g)
print(MTOW_jet/g, OEW_jet/g, W_fuel_jet/g)
