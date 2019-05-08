# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from fuel_fraction import fuel_fraction
from class1_sizing

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
