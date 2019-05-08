# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from fuel_fraction import fuel_fraction
from class1_sizing

# Initial mass and fractions
M_payload = 10
M_crew = 10
f_trapped_fuel = 0.003      # Range 0.001-0.005

## Initial jet and tbp aircraft parameters
C_fe
S
S_wet = 5 * S

# Jet
A_jet =
e_jet =
LD_loiter_jet = 16          # (14-18) [-]
LD_cruise_jet = 12          # [-]
cj_loiter_jet = 0.5         # (0.4-0.6) [lbs/lbs/hr]
cj_cruise_jet = 0.6         # (0.5-0.9) [lbs/lbs/hr]

# Tbp
A_tbp =
e_tbp = 
eff_cruise_tbp = 0.85       # [-]
eff_loiter_tbp = 0.77       # [-]
cp_cruise_tbp = 0.5         # (0.4-0.6) [lbs/hp/hr]
cp_loiter_tbp = 0.6         # (0.5-0.7) [lbs/hp/hr]
LD_loiter_tbp = 15          # (14-16) [-]
LD_cruise_tbp = 14          # [-]

for iter in range(1):

    # Fuel fraction
    f_fuel_jet = fuel_fraction(LD_cruise_jet = LD_cruise_jet, cj_cruise_jet = cj_cruise_jet, cj_loiter_jet = cj_loiter_jet, LD_loiter_jet = LD_loiter_jet, jet = True)
    f_fuel_tbp = fuel_fraction(LD_cruise_tbp = LD_cruise_tbp, eff_cruise_tbp = eff_cruise_tbp, eff_loiter_tbp = eff_loiter_tbp, cp_cruise_tbp = cp_cruise_tbp, cp_loiter_tbp = cp_loiter_tbp, tbp = True)
