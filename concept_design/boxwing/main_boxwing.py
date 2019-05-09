# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from class1_boxwing import Weights_Class_I

# Gravitional constant
g = 9.8065

# Passengers and crew
n_passenger = 60
M_passenger = 105           # (Including luggage)
n_crew = 4
M_crew_member = 100

#M_OEM = 17133.6
#M_MTOM = 25584.5


# Initial mass and fractions
M_payload = n_passenger * M_passenger
M_crew = n_crew * M_crew_member
f_trapped_fuel = 0.003      # Range 0.001-0.005
#M_empty_tbp = M_OEM-M_crew-f_trapped_fuel*M_MTOM
M_empty_jet = 13874.75

#M_empty_jet = M_OEM-M_crew-f_trapped_fuel*M_MTOM


# Convert to weights
W_payload = M_payload * g
W_crew = M_crew * g

W_empty_jet = M_empty_jet * g

## Initial jet and tbp aircraft parameters
C_fe = 0.003
S = 1
S_wet = 5 * S

# Jet
A_jet = 10
e_jet = 1.2
cj_loiter_jet = 0.5         # (0.4-0.6) [lbs/lbs/hr]
cj_cruise_jet = 0.5        # (0.5-0.9) [lbs/lbs/hr]



    # CLASS I ESTIMATION
for iter in range(3):
    MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet =  Weights_Class_I(W_empty_jet, W_payload, W_crew, C_fe, S, S_wet, A_jet, e_jet, cj_loiter_jet, cj_cruise_jet, f_trapped_fuel, jet = True)


    W_empty_jet = (OEW_jet-W_crew)-f_trapped_fuel*MTOW_jet

print(MTOW_jet/g, OEW_jet/g, W_fuel_jet/g)
