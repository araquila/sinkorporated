# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from fuel_fraction import fuel_fraction

# Weights
W_payload = 10
W_crew = 10
MTOW_jet = 10
MTOW_tbp = 10

# Fuel fractions
ff_jet = fuel_fraction(jet = True)
ff_tbp = fuel_fraction(tbp = True)
