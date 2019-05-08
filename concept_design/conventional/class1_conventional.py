
# Trapped fuel
# Range for trapped fuel fractions: 0.001 - 0.005
fraction_trapped_fuel = 0.003
W_trapped_fuel_jet = fraction_trapped_fuel * MTOW_jet

# Useful fuel
W_fuel_jet = W_used_fuel_jet 

# Determine take-off weight for a jet
MTOW_jet = OEW_jet + W_fuel_jet + W_payload
# Determine the operative empty weight for a jet
OEW_jet = W_empty_jet + W_trapped_fuel_jet + W_crew




W_trapped_fuel_tbp = fraction_trapped_fuel * MTOW_tbp

# Determine take-off weight for a turboprop
MTOW_tbp = OEW_tbp + W_fuel_tbp + W_payload
# Determine the operative empty weight for a turboprop
OEW_tbp = W_empty_tbp + W_trapped_fuel_tbp + W_crew
