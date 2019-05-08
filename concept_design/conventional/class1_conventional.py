# IMPORTS


# THIS SECTION IS USED FOR JET AIRCRAFT; FOR TURBOPROP AIRCRAFT, SCROLL DOWN

# Trapped fuel
W_trapped_fuel_jet = fraction_trapped_fuel * MTOW_jet

# Useful fuel
# Reserve fuel
W_reserve_fuel_jet = f_reserve_fuel * W_used_fuel_jet
# Used fuel
f_used_fuel_jet = 1 - f_fuel_jet
W_used_fuel_jet = f_used_fuel_jet * MTOW_jet
# Total useful fuel
W_fuel_jet = W_used_fuel_jet + W_reserve_fuel_jet

# Determine take-off weight for a jet
MTOW_jet = OEW_jet + W_fuel_jet + W_payload
# Determine the operative empty weight for a jet
OEW_jet = W_empty_jet + W_trapped_fuel_jet + W_crew


# Determine range 





















W_trapped_fuel_tbp = fraction_trapped_fuel * MTOW_tbp

# Determine take-off weight for a turboprop
MTOW_tbp = OEW_tbp + W_fuel_tbp + W_payload
# Determine the operative empty weight for a turboprop
OEW_tbp = W_empty_tbp + W_trapped_fuel_tbp + W_crew
