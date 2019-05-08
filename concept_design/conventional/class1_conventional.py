# IMPORTS


# THIS SECTION IS USED FOR JET AIRCRAFT; FOR TURBOPROP AIRCRAFT, SCROLL DOWN

def Weights(jet = False, tbp = False)
    # Calculate drag polar
    # Find C_D_0
    C_D_0 = C_fe * S_wet/S
    # Find C_D
    # Assume that e is in the range of 0.80 and 0.85 for turboprops and 0.75 - 0.85 for jets
    C_D = C_D_0 + C_L**2 / (pi * A * e)

    if jet:
    # Determine the maximum L/D for jets
        LD_cruise_jet = 3/4 * sqrt((pi * A * e) / (3 * C_D_0))

        # Requirements for loiter
        # CL/CD during loiter for a jet
        C_L_loiter_jet = sqrt(C_D_0 * pi * A * e)
        C_D_loiter_jet = 2 * C_D_0
        LD_loiter_jet = C_L_loiter_jet / C_D_loiter_jet

        # Formula for the maximum take-off weight
        MTOW_jet = (W_empty_jet + W_crew + W_payload) / (1 - f_trapped_fuel + f_used_fuel_jet * (1 + f_reserve_fuel))

        # Trapped fuel
        W_trapped_fuel_jet = fraction_trapped_fuel * MTOW_jet

        # Useful fuel
        # Used fuel
        f_used_fuel_jet = 1 - f_fuel_jet
        W_used_fuel_jet = f_used_fuel_jet * MTOW_jet
        # Reserve fuel
        W_reserve_fuel_jet = f_reserve_fuel * W_used_fuel_jet
        # Total useful fuel
        W_fuel_jet = W_used_fuel_jet + W_reserve_fuel_jet

        # Determine take-off weight for a jet
        # MTOW_jet = OEW_jet + W_fuel_jet + W_payload
        # Determine the operative empty weight for a jet
        OEW_jet = W_empty_jet + W_trapped_fuel_jet + W_crew

    if tbp:

W_trapped_fuel_tbp = fraction_trapped_fuel * MTOW_tbp

# Determine take-off weight for a turboprop
MTOW_tbp = OEW_tbp + W_fuel_tbp + W_payload
# Determine the operative empty weight for a turboprop
OEW_tbp = W_empty_tbp + W_trapped_fuel_tbp + W_crew
