# IMPORTS
import numpy as np
from fuel_fraction import fuel_fraction

# THIS SECTION IS USED FOR JET AIRCRAFT; FOR TURBOPROP AIRCRAFT, SCROLL DOWN

def Weights(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, jet = False, tbp = False):
    # Calculate drag polar
    # Find C_D_0
    C_D_0 = C_fe * S_wet/S
    # Find C_D
    # Assume that e is in the range of 0.80 and 0.85 for turboprops and 0.75 - 0.85 for jets
    # C_D = C_D_0 + C_L**2 / (np.pi * A * e)

    if jet:
    # Determine the maximum L/D for jets
        LD_cruise_jet = 3/4 * np.sqrt((np.pi * A_jet * e_jet) / (3 * C_D_0))

        # Requirements for loiter
        # CL/CD during loiter for a jet
        C_L_loiter_jet = np.sqrt(C_D_0 * np.pi * A_jet * e_jet)
        C_D_loiter_jet = 2 * C_D_0
        LD_loiter_jet = C_L_loiter_jet / C_D_loiter_jet

        f_fuel_jet, f_reserve_jet = fuel_fraction(LD_cruise_jet = LD_cruise_jet, cj_cruise_jet = cj_cruise_jet, cj_loiter_jet = cj_loiter_jet, LD_loiter_jet = LD_loiter_jet, jet = True)

        # Used fuel
        f_used_fuel_jet = 1 - f_fuel_jet

        # Formula for the maximum take-off weight
        MTOW_jet = (W_empty_jet + W_crew + W_payload) / (1 - (f_trapped_fuel + f_used_fuel_jet * (1 + f_reserve_jet)))

        # Trapped fuel
        W_trapped_fuel_jet = f_trapped_fuel * MTOW_jet

        # Useful fuel
        W_used_fuel_jet = f_used_fuel_jet * MTOW_jet
        # Reserve fuel
        W_reserve_fuel_jet = f_reserve_jet * W_used_fuel_jet
        # Total useful fuel
        W_fuel_jet = W_used_fuel_jet + W_reserve_fuel_jet

        # Determine take-off weight for a jet
        # MTOW_jet = OEW_jet + W_fuel_jet + W_payload
        # Determine the operative empty weight for a jet
        OEW_jet = W_empty_jet + W_trapped_fuel_jet + W_crew

        return MTOW_jet, OEW_jet, W_fuel_jet

    if tbp:
    # Determine the maximum L/D for a turboprop
        LD_cruise_tbp = np.sqrt((np.pi * A_tbp * e_tbp) / (4 * C_D_0))

        # Requirements for loiter
        # CL/CD during loiter for a turboprop
        C_L_loiter_tbp = np.sqrt(3 * C_D_0 * np.pi * A_tbp * e_tbp)
        C_D_loiter_tbp = 4 * C_D_0
        LD_loiter_tbp = C_L_loiter_tbp / C_D_loiter_tbp

        # Determine fuel fractions for a turboprop
        f_fuel_tbp, f_reserve_tbp = fuel_fraction(LD_cruise_tbp = LD_cruise_tbp, eff_cruise_tbp = eff_cruise_tbp, eff_loiter_tbp = eff_loiter_tbp, cp_cruise_tbp = cp_cruise_tbp, cp_loiter_tbp = cp_loiter_tbp, tbp = True)

        # Used fuel
        f_used_fuel_tbp = 1 - f_fuel_tbp

        # Formula for the maximum take-off weight
        MTOW_tbp = (W_empty_tbp + W_crew + W_payload) / (1 - (f_trapped_fuel + f_used_fuel_tbp * (1 + f_reserve_tbp)))

        # Trapped fuel
        W_trapped_fuel_tbp = f_trapped_fuel * MTOW_tbp

        # Useful fuel

        W_used_fuel_tbp = f_used_fuel_tbp * MTOW_tbp
        # Reserve fuel
        W_reserve_fuel_tbp = f_reserve_tbp * W_used_fuel_tbp
        # Total useful fuel
        W_fuel_tbp = W_used_fuel_tbp + W_reserve_fuel_tbp

        # Determine take-off weight for a turboprop
        # MTOW_tbp = OEW_tbp + W_fuel_tbp + W_payload
        # Determine the operative empty weight for a turboprop
        OEW_tbp = W_empty_tbp + W_trapped_fuel_tbp + W_crew

        return MTOW_tbp, OEW_tbp, W_fuel_tbp