# IMPORTS
import numpy as np
from fuel_fraction import fuel_fraction

# FUNCTION
def Weights_Class_I_jet(W_empty_jet, W_payload, W_crew, C_fe_jet, S_jet, S_wet_jet, A_jet, e_jet, cj_loiter_jet, cj_cruise_jet, f_trapped_fuel):
# THIS SECTION IS USED FOR JET AIRCRAFT; FOR TURBOPROP AIRCRAFT, SCROLL DOWN
    # Find C_D_0
    C_D_0_jet = C_fe_jet * S_wet_jet/S_jet

    # Improved efficiency compared to regular desicrete wing and body aircraft
    eff_liftingbody = 1.19

# Determine the maximum L/D for jets
    LD_cruise_jet = eff_liftingbody * 3/4 * np.sqrt((np.pi * A_jet * e_jet) / (3 * C_D_0_jet))
    print(LD_cruise_jet)
    # Requirements for loiter
    # CL/CD during loiter for a jet
    C_L_loiter_jet = np.sqrt(C_D_0_jet * np.pi * A_jet * e_jet)
    C_D_loiter_jet = 2 * C_D_0_jet
    LD_loiter_jet = eff_liftingbody *  C_L_loiter_jet / C_D_loiter_jet

    f_fuel_jet, f_reserve_jet, f_cruise_start_jet, f_cruise_end_jet = fuel_fraction(LD_cruise_jet = LD_cruise_jet, cj_cruise_jet = cj_cruise_jet, cj_loiter_jet = cj_loiter_jet, LD_loiter_jet = LD_loiter_jet, jet = True)

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

    # Determine the operative empty weight for a jet
    OEW_jet = W_empty_jet + W_trapped_fuel_jet + W_crew

    return MTOW_jet, OEW_jet, W_fuel_jet, C_D_0_jet, f_cruise_start_jet, f_cruise_end_jet


def Weights_Class_I_tbp(W_empty_tbp, W_payload, W_crew, C_fe_tbp, S_tbp, S_wet_tbp, A_tbp, e_tbp,  eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel):
# THIS SECTION IS USED FOR TURBOPROP AIRCRAFT
    # Find C_D_0
    C_D_0_tbp = C_fe_tbp * S_wet_tbp/S_tbp

    # Improved efficiency compared to regular desicrete wing and body aircraft
    eff_liftingbody = 1.19
# Determine the maximum L/D for a turboprop
    LD_cruise_tbp = eff_liftingbody * np.sqrt((np.pi * A_tbp * e_tbp) / (4 * C_D_0_tbp))
    print(LD_cruise_tbp)
    # Requirements for loiter
    # CL/CD during loiter for a turboprop
    C_L_loiter_tbp = np.sqrt(3 * C_D_0_tbp * np.pi * A_tbp * e_tbp)
    C_D_loiter_tbp = 4 * C_D_0_tbp
    LD_loiter_tbp = eff_liftingbody * C_L_loiter_tbp / C_D_loiter_tbp

    # Determine fuel fractions for a turboprop
    f_fuel_tbp, f_reserve_tbp, f_cruise_start_tbp, f_cruise_end_tbp = fuel_fraction(LD_cruise_tbp = LD_cruise_tbp, eff_cruise_tbp = eff_cruise_tbp, eff_loiter_tbp = eff_loiter_tbp, cp_cruise_tbp = cp_cruise_tbp, cp_loiter_tbp = cp_loiter_tbp, tbp = True)

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

    # Determine the operative empty weight for a turboprop
    OEW_tbp = W_empty_tbp + W_trapped_fuel_tbp + W_crew

    return MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp, f_cruise_start_tbp, f_cruise_end_tbp
