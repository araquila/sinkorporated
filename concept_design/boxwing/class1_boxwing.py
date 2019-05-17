# IMPORTS
import numpy as np
from fuel_fraction import fuel_fraction
from constant_variables import *

# FUNCTION
def Weights_Class_I(W_empty_jet, W_payload, W_crew, C_fe, S, S_wet, A_jet, e_jet, cj_loiter_jet, cj_cruise_jet, f_trapped_fuel,  V_cruise_jet, range_cruise_jet, endurance_loiter_jet, jet = False, tbp = False):
    # Find C_D_0
    C_D_0 = C_fe * S_wet/S
    print(C_D_0,'cd0')
# THIS SECTION IS USED FOR JET AIRCRAFT; FOR TURBOPROP AIRCRAFT, SCROLL DOWN
    if jet:
    # Determine the maximum L/D for jets
        LD_cruise_jet = 3/4 * np.sqrt((np.pi * A_jet * e_jet) / (3 * C_D_0))
        print(LD_cruise_jet,"L/D")
        # Requirements for loiter
        # CL/CD during loiter for a jet
        C_L_loiter_jet = np.sqrt(C_D_0 * np.pi * A_jet * e_jet)
        C_D_loiter_jet = 2 * C_D_0
        LD_loiter_jet = C_L_loiter_jet / C_D_loiter_jet

        f_fuel_jet, f_reserve_jet, f_cruise_start_jet, f_cruise_end_jet = fuel_fraction(LD_cruise_jet = LD_cruise_jet, cj_cruise_jet = cj_cruise_jet, cj_loiter_jet = cj_loiter_jet, LD_loiter_jet = LD_loiter_jet, V_cruise_jet = V_cruise_jet, range_cruise_jet = range_cruise_jet, endurance_loiter_jet = endurance_loiter_jet, jet = True)

        # Used fuel
        f_used_fuel_jet = 1 - f_fuel_jet

        # Formula for the maximum take-off weight
        MTOW_jet = (W_empty_jet + W_crew + W_payload) / (1 - (f_trapped_fuel + f_used_fuel_jet))

        # Trapped fuel
        W_trapped_fuel_jet = f_trapped_fuel * MTOW_jet

        # Useful fuel
        W_used_fuel_jet = f_used_fuel_jet * MTOW_jet
        # Reserve fuel
        W_reserve_fuel_jet = f_reserve_jet * W_used_fuel_jet
        # Total useful fuel
        W_fuel_jet = W_used_fuel_jet #+ W_reserve_fuel_jet


        # Determine the operative empty weight for a jet
        OEW_jet = W_empty_jet + W_trapped_fuel_jet + W_crew
        print('fuel jet',W_fuel_jet/9.80665)
        return MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet, W_used_fuel_jet, f_cruise_start_jet, f_cruise_end_jet
