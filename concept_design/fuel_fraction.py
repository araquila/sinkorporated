# Fuel fraction Calculation
# Import packages
import numpy as np
import variables

def fuel_fraction(jet = False, tbp = False):

    # Unit conversions cruise
    range_cruise_jet = range_cruise_jet * 0.621371192       # [km -> miles]
    range_cruise_tbp = range_cruise_tbp * 0.621371192       # [km -> miles]
    V_cruise_jet = V_cruise_jet * 0.621371192               # [km -> miles/hr]

    # Unit conversions loiter
    endurance_loiter_jet = endurance_loiter_jet / 60        # [minutes -> hr]
    endurance_loiter_tbp = endurance_loiter_tbp / 60        # [minutes -> hr]
    V_loiter_tbp = V_loiter_tbp * 0.621371192               # [km/hr -> miles/hr]

    if tbp:
        # Fuel fractions from Roskam for regional tbp
        f1_tbp = 0.990      # W_1 / W_TO (Engine start, warm-up)
        f2_tbp = 0.995      # W_2 / W_1 (Taxi)
        f3_tbp = 0.995      # W_3 / W_2 (Take-off)
        f4_tbp = 0.985      # W_4 / W_3 (Climb)
        f5_tbp = None       # W_5 / W_4 (Cruise)
        f6_tbp = None       # W_6 / W_5 (Loiter)
        f7_tbp = 0.985      # W_7 / W_6 (Descent)
        f8_tbp = 0.995      # W_8 / W_7 (Landing, taxi, shutdown)

        # Calculation of cruise fuel fraction
        f5_tbp = 1/np.exp(range_cruise_tbp/(375*(eff_cruise_tbp/cp_cruise_tbp)*LD_cruise_tbp))

        # Calculation of loiter fuel fraction
        f6_tbp = 1/np.exp(endurance_loiter_tbp/((375/V_loiter_tbp)*(eff_loiter_tbp/cp_loiter_tbp)*LD_loiter_tbp))

        # Insert fuel fractions in list
        ff_tbp = [f1_tbp, f2_tbp, f3_tbp, f4_tbp, f5_tbp, f6_tbp, f7_tbp, f8_tbp]

        return ff_tbp

    if jet:
        # Fuel fractions from Roskam for jet
        f1_jet = 0.990      # W_1 / W_TO (Engine start, warm-up)
        f2_jet = 0.990      # W_2 / W_1 (Taxi)
        f3_jet = 0.995      # W_3 / W_2 (Take-off)
        f4_jet = 0.980      # W_4 / W_3 (Climb)
        f5_jet = None       # W_5 / W_4 (Cruise)
        f6_jet = None       # W_6 / W_5 (Loiter)
        f7_jet = 0.990      # W_7 / W_6 (Descent)
        f8_jet = 0.992      # W_8 / W_7 (Landing, taxi, shutdown)

        # Calculation of cruise fuel fraction
        f5_jet = 1/np.exp(range_cruise_jet/(((V_cruise_jet/cj_cruise_jet)*LD_cruise_jet)))

        # Calculation of loiter fuel fraction
        f6_jet = 1/np.exp(endurance_loiter_jet/((1/cj_loiter_jet)*LD_loiter_jet))

        # Insert fuel fractions in list
        ff_jet = [f1_jet, f2_jet, f3_jet, f4_jet, f5_jet, f6_jet, f7_jet, f8_jet]

        return ff_jet

    else:
        message = 'Specify type of aircraft (tbp or jet).'
        return message

ff = fuel_fraction(500,)
