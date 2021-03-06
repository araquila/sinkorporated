# FUEL FRACTIONS
# Inputs are: LD_cruise_jet or LD_cruise_tbp and V_cruise_jet or V_loiter_tbp and jet or tbp = True
#test
# Import modules
import numpy as np

g = 9.80665

def fuel_fraction(eff_cruise_tbp = 0.85, eff_loiter_tbp = 0.77, LD_cruise_jet = 12, cp_cruise_tbp = 0.5, LD_cruise_tbp = 14, cj_cruise_jet = 0.6, cp_loiter_tbp = 0.6, cj_loiter_jet = 0.5, LD_loiter_tbp = 15, LD_loiter_jet = 16 ,V_cruise_jet = 5, V_loiter_tbp = 5, range_cruise_jet = 1850000, range_cruise_tbp = 1850000, endurance_loiter_jet = 2700, endurance_loiter_tbp = 2700, jet = False, tbp = False):

    # CONVERSIONS
    # Unit conversions cruise
#    range_cruise_jet = range_cruise_jet * 0.001              # [m -> km]
#    range_cruise_tbp = range_cruise_tbp * 0.001              # [m -> km]
#    V_cruise_jet = V_cruise_jet * 2.23693629                # [m/s -> miles/hr]

    # Unit conversions loiter
#    endurance_loiter_jet = endurance_loiter_jet             # [sec -> min]
#    endurance_loiter_tbp = endurance_loiter_tbp             # [sec -> min]
#    V_loiter_tbp = V_loiter_tbp * 2.23693629                # [m/s -> miles/hr]

    # FUEL FRACTION FOR TBP
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
        f5_tbp = 1/np.exp(range_cruise_tbp/((eff_cruise_tbp/(g*cp_cruise_tbp))*LD_cruise_tbp))
#        print('tbp cruise: ' + str(f5_tbp))

        # Calculation of loiter fuel fraction
        f6_tbp = 1/np.exp(endurance_loiter_tbp/((eff_loiter_tbp/(V_loiter_tbp*cp_loiter_tbp))*LD_loiter_tbp))
#        print('tbp loiter: ' + str(f6_tbp))

        # Find fuel fraction tbp
        f_fuel_tbp = f1_tbp * f2_tbp * f3_tbp * f4_tbp * f5_tbp * f6_tbp * f7_tbp * f8_tbp

        f_reserve_tbp = f6_tbp
        f_cruise_start_tbp = f1_tbp*f2_tbp*f3_tbp*f4_tbp
        f_cruise_end_tbp = f_cruise_start_tbp*f5_tbp

        return f_fuel_tbp, f_reserve_tbp,f_cruise_start_tbp, f_cruise_end_tbp

    # FUEL FRACTIONS FOR JET
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
        f5_jet = 1/np.exp(range_cruise_jet/(((V_cruise_jet/(g*cj_cruise_jet))*LD_cruise_jet)))
#        print('jet cruise: ' + str(f5_jet))

        # Calculation of loiter fuel fraction
        f6_jet = 1/np.exp(endurance_loiter_jet/((1/(g*cj_loiter_jet))*LD_loiter_jet))
#        print('jet loiter: ' + str(f6_jet))

        # Find fuel fraction jet
        f_fuel_jet = f1_jet * f2_jet * f3_jet * f4_jet * f5_jet * f6_jet * f7_jet * f8_jet
        f_reserve_jet = f6_jet

        f_cruise_start_jet = f1_jet*f2_jet*f3_jet*f4_jet
        f_cruise_end_jet = f_cruise_start_jet*f5_jet

        return f_fuel_jet, f_reserve_jet, f_cruise_start_jet, f_cruise_end_jet

    # NO AIRCRAFT TYPE SPECIFIED
    else:
        message = 'Specify type of aircraft (tbp or jet).'
        return message
