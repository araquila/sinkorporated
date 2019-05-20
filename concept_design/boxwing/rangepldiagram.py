from math import *
from matplotlib import pyplot as plt
#import the following from your own directory: MTOW, OEW, Fuel weight, L/D in cruise
from main_boxwing import MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet

def payloadrange(MTOWinit, OEWinit, W_fuel_init, LD_cruise_jet, LD_cruise_tbp, A_jet, A_tbp, eff_cruise_tbp, eff_loiter_tbp, e_jet, e_tbp, V_cruise_jet, V_cruise_tbp, V_loiter_tbp, jet = False, tbp = False):

    # Gravitional constant
    g = 9.8065

    # Passengers and crew
    n_passenger = 60
    M_passenger = 105           # (Including luggage)

    # Initial mass and fractions
    M_payload = n_passenger * M_passenger

    # Convert to weights
    W_payload = M_payload * g

    ## Initial jet and tbp aircraft parameters
    C_fe = 0.003
    S = 1
    S_wet = 5 * S

    if jet:
        cj_loiter_jet = 19e-6/1.16       # (0.4-0.6) [lbs/lbs/hr]
        cj_cruise_jet = 19e-6/1.16      # (0.5-0.9) [lbs/lbs/hr]
        range_cruise_jet = 1850000
        endurance_loiter_jet = 2700

        f1_jet = 0.990      # W_1 / W_TO (Engine start, warm-up)
        f2_jet = 0.990      # W_2 / W_1 (Taxi)
        f3_jet = 0.995      # W_3 / W_2 (Take-off)
        f4_jet = 0.980      # W_4 / W_3 (Climb)
        f5_jet = None       # W_5 / W_4 (Cruise)
        f6_jet = None       # W_6 / W_5 (Loiter)
        f7_jet = 0.990      # W_7 / W_6 (Descent)
        f8_jet = 0.992      # W_8 / W_7 (Landing, taxi, shutdown)

        C_D_0 = C_fe * S_wet/S
        C_L_loiter_jet = sqrt(C_D_0 * pi * A_jet * e_jet)
        C_D_loiter_jet = 2 * C_D_0
        LD_loiter_jet = C_L_loiter_jet / C_D_loiter_jet
        f6_jet = 1/exp(endurance_loiter_jet/((1/(g*cj_loiter_jet))*LD_loiter_jet))

        #point A
        RA = 0
        PA = W_payload

        #point B
        Wf = MTOWinit - OEWinit - W_payload
        Mff = 1 - (Wf/MTOWinit)

        W54 = Mff/(f1_jet*f2_jet*f3_jet*f4_jet*f6_jet*f7_jet*f8_jet)
        W45 = 1/W54
        R45 = (V_cruise_jet/(g*cj_cruise_jet))*LD_cruise_jet*log(W45)
        RB = R45
        PB = W_payload

        #point C
        Wpay = MTOWinit - OEWinit - W_fuel_init

        Mff = 1 - (W_fuel_init/MTOWinit)
        W54 = Mff/(f1_jet*f2_jet*f3_jet*f4_jet*f6_jet*f7_jet*f8_jet)
        W45 = 1/W54
        R45 = (V_cruise_jet/(g*cj_cruise_jet))*LD_cruise_jet*log(W45)
        RC = R45
        PC = Wpay

        #point D
        MTOWnew = W_fuel_init + OEWinit
        Wf = W_fuel_init
        Mff = 1 - (Wf/MTOWnew)

        W54 = Mff/(f1_jet*f2_jet*f3_jet*f4_jet*f6_jet*f7_jet*f8_jet)
        W45 = 1/W54
        R45 = (V_cruise_jet/(g*cj_cruise_jet))*LD_cruise_jet*log(W45)
        RD = R45
        PD = 0

        Rlist = [RA/1000, RB/1000, RC/1000, RD/1000]
        Plist = [PA/g, PB/g, PC/g, PD/g]
        print('JET')
        return Rlist, Plist, M_payload

    elif tbp:
        cp_loiter_tbp = 90e-9       # (0.4-0.6) [lbs/lbs/hr]
        cp_cruise_tbp = 90e-9      # (0.5-0.9) [lbs/lbs/hr]
        range_cruise_tbp = 1850000
        endurance_loiter_tbp = 2700

        f1_tbp = 0.990      # W_1 / W_TO (Engine start, warm-up)
        f2_tbp = 0.995      # W_2 / W_1 (Taxi)
        f3_tbp = 0.995      # W_3 / W_2 (Take-off)
        f4_tbp = 0.985      # W_4 / W_3 (Climb)
        f5_tbp = None       # W_5 / W_4 (Cruise)
        f6_tbp = None      # W_6 / W_5 (Loiter)
        f7_tbp = 0.985      # W_7 / W_6 (Descent)
        f8_tbp = 0.995      # W_8 / W_7 (Landing, taxi, shutdown)

        C_D_0 = C_fe * S_wet/S
        #LD_cruise_tbp = sqrt((pi * A_tbp * e_tbp) / (4 * C_D_0))
        C_L_loiter_tbp = sqrt(3 * C_D_0 * pi * A_tbp * e_tbp)
        C_D_loiter_tbp = 4 * C_D_0
        LD_loiter_tbp = C_L_loiter_tbp / C_D_loiter_tbp
        f6_tbp = 1/exp(endurance_loiter_tbp/((eff_loiter_tbp/(V_loiter_tbp*cp_loiter_tbp))*LD_loiter_tbp))

        #point A
        RA = 0
        PA = W_payload

        #point B
        Wfb = MTOWinit - OEWinit - W_payload
        Mff = 1 - (Wfb/MTOWinit)

        W54 = Mff/(f1_tbp*f2_tbp*f3_tbp*f4_tbp*f6_tbp*f7_tbp*f8_tbp)
        W45 = 1/W54
        R45 = (eff_cruise_tbp/(g*cp_cruise_tbp))*LD_cruise_tbp*log(W45)
        RB = R45
        PB = W_payload

        #point C
        Wpay = MTOWinit - OEWinit - W_fuel_init

        Mff = 1 - (W_fuel_init/MTOWinit)
        W54 = Mff/(f1_tbp*f2_tbp*f3_tbp*f4_tbp*f6_tbp*f7_tbp*f8_tbp)
        W45 = 1/W54
        R45 = (eff_cruise_tbp/(g*cp_cruise_tbp))*LD_cruise_tbp*log(W45)
        RC = R45
        PC = Wpay

        #point D
        MTOWnew = W_fuel_init + OEWinit
        Wfd = W_fuel_init
        Mff = 1 - (Wfd/MTOWnew)

        W54 = Mff/(f1_tbp*f2_tbp*f3_tbp*f4_tbp*f6_tbp*f7_tbp*f8_tbp)
        W45 = 1/W54
        R45 = (eff_cruise_tbp/(g*cp_cruise_tbp))*LD_cruise_tbp*log(W45)
        RD = R45
        PD = 0

        Rlist = [RA/1000, RB/1000, RC/1000, RD/1000]
        Plist = [PA/g, PB/g, PC/g, PD/g]
        print('TURBOPROP')
        return Rlist, Plist, M_payload
    else:
        print('wrong entry type')
        return 0, 0, 0


Rlist, Plist, M_payload = payloadrange(MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet, 16, 12, 12,0.85, 0.77, 1.2, 1.2, 229, 229, 100, jet = True, tbp = False)
#Rlist, Plist, M_payload = payloadrange(17832*9.80665,9434*9.80665 , 2097*9.80665, 16, 28.3, 12, 18,0.85, 0.77, 1.2, 0.85, 229, 180, 150, jet = False, tbp = True)


plt.plot(Rlist, Plist)
plt.axis([0,5000, 0, 7000])
plt.ylabel('Payload Mass [kg]', fontsize = 13)
plt.xlabel('Range [km]', fontsize = 13)
plt.show()
