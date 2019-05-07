# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:31:19 2019

@author: Matthijs Torsij
"""
# Import packages
import numpy as np

# MTOW input variables
MTOW_tbp = 10               # [N]
MTOW_jet = 8                # [N]

# Cruise input variables
range_cruise_tbp = 500      # [km]
range_cruise_jet = 600      # [km]
V_cruise_jet = 100          # [km/hr]

# Loiter input variables
endurance_loiter_tbp = 1    # [minutes]
endurance_loiter_jet = 1    # [minutes]
V_loiter_tbp = 100          # [km/u]

# Unit conversions cruise
range_cruise_jet = range_cruise_jet * 0.621371192       # [-> miles]
range_cruise_tbp = range_cruise_tbp * 0.621371192       # [-> miles]
V_cruise_jet = V_cruise_jet * 0.621371192               # [-> miles/hr]

# Unit conversions loiter
endurance_loiter_jet = endurance_loiter_jet / 60        # [-> hr]
endurance_loiter_tbp = endurance_loiter_tbp / 60        # [-> hr]
V_loiter_tbp = V_loiter_tbp * 0.621371192               # [-> miles/hr]

# Jet and tbp characteristics (range)[unit]
eff_cruise_tbp = 0.85       # [-]
eff_loiter_tbp = 0.77       # [-]
cp_cruise_tbp = 0.5         # (0.4-0.6) [lbs/hp/hr]
cj_cruise_jet = 0.6         # (0.5-0.9) [lbs/lbs/hr]
cp_loiter_tbp = 0.6         # (0.5-0.7) [lbs/hp/hr]
cj_loiter_jet = 0.5         # (0.4-0.6) [lbs/lbs/hr]
LD_cruise_tbp = 12          # (11-13) [-]
LD_cruise_jet = 14          # (13-15) [-]
LD_loiter_tbp = 15          # (14-16) [-]
LD_loiter_jet = 16          # (14-18) [-]

# Fuel fractions from Roskam for regional tbp
f1_tbp = 0.990      # W_1 / W_TO (Engine start, warm-up)
f2_tbp = 0.995      # W_2 / W_1 (Taxi)
f3_tbp = 0.995      # W_3 / W_2 (Take-off)
f4_tbp = 0.985      # W_4 / W_3 (Climb)
f5_tbp = None       # W_5 / W_4 (Cruise)
f6_tbp = None       # W_6 / W_5 (Loiter)
f7_tbp = 0.985      # W_7 / W_6 (Descent)
f8_tbp = 0.995      # W_8 / W_7 (Landing, taxi, shutdown)

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
f5_tbp = 1/(np.exp(range_cruise_tbp/(375*(eff_cruise_tbp/cp_cruise_tbp)*LD_cruise_tbp)))
f5_jet = 1/(np.exp(range_cruise_jet/((V_cruise_jet/cj_cruise_jet)*LD_cruise_jet)))

# Calculation of loiter fuel fraction
f6_tbp = 1/()
f6_jet = 1/()
