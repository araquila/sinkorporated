# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM
import os
import sys
sys.path.append(os.getcwd())

# Import modules
from class1_conventional import Weights_Class_I
from power_wingloading_conventional import wingloading_jet, wingloading_tbp
from class1sizing_strutwing import *
from conversion_formulas import *
from class2_strutwing import *
import numpy as np

# Gravitional constant
g = 9.8065

# Passengers and crew
n_passenger = 60
M_passenger = 102                   #(including luggage)
n_crew = 4
M_crew_member = 100

# Initial mass and fractions
M_payload = n_passenger * M_passenger
M_crew = n_crew * M_crew_member
f_trapped_fuel = 0.003              # Range 0.001-0.005
M_empty_tbp = 16000
M_empty_jet = 16300

# Convert to weights
W_payload = M_payload * g
W_crew = M_crew * g
W_empty_tbp = M_empty_tbp * g
W_empty_jet = M_empty_jet * g

# Initial jet and tbp aircraft parameters
C_fe = 0.003
S = 1
S_wet = 5 * S

# Jet
A_jet = 19.5
e_jet = 0.8                         # Adjust per concept
cj_loiter_jet = 19e-6               # (0.4-0.6) [g/j] Propfan: 0.441
cj_cruise_jet = 19e-6               # (0.5-0.9) [g/j] Propfan: 0.441
V_cruise_jet =  200                 # [m/s]
S_jet = 61

# Tbp
A_tbp = 18
e_tbp = 0.85                        # Adjust per concept
eff_cruise_tbp = 0.85               # [-]
eff_loiter_tbp = 0.77               # [-]
cp_cruise_tbp = 90e-9               # (0.4-0.6) [kg/ns]
cp_loiter_tbp = 90e-9               # (0.5-0.7) [kg/ns]
V_cruise_tbp = 180                  # [m/s]
M_cruise_tbp = 0.72                 # [-]
C_L_cruise = 0.8                    # [-]
S_tbp = 66                          # [m^2]
V_stall_tbp = 46.3                  # [m/s]

S_control = 3

# Engine
n_engines = 2                       # [-]
P_TO_tbp = 5.8e6                    # [W]

# Empennage
V_h = 1.57                          # [-]
V_v = 0.07                          # [-]
l_h = 11                            # [m]
l_v = 11                            # [m]

# Fuselage
n_seats_abreast = 4
n_aisles = 1
main_landing_pos = 11               # [m]
nose_landing_pos = 3                # [m]

# Iterator
for iter in range(50):
    MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp, LD_cruise_tbp = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, tbp = True)
    length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside, diameter_fuselage_inside = fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles)
    sweepqc = det_quarter_chord_sweep(M_cruise_tbp)
    dihedral_angle = det_dihedral_angle(sweepqc, high=True)
    b, taper, root_chord, tip_chord, t_c_ratio = det_planform(S_tbp, A_tbp, M_cruise_tbp, C_L_cruise, sweepqc)
    c = 0.5*root_chord + 0.5*tip_chord
    diameter_engine, length_engine, diameter_propeller = enginedimensions(n_engines, P_TO_tbp)
    AR_h, AR_v, S_h, span_h, root_chord_h, tip_chord_h, sweepqc_h, sweepLE_h, S_v, span_v, root_chord_v, tip_chord_v, sweepLE_v = empennage(V_h, V_v, l_h, l_v, S_tbp, b, c)
    wheel_height, lateral_position = undercarriage(main_landing_pos, nose_landing_pos, length_fuselage, length_tail, diameter_fuselage_outside)

#    wingloading_jet(MTOW_jet,OEW_jet,V_cruise_jet,e_jet,C_D_0_jet,A_jet,S_jet)
    #wingloading_tbp(MTOW_tbp, OEW_tbp, S_tbp, A_tbp, V_cruise_tbp, e_tbp, eff_cruise_tbp, C_D_0_tbp)

print(det_wing_weight(kg_to_pounds(MTOW_tbp), ult_load_factor(kg_to_pounds(MTOW_tbp)), metersquared_to_feetsquared(S_tbp), A_tbp, t_c_ratio, taper, sweepqc, metersquared_to_feetsquared(S_control)))
print(det_hor_tail_weight(meter_to_feet(diameter_fuselage_outside), meter_to_feet(span_h), kg_to_pounds(MTOW_tbp),  1.5*ult_load_factor(kg_to_pounds(MTOW_tbp)), meter_to_feet(S_h), meter_to_feet(l_h), np.radians(sweepqc_h), AR_h, meter_to_feet(0.3*S_h)))
print(det_vert_tail_weight(meter_to_feet(span_v), meter_to_feet(span_v), kg_to_pounds(MTOW_tbp), 1.5*ult_load_factor(kg_to_pounds(MTOW_tbp)), l_v, S_v, np.radians(sweepLE_v), AR_v, t_c_ratio))
print(det_fuselage_weight(kg_to_pounds(MTOW_tbp), 1.5*ult_load_factor(kg_to_pounds(MTOW_tbp)), meter_to_feet(length_cabin), metersquared_to_feetsquared(np.pi*diameter_fuselage_outside*length_fuselage), taper, b, sweepqc, LD_cruise_tbp, fuselage_mounted_lg=True))
print(det_main_lg_weight(kg_to_pounds(MTOW_tbp), 4.5, meter_to_inch(wheel_height), 4, 2, ms_to_knots(V_stall_tbp)))
print(det_nose_lg_weight(kg_to_pounds(MTOW_tbp), 4.5, meter_to_inch(wheel_height), 2))



instruments_weight = det_instruments_weight(n_pilots, n_engines, length_fuselage, wingspan, turboprop = True)
hydraulics_weight = det_hydraulics_weight(meter_to_feet(length_fuselage), wingspan)
electrical_weight = det_electrical_weight(meter_to_feet(length_fuselage), N_gen = 3)
avionics_weight = det_avionics_weight()
det_furnishings_weight(n_pilots, 13.608*60, metersquared_to_feetsquared(np.pi * diameter_fuselage_outside * length_fuselage))
pres_vol = np.pi / 4 * diameter_fuselage_inside**2 * (length_nose + length_nose)
aircond_weight = det_aircond_weight(n_passenger, metercubed_to_feetcubed(pres_vol))
anti_ice_weight = det_anti_ice_weight(kg_to_pounds(MTOW_tbp))
handling_gear_weight = det_handling_gear_weight(kg_to_pounds(MTOW_tbp))


MTOM_tbp = MTOW_tbp / g

print('Tbp: ' + str(MTOM_tbp) , str(OEW_tbp/g) , str(W_fuel_tbp/g))
