# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM
import os
import sys
sys.path.append(os.getcwd())

# Import modules
from class1_conventional import Weights_Class_I
from power_wingloading_conventional_redone import wingloading_jet, wingloading_tbp
import sustainability_functions as sf
from wingloadingfunctions import T_W_calc, W_P_climb_calc
from class1sizing_strutwing import *
from conversion_formulas import *
from class2_strutwing import *
import numpy as np
import fuel_fraction as ff
import matplotlib.pyplot as plt
import rangepldiagram as pld
import cost_equations as ceq

#Do you want a pie chart?
piechart = False
print_payloadrange = False


# Gravitional constant
g = 9.8065

# Passengers and crew
n_passenger = 60
M_passenger = 105                   #(including luggage)
n_crew = 4
n_pilots = 2
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
V_loiter_jet = 150
S_jet = 61

# Tbp
A_tbp = 18
e_tbp = 0.85                        # Adjust per concept
eff_cruise_tbp = 0.85               # [-]
eff_loiter_tbp = 0.77               # [-]
cp_cruise_tbp = 90e-9               # (0.4-0.6) [kg/ns]
cp_loiter_tbp = 90e-9               # (0.5-0.7) [kg/ns]
V_cruise_tbp = 180                  # [m/s]
V_loiter_tbp = 80                   # [m/s]
M_cruise_tbp = 0.72                 # [-]
C_L_cruise = 0.8                    # [-]
S_tbp = 66                          # [m^2]
V_stall_tbp = 46.3                  # [m/s]
C_L_max_land_tbp = 2.4              # [-]
C_L_max_TO_tbp = 1.4                # [-]

range_cruise_tbp = 1850000          # [m]
endurance_loiter_tbp = 2700         # [s]


# Engine
n_engines = 2                       # [-]
P_TO_tbp = 5.8e6                    # [W]
pos_engine = 10                     # [m]
mass_engine = 450                   # [kg]
n_fueltanks = 2                     # [-]

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

class1 = {"MTOW": [], "OEW": [], "Fuel": []}
class1sizing_fuselage = {"Length Fuselage": [], "Diameter Fuselage": []}
class1sizing_wing = {"Wing Span": [], "Quarter Chord Sweep": [], "Dihedral": [], "Taper": [], "Root Chord": [], "T/C Ratio": [], "Aspect Ratio": []}
class1sizing_engine = {"Engine Diameter": [], "Engine Length": [], "Propeller Diameter": []}
class1sizing_htail = {"Horizontal Tail Span": [], "Quarter Chord Sweep": [], "Taper": [], "Root Chord": [], "T/C Ratio": [], "Aspect Ratio": []}
class1sizing_vtail = {"Vertical Tail Span": [], "Leading Edge Sweep": [], "Taper": [], "Root Chord": [], "T/C Ratio": [], "Aspect Ratio": []}
class1sizing_gear = {"Gear Vertical Height": [], "Gear Lateral Position": []}

class2 = {"wing weight": [], "horizontal tail weight": [], "vertical tail weight": [], "fuselage weight": [], "main landing gear weight": [], "nose landing gear weight": [], "nacelle group weight": [], "engine weight": [], "engine controls weight": [], "starter weight": [], "fuel system weight" : []}
class2["flight controls weight"] = []
class2["APU weight"] = []
class2["instruments weight"] = []
class2["hydraulics weight"] = []
class2["electrical system weight"] = []
class2["avionics weight"] = []
class2["furnishings weight"] = []
class2["airconditioning weight"] = []
class2["anti icing system weight"] = []
class2["handling gear weight"] = []


# Iterator
for iter in range(5):
    MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp, f_cruise_start_tbp, f_cruise_end_tbp, LD_cruise_tbp, L_D_loiter_tbp = Weights_Class_I(0, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, 0, A_tbp, 0, e_tbp, 0, 0, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, 0, V_loiter_tbp, 0, range_cruise_tbp, 0, endurance_loiter_tbp, tbp = True, jet = False)

    class1["MTOW"].append(MTOW_tbp)
    class1["OEW"].append(OEW_tbp)
    class1["Fuel"].append(W_fuel_tbp)

    #Wing and Power Loading
    W_S_landing_tbp, W_P_critical = wingloading_tbp(MTOW_tbp, OEW_tbp, S_tbp, A_tbp, V_cruise_tbp, e_tbp, eff_cruise_tbp, C_D_0_tbp, C_L_max_land_tbp, C_L_max_TO_tbp)

    # S for jet and tbp
    S_tbp = MTOW_tbp/W_S_landing_tbp
    P_TO_tbp = MTOW_tbp / W_P_critical
    MTOM_tbp = MTOW_tbp / g

    length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside, diameter_fuselage_inside = fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles)
    class1sizing_fuselage["Length Fuselage"].append(length_fuselage)
    class1sizing_fuselage["Diameter Fuselage"].append(diameter_fuselage_outside)

    sweepqc = det_quarter_chord_sweep(M_cruise_tbp)
    dihedral_angle = det_dihedral_angle(sweepqc, high=True)
    b, taper, root_chord, tip_chord, t_c_ratio = det_planform(S_tbp, A_tbp, M_cruise_tbp, C_L_cruise, sweepqc)
    c = 0.5*root_chord + 0.5*tip_chord
    class1sizing_wing["Wing Span"].append(b)
    class1sizing_wing["Quarter Chord Sweep"].append(sweepqc)
    class1sizing_wing["Dihedral"].append(dihedral_angle)
    class1sizing_wing["Taper"].append(taper)
    class1sizing_wing["Root Chord"].append(root_chord)
    class1sizing_wing["T/C Ratio"].append(t_c_ratio)
    class1sizing_wing["Aspect Ratio"].append(A_tbp)

    diameter_engine, length_engine, diameter_propeller = enginedimensions(n_engines, P_TO_tbp)
    class1sizing_engine["Engine Diameter"].append(diameter_engine)
    class1sizing_engine["Engine Length"].append(length_engine)
    class1sizing_engine["Propeller Diameter"].append(diameter_propeller)

    AR_h, AR_v, S_h, span_h, root_chord_h, tip_chord_h, sweepqc_h, sweepLE_h, S_v, span_v, root_chord_v, tip_chord_v, sweepLE_v = empennage(V_h, V_v, l_h, l_v, S_tbp, b, c)
    class1sizing_htail["Horizontal Tail Span"].append(span_h)
    class1sizing_htail["Quarter Chord Sweep"].append(sweepqc_h)
    class1sizing_htail["Taper"].append(tip_chord_h/root_chord_h)
    class1sizing_htail["Root Chord"].append(root_chord_h)
    class1sizing_htail["T/C Ratio"].append(0)
    class1sizing_htail["Aspect Ratio"].append(AR_h)

    class1sizing_vtail["Vertical Tail Span"].append(span_v)
    class1sizing_vtail["Leading Edge Sweep"].append(sweepLE_v)
    class1sizing_vtail["Taper"].append(tip_chord_v/root_chord_v)
    class1sizing_vtail["Root Chord"].append(root_chord_v)
    class1sizing_vtail["T/C Ratio"].append(0)
    class1sizing_vtail["Aspect Ratio"].append(AR_v)

    wheel_height, lateral_position = undercarriage(main_landing_pos, nose_landing_pos, length_fuselage, length_tail, diameter_fuselage_outside)
    class1sizing_gear["Gear Vertical Height"].append(wheel_height)
    class1sizing_gear["Gear Lateral Position"].append(lateral_position)


    #W_wing = pounds_to_kg(0.82 * det_wing_weight(kg_to_pounds(MTOM_tbp), 1.5*ult_load_factor(kg_to_pounds(MTOM_tbp)), metersquared_to_feetsquared(S_tbp), A_tbp, t_c_ratio, taper, np.radians(sweepqc), metersquared_to_feetsquared(0.05*S_tbp)))
    W_wing = det_wing_weight_new(b, S_tbp, sweepqc, taper, MTOM_tbp, V_cruise_tbp, t_c_ratio)
    W_h = pounds_to_kg(det_hor_tail_weight(meter_to_feet(diameter_fuselage_outside), meter_to_feet(span_h), kg_to_pounds(MTOM_tbp),  1.5*ult_load_factor(kg_to_pounds(MTOM_tbp)), metersquared_to_feetsquared(S_h), meter_to_feet(l_h), np.radians(sweepqc_h), AR_h, metersquared_to_feetsquared(0.3*S_h)))
    W_v = pounds_to_kg(det_vert_tail_weight(meter_to_feet(span_v), meter_to_feet(span_v), kg_to_pounds(MTOM_tbp), 1.5*ult_load_factor(kg_to_pounds(MTOM_tbp)), meter_to_feet(l_v), metersquared_to_feetsquared(S_v), np.radians(sweepLE_v), AR_v, t_c_ratio))
    W_fus = pounds_to_kg(det_fuselage_weight(kg_to_pounds(MTOM_tbp), 1.5*ult_load_factor(kg_to_pounds(MTOM_tbp)), meter_to_feet(length_fuselage), metersquared_to_feetsquared(np.pi*diameter_fuselage_outside*length_fuselage), taper, meter_to_feet(b), sweepqc, LD_cruise_tbp, fuselage_mounted_lg=True))
    W_ml = pounds_to_kg(det_main_lg_weight(kg_to_pounds(MTOM_tbp), 4.5, meter_to_inch(wheel_height), 4, 2, ms_to_knots(V_stall_tbp)))
    W_nl = pounds_to_kg(det_nose_lg_weight(kg_to_pounds(MTOM_tbp), 4.5, meter_to_inch(wheel_height), 2))
    W_nacelle = pounds_to_kg(det_nacelle_group_weight(meter_to_feet(length_engine), meter_to_feet(diameter_engine), 1.5*ult_load_factor(kg_to_pounds(MTOM_tbp)), n_engines, metersquared_to_feetsquared(np.pi * diameter_engine * length_engine), W_ec = kg_to_pounds(mass_engine), propeller=True, thrust_reverser=True))
    W_engine = 2 * mass_engine
    W_engine_controls = pounds_to_kg(det_engine_controls_weight(n_engines, n_engines*meter_to_feet(pos_engine)))
    W_starter = pounds_to_kg(det_starter_weight(n_engines, kg_to_pounds(mass_engine)))
    W_fuel_system = pounds_to_kg(det_fuel_system_weight(kg_to_pounds(W_fuel_tbp/g)/6.67632, kg_to_pounds(W_fuel_tbp/g)/6.67632, 0, n_fueltanks))
    W_flight_control = pounds_to_kg(det_flight_controls_weight(metersquared_to_feetsquared(0.3*S_h+0.05*S_tbp), (meter_to_feet(length_fuselage)**2 * kg_to_pounds(MTOM_tbp)*0.34**2)/(4*32.19)))
    APU_weight = det_APU_weight(200)
    W_instruments = pounds_to_kg(det_instruments_weight(n_pilots, n_engines, meter_to_feet(length_fuselage), meter_to_feet(b), turboprop = True))
    W_hydraulics = pounds_to_kg(det_hydraulics_weight(meter_to_feet(length_fuselage), meter_to_feet(b)))
    W_electrical = pounds_to_kg(det_electrical_weight(meter_to_feet(length_fuselage), N_gen = 3))
    W_avionics = pounds_to_kg(det_avionics_weight())
    W_furnishings = pounds_to_kg(det_furnishings_weight(n_pilots, 30*60, metersquared_to_feetsquared(np.pi * diameter_fuselage_outside * length_fuselage)))
    pres_vol = np.pi / 4 * diameter_fuselage_inside**2 * (length_nose + length_cabin)
    W_airco = pounds_to_kg(det_aircond_weight(n_passenger+n_crew, metercubed_to_feetcubed(pres_vol)))
    W_anti_ice = pounds_to_kg(det_anti_ice_weight(kg_to_pounds(MTOM_tbp)))
    W_handling_gear = pounds_to_kg(det_handling_gear_weight(kg_to_pounds(MTOM_tbp)))
    class2["wing weight"].append(W_wing)
    class2["horizontal tail weight"].append(W_h)
    class2["vertical tail weight"].append(W_v)
    class2["fuselage weight"].append(W_fus)
    class2["main landing gear weight"].append(W_ml)
    class2["nose landing gear weight"].append(W_nl)
    class2["nacelle group weight"].append(W_nacelle)
    class2["engine weight"].append(W_engine)
    class2["engine controls weight"].append(W_engine_controls)
    class2["starter weight"].append(W_starter)
    class2["fuel system weight"].append(W_fuel_system)
    class2["flight controls weight"].append(W_flight_control)
    class2["APU weight"].append(APU_weight)
    class2["instruments weight"].append(W_instruments)
    class2["hydraulics weight"].append(W_hydraulics)
    class2["electrical system weight"].append(W_electrical)
    class2["avionics weight"].append(W_avionics)
    class2["furnishings weight"].append(W_furnishings)
    class2["airconditioning weight"].append(W_airco)
    class2["anti icing system weight"].append(W_anti_ice)
    class2["handling gear weight"].append(W_handling_gear)

    total = 0
    weight_fractions = []
    weight_label = []
    for item in class2:
        total += class2[item][-1]
        weight_fractions.append(class2[item][-1])
        weight_label.append(item)

    W_empty_tbp = total * g

print('MTOW:', MTOM_tbp, 'kg')
print('OEW:', OEW_tbp/g, 'kg')
print('FUEL:', W_fuel_tbp/g, 'kg')

#Design Cruise CL
q_tbp = 0.5 * 0.525168 * V_cruise_tbp**2
C_L_cruise_tbp = C_L_des(q_tbp,f_cruise_start_tbp*MTOW_tbp/S_tbp,f_cruise_end_tbp*MTOW_tbp/S_tbp)
C_D_cruise_tbp = C_D_0_tbp + (1/(A_tbp*e_tbp*np.pi)) * C_L_cruise_tbp**2

print('Design CL:', C_L_cruise_tbp)
print('Design CD:', C_D_cruise_tbp)
print('Resulting CL/CD:', C_L_cruise_tbp/C_D_cruise_tbp)

f_fuel_tbp_1000, f_reserve_tbp_1000, f_cruise_start_tbp_1000, f_cruise_end_tbp_1000 = ff.fuel_fraction(LD_loiter_tbp = L_D_loiter_tbp, LD_cruise_tbp = LD_cruise_tbp, eff_cruise_tbp = eff_cruise_tbp, eff_loiter_tbp = eff_loiter_tbp, cp_cruise_tbp = cp_cruise_tbp, cp_loiter_tbp = cp_loiter_tbp, V_loiter_tbp = V_loiter_tbp, range_cruise_tbp = 1000000, endurance_loiter_tbp = endurance_loiter_tbp, tbp = True)
W_fuel_tbp_1000 = (1 - f_fuel_tbp_1000) * MTOW_tbp
fuel_per_passenger_tbp_1000 = (W_fuel_tbp_1000/n_passenger)/g

print('Fuel burn 1000km trip:',fuel_per_passenger_tbp_1000,'kg')


#Emission Calculations
CO2_emission = sf.CO2_calc(W_fuel_tbp/(g*n_passenger), 43)
NOX_emission = sf.NOX_calc(W_fuel_tbp/(g*n_passenger), 43)

#Noise Calculations
#noise_prop = sf.prop_noise(diameter_propeller, 4, n_p, P_TO_tbp/n_engines, N_p, n_engines)
#noise_airframe = sf.airframe_noise(V_cruise_tbp, MTOW_tbp)
#noise_total = sf.total_noise(noise_prop, noise_airframe)

#Payload Range Diagram
range_list, payload_list, M_payload = pld.payloadrange(MTOW_tbp, OEW_tbp, W_fuel_tbp, 0, LD_cruise_tbp, 0, A_tbp, eff_cruise_tbp, eff_loiter_tbp, 0, e_tbp, 0, V_cruise_tbp, V_loiter_tbp, jet = False, tbp = True)

if print_payloadrange:
    plt.plot(range_list, payload_list)
    plt.xlim([0,4000])
    plt.ylim([0,7000])
    plt.title('Payload Range Diagram')
    plt.ylabel('Payload Mass [kg]')
    plt.xlabel('Range [km]')
    plt.show()

if piechart:
    patches, texts, values = plt.pie(weight_fractions, counterclock = False, startangle=90, autopct='%1.11f%%')
    plt.legend(patches, weight_label, loc="best", fontsize = 'x-large')
    plt.axis('equal')
    plt.show()

print('Development cost :', ceq.non_recurring_cost(W_wing,W_h+W_v,W_fus,W_nl+W_ml,W_engine, W_engine_controls + W_starter + W_fuel_system + W_flight_control + W_instruments + W_hydraulics + W_electrical + W_avionics + W_furnishings + W_airco + W_anti_ice + W_handling_gear, W_payload/g),'Million USD (2019)')
print('Production cost per unit :', ceq.recurring_cost(500,W_wing,W_h+W_v,W_fus,W_nl+W_ml,W_engine, W_engine_controls + W_starter + W_fuel_system + W_flight_control + W_instruments + W_hydraulics + W_electrical + W_avionics + W_furnishings + W_airco + W_anti_ice + W_handling_gear, W_payload/g,W_empty_tbp/g)/500,'Million USD (2019)')
print('Total cost per unit', ceq.total_cost(500,W_wing,W_h+W_v,W_fus,W_nl+W_ml,W_engine, W_engine_controls + W_starter + W_fuel_system + W_flight_control + W_instruments + W_hydraulics + W_electrical + W_avionics + W_furnishings + W_airco + W_anti_ice + W_handling_gear, W_payload/g,W_empty_tbp/g),'Million USD (2019)')
