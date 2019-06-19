# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM
#import os
#import sys
#sys.path.append(os.getcwd())

# Import modules
from concept_design.strutwing.class1_conventional import Weights_Class_I
from concept_design.strutwing.power_wingloading_conventional_redone import wingloading_jet, wingloading_tbp
import concept_design.strutwing.sustainability_functions as sf
from concept_design.strutwing.wingloadingfunctions import T_W_calc, W_P_climb_calc
from concept_design.strutwing.class1sizing_strutwing import *
from concept_design.strutwing.conversion_formulas import *
from concept_design.strutwing.class2_strutwing import *
from concept_design.strutwing.atmosphere import atmosphere_calc
import numpy as np
import concept_design.strutwing.fuel_fraction as ff
import matplotlib.pyplot as plt
import concept_design.strutwing.rangepldiagram as pld
import concept_design.strutwing.cost_equations as ceq
from concept_design.strutwing.LNG import *
import detailed_design.parameters as p
from concept_design.strutwing.class1sizing_strutwing import fuselage
from concept_design.strutwing.class1sizing_strutwing import det_planform

#Do you want a pie chart?
piechart = False
print_payloadrange = True

altitude = 8000
temperature0 = 288.15
temperature_gradient = -0.0065
gamma = 1.4
rho0 = 1.225
g = 9.80665
R = 287

temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, temperature0, temperature_gradient, g, R, gamma)

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
M_empty_tbp = p.W_empty / p.g
#M_empty_jet = 16300

# Convert to weights
W_payload = M_payload * g
W_crew = M_crew * g
W_empty_tbp = M_empty_tbp  * g
#W_empty_jet = p.M_empty_jet * g

# Initial jet and tbp aircraft parameters
C_fe = 0.003
S = p.S
S_wet = 5.5 * S


energy_density_LNG = 53.6 #[MJ/kg]
energy_density_kerosene = 43 #[MJ/kg]
chosen_fuel_energy_density = energy_density_LNG
fuel_efficiency_factor = energy_density_kerosene/chosen_fuel_energy_density

## Jet
#A_jet = 19.5
#e_jet = 0.8                         # Adjust per concept
#cj_loiter_jet = fuel_efficiency_factor*12.5e-6 # oude waarde 19e-6               # (0.4-0.6) [g/j] Propfan: 0.441
#cj_cruise_jet = fuel_efficiency_factor*12.5e-6 # oude waarde 19e-6               # (0.5-0.9) [g/j] Propfan: 0.441
#V_cruise_jet =  200                 # [m/s]
#V_loiter_jet = 150
#S_jet = 61

# Tbp
A_tbp = p.A
e_tbp = 0.85                        # Adjust per concept
eff_cruise = 0.85               # [-]
eff_loiter = 0.77               # [-]
cp_cruise = p.cp_cruise # oude waarde 90e-9              # (0.4-0.6) [kg/ns]
cp_loiter = p.cp_loiter # oude waarde 90e-9               # (0.5-0.7) [kg/ns]
M_cruise_tbp = 0.55                  # [-]
V_cruise_tbp = M_cruise_tbp*speed_of_sound     # [m/s]
V_loiter_tbp = 80                   # [m/s]
C_L_cruise = p.C_L_cruise                    # [-]
S_tbp = p.S                       # [m^2]
V_stall_tbp = p.V_stall                 # [m/s]
C_L_max_land_tbp = p.C_L_max_land              # [-]
C_L_max_TO_tbp = p.C_L_max_TO                # [-]

range_cruise_tbp = 1850000          # [m]
endurance_loiter_tbp = 2700         # [s]


# Engine
n_engines = 2                       # [-]
P_TO_tbp = p.P_TO                   # [W]
pos_engine = p.y_engine                     # [m]
mass_engine = p.M_engine                   # [kg]
n_fueltanks = 2                     # [-]
n_blades = 6                        # [-]

# Empennage
V_h = 1.87                        # [-]
V_v = 0.07                          # [-]
l_h = p.l_h                          # [m]
l_v = p.l_v                           # [m]

# Fuselage
n_seats_abreast = 4
n_aisles = 1
main_landing_pos = p.main_landing_pos               # [m]
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
iterlist = []
emptylist = []
maxlist =[]
# Iterator
for iter in range(10):
#    MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp, f_cruise_start_tbp, f_cruise_end_tbp, LD_cruise_tbp, L_D_loiter_tbp = Weights_Class_I(0, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, 0, A_tbp, 0, e_tbp, 0, 0, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, 0, V_loiter_tbp, 0, range_cruise_tbp, 0, endurance_loiter_tbp, tbp = True, jet = False)

    MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp, f_cruise_start_tbp, f_cruise_end_tbp, LD_cruise_tbp, L_D_loiter_tbp = Weights_Class_I(0, W_empty_tbp, W_payload, W_crew, C_fe, p.S, S_wet, 0, p.A, 0, p.e, 0, 0, eff_loiter, eff_cruise, cp_loiter, cp_cruise, f_trapped_fuel, 0, V_loiter_tbp, 0, range_cruise_tbp, 0, endurance_loiter_tbp, tbp = True, jet = False)

#    MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp, f_cruise_start_tbp, f_cruise_end_tbp, LD_cruise_tbp, L_D_loiter_tbp = p.MTOW, p.OEW, p.W_fuel,p.Cd0,  p.f_cruise_start, p.f_cruise_end, p.LD_ratio, p.LD_loiter

    class1["MTOW"].append(MTOW_tbp)
    class1["OEW"].append(OEW_tbp)
    class1["Fuel"].append(W_fuel_tbp)

    #Wing and Power Loading
    W_S_landing_tbp, W_P_critical = wingloading_tbp(MTOW_tbp, OEW_tbp, S_tbp, A_tbp, V_cruise_tbp, e_tbp, eff_cruise, C_D_0_tbp, C_L_max_land_tbp, C_L_max_TO_tbp)

    # S for jet and tbp
    S_tbp = p.S
    P_TO_tbp = p.P_TO
    MTOM_tbp = MTOW_tbp / p.g

    length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside, diameter_fuselage_inside = fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles)
    class1sizing_fuselage["Length Fuselage"].append(length_fuselage)
    class1sizing_fuselage["Diameter Fuselage"].append(diameter_fuselage_outside)

#    sweepqc = det_quarter_chord_sweep(M_cruise_tbp)
    sweepqc = p.sweep_qc
#    dihedral_angle = det_dihedral_angle(sweepqc, high=True)
    dihedral_angle = p.dihedral
#    b, taper, root_chord, tip_chord, t_c_ratio = det_planform(S_tbp, A_tbp, M_cruise_tbp, C_L_cruise, sweepqc)
    b, taper, root_chord, tip_chord, t_c_ratio = p.b, p.taper, p.root_chord, p.tip_chord, p.tc_ratio_root
#    c = 0.5*root_chord + 0.5*tip_chord
    c = p.MAC
    class1sizing_wing["Wing Span"].append(b)
    class1sizing_wing["Quarter Chord Sweep"].append(sweepqc)
    class1sizing_wing["Dihedral"].append(dihedral_angle)
    class1sizing_wing["Taper"].append(taper)
    class1sizing_wing["Root Chord"].append(root_chord)
    class1sizing_wing["T/C Ratio"].append(t_c_ratio)
    class1sizing_wing["Aspect Ratio"].append(A_tbp)

    diameter_engine, length_engine, diameter_propeller = 0.7, p.l_engine, p.d_prop
    class1sizing_engine["Engine Diameter"].append(diameter_engine)
    class1sizing_engine["Engine Length"].append(length_engine)
    class1sizing_engine["Propeller Diameter"].append(diameter_propeller)

#    AR_h, AR_v, S_h, span_h, root_chord_h, tip_chord_h, sweepqc_h, sweepLE_h, S_v, span_v, root_chord_v, tip_chord_v, sweepLE_v = detailed_design.empennage.empennage(V_h, V_v, l_h, l_v, S_tbp, b, c)
    
    AR_h, AR_v, S_h, span_h, root_chord_h, tip_chord_h, sweepqc_h, sweepLE_h, S_v, span_v, root_chord_v, tip_chord_v, sweepLE_v = p.A_h, p.A_v, p.S_h, p.b_h, p.root_chord_h, p.tip_chord_h, p.sweep_h, p.sweep_h, p.S_v, p.b_v, p.root_chord_v, p.tip_chord_v, p.sweep_v
    
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

#    wheel_height, lateral_position = undercarriage(main_landing_pos, nose_landing_pos, length_fuselage, length_tail, diameter_fuselage_outside)
    wheel_height, lateral_position = p.h_wheel, p.lateral_pos
    
    class1sizing_gear["Gear Vertical Height"].append(wheel_height)  
    class1sizing_gear["Gear Lateral Position"].append(lateral_position)


    #Methane tank mass
#    W_tank = LNG_system(2, LNG_volume(W_fuel_tbp/g), podded = True, shell_length_ratio =  0.5, tank_circular_ratio = 1, tank_head_ratio = 3)
    W_tank = p.W_fuel_tanks / g

    #W_wing = pounds_to_kg(0.82 * det_wing_weight(kg_to_pounds(MTOM_tbp), 1.5*ult_load_factor(kg_to_pounds(MTOM_tbp)), metersquared_to_feetsquared(S_tbp), A_tbp, t_c_ratio, taper, np.radians(sweepqc), metersquared_to_feetsquared(0.05*S_tbp)))
#    W_wing = det_wing_weight_new(b, S_tbp, sweepqc, taper, MTOM_tbp, V_cruise_tbp, t_c_ratio)
    W_wing = p.W_wing / p.g
    W_h = pounds_to_kg(det_hor_tail_weight(meter_to_feet(diameter_fuselage_outside), meter_to_feet(span_h), kg_to_pounds(MTOM_tbp),  1.5*ult_load_factor(kg_to_pounds(MTOM_tbp)), metersquared_to_feetsquared(S_h), meter_to_feet(l_h), np.radians(sweepqc_h), AR_h, metersquared_to_feetsquared(0.3*S_h)))
    W_v = pounds_to_kg(det_vert_tail_weight(meter_to_feet(span_v), meter_to_feet(span_v), kg_to_pounds(MTOM_tbp), 1.5*ult_load_factor(kg_to_pounds(MTOM_tbp)), meter_to_feet(l_v), metersquared_to_feetsquared(S_v), np.radians(sweepLE_v), AR_v, t_c_ratio))
    W_fus = pounds_to_kg(det_fuselage_weight(kg_to_pounds(MTOM_tbp), 1.5*ult_load_factor(kg_to_pounds(MTOM_tbp)), meter_to_feet(length_fuselage), metersquared_to_feetsquared(np.pi*diameter_fuselage_outside*length_fuselage), taper, meter_to_feet(b), sweepqc, LD_cruise_tbp, fuselage_mounted_lg=True))
    W_ml = pounds_to_kg(det_main_lg_weight(kg_to_pounds(MTOM_tbp), 4.5, meter_to_inch(wheel_height), 4, 2, ms_to_knots(V_stall_tbp)))
    W_nl = pounds_to_kg(det_nose_lg_weight(kg_to_pounds(MTOM_tbp), 4.5, meter_to_inch(wheel_height), 2))
#    W_nacelle = pounds_to_kg(det_nacelle_group_weight(meter_to_feet(length_engine), meter_to_feet(diameter_engine), 1.5*ult_load_factor(kg_to_pounds(MTOM_tbp)), n_engines, metersquared_to_feetsquared(np.pi * diameter_engine * length_engine), W_ec = kg_to_pounds(mass_engine), propeller=True, thrust_reverser=True))
    W_nacelle = 357.38
    W_engine = 2 * mass_engine + 330
#    W_engine_controls = pounds_to_kg(det_engine_controls_weight(n_engines, n_engines*meter_to_feet(pos_engine)))
    W_engine_controls = 0
#    W_starter = pounds_to_kg(det_starter_weight(n_engines, kg_to_pounds(mass_engine)))
    W_starter = 0
#    W_fuel_system = pounds_to_kg(det_fuel_system_weight(kg_to_pounds(W_fuel_tbp/g)/6.67632, kg_to_pounds(W_fuel_tbp/g)/6.67632, 0, n_fueltanks)) + W_tank
    W_fuel_system = 358.7959 + W_tank
    W_flight_control = pounds_to_kg(det_flight_controls_weight(metersquared_to_feetsquared(0.3*S_h+0.05*S_tbp), (meter_to_feet(length_fuselage)**2 * kg_to_pounds(MTOM_tbp)*0.34**2)/(4*32.19)))
#    APU_weight = det_APU_weight(200)
    APU_weight = 61.23
    
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
    iterlist.append(iter)
    emptylist.append(total)
    maxlist.append(MTOM_tbp)
print('MTOW:', MTOM_tbp, 'kg')
print('OEW:', OEW_tbp/g, 'kg')
print('FUEL:', W_fuel_tbp/g, 'kg')
print('Wing area:', S_tbp, 'm^2')
plt.figure(1)
plt.plot(iterlist, emptylist)
plt.figure(2)
plt.plot(iterlist, maxlist)
#Design Cruise CL
q_tbp = 0.5 * 0.525168 * V_cruise_tbp**2
#C_L_cruise_tbp = C_L_des(q_tbp,f_cruise_start_tbp*MTOW_tbp/S_tbp,f_cruise_end_tbp*MTOW_tbp/S_tbp)
#C_D_cruise_tbp = C_D_0_tbp + (1/(A_tbp*e_tbp*np.pi)) * C_L_cruise_tbp**2

#print('Design CL:', C_L_cruise_tbp)
#print('Design CD:', C_D_cruise_tbp)
#print('Resulting CL/CD:', C_L_cruise_tbp/C_D_cruise_tbp)

f_fuel_tbp_1000, f_reserve_tbp_1000, f_cruise_start_tbp_1000, f_cruise_end_tbp_1000 = ff.fuel_fraction(LD_loiter_tbp = p.LD_loiter, LD_cruise_tbp = p.LD_ratio, eff_cruise_tbp = eff_cruise, eff_loiter_tbp = eff_loiter, cp_cruise_tbp = cp_cruise, cp_loiter_tbp = cp_loiter, V_loiter_tbp = V_loiter_tbp, range_cruise_tbp = 1000000, endurance_loiter_tbp = endurance_loiter_tbp, tbp = True)
W_fuel_tbp_1000 = (1 - f_fuel_tbp_1000) * MTOW_tbp
fuel_per_passenger_tbp_1000 = (W_fuel_tbp_1000/n_passenger)/g

print('Fuel burn 1000km trip:',fuel_per_passenger_tbp_1000,'kg')


#Emission Calculations
CO2_emission = sf.CO2_calc(fuel_per_passenger_tbp_1000, 53.6)
NOX_emission = sf.NOX_calc(fuel_per_passenger_tbp_1000, 53.6)

print('CO2 emission 1000km trip:',CO2_emission,'kg')
print('NOX emission 1000km trip:',NOX_emission,'kg')

#Noise Calculations
noise_prop = sf.prop_noise(diameter_propeller, n_blades, 1200, P_TO_tbp/n_engines, n_engines, 308.063)
noise_airframe = sf.airframe_noise(V_cruise_tbp, MTOW_tbp)
noise_total = sf.total_noise(noise_prop, noise_airframe)
noise_at_distance = sf.noise_distance(noise_total,1,2500)

print('Total noise production:', noise_at_distance, 'dB')

#Payload Range Diagram
range_list, payload_list, M_payload = pld.payloadrange(MTOW_tbp, OEW_tbp, W_fuel_tbp, 0, LD_cruise_tbp, 0, A_tbp, eff_cruise, eff_loiter, 0, e_tbp, 0, V_cruise_tbp, V_loiter_tbp, jet = False, tbp = True)

#if print_payloadrange:
#    plt.plot(range_list, payload_list)
#    plt.axis([0,5000, 0, 7000])
#    plt.ylabel('Payload Mass [kg]', fontsize = 13)
#    plt.xlabel('Range [km]', fontsize = 13)
#    plt.show()

#if piechart:
#    patches, texts, values = plt.pie(weight_fractions, counterclock = False, startangle=90, autopct='%1.11f%%')
#    plt.legend(patches, weight_label, loc="best", fontsize = 'x-large')
#    plt.axis('equal')
#    plt.show()

print('Development cost :', ceq.non_recurring_cost(W_wing,W_h+W_v,W_fus,W_nl+W_ml,W_engine, W_engine_controls + W_starter + W_fuel_system + W_flight_control + W_instruments + W_hydraulics + W_electrical + W_avionics + W_furnishings + W_airco + W_anti_ice + W_handling_gear, W_payload/g),'Million USD (2019)')
print('Production cost per unit :', ceq.recurring_cost(500,W_wing,W_h+W_v,W_fus,W_nl+W_ml,W_engine, W_engine_controls + W_starter + W_fuel_system + W_flight_control + W_instruments + W_hydraulics + W_electrical + W_avionics + W_furnishings + W_airco + W_anti_ice + W_handling_gear, W_payload/g,W_empty_tbp/g)/500,'Million USD (2019)')
print('Total cost per unit', ceq.total_cost(500,W_wing,W_h+W_v,W_fus,W_nl+W_ml,W_engine, W_engine_controls + W_starter + W_fuel_system + W_flight_control + W_instruments + W_hydraulics + W_electrical + W_avionics + W_furnishings + W_airco + W_anti_ice + W_handling_gear, W_payload/g,W_empty_tbp/g),'Million USD (2019)')

print()
print('-----------------  Sizes for CATIA   ------------------')
print("Nose length =",length_nose)
print("Cabin length =", length_cabin)
print("Tail length =", length_tail)
print("Wing area =", S_tbp)
print(lateral_position)
