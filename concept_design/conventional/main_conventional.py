# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from constant_variables import *
from class1_conventional import Weights_Class_I
from power_wingloading_conventional_redone import wingloading_jet, wingloading_tbp
from wingloadingfunctions import T_W_calc, W_P_climb_calc
from class1sizing_conventional import fuselage, det_quarter_chord_sweep, det_planform, det_dihedral_angle, enginedimensions, MAC, empennage, undercarriage, tiresizing
from atmosphere import atmosphere_calc
from cg_determination import x_lemac_tbp_calc, x_lemac_jet_calc
from fuel_fraction import fuel_fraction
from conversion_formulas import *
import class2_conventional as class2
from sustainability_functions import CO2_calc, NOX_calc, prop_noise, airframe_noise, total_noise, noise_distance, turbofan_noise
import numpy as np

## INPUTS AND CONSTANTS
# fuel efficiency
chosen_fuel_energy_density = energy_density_kerosene
fuel_efficiency_factor = energy_density_kerosene/chosen_fuel_energy_density

# Flight parameters
s_landing = 1400                    #[m]
altitude = 8000
V_landing = 48.93                   #[m/s] maximum landing speed that is allowed on a runway of 1400 m this is set for all aircraft

# Atmospherical parameters at cruise altitude
temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, temperature0, temperature_gradient, g, R, gamma)

# Initial jet and tbp aircraft parameters
C_fe = 0.003
S = 1
S_wet = 5 * S
c = 10                                  #[m/s]

# Other jet parameters
A_jet = 10
e_jet = 0.8
V_cruise_jet =  236.11                  # [m/s]
S_jet = 61
TOP_jet = 6698
M_cruise_jet = V_cruise_jet/speed_of_sound
C_L_cruise_jet = 0.4
C_L_max_jet = 2.3
C_L_max_land_jet = 2.3
C_L_max_TO_jet = 1.7
range_cruise_jet = 1850000               # [m]
endurance_loiter_jet = 2700              # [s]

# Empennage jet
V_h_jet = 1.07                           # [-]
V_v_jet = 0.085                          # [-]
nose_landing_pos_jet = 3                 # [m]

# Other tbp parameters
A_tbp = 16
e_tbp = 0.85                             # Adjust per concept
V_loiter_tbp = 80                        # [m/s]
V_cruise_tbp = 180                       # [m/s]
M_cruise_tbp = V_cruise_tbp/speed_of_sound
C_L_cruise_tbp = 0.8
S_tbp = 76
TOP_tbp = 139
C_L_max_tbp = 2.6
C_L_max_land_tbp = 2.6
C_L_max_TO_tbp = 1.6
range_cruise_tbp = 1850000               # [m]
endurance_loiter_tbp = 2700              # [s]

# Empennage tbp
V_h_tbp = 1.57                           # [-]
V_v_tbp = 0.07                           # [-]
nose_landing_pos_tbp = 3                 # [m]

# Dynamic pressure
q_jet = 0.5*rho*V_cruise_jet**2          # [n/m2]
q_tbp = 0.5*rho*V_cruise_tbp**2          # [n/m2]

# Engine characteristics
thrust_to_weight_jet = 73.21         # [N/kg] #add 2/3 if propfan is used
cj_loiter_jet = fuel_efficiency_factor*19.8e-6                  # (0.4-0.6) [g/j] Propfan: 0.441
cj_cruise_jet = fuel_efficiency_factor*19.8e-6                  # (0.5-0.9) [g/j] Propfan: 0.441

power_to_weight_tbp = 4000               # [W/kg]
eff_cruise_tbp = 0.85                    # [-]
eff_loiter_tbp = 0.77                    # [-]
cp_cruise_tbp = 0.8*fuel_efficiency_factor * 90e-9              # (0.4-0.6) [kg/ns]
cp_loiter_tbp = 0.8*fuel_efficiency_factor * 90e-9              # (0.5-0.7) [kg/ns]

# Iterative sizing process
for iter in range(10):

    ## CLASS I
    MTOW_jet, OEW_jet, W_fuel_jet, C_D_0_jet, f_cruise_start_jet, f_cruise_end_jet, L_D_jet = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, V_cruise_jet, V_loiter_tbp, range_cruise_jet, range_cruise_tbp, endurance_loiter_jet, endurance_loiter_tbp, jet = True, tbp = False)
    MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp, f_cruise_start_tbp, f_cruise_end_tbp, L_D_tbp = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, V_cruise_jet, V_loiter_tbp, range_cruise_jet, range_cruise_tbp, endurance_loiter_jet, endurance_loiter_tbp, tbp = True, jet = False)

    MTOM_jet = MTOW_jet/g
    MTOM_tbp = MTOW_tbp/g

    #--------------------------------------------------------------------------

    ## WING LOADING AND POWER LOADING
    W_S_landing_jet = wingloading_jet(MTOW_jet,OEW_jet,V_cruise_jet,e_jet,C_D_0_jet,A_jet,S_jet,C_L_max_land_jet,C_L_max_TO_jet)
    W_S_landing_tbp, W_P_critical = wingloading_tbp(MTOW_tbp, OEW_tbp, S_tbp, A_tbp, V_cruise_tbp, e_tbp, eff_cruise_tbp, C_D_0_tbp, C_L_max_land_tbp,C_L_max_TO_tbp)

    # T/W for jet
    T_W_jet_range = T_W_calc(W_S_landing_jet,TOP_jet,1.32)

    # S for jet and tbp
    S_jet = MTOW_jet/W_S_landing_jet
    S_tbp = MTOW_tbp/W_S_landing_tbp

    #--------------------------------------------------------------------------

    ## SIZING
    # Fuselage
    length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside, diameter_fuselage_inside = fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles)

    # Wing tbp
    sweep_tbp = det_quarter_chord_sweep(M_cruise_tbp)
    b_tbp, taper_tbp, root_chord_tbp, tip_chord_tbp, t_c_ratio_tbp = det_planform(S_tbp, A_tbp, M_cruise_tbp, C_L_cruise_tbp, sweep_tbp)
    dihedral_angle_tbp = det_dihedral_angle(sweep_tbp, high=True)
    MAC_tbp = MAC(root_chord_tbp, t_c_ratio_tbp)

    # Wing jet
    sweep_jet = det_quarter_chord_sweep(M_cruise_jet)
    b_jet, taper_jet, root_chord_jet, tip_chord_jet, t_c_ratio_jet = det_planform(S_jet, A_jet, M_cruise_jet, C_L_cruise_jet, sweep_jet)
    dihedral_angle_jet = det_dihedral_angle(sweep_tbp, low=True)
    MAC_jet = MAC(root_chord_jet, t_c_ratio_jet)

    # Engines for jet and tbp
    P_TO_tbp = MTOW_tbp / W_P_critical                      # Take-off power tbp [W]
    T_TO_jet = T_W_jet_range * MTOW_jet              # Take-off thrust jet [N]
    diameter_engine_tbp, length_engine_tbp, diameter_propeller_tbp = enginedimensions(rho0,n_engines_tbp, P_TO_tbp, T_TO_jet, tbp=True)
    length_nacelle_jet, length_fan_cowling_jet, diameter_highlight_jet, diameter_exit_fan_jet, diameter_gas_generator_jet, diameter_nacelle_jet = enginedimensions(rho0,n_engines_jet, P_TO_tbp, T_TO_jet, jettypeB=True)

    #CG and undercarriage
    x_lemac_tbp = x_lemac_tbp_calc(0.4*length_fuselage,MAC_tbp)
    x_lemac_jet = x_lemac_jet_calc(0.4*length_fuselage,MAC_jet)

    l_h_jet, l_v_jet = 0.9*length_fuselage-(x_lemac_jet+0.25*MAC_jet), 0.9*length_fuselage-(x_lemac_jet+0.25*MAC_jet)                           # [m]
    l_h_tbp, l_v_tbp = 0.9*length_fuselage-(x_lemac_tbp+0.25*MAC_tbp), 0.9*length_fuselage-(x_lemac_tbp+0.25*MAC_tbp)

    main_landing_pos_jet = x_lemac_jet+0.4*MAC_jet      # [m]
    main_landing_pos_tbp = x_lemac_tbp+0.4*MAC_tbp      # [m]

    AR_h_jet, AR_v_jet, S_h_jet, span_h_jet, root_chord_h_jet, tip_chord_h_jet, sweepqc_h_jet, sweepLE_h_jet, S_v_jet, span_v_jet, root_chord_v_jet, tip_chord_v_jet, sweepLE_v_jet = empennage(V_h_jet, V_v_jet, l_h_jet, l_v_jet, S_jet, b_jet, MAC_jet)
    wheel_height_jet, lateral_position_jet = undercarriage(main_landing_pos_jet, nose_landing_pos_jet, length_fuselage, length_tail, diameter_fuselage_outside)
    AR_h_tbp, AR_v_tbp, S_h_tbp, span_h_tbp, root_chord_h_tbp, tip_chord_h_tbp, sweepqc_h_tbp, sweepLE_h_tbp, S_v_tbp, span_v_tbp, root_chord_v_tbp, tip_chord_v_tbp, sweepLE_v_tbp = empennage(V_h_tbp, V_v_tbp, l_h_tbp, l_v_tbp, S_tbp, b_tbp, MAC_tbp)
    wheel_height_tbp, lateral_position_tbp = undercarriage(main_landing_pos_tbp, nose_landing_pos_tbp, length_fuselage, length_tail, diameter_fuselage_outside)

    #--------------------------------------------------------------------------

    ## CLASS II
    # C_L and C_l des
    C_L_cruise_jet = class2.C_L_des(q_jet,f_cruise_start_jet*MTOW_jet/S_jet,f_cruise_end_jet*MTOW_jet/S_jet)
    C_l_des_jet = class2.C_l_des(C_L_cruise_jet,sweep_jet)
    C_L_cruise_tbp = class2.C_L_des(q_tbp,f_cruise_start_tbp*MTOW_tbp/S_tbp,f_cruise_end_tbp*MTOW_tbp/S_tbp)
    C_l_des_tbp = class2.C_l_des(C_L_cruise_tbp,sweep_tbp)

    # Wing weight
    n_max_jet = class2.ult_load_factor(kg_to_pounds(MTOM_jet))
    n_max_tbp = class2.ult_load_factor(kg_to_pounds(MTOM_tbp))

    qc_sweep_tbp = det_quarter_chord_sweep(M_cruise_tbp)
    qc_sweep_jet = det_quarter_chord_sweep(M_cruise_jet)

    wing_weight_jet = class2.det_wing_weight_revised(0.034,b_jet,S_jet,qc_sweep_jet,taper_jet,MTOM_jet,1.5*n_max_jet,1.4*V_cruise_jet,t_c_ratio_jet)
    wing_weight_tbp = class2.det_wing_weight_revised(0.034,b_tbp,S_tbp,qc_sweep_tbp,taper_tbp,MTOM_tbp,1.5*n_max_tbp,1.4*V_cruise_tbp,t_c_ratio_tbp)

#    wing_weight_jet = pounds_to_kg(class2.det_wing_weight(kg_to_pounds(MTOM_jet),(n_max_jet*1.5),metersquared_to_feetsquared(S_jet),A_jet,t_c_ratio_jet,taper_jet,qc_sweep_jet,metersquared_to_feetsquared(0.05*S_jet)))
#    wing_weight_tbp = pounds_to_kg(class2.det_wing_weight(kg_to_pounds(MTOM_tbp),(n_max_tbp*1.5),metersquared_to_feetsquared(S_tbp),A_tbp,t_c_ratio_tbp,taper_tbp,qc_sweep_tbp,metersquared_to_feetsquared(0.05*S_tbp)))

    # Horizontal tail
    hor_tail_weight_jet = pounds_to_kg(class2.det_hor_tail_weight(meter_to_feet(diameter_fuselage_outside),meter_to_feet(span_h_jet),kg_to_pounds(MTOW_jet/g),n_max_jet*1.5,metersquared_to_feetsquared(S_h_jet),meter_to_feet(l_h_jet),np.radians(sweepqc_h_jet),AR_h_jet,metersquared_to_feetsquared(0.3*S_h_jet)))
    hor_tail_weight_tbp = pounds_to_kg(class2.det_hor_tail_weight(meter_to_feet(diameter_fuselage_outside),meter_to_feet(span_h_tbp),kg_to_pounds(MTOW_tbp/g),n_max_tbp*1.5,metersquared_to_feetsquared(S_h_tbp),meter_to_feet(l_h_tbp),np.radians(sweepqc_h_tbp),AR_h_tbp,metersquared_to_feetsquared(0.3*S_h_tbp)))

    # Vertical tail
    ver_tail_weight_jet = pounds_to_kg(class2.det_vert_tail_weight(0,1,kg_to_pounds(MTOW_jet/g),n_max_jet*1.5,meter_to_feet(l_h_jet),metersquared_to_feetsquared(S_v_jet),np.radians(sweepLE_v_jet),AR_v_jet,t_c_ratio_jet))
    ver_tail_weight_tbp = pounds_to_kg(class2.det_vert_tail_weight(0,1,kg_to_pounds(MTOW_tbp/g),n_max_tbp*1.5,meter_to_feet(l_h_tbp),metersquared_to_feetsquared(S_v_tbp),np.radians(sweepLE_v_tbp),AR_v_tbp,t_c_ratio_tbp))

    # Fuselage
    fuselage_weight_jet = pounds_to_kg(class2.det_fuselage_weight(kg_to_pounds(MTOW_jet/g), 1.5*n_max_jet, meter_to_feet(length_fuselage), metersquared_to_feetsquared(np.pi*diameter_fuselage_outside*length_fuselage), taper_jet, meter_to_feet(b_jet), sweep_jet, L_D_jet, cargo_doors = 1, fuselage_mounted_lg = False))
    fuselage_weight_tbp = pounds_to_kg(class2.det_fuselage_weight(kg_to_pounds(MTOW_tbp/g), 1.5*n_max_tbp, meter_to_feet(length_fuselage), metersquared_to_feetsquared(np.pi*diameter_fuselage_outside*length_fuselage), taper_tbp, meter_to_feet(b_tbp), sweep_tbp, L_D_tbp, cargo_doors = 1, fuselage_mounted_lg = False))

    # Main landing gear
    V_stall_jet = np.sqrt(2*MTOW_jet/(rho0*S_jet*C_L_max_jet))
    V_stall_tbp = np.sqrt(2*MTOW_tbp/(rho0*S_tbp*C_L_max_tbp))
    main_lg_weight_jet = pounds_to_kg(class2.det_main_lg_weight(kg_to_pounds(MTOW_jet/g),1.5*n_max_jet,meter_to_inch(wheel_height_jet),4,2,ms_to_knots(V_stall_jet),kneeling_main_lg = False))
    main_lg_weight_tbp = pounds_to_kg(class2.det_main_lg_weight(kg_to_pounds(MTOW_tbp/g),1.5*n_max_tbp,meter_to_inch(wheel_height_tbp),4,2,ms_to_knots(V_stall_tbp),kneeling_main_lg = False))

    # Nose landing gear
    nose_lg_weight_jet = pounds_to_kg(class2.det_nose_lg_weight(kg_to_pounds(MTOW_jet/g), 1.5*n_max_jet, meter_to_inch(wheel_height_jet), 2, kneeling_nose_lg = False))
    nose_lg_weight_tbp = pounds_to_kg(class2.det_nose_lg_weight(kg_to_pounds(MTOW_tbp/g), 1.5*n_max_tbp, meter_to_inch(wheel_height_tbp), 2, kneeling_nose_lg = False))

    # Engine weight
    M_engine_tbp = P_TO_tbp / power_to_weight_tbp
    M_engine_jet = T_TO_jet / thrust_to_weight_jet

    # Nacelle
    nacelle_group_weight_jet = pounds_to_kg(class2.det_nacelle_group_weight(meter_to_feet(length_nacelle_jet), meter_to_feet(diameter_nacelle_jet), 1.5*n_max_jet, 2, metersquared_to_feetsquared(np.pi * diameter_nacelle_jet * length_nacelle_jet), pylon_mounted = True, W_ec = 0, W_engine = kg_to_pounds(M_engine_jet/2), propeller = False, thrust_reverser = False))
    nacelle_group_weight_tbp = pounds_to_kg(class2.det_nacelle_group_weight(meter_to_feet(length_engine_tbp), meter_to_feet(diameter_engine_tbp), 1.5*n_max_tbp, 2, metersquared_to_feetsquared(np.pi * diameter_engine_tbp * length_engine_tbp), pylon_mounted = True, W_ec = 0, W_engine = kg_to_pounds(M_engine_jet/2), propeller = True, thrust_reverser = False))

    # Engine controls weight
    engine_controls_weight_jet = pounds_to_kg(class2.det_engine_controls_weight(2,40))
    engine_controls_weight_tbp = pounds_to_kg(class2.det_engine_controls_weight(2,40))

    # Starter weight
    starter_weight_jet = pounds_to_kg(class2.det_starter_weight(2, kg_to_pounds(M_engine_jet/2)))
    starter_weight_tbp = pounds_to_kg(class2.det_starter_weight(2, kg_to_pounds(M_engine_tbp/2)))

    # Fuel system
    n_fueltanks = 2
    W_fuel_system_jet = pounds_to_kg(class2.det_fuel_system_weight(kg_to_pounds(W_fuel_jet/g)/6.67632, kg_to_pounds(W_fuel_jet/g)/6.67632, 0, n_fueltanks))
    W_fuel_system_tbp = pounds_to_kg(class2.det_fuel_system_weight(kg_to_pounds(W_fuel_tbp/g)/6.67632, kg_to_pounds(W_fuel_tbp/g)/6.67632, 0, n_fueltanks))

    # Flight controls
    flight_controls_weight_jet = pounds_to_kg(class2.det_flight_controls_weight(metersquared_to_feetsquared(0.3*S_h_jet+0.05*S_jet), (meter_to_feet(length_fuselage)**2*kg_to_pounds(MTOM_jet)*0.34**2)/(4*32.19), N_f = 6, N_m = 1))
    flight_controls_weight_tbp = pounds_to_kg(class2.det_flight_controls_weight(metersquared_to_feetsquared(0.3*S_h_tbp+0.05*S_tbp), (meter_to_feet(length_fuselage)**2*kg_to_pounds(MTOM_tbp)*0.34**2)/(4*32.19), N_f = 6, N_m = 1))

    # Instruments
    instruments_weight_jet = pounds_to_kg(class2.det_instruments_weight(2, 2, meter_to_feet(length_fuselage), meter_to_feet(b_jet), reciprocating = False, turboprop = True))
    instruments_weight_tbp = pounds_to_kg(class2.det_instruments_weight(2, 2, meter_to_feet(length_fuselage), meter_to_feet(b_tbp), reciprocating = False, turboprop = True))

    # Hydraulics
    hydraulics_weight_jet = pounds_to_kg(class2.det_hydraulics_weight(meter_to_feet(length_fuselage), meter_to_feet(b_jet), N_f = 6))
    hydraulics_weight_tbp = pounds_to_kg(class2.det_hydraulics_weight(meter_to_feet(length_fuselage), meter_to_feet(b_tbp), N_f = 6))

    # Electrical
    electrical_weight_jet = pounds_to_kg(class2.det_electrical_weight(200, R_kva = 50, N_gen = 0, N_en = 2)) #plug in better numbers
    electrical_weight_tbp = pounds_to_kg(class2.det_electrical_weight(200, R_kva = 50, N_gen = 0, N_en = 2)) #plug in better numbers

    # Avionics
    avionics_weight_jet = pounds_to_kg(class2.det_avionics_weight(W_uav = 1100))
    avionics_weight_tbp = pounds_to_kg(class2.det_avionics_weight(W_uav = 1100))

    # Furnishings
    furnishings_weight_jet = pounds_to_kg(class2.det_furnishings_weight(2, kg_to_pounds(13.608*60), metersquared_to_feetsquared(np.pi * diameter_fuselage_outside * length_fuselage)))
    furnishings_weight_tbp = pounds_to_kg(class2.det_furnishings_weight(2, kg_to_pounds(13.608*60), metersquared_to_feetsquared(np.pi * diameter_fuselage_outside * length_fuselage)))

    # Airconditioning
    pres_vol = np.pi / 4 * diameter_fuselage_inside**2 * (length_nose + length_nose)
    aircond_weight_jet = pounds_to_kg(class2.det_aircond_weight(n_passenger+n_crew, metercubed_to_feetcubed(pres_vol), W_uav = 1100))
    aircond_weight_tbp = pounds_to_kg(class2.det_aircond_weight(n_passenger+n_crew, metercubed_to_feetcubed(pres_vol), W_uav = 1100))

    # Anti-ice
    anti_ice_weight_jet = pounds_to_kg(class2.det_anti_ice_weight(kg_to_pounds(MTOW_jet)))
    anti_ice_weight_tbp = pounds_to_kg(class2.det_anti_ice_weight(kg_to_pounds(MTOW_tbp)))

    # Handling gear
    handling_gear_weight_jet = pounds_to_kg(class2.det_handling_gear_weight(kg_to_pounds(MTOW_jet)))
    handling_gear_weight_tbp = pounds_to_kg(class2.det_handling_gear_weight(kg_to_pounds(MTOW_tbp)))

    # Total weight
    M_empty_jet = wing_weight_jet + hor_tail_weight_jet + ver_tail_weight_jet + fuselage_weight_jet + main_lg_weight_jet + nose_lg_weight_jet + nacelle_group_weight_jet+ engine_controls_weight_jet + starter_weight_jet + W_fuel_system_jet+ flight_controls_weight_jet + instruments_weight_jet + hydraulics_weight_jet + electrical_weight_jet + avionics_weight_jet + furnishings_weight_jet + aircond_weight_jet + anti_ice_weight_jet + handling_gear_weight_jet + M_engine_jet
    M_empty_tbp = wing_weight_tbp + hor_tail_weight_tbp + ver_tail_weight_tbp + fuselage_weight_tbp + main_lg_weight_tbp + nose_lg_weight_tbp + nacelle_group_weight_tbp + engine_controls_weight_tbp + starter_weight_tbp + W_fuel_system_tbp+ flight_controls_weight_tbp + instruments_weight_tbp + hydraulics_weight_tbp + electrical_weight_tbp + avionics_weight_tbp + furnishings_weight_jet + aircond_weight_tbp + anti_ice_weight_tbp + handling_gear_weight_tbp + M_engine_tbp

    W_empty_jet =  M_empty_jet * g
    W_empty_tbp =  M_empty_tbp * g

    fuel_per_passenger_jet = (W_fuel_jet/n_passenger)/g
    fuel_per_passenger_tbp = (W_fuel_tbp/n_passenger)/g

def print_list(items):
    print(items[0])
    for item in items[1:]:
        print(item[0] + ': ' + str(item[1]))
    print()
    print('------------------------------------------------------')
    print()

def print_mass_data():
    mass_data_jet = []
    mass_data_tbp = []

    mass_data_jet.append('### Masses of components for jet in [kg] ###')
    mass_data_tbp.append('### Masses of components for tbp in [kg] ###')
    mass_data_jet.append(('MTOM_jet',MTOM_jet))
    mass_data_tbp.append(('MTOM_tbp',MTOM_tbp))
    mass_data_tbp.append(('M_empty_tbp', M_empty_tbp))
    mass_data_jet.append(('M_empty_jet', M_empty_jet))
    mass_data_tbp.append(('wing_weight_tbp', wing_weight_tbp))
    mass_data_jet.append(('wing_weight_jet', wing_weight_jet))
    mass_data_tbp.append(('hor_tail_weight_tbp', hor_tail_weight_tbp))
    mass_data_jet.append(('hor_tail_weight_jet', hor_tail_weight_jet))
    mass_data_tbp.append(('ver_tail_weight_tbp', ver_tail_weight_tbp))
    mass_data_jet.append(('ver_tail_weight_jet', ver_tail_weight_jet))
    mass_data_tbp.append(('fuselage_weight_tbp', fuselage_weight_tbp))
    mass_data_jet.append(('fuselage_weight_jet', fuselage_weight_jet))
    mass_data_tbp.append(('main_lg_weight_tbp', main_lg_weight_tbp))
    mass_data_jet.append(('main_lg_weigth_jet', main_lg_weight_jet))
    mass_data_tbp.append(('nose_lg_weight_tbp', nose_lg_weight_tbp))
    mass_data_jet.append(('nose_lg_weight_jet', nose_lg_weight_jet))
    mass_data_tbp.append(('nacelle_group_weight_tbp', nacelle_group_weight_tbp))
    mass_data_jet.append(('nacelle_group_weight_jet', nacelle_group_weight_jet))
    mass_data_tbp.append(('engine_controls_weight_tbp', engine_controls_weight_tbp))
    mass_data_tbp.append(('engine_controls_weight_jet', engine_controls_weight_jet))
    mass_data_tbp.append(('starter_weight_tbp', starter_weight_tbp))
    mass_data_jet.append(('starter_weight_jet', starter_weight_jet))
    mass_data_tbp.append(('W_fuel_system_tbp', W_fuel_system_tbp))
    mass_data_jet.append(('W_fuel_system_jet', W_fuel_system_jet))
    mass_data_tbp.append(('flight_controls_weight_tbp', flight_controls_weight_tbp))
    mass_data_jet.append(('flight_controls_weight_jet', flight_controls_weight_jet))
    mass_data_tbp.append(('instruments_weight_tbp', instruments_weight_tbp))
    mass_data_jet.append(('instruments_weight_jet', instruments_weight_jet))
    mass_data_tbp.append(('hydraulics_weight_tbp', hydraulics_weight_tbp))
    mass_data_jet.append(('hydraulics_weight_jet', hydraulics_weight_jet))
    mass_data_tbp.append(('electrical_weight_tbp', electrical_weight_tbp))
    mass_data_jet.append(('electrical_weight_jet', electrical_weight_jet))
    mass_data_tbp.append(('avionics_weight_tbp', avionics_weight_tbp))
    mass_data_jet.append(('avionics_weight_jet', avionics_weight_jet))
    mass_data_tbp.append(('furnishings_weight_tbp', furnishings_weight_tbp))
    mass_data_jet.append(('furnishings_weight_jet', furnishings_weight_jet))
    mass_data_tbp.append(('aircond_weight_tbp', aircond_weight_tbp))
    mass_data_jet.append(('aircond_weight_jet', aircond_weight_jet))
    mass_data_tbp.append(('anti_ice_weight_tbp', anti_ice_weight_tbp))
    mass_data_jet.append(('anti_ice_weight_jet', anti_ice_weight_jet))
    mass_data_tbp.append(('handling_gear_weight_tbp', handling_gear_weight_tbp))
    mass_data_jet.append(('handling_gear_weight_jet', handling_gear_weight_jet))

    print_list(mass_data_jet)
    print_list(mass_data_tbp)

def print_size_data():
    size_data_jet = []
    size_data_tbp = []

    size_data_jet.append('### Dimensions and areas for jet in [m] and [m^2] ###')
    size_data_tbp.append('### Dimensions and areas for tbp in [m] and [m^2] ###')
    size_data_jet.append(('L_fuselage_jet ', length_fuselage))
    size_data_tbp.append(('L_fuselage_tbp ', length_fuselage))
    size_data_tbp.append(('diameter_propeller_tbp ', diameter_propeller_tbp))
    size_data_jet.append(('diameter_highlight_jet ', diameter_highlight_jet))
    size_data_jet.append(('b_jet ', b_jet))
    size_data_tbp.append(('b_tbp ', b_tbp))
    size_data_jet.append(('S_jet', S_jet))
    size_data_tbp.append(('S_tbp', S_tbp))
    size_data_tbp.append(('S_h_tbp', S_h_tbp))
    size_data_jet.append(('S_h_jet', S_h_jet))
    size_data_tbp.append(('S_v_tbp', S_v_tbp))
    size_data_jet.append(('S_v_jet', S_v_jet))

    print_list(size_data_jet)
    print_list(size_data_tbp)

def print_flight_char_data():
    flight_char_data_jet = []
    flight_char_data_tbp = []

    flight_char_data_jet.append('### Flight char data for jet ###')
    flight_char_data_tbp.append('### Flight char data for tbp ###')
    flight_char_data_tbp.append(('C_l_des_tbp', C_l_des_tbp))
    flight_char_data_jet.append(('C_l_des_jet', C_l_des_jet))
    flight_char_data_jet.append(('fuel_per_passenger',fuel_per_passenger_jet))
    flight_char_data_tbp.append(('fuel_per_passenger',fuel_per_passenger_tbp))

    print_list(flight_char_data_jet)
    print_list(flight_char_data_tbp)

#Calculate performance for 1000 km trip
MTOW_jet_1000, OEW_jet_1000, W_fuel_jet_1000, C_D_0, f_cruise_start_jet, f_cruise_end_jet, LD_cruise_jet = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, V_cruise_jet, V_loiter_tbp, 1000*1000, 1000*1000, 2700, 2700, jet = True, tbp = False)
MTOW_tbp_1000, OEW_tbp_1000, W_fuel_tbp_1000, C_D_0, f_cruise_start_jet, f_cruise_end_jet, LD_cruise_jet = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, V_cruise_jet, V_loiter_tbp, 1000*1000, 1000*1000, 2700, 2700, tbp = True, jet = False)

fuel_per_passenger_jet_1000 = (W_fuel_jet_1000/n_passenger)/g
fuel_per_passenger_tbp_1000 = (W_fuel_tbp_1000/n_passenger)/g

CO2_tbp = CO2_calc(fuel_per_passenger_tbp_1000,chosen_fuel_energy_density)
CO2_jet = CO2_calc(fuel_per_passenger_jet_1000,chosen_fuel_energy_density)

NOX_tbp = NOX_calc(fuel_per_passenger_tbp_1000,chosen_fuel_energy_density)
NOX_jet = NOX_calc(fuel_per_passenger_jet_1000,chosen_fuel_energy_density)

range_cruise_jet_time = 1000000
range_cruise_tbp_time = 1000000
t_climb = altitude/c
d_horizontal_climb_jet = altitude/0.2
d_horizontal_climb_tbp = altitude/0.083
t_cruise_jet = (range_cruise_jet_time-d_horizontal_climb_jet)/V_cruise_jet
t_cruise_tbp = (range_cruise_tbp_time-d_horizontal_climb_tbp)/V_cruise_tbp
t_descent_jet = altitude/7.112 #descent of 1400 feet per minute
t_descent_tbp = altitude/7.112 #descent of 1400 feet per minute

t_jet = (t_climb+t_cruise_jet+t_descent_jet)/3600 #hours
t_tbp = (t_climb+t_cruise_tbp+t_descent_tbp)/3600 #hours

print('MTOM tbp: ' + str(round(MTOM_tbp,2)) + ' [kg]')
print('Fueltype: HHV')
print('Fuel per passenger per 1000 km tbp: ' + str(round(fuel_per_passenger_tbp_1000,2)) + ' [kg]')
print('CO2 per passanger per 1000 km tbp: ' + str(CO2_tbp) + ' [kg]')
print('NOX per passenger per 1000 km tbp: ' + str(NOX_tbp) + ' [kg]')
print('Time for a ' + str(range_cruise_tbp_time/1000) + 'km trip is ' + str(round(t_tbp,2)) + '[h]')
print()
print('MTOM propfan: ' + str(round(MTOM_jet,2)) + ' [kg]')
print('Fueltype: HHV')
print('Fuel per passenger per 1000 km propfan: ' + str(round(fuel_per_passenger_jet_1000,2)) + ' [kg]')
print('CO2 per passenger per 1000 km propfan: ' + str(CO2_jet) + ' [kg]')
print('NOX per passenger per 1000 km propfan: ' + str(NOX_jet) + ' [kg]')
print('Time for a ' + str(range_cruise_jet_time/1000) + 'km trip is ' + str(round(t_jet,2)) + '[h]')

SPL_engine_prop = prop_noise(diameter_propeller_tbp,5,1400,1850,1,speed_of_sound)
SPL_turbofan = turbofan_noise()
SPL_airframe = airframe_noise(V_cruise_jet,MTOW_jet)
SPL_total = total_noise(SPL_turbofan,SPL_airframe)
SPL_distance = noise_distance(SPL_total,450,2500)
print(SPL_turbofan,SPL_airframe,SPL_total,SPL_distance)
