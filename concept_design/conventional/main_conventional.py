# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from constant_variables import *
from class1_conventional import Weights_Class_I
from power_wingloading_conventional import wingloading_jet, wingloading_tbp
from wingloadingfunctions import T_W_calc, W_P_climb_calc
from class1sizing_conventional import fuselage, det_quarter_chord_sweep, det_planform, det_dihedral_angle, enginedimensions, MAC, empennage, undercarriage, tiresizing
from atmosphere import atmosphere_calc
from cg_determination import x_lemac_tbp_calc, x_lemac_jet_calc
from fuel_fraction import fuel_fraction
from conversion_formulas import *
import class2_conventional as class2
import numpy as np

## INPUTS AND CONSTANTS



# Flight parameters
s_landing = 1400                    #[m]
altitude = 8000
V_landing = 48.93                   #[m/s] maximum landing speed that is allowed on a runway of 1400 m this is set for all aircraft

# Atmospherical parameters at cruise altitude
temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, temperature0, temperature_gradient, g, R, gamma)
c = 10                              #[m/s] climb rate THIS IS INPUT








# Initial jet and tbp aircraft parameters
C_fe = 0.003
S = 1
S_wet = 5 * S

# Other jet parameters
A_jet = 10
e_jet = 0.8                         # Adjust per concept
cj_loiter_jet = 19e-6               # (0.4-0.6) [g/j] Propfan: 0.441
cj_cruise_jet = 19e-6               # (0.5-0.9) [g/j] Propfan: 0.441
V_cruise_jet =  236.11                 # [m/s]
S_jet = 61
TOP_jet = 6698
M_cruise_jet = V_cruise_jet/speed_of_sound
C_L_cruise_jet = 0.4
n_engines_jet = 2
C_L_max_jet = 2.3

# Empennage jet
V_h_jet = 1.07                         # [-]
V_v_jet = 0.085                          # [-]
nose_landing_pos_jet = 3                # [m]

# Other tbp parameters
A_tbp = 12
e_tbp = 0.85                        # Adjust per concept
eff_cruise_tbp = 0.85               # [-]
eff_loiter_tbp = 0.77               # [-]
cp_cruise_tbp = 90e-9               # (0.4-0.6) [kg/ns]
cp_loiter_tbp = 90e-9               # (0.5-0.7) [kg/ns]
V_loiter_tbp = 80                   # [m/s]
V_cruise_tbp = 150                  # [m/s]
M_cruise_tbp = V_cruise_tbp/speed_of_sound
C_L_cruise_tbp = 0.4
S_tbp = 76
TOP_tbp = 139
n_engines_tbp = 2
C_L_max_tbp = 2.6

# Empennage tbp
V_h_tbp = 1.57                          # [-]
V_v_tbp = 0.07                          # [-]
nose_landing_pos_tbp = 3                # [m]

# Dynamic pressure
q_jet = 0.5*rho*V_cruise_jet**2     # [n/m2]
q_tbp = 0.5*rho*V_cruise_tbp**2     # [n/m2]

# Iterative sizing process
for iter in range(1):
    
    jet_data_list = []
    tbp_data_list = []
    
    ## CLASS I
    MTOW_jet, OEW_jet, W_fuel_jet, C_D_0_jet, f_cruise_start_jet, f_cruise_end_jet, L_D_jet = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, V_cruise_jet, V_loiter_tbp, jet = True, tbp = False)
    MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp, f_cruise_start_tbp, f_cruise_end_tbp, L_D_tbp = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, V_cruise_jet, V_loiter_tbp, tbp = True, jet = False)

    jet_data_list.append(('MTOW_jet',MTOW_jet))
    tbp_data_list.append(('MTOW_tbp',MTOW_tbp))

    MTOM_jet = MTOW_jet/g
    MTOM_tbp = MTOW_jet/g

    ## WING LOADING AND POWER LOADING
    W_S_landing_jet = wingloading_jet(MTOW_jet,OEW_jet,V_cruise_jet,e_jet,C_D_0_jet,A_jet,S_jet)
    W_S_landing_tbp = wingloading_tbp(MTOW_tbp, OEW_tbp, S_tbp, A_tbp, V_cruise_tbp, e_tbp, eff_cruise_tbp, C_D_0_tbp)

    # T/W for jet
    T_W_jet_range = np.zeros(len(W_S_landing_jet))
    for i in range(len(W_S_landing_jet)):
        T_W_jet_range[i] = T_W_calc(W_S_landing_jet[i],TOP_jet,1.32)

    # W/P for tbp
    W_P_tbp = W_P_climb_calc(eff_cruise_tbp,c,W_S_landing_tbp[0][0],rho,A_tbp,e_tbp,C_D_0_tbp)

    # S for jet and tbp
    S_jet = MTOW_jet/W_S_landing_jet[0]
    S_tbp = MTOW_tbp/W_S_landing_tbp[0][0]

    # Append to data list
    jet_data_list.append(('S_jet', S_jet))
    tbp_data_list.append(('S_tbp', S_tbp))

    ## SIZING
    # Fuselage
    length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside, diameter_fuselage_inside = fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles)

    # Append to data list
    jet_data_list.append(('L_fuselage_jet ', length_fuselage))
    tbp_data_list.append(('L_fuselage_tbp ', length_fuselage))

    # Wing tbp
    sweep_tbp = det_quarter_chord_sweep(M_cruise_tbp)
    b_tbp, taper_tbp, root_chord_tbp, tip_chord_tbp, t_c_ratio_tbp = det_planform(S_tbp, A_tbp, M_cruise_tbp, C_L_cruise_tbp, sweep_tbp)
    dihedral_angle_tbp = det_dihedral_angle(sweep_tbp, high=True)
    MAC_tbp = MAC(root_chord_tbp, t_c_ratio_tbp)

    # Append to data list
    tbp_data_list.append(('b_tbp ', b_tbp))

    # Wing jet
    sweep_jet = det_quarter_chord_sweep(M_cruise_jet)
    b_jet, taper_jet, root_chord_jet, tip_chord_jet, t_c_ratio_jet = det_planform(S_jet, A_jet, M_cruise_jet, C_L_cruise_jet, sweep_jet)
    dihedral_angle_jet = det_dihedral_angle(sweep_tbp, low=True)
    MAC_jet = MAC(root_chord_jet, t_c_ratio_jet)

    # Append to data list
    jet_data_list.append(('b_jet ', b_jet))

    # Engines for jet and tbp
    P_TO_tbp = MTOW_tbp / W_P_tbp                       # Take-off power tbp [W]
    T_TO_jet = T_W_jet_range[0] * MTOW_jet              # Take-off thrust jet [N]
    diameter_engine_tbp, length_engine_tbp, diameter_propeller_tbp = enginedimensions(rho0,n_engines_tbp, P_TO_tbp, T_TO_jet, tbp=True)
    length_nacelle_jet, length_fan_cowling_jet, diameter_highlight_jet, diameter_exit_fan_jet, diameter_gas_generator_jet, diameter_nacelle_jet = enginedimensions(rho0,n_engines_jet, P_TO_tbp, T_TO_jet, jettypeB=True)

    # Append to data list
    tbp_data_list.append(('diameter_propeller_tbp ', diameter_propeller_tbp))
    jet_data_list.append(('diameter_highlight_jet ', diameter_highlight_jet))

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

    # Append to data list
    tbp_data_list.append(('S_h_tbp', S_h_tbp))
    jet_data_list.append(('S_h_jet', S_h_jet))
    tbp_data_list.append(('S_v_tbp', S_v_tbp))
    jet_data_list.append(('S_v_jet', S_v_jet))

    ## CLASS II
    # C_L and C_l des
    C_L_des_jet = class2.C_L_des(q_jet,f_cruise_start_jet*MTOW_jet/S_jet,f_cruise_end_jet*MTOW_jet/S_jet)
    C_l_des_jet = class2.C_l_des(C_L_des_jet,sweep_jet)
    C_L_des_tbp = class2.C_L_des(q_tbp,f_cruise_start_tbp*MTOW_tbp/S_tbp,f_cruise_end_tbp*MTOW_tbp/S_tbp)
    C_l_des_tbp = class2.C_l_des(C_L_des_tbp,sweep_tbp)

    # Append to data list
    tbp_data_list.append(('C_l_des_tbp', C_l_des_tbp))
    jet_data_list.append(('C_l_des_jet', C_l_des_jet))

    # Wing weight
    n_max_jet = class2.ult_load_factor(kg_to_pounds(MTOW_jet/g))
    n_max_tbp = class2.ult_load_factor(kg_to_pounds(MTOW_tbp/g))
    
    qc_sweep_tbp = det_quarter_chord_sweep(M_cruise_tbp)
    qc_sweep_jet = det_quarter_chord_sweep(M_cruise_jet)

    wing_weight_jet = pounds_to_kg(class2.det_wing_weight(kg_to_pounds(MTOW_jet/g),n_max_jet*1.5,metersquared_to_feetsquared(S_jet),A_jet,t_c_ratio_jet,taper_jet,qc_sweep_jet,metersquared_to_feetsquared(0.05*S_jet))) #DIT NOG FF CHECKEN
    wing_weight_tbp = pounds_to_kg(class2.det_wing_weight(kg_to_pounds(MTOW_tbp/g),n_max_tbp*1.5,metersquared_to_feetsquared(S_tbp),A_tbp,t_c_ratio_tbp,taper_tbp,qc_sweep_tbp,metersquared_to_feetsquared(0.05*S_tbp))) #DIT NOG FF CHECKEN

    # Append to data list
    tbp_data_list.append(('wing_weight_tbp', wing_weight_tbp))
    jet_data_list.append(('wing_weight_jet', wing_weight_jet))

    # Horizontal tail
    hor_tail_weight_jet = pounds_to_kg(class2.det_hor_tail_weight(meter_to_feet(diameter_fuselage_outside),meter_to_feet(span_h_jet),kg_to_pounds(MTOW_jet/g),n_max_jet*1.5,metersquared_to_feetsquared(S_h_jet),meter_to_feet(l_h_jet),np.radians(sweepqc_h_jet),AR_h_jet,metersquared_to_feetsquared(0.3*S_h_jet)))
    hor_tail_weight_tbp = pounds_to_kg(class2.det_hor_tail_weight(meter_to_feet(diameter_fuselage_outside),meter_to_feet(span_h_tbp),kg_to_pounds(MTOW_tbp/g),n_max_tbp*1.5,metersquared_to_feetsquared(S_h_tbp),meter_to_feet(l_h_tbp),np.radians(sweepqc_h_tbp),AR_h_tbp,metersquared_to_feetsquared(0.3*S_h_tbp)))

    # Append to data list
    tbp_data_list.append(('hor_tail_weight_tbp', hor_tail_weight_tbp))
    jet_data_list.append(('hor_tail_weight_jet', hor_tail_weight_jet))

    # Vertical tail

    ver_tail_weight_jet = pounds_to_kg(class2.det_vert_tail_weight(0,1,kg_to_pounds(MTOW_jet/g),n_max_jet*1.5,meter_to_feet(l_h_jet),metersquared_to_feetsquared(S_v_jet),np.radians(sweepLE_v_jet),AR_v_jet,t_c_ratio_jet))
    ver_tail_weight_tbp = pounds_to_kg(class2.det_vert_tail_weight(0,1,kg_to_pounds(MTOW_tbp/g),n_max_tbp*1.5,meter_to_feet(l_h_tbp),metersquared_to_feetsquared(S_v_tbp),np.radians(sweepLE_v_tbp),AR_v_tbp,t_c_ratio_tbp))

    # Append to data list
    tbp_data_list.append(('ver_tail_weight_tbp', ver_tail_weight_tbp))
    jet_data_list.append(('ver_tail_weight_jet', ver_tail_weight_jet))

    # Fuselage
    fuselage_weight_jet = pounds_to_kg(class2.det_fuselage_weight(kg_to_pounds(MTOW_jet/g), 1.5*n_max_jet, meter_to_feet(length_fuselage), metersquared_to_feetsquared(np.pi*diameter_fuselage_outside*length_fuselage), taper_jet, meter_to_feet(b_jet), sweep_jet, L_D_jet, cargo_doors = 1, fuselage_mounted_lg = False))
    fuselage_weight_tbp = pounds_to_kg(class2.det_fuselage_weight(kg_to_pounds(MTOW_tbp/g), 1.5*n_max_tbp, meter_to_feet(length_fuselage), metersquared_to_feetsquared(np.pi*diameter_fuselage_outside*length_fuselage), taper_tbp, meter_to_feet(b_tbp), sweep_tbp, L_D_tbp, cargo_doors = 1, fuselage_mounted_lg = False))

    # Append to data list
    tbp_data_list.append(('fuselage_weight_tbp', fuselage_weight_tbp))
    jet_data_list.append(('fuselage_weight_jet', fuselage_weight_jet))

    # Main landing gear
    V_stall_jet = np.sqrt(2*MTOW_jet/(rho0*S_jet*C_L_max_jet))
    V_stall_tbp = np.sqrt(2*MTOW_tbp/(rho0*S_tbp*C_L_max_tbp))
    main_lg_weight_jet = pounds_to_kg(class2.det_main_lg_weight(kg_to_pounds(MTOW_jet/g),1.5*n_max_jet,meter_to_inch(wheel_height_jet),4,2,ms_to_knots(V_stall_jet),kneeling_main_lg = False))
    main_lg_weight_tbp = pounds_to_kg(class2.det_main_lg_weight(kg_to_pounds(MTOW_tbp/g),1.5*n_max_tbp,meter_to_inch(wheel_height_tbp),4,2,ms_to_knots(V_stall_tbp),kneeling_main_lg = False))

    # Append to data list
    tbp_data_list.append(('main_lg_weight_tbp', main_lg_weight_tbp))
    jet_data_list.append(('main_lg_weigth_jet', main_lg_weight_jet))

    #Nose landing gear
    nose_lg_weight_jet = pounds_to_kg(class2.det_nose_lg_weight(kg_to_pounds(MTOW_jet/g), 1.5*n_max_jet, meter_to_inch(wheel_height_jet), 2, kneeling_nose_lg = False))
    nose_lg_weight_tbp = pounds_to_kg(class2.det_nose_lg_weight(kg_to_pounds(MTOW_tbp/g), 1.5*n_max_tbp, meter_to_inch(wheel_height_tbp), 2, kneeling_nose_lg = False))

    # Append to data list
    tbp_data_list.append(('nose_lg_weight_tbp', nose_lg_weight_tbp))
    jet_data_list.append(('nose_lg_weight_jet', nose_lg_weight_jet))

    # Nacelle
    nacelle_group_weight_jet = pounds_to_kg(class2.det_nacelle_group_weight(meter_to_feet(length_nacelle_jet), meter_to_feet(diameter_nacelle_jet), 1.5*n_max_jet, 2, metersquared_to_feetsquared(np.pi * diameter_nacelle_jet * length_nacelle_jet), pylon_mounted = True, W_ec = 0, W_engine = kg_to_pounds(820), propeller = False, thrust_reverser = False))
    nacelle_group_weight_tbp = pounds_to_kg(class2.det_nacelle_group_weight(meter_to_feet(length_engine_tbp), meter_to_feet(diameter_engine_tbp), 1.5*n_max_tbp, 2, metersquared_to_feetsquared(np.pi * diameter_engine_tbp * length_engine_tbp), pylon_mounted = True, W_ec = 0, W_engine = kg_to_pounds(700), propeller = True, thrust_reverser = False))

    #Append to data list
    tbp_data_list.append(('nacelle_group_weight_tbp', nacelle_group_weight_tbp))
    jet_data_list.append(('nacelle_group_weight_jet', nacelle_group_weight_jet))


    # Engine controls weight
    engine_controls_weight_jet = pounds_to_kg(class2.det_engine_controls_weight(2,40))
    engine_controls_weight_tbp = pounds_to_kg(class2.det_engine_controls_weight(2,40))

    # Append to data list
    tbp_data_list.append(('engine_controls_weight_tbp', engine_controls_weight_tbp))
    jet_data_list.append(('engine_controls_weight_jet', engine_controls_weight_jet))

    # Starter weight
    starter_weight_jet = pounds_to_kg(class2.det_starter_weight(2, kg_to_pounds(800)))
    starter_weight_tbp = pounds_to_kg(class2.det_starter_weight(2, kg_to_pounds(600)))

    # Append to data list
    tbp_data_list.append(('starter_weight_tbp', starter_weight_tbp))
    jet_data_list.append(('starter_weight_jet', starter_weight_jet))

    # Fuel system
    flight_controls_weight_jet = pounds_to_kg(class2.det_flight_controls_weight(metersquared_to_feetsquared(0.3*S_h_jet+0.05*S_jet), (meter_to_feet(length_fuselage)**2*kg_to_pounds(MTOM_jet)*0.34**2)/(4*32.19), N_f = 6, N_m = 1))
    flight_controls_weight_tbp = pounds_to_kg(class2.det_flight_controls_weight(metersquared_to_feetsquared(0.3*S_h_tbp+0.05*S_tbp), (meter_to_feet(length_fuselage)**2*kg_to_pounds(MTOM_tbp)*0.34**2)/(4*32.19), N_f = 6, N_m = 1))

    # Append to data list
    tbp_data_list.append(('flight_controls_weight_tbp', flight_controls_weight_tbp))
    jet_data_list.append(('flight_controls_weight_jet', flight_controls_weight_jet))

    # Instruments
    instruments_weight_jet = pounds_to_kg(class2.det_instruments_weight(2, 2, meter_to_feet(length_fuselage), meter_to_feet(b_jet), reciprocating = False, turboprop = True))
    instruments_weight_tbp = pounds_to_kg(class2.det_instruments_weight(2, 2, meter_to_feet(length_fuselage), meter_to_feet(b_tbp), reciprocating = False, turboprop = True))

    # Append to data list
    tbp_data_list.append(('instruments_weight_tbp', instruments_weight_tbp))
    jet_data_list.append(('instruments_weight_jet', instruments_weight_jet))

    # Hydraulics

    hydraulics_weight_jet = pounds_to_kg(class2.det_hydraulics_weight(meter_to_feet(length_fuselage), meter_to_feet(b_jet), N_f = 6))
    hydraulics_weight_tbp = pounds_to_kg(class2.det_hydraulics_weight(meter_to_feet(length_fuselage), meter_to_feet(b_tbp), N_f = 6))

    # Append to data list
    tbp_data_list.append(('hydraulics_weight_tbp', hydraulics_weight_tbp))
    jet_data_list.append(('hydraulics_weight_jet', hydraulics_weight_jet))

    # Electrical
    electrical_weight_jet = pounds_to_kg(class2.det_electrical_weight(200, R_kva = 50, N_gen = 0, N_en = 2)) #plug in better numbers
    electrical_weight_tbp = pounds_to_kg(class2.det_electrical_weight(200, R_kva = 50, N_gen = 0, N_en = 2)) #plug in better numbers

    # Append to data list
    tbp_data_list.append(('electrical_weight_tbp', electrical_weight_tbp))
    jet_data_list.append(('electrical_weight_jet', electrical_weight_jet))

    # Avionics

    avionics_weight_jet = pounds_to_kg(class2.det_avionics_weight(W_uav = 1100))
    avionics_weight_tbp = pounds_to_kg(class2.det_avionics_weight(W_uav = 1100))

    # Append to data list
    tbp_data_list.append(('avionics_weight_tbp', avionics_weight_tbp))
    jet_data_list.append(('avionics_weight_jet', avionics_weight_jet))

    # Furnishings
    furnishings_weight_jet = pounds_to_kg(class2.det_furnishings_weight(2, kg_to_pounds(13.608*60), metersquared_to_feetsquared(np.pi * diameter_fuselage_outside * length_fuselage)))
    furnishings_weight_tbp = pounds_to_kg(class2.det_furnishings_weight(2, kg_to_pounds(13.608*60), metersquared_to_feetsquared(np.pi * diameter_fuselage_outside * length_fuselage)))

    # Append to data list
    tbp_data_list.append(('furnishings_weight_tbp', furnishings_weight_tbp))
    jet_data_list.append(('furnishings_weight_jet', furnishings_weight_jet))

    # Airconditioning
    pres_vol = np.pi / 4 * diameter_fuselage_inside**2 * (length_nose + length_nose)
    aircond_weight_jet = pounds_to_kg(class2.det_aircond_weight(n_passenger+n_crew, metercubed_to_feetcubed(pres_vol), W_uav = 1100))
    aircond_weight_tbp = pounds_to_kg(class2.det_aircond_weight(n_passenger+n_crew, metercubed_to_feetcubed(pres_vol), W_uav = 1100))

    # Append to data list
    tbp_data_list.append(('aircond_weight_tbp', aircond_weight_tbp))
    jet_data_list.append(('aircond_weight_jet', aircond_weight_jet))

    # Anti-ice
    anti_ice_weight_jet = pounds_to_kg(class2.det_anti_ice_weight(kg_to_pounds(MTOW_jet)))
    anti_ice_weight_tbp = pounds_to_kg(class2.det_anti_ice_weight(kg_to_pounds(MTOW_tbp)))

    # Append to data list
    tbp_data_list.append(('anti_ice_weight_tbp', anti_ice_weight_tbp))
    jet_data_list.append(('anti_ice_weight_jet', anti_ice_weight_jet))

    # Handling gear
    handling_gear_weight_jet = pounds_to_kg(class2.det_handling_gear_weight(kg_to_pounds(MTOW_jet)))
    handling_gear_weight_tbp = pounds_to_kg(class2.det_handling_gear_weight(kg_to_pounds(MTOW_tbp)))

    # Append to data list
    tbp_data_list.append(('handling_gear_weight_tbp', handling_gear_weight_tbp))
    jet_data_list.append(('handling_gear_weight_jet', handling_gear_weight_jet))

    # Total weight
    W_empty_jet = wing_weight_jet + hor_tail_weight_jet + ver_tail_weight_jet + fuselage_weight_jet + main_lg_weight_jet + nose_lg_weight_jet + engine_controls_weight_jet + starter_weight_jet + flight_controls_weight_jet + instruments_weight_jet + hydraulics_weight_jet + electrical_weight_jet + avionics_weight_jet + furnishings_weight_jet + aircond_weight_jet + anti_ice_weight_jet + handling_gear_weight_jet
    W_empty_tbp = wing_weight_tbp + hor_tail_weight_tbp + ver_tail_weight_tbp + fuselage_weight_tbp + main_lg_weight_tbp + nose_lg_weight_tbp + nacelle_group_weight_tbp + engine_controls_weight_tbp + starter_weight_tbp + flight_controls_weight_tbp + instruments_weight_tbp + hydraulics_weight_tbp + electrical_weight_tbp + avionics_weight_tbp + furnishings_weight_jet + aircond_weight_tbp + anti_ice_weight_tbp + handling_gear_weight_tbp

    # Append to data list
    tbp_data_list.append(('W_empty_jet', W_empty_jet))
    jet_data_list.append(('W_empty_tbp', W_empty_tbp))

    ## PRINT RELEVANT DATA
    print('Iteration: ' + str(iter+1))
    print()
    print('### JET VALUES ###')
    for value in jet_data_list:
        print(value[0] + ': ' + str(value[1]))
    print()
    print('### TBP VALUES ###')
    for value in tbp_data_list:
        print(value[0] + ': ' + str(value[1]))
    print('----------------------------------------------------')
