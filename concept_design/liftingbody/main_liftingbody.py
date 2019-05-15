# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from constant_variables import *
from design_parameters import *
from class1_liftingbody import Weights_Class_I
from power_wingloading_liftingbody import wingloading_jet, wingloading_tbp
from wingloadingfunctions import T_W_calc, W_P_climb_calc
from class1sizing_liftingbody import fuselage, det_quarter_chord_sweep, det_planform, det_dihedral_angle, enginedimensions_jet, enginedimensions_tbp, MAC, empennage, undercarriage, tiresizing
from atmosphere import atmosphere_calc
from cg_determination import x_lemac_tbp_calc, x_lemac_jet_calc
from fuel_fraction import fuel_fraction
from conversion_formulas import *
import class2_liftingbody as class2
import numpy as np


# Atmospherical parameters at cruise altitude
temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, temperature0, temperature_gradient, g, R, gamma)
c = 12                              #[m/s] climb rate THIS IS INPUT

# Dynamic pressure
q_jet = 0.5*rho*V_cruise_jet**2     # [n/m2]
q_tbp = 0.5*rho*V_cruise_tbp**2     # [n/m2]


# Weight estimation and wing loading----------------------------------------
for iter in range(10):

    jet_data_list = []
    tbp_data_list = []
    ##################### CLASS I ############################
    # Weight estimations
    MTOW_jet, OEW_jet, W_fuel_jet, C_D_0_jet, f_cruise_start_jet, f_cruise_end_jet, LD_cruise_jet = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe_jet, C_fe_tbp, S_jet, S_tbp, S_wet_jet, S_wet_tbp, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel,  V_cruise_jet, V_loiter_tbp, range_cruise_jet, range_cruise_tbp, endurance_loiter_jet, endurance_loiter_tbp, jet = True, tbp = False)
    MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp, f_cruise_start_tbp, f_cruise_end_tbp, LD_cruise_tbp = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe_jet, C_fe_tbp, S_jet, S_tbp, S_wet_jet, S_wet_tbp, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel,  V_cruise_jet, V_loiter_tbp, range_cruise_jet, range_cruise_tbp, endurance_loiter_jet, endurance_loiter_tbp, jet = False, tbp = True)

    MTOM_jet = MTOW_jet/g
    MTOM_tbp = MTOW_tbp/g

    # Wing loading
    W_S_landing_jet = wingloading_jet(MTOW_jet,OEW_jet,V_cruise_jet,e_jet,C_D_0_jet,A_jet,S_jet,C_L_max_land_jet,C_L_max_TO_jet)
    W_S_landing_tbp, W_P_critical = wingloading_tbp(MTOW_tbp, OEW_tbp, S_tbp, A_tbp, V_cruise_tbp, e_tbp, eff_cruise_tbp, C_D_0_tbp, C_L_max_land_tbp,C_L_max_TO_tbp)

    # T/W for jet
    T_W_jet_range = T_W_calc(W_S_landing_jet,TOP_jet,1.32)

    ################### Sizing #################
    # Fuselage sizing
    length_nose, length_nosecone, length_cabin, length_tail, length_tailcone, length_fuselage, diameter_fuselage_inside, diameter_fuselage_outside, width_fuselage_outside = fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles)

    # Wing sizing
    # Required inputs
    # S for jet and tbp
    S_jet = MTOW_jet/W_S_landing_jet
    S_tbp = MTOW_tbp/W_S_landing_tbp

              # Depends on loading diagrams!
    S_wet = 4.5 * S_jet
    C_L_cruise = MTOW_jet / (0.5 * rho * V_cruise_jet**2 * S_jet)

    # Determine lift generated by the fuselage
    L_fuselage_jet = W_S_landing_jet*C_L_fuselage*(length_cabin*width_fuselage_outside)
    L_fuselage_tbp = W_S_landing_tbp*C_L_fuselage*(length_cabin*width_fuselage_outside)

    # Determine required surface area of the wing
    S_wing_jet = (MTOW_jet - L_fuselage_jet)/W_S_landing_jet
    S_wing_tbp = (MTOW_tbp - L_fuselage_tbp)/W_S_landing_tbp

    # Wing jet
    sweep_jet = det_quarter_chord_sweep(M_cruise_jet)
    b_jet, taper_jet, root_chord_jet, tip_chord_jet, t_c_ratio_jet = det_planform(S_wing_jet, A_jet, M_cruise_jet, C_L_cruise_jet, sweep_jet)
    dihedral_angle_jet = det_dihedral_angle(sweep_jet, low=True)
    MAC_jet = MAC(root_chord_jet, t_c_ratio_jet)

    # Wing tbp
    sweep_tbp = det_quarter_chord_sweep(M_cruise_tbp)
    b_tbp, taper_tbp, root_chord_tbp, tip_chord_tbp, t_c_ratio_tbp = det_planform(S_wing_tbp, A_tbp, M_cruise_tbp, C_L_cruise_tbp, sweep_tbp)
    dihedral_angle_tbp = det_dihedral_angle(sweep_tbp, high=True)
    MAC_tbp = MAC(root_chord_tbp, t_c_ratio_tbp)

    # Engine sizing
    T_TO_jet = T_W_jet_range * MTOW_jet              # Take-off thrust jet [N]
    P_TO_tbp = MTOW_tbp / W_P_critical                       # Take-off power tbp [W]
    length_nacelle_jet, length_fan_cowling_jet, diameter_highlight_jet, diameter_exit_fan_jet, diameter_gas_generator_jet, diameter_nacelle_jet = enginedimensions_jet(rho0, n_engines, T_TO_jet, jettypeC=True)
    diameter_engine_tbp, length_engine_tbp, diameter_propeller_tbp = enginedimensions_tbp(rho0,n_engines_tbp, P_TO_tbp)

    #CG determination
    x_lemac_jet = x_lemac_jet_calc(0.4*length_fuselage,MAC_jet)
    x_lemac_tbp = x_lemac_tbp_calc(0.4*length_fuselage,MAC_tbp)

    main_landing_pos_jet = x_lemac_jet+0.4*MAC_jet      # [m]
    main_landing_pos_tbp = x_lemac_tbp+0.4*MAC_tbp      # [m]

    # Empennage
    l_h_jet, l_v_jet = 0.9*length_fuselage-(x_lemac_jet+0.25*MAC_jet), 0.9*length_fuselage-(x_lemac_jet+0.25*MAC_jet)                           # [m]
    l_h_tbp, l_v_tbp = 0.9*length_fuselage-(x_lemac_tbp+0.25*MAC_tbp), 0.9*length_fuselage-(x_lemac_tbp+0.25*MAC_tbp)
    AR_h_jet, AR_v_jet, S_h_jet, span_h_jet, root_chord_h_jet, tip_chord_h_jet, sweepqc_h_jet, sweepLE_h_jet, S_v_jet, span_v_jet, root_chord_v_jet, tip_chord_v_jet, sweepLE_v_jet = empennage(V_h_jet, V_v_jet, l_h_jet, l_v_jet, S_jet, b_jet, MAC_jet)
    AR_h_tbp, AR_v_tbp, S_h_tbp, span_h_tbp, root_chord_h_tbp, tip_chord_h_tbp, sweepqc_h_tbp, sweepLE_h_tbp, S_v_tbp, span_v_tbp, root_chord_v_tbp, tip_chord_v_tbp, sweepLE_v_tbp = empennage(V_h_tbp, V_v_tbp, l_h_tbp, l_v_tbp, S_tbp, b_tbp, MAC_tbp)

    # Undercarriage sizing
    wheel_height_jet, lateral_position_jet = undercarriage(main_landing_pos_jet, nose_landing_pos_jet, length_fuselage, length_tail, diameter_fuselage_outside)
    wheel_height_tbp, lateral_position_tbp = undercarriage(main_landing_pos_tbp, nose_landing_pos_tbp, length_fuselage, length_tail, diameter_fuselage_outside)

    # Tire sizing
    tire_pressure_jet, P_mw_jet, P_nw_jet = tiresizing(MTOW_jet, 25)
    tire_pressure_tbp, P_mw_tbp, P_nw_tbo = tiresizing(MTOW_tbp, 25)


    ##################### CLASS II ############################
    # C_L and C_l des
    C_L_des_jet = class2.C_L_des(q_jet,f_cruise_start_jet*MTOW_jet/S_jet,f_cruise_end_jet*MTOW_jet/S_jet)
    C_l_des_jet = class2.C_l_des(C_L_des_jet,sweep_jet)
    C_L_des_tbp = class2.C_L_des(q_tbp,f_cruise_start_tbp*MTOW_tbp/S_tbp,f_cruise_end_tbp*MTOW_tbp/S_tbp)
    C_l_des_tbp = class2.C_l_des(C_L_des_tbp,sweep_tbp)

    # Wing weight
    n_max_jet = class2.ult_load_factor(kg_to_pounds(MTOM_jet))
    n_max_tbp = class2.ult_load_factor(kg_to_pounds(MTOM_tbp))

    qc_sweep_jet = det_quarter_chord_sweep(M_cruise_jet)
    qc_sweep_tbp = det_quarter_chord_sweep(M_cruise_tbp)

    wing_weight_jet = pounds_to_kg(class2.det_wing_weight(kg_to_pounds(MTOM_jet),(n_max_jet*1.5),metersquared_to_feetsquared(S_wing_jet),A_jet,t_c_ratio_jet,taper_jet,qc_sweep_jet,metersquared_to_feetsquared(0.05*S_wet_jet)))
    wing_weight_tbp = pounds_to_kg(class2.det_wing_weight(kg_to_pounds(MTOM_tbp),(n_max_tbp*1.5),metersquared_to_feetsquared(S_wing_tbp),A_tbp,t_c_ratio_tbp,taper_tbp,qc_sweep_tbp,metersquared_to_feetsquared(0.05*S_wet_tbp)))

    # Horizontal tail
    hor_tail_weight_jet = pounds_to_kg(class2.det_hor_tail_weight(meter_to_feet(width_fuselage_outside),meter_to_feet(span_h_jet),kg_to_pounds(MTOW_jet/g),n_max_jet*1.5,metersquared_to_feetsquared(S_h_jet),meter_to_feet(l_h_jet),np.radians(sweepqc_h_jet),AR_h_jet,metersquared_to_feetsquared(0.3*S_h_jet)))
    hor_tail_weight_tbp = pounds_to_kg(class2.det_hor_tail_weight(meter_to_feet(width_fuselage_outside),meter_to_feet(span_h_tbp),kg_to_pounds(MTOW_tbp/g),n_max_tbp*1.5,metersquared_to_feetsquared(S_h_tbp),meter_to_feet(l_h_tbp),np.radians(sweepqc_h_tbp),AR_h_tbp,metersquared_to_feetsquared(0.3*S_h_tbp)))

    # Vertical tail
    ver_tail_weight_jet = pounds_to_kg(class2.det_vert_tail_weight(span_v_jet,kg_to_pounds(MTOW_jet/g),n_max_jet*1.5,meter_to_feet(l_h_jet),metersquared_to_feetsquared(S_v_jet),np.radians(sweepLE_v_jet),AR_v_jet,t_c_ratio_jet))
    ver_tail_weight_tbp = pounds_to_kg(class2.det_vert_tail_weight(span_v_tbp,kg_to_pounds(MTOW_tbp/g),n_max_tbp*1.5,meter_to_feet(l_h_tbp),metersquared_to_feetsquared(S_v_tbp),np.radians(sweepLE_v_tbp),AR_v_tbp,t_c_ratio_tbp))

    # Fuselage
    fuselage_weight_jet = pounds_to_kg(class2.det_fuselage_weight(kg_to_pounds(MTOW_jet/g), 1.5*n_max_jet, meter_to_feet(length_fuselage), metersquared_to_feetsquared(length_fuselage*(diameter_fuselage_outside+2.4)*0.9), taper_jet, meter_to_feet(b_jet), sweep_jet, LD_cruise_jet, cargo_doors = 1, fuselage_mounted_lg = False))
    fuselage_weight_tbp = pounds_to_kg(class2.det_fuselage_weight(kg_to_pounds(MTOW_tbp/g), 1.5*n_max_tbp, meter_to_feet(length_fuselage), metersquared_to_feetsquared(length_fuselage*(diameter_fuselage_outside+2.4)*0.9), taper_tbp, meter_to_feet(b_tbp), sweep_tbp, LD_cruise_tbp, cargo_doors = 1, fuselage_mounted_lg = False))

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
    M_engine_jet = T_TO_jet / thrust_to_weight_jet + engine_gear_mass

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
    instruments_weight_jet = pounds_to_kg(class2.det_instruments_weight(2, 2, meter_to_feet(length_fuselage), meter_to_feet(b_jet), reciprocating = False, turboprop = False))
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
    furnishings_weight_jet = pounds_to_kg(class2.det_furnishings_weight(2, kg_to_pounds(n_passenger*23), metersquared_to_feetsquared(length_fuselage*(diameter_fuselage_outside+2.4)*0.9)))
    furnishings_weight_tbp = pounds_to_kg(class2.det_furnishings_weight(2, kg_to_pounds(n_passenger*23), metersquared_to_feetsquared(length_fuselage*(diameter_fuselage_outside+2.4)*0.9)))

    # Airconditioning
    pres_vol = 1.2 * np.pi / 4 * diameter_fuselage_inside**2 * (length_nose + length_nose)
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

#Calculate performance for 1000 km trip
MTOW_jet_1000, OEW_jet_1000, W_fuel_jet_1000, C_D_0, f_cruise_start_jet, f_cruise_end_jet, LD_cruise_jet = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe_jet, C_fe_tbp, S_wing_jet, S_wing_tbp, S_wet_jet, S_wet_tbp, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, V_cruise_jet, V_loiter_tbp, 1000*1000, 1000*1000, 2700, 2700, jet = True, tbp = False)
MTOW_tbp_1000, OEW_tbp_1000, W_fuel_tbp_1000, C_D_0, f_cruise_start_jet, f_cruise_end_jet, LD_cruise_jet = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe_jet, C_fe_tbp, S_wing_jet, S_wing_tbp, S_wet_jet, S_wet_tbp, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, V_cruise_jet, V_loiter_tbp, 1000*1000, 1000*1000, 2700, 2700, tbp = True, jet = False)

fuel_per_passenger_jet_1000 = (W_fuel_jet_1000/n_passenger)/g
fuel_per_passenger_tbp_1000 = (W_fuel_tbp_1000/n_passenger)/g


## PRINT RELEVANT DATA
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
    mass_data_jet.append(('engine_controls_weight_jet', engine_controls_weight_jet))
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

print('MTOM tbp: ' + str(MTOM_tbp))
print('MTOM jet: ' + str(MTOM_jet))


# Print data
print_mass_data()

print_size_data()

print_flight_char_data()
