# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from class1_conventional import Weights_Class_I
from power_wingloading_conventional import wingloading_jet, wingloading_tbp
from wingloadingfunctions import T_W_calc, W_P_climb_calc
from class1sizing_conventional import fuselage, det_quarter_chord_sweep, det_planform, det_dihedral_angle, enginedimensions, MAC, empennage, undercarriage, tiresizing
from atmosphere import atmosphere_calc
from cg_determination import x_lemac_tbp, x_lemac_jet
from fuel_fraction import fuel_fraction
from class2_conventional import C_L_des, C_l_des
import numpy as np

jet_data_list = []
tbp_data_list = []

## INPUTS AND CONSTANTS

# Gravitional constant
g = 9.80665
R = 287

# Atmospherical parameters
temperature0 = 288.15
temperature_gradient = -0.0065
gamma = 1.4
rho0 = 1.225                        #[kg/m3]

# Flight parameters
s_landing = 1400                    #[m]
altitude = 8000
V_landing = 48.93                   #[m/s] maximum landing speed that is allowed on a runway of 1400 m this is set for all aircraft

# Atmospherical parameters at cruise altitude
temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, temperature0, temperature_gradient, g, R, gamma)
c = 10                              #[m/s] climb rate THIS IS INPUT

# Passengers and crew
n_passenger = 60
M_passenger = 102                   #(including luggage)
n_crew = 4
M_crew_member = 100

# Cabin layout
n_seats_abreast = 4
n_aisles = 1

# Initial mass and fractions
M_payload = n_passenger * M_passenger
M_crew = n_crew * M_crew_member
f_trapped_fuel = 0.003              # Range 0.001-0.005
M_empty_tbp = 14400
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

# Other jet parameters
A_jet = 10
e_jet = 0.8                         # Adjust per concept
cj_loiter_jet = 19e-6               # (0.4-0.6) [g/j] Propfan: 0.441
cj_cruise_jet = 19e-6               # (0.5-0.9) [g/j] Propfan: 0.441
V_cruise_jet =  200                 # [m/s]
S_jet = 61
TOP_jet = 6698
M_cruise_jet = V_cruise_jet/speed_of_sound
C_L_cruise_jet = 0.4
n_engines_jet = 2

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
V_cruise_tbp = 150                  # [m/s]
M_cruise_tbp = V_cruise_tbp/speed_of_sound
C_L_cruise_tbp = 0.4
S_tbp = 76
TOP_tbp = 139
n_engines_tbp = 2

# Empennage tbp
V_h_tbp = 1.57                          # [-]
V_v_tbp = 0.07                          # [-]
nose_landing_pos_tbp = 3                # [m]


# Iterative sizing process
for iter in range(1):

    ## CLASS I
    MTOW_jet, OEW_jet, W_fuel_jet, C_D_0_jet, f_cruise_start_jet, f_cruise_end_jet = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, jet = True)
    MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp, f_cruise_start_tbp, f_cruise_end_tbp = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, tbp = True)

    jet_data_list.append(('MTOW_jet',MTOW_jet))
    tbp_data_list.append(('MTOW_tbp',MTOW_tbp))

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

    jet_data_list.append(('S_jet', S_jet))
    tbp_data_list.append(('S_tbp', S_tbp))

    ## SIZING
    # Fuselage
    length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside = fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles)
    jet_data_list.append(('L_fuselage_jet ', length_fuselage))
    tbp_data_list.append(('L_fuselage_tbp ', length_fuselage))

    # Wing tbp
    sweep_tbp = det_quarter_chord_sweep(M_cruise_tbp)
    b_tbp, taper_tbp, root_chord_tbp, tip_chord_tbp, t_c_ratio_tbp = det_planform(S_tbp, A_tbp, M_cruise_tbp, C_L_cruise_tbp, sweep_tbp)
    dihedral_angle_tbp = det_dihedral_angle(sweep_tbp, high=True)
    MAC_tbp = MAC(root_chord_tbp, t_c_ratio_tbp)
    tbp_data_list.append(('b_tbp ', b_tbp))

    # Wing jet
    sweep_jet = det_quarter_chord_sweep(M_cruise_jet)
    b_jet, taper_jet, root_chord_jet, tip_chord_jet, t_c_ratio_jet = det_planform(S_jet, A_jet, M_cruise_jet, C_L_cruise_jet, sweep_jet)
    dihedral_angle_jet = det_dihedral_angle(sweep_tbp, low=True)
    MAC_jet = MAC(root_chord_jet, t_c_ratio_jet)
    jet_data_list.append(('b_jet ', b_jet))

    # Engines for jet and tbp
    P_TO_tbp = MTOW_tbp / W_P_tbp                    # Take-off power tbp [W]
    T_TO_jet = T_W_jet_range[0] * MTOW_jet              # Take-off thrust jet [N]
    diameter_engine_tbp, length_engine_tbp, diameter_propeller_tbp = enginedimensions(rho0,n_engines_tbp, P_TO_tbp, T_TO_jet, tbp=True)
    length_nacelle_jet, length_fan_cowling_jet, diameter_highlight_jet, diameter_exit_fan_jet, diameter_gas_generator_jet = enginedimensions(rho0,n_engines_jet, P_TO_tbp, T_TO_jet, jettypeB=True)

    tbp_data_list.append(('diameter_propeller_tbp ', diameter_propeller_tbp))
    jet_data_list.append(('diameter_highlight_jet ', diameter_highlight_jet))

    #CG and undercarriage
    x_lemac_tbp = x_lemac_tbp(0.4*length_fuselage,MAC_tbp)
    x_lemac_jet = x_lemac_jet(0.4*length_fuselage,MAC_jet)

    l_h_jet, l_v_jet = 0.9*length_fuselage-(x_lemac_jet+0.25*MAC_jet), 0.9*length_fuselage-(x_lemac_jet+0.25*MAC_jet)                           # [m]
    l_h_tbp, l_v_tbp = 0.9*length_fuselage-(x_lemac_tbp+0.25*MAC_tbp), 0.9*length_fuselage-(x_lemac_tbp+0.25*MAC_tbp)

    main_landing_pos_jet = x_lemac_jet+0.4*MAC_jet               # [m]
    main_landing_pos_tbp = x_lemac_tbp+0.4*MAC_tbp             # [m]

    AR_h_jet, AR_v_jet, S_h_jet, span_h_jet, root_chord_h_jet, tip_chord_h_jet, sweepqc_h_jet, sweepLE_h_jet, S_v_jet, span_v_jet, root_chord_v_jet, tip_chord_v_jet, sweepLE_v_jet = empennage(V_h_jet, V_v_jet, l_h_jet, l_v_jet, S_jet, b_jet, MAC_jet)
    wheel_height_jet, lateral_position_jet = undercarriage(main_landing_pos_jet, nose_landing_pos_jet, length_fuselage, length_tail, diameter_fuselage_outside)
    AR_h_tbp, AR_v_tbp, S_h_tbp, span_h_tbp, root_chord_h_tbp, tip_chord_h_tbp, sweepqc_h_tbp, sweepLE_h_tbp, S_v_tbp, span_v_tbp, root_chord_v_tbp, tip_chord_v_tbp, sweepLE_v_tbp = empennage(V_h_tbp, V_v_tbp, l_h_tbp, l_v_tbp, S_tbp, b_tbp, MAC_tbp)
    wheel_height_tbp, lateral_position_tbp = undercarriage(main_landing_pos_tbp, nose_landing_pos_tbp, length_fuselage, length_tail, diameter_fuselage_outside)

    tbp_data_list.append(('S_h_tbp', S_h_tbp))
    jet_data_list.append(('S_h_jet', S_h_jet))
    tbp_data_list.append(('S_v_tbp', S_v_tbp))
    jet_data_list.append(('S_v_jet', S_v_jet))

    ## CLASS II
    q_jet = 0.5*rho*V_cruise_jet**2
    q_tbp = 0.5*rho*V_cruise_tbp**2
    C_L_des_jet = C_L_des(q_jet,f_cruise_start_jet*MTOW_jet/S_jet,f_cruise_end_jet*MTOW_jet/S_jet)
    C_l_des_jet = C_l_des(C_L_des_jet,sweep_jet)
    C_L_des_tbp = C_L_des(q_tbp,f_cruise_start_tbp*MTOW_tbp/S_tbp,f_cruise_end_tbp*MTOW_tbp/S_tbp)
    C_l_des_tbp = C_l_des(C_L_des_tbp,sweep_tbp)

    print(f_cruise_start_jet, f_cruise_start_tbp)
    ## PRINT RELEVANT DATA
    print('### JET VALUES ###')
    for value in jet_data_list:
        print(value[0] + ': ' + str(value[1]))
    print()
    print('### TBP VALUES ###')
    for value in tbp_data_list:
        print(value[0] + ': ' + str(value[1]))
