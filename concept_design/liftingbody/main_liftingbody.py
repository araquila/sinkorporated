# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from class1_liftingbody import *
from class2_liftingbody import *
from class1sizing_liftingbody import *
from power_wingloading_liftingbody import *
from atmosphere import *
from conversion_formulas import *

# Decide whether you like jet or turboprop:
jet = True
tbp = False

# Constants
g = 9.8065
R = 287
gamma = 1.4

# Atmospherical parameters
temperature0 = 288.15
temperature_gradient = -0.0065
altitude = 25000 * 0.3048 # altitude in meters; this can be adjusted!
temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, temperature0, temperature_gradient, g, R, gamma)

# Passengers and crew
n_passenger = 60
M_passenger = 102                   #(including luggage)
n_crew = 4
M_crew_member = 100
n_seats_abreast = 6
n_aisles = 2

# Initial mass and fractions
M_payload = n_passenger * M_passenger
M_crew = n_crew * M_crew_member
f_trapped_fuel = 0.003              # Range 0.001-0.005
M_empty_tbp = 14827                 # Adjust per concept
M_empty_jet = 18783                 # Adjust per concept

# Convert to weights
W_payload = M_payload * g
W_crew = M_crew * g
W_empty_tbp = M_empty_tbp * g
W_empty_jet = M_empty_jet * g

# More general data
C_L_fuselage = 0.07
n_engines = 2
n_max_flap = 2
n_max_clean = 2.5
n_min = -1
n_max_man = 2.5
s_landing = 1400 #[m]
rho0 = 1.225 #[kg/m3]
V_landing = 48.93 #[m/s] maximum landing speed that is allowed on a runway of 1400 m this is set for all aircraft
#INPUTS!!!!!!!!!!!!###
rho = 0.4 #[kg/m3] altitude for cruise flight THIS IS INPUT
c = 10 #[m/s] climb rate THIS IS INPUT
weight_fraction = 0.8 #weight fraction of aircraft in cruise wrt MTOW THIS IS INPUT

#graph data
wing_loading_x = np.linspace(0.1,6000,200)

# Initial jet and tbp aircraft parameters
C_L = 0.4                            # during cruise

################################## Jet #######################################
C_fe_jet = 0.003
S_jet = 60                               # Adjust per concept
S_wet_jet = 4.5 * S_jet                  # Adjust per concept
A_jet = 12                          # Adjust per concept
e_jet = 0.85                        # Adjust per concept
cj_loiter_jet = 17/1e6         # [g/Ns]
cj_cruise_jet = 19/1e6         # [g/Ns]
Mach_cruise_jet = 0.8
W_S_jet = 3585

#Coefficients
C_L_max_jet_clean_min = 1.2
C_L_max_jet_clean_max = 1.8
C_L_max_jet_take_min = 1.6
C_L_max_jet_take_max = 2.2
C_L_max_jet_land_min = 2.4
C_L_max_jet_land_max = 2.8

#take off parameter jet
TOP_aquila_jet = 6698 #Take from statistics THIS IS INPUT
V_cruise_jet = Mach_cruise_jet * speed_of_sound #[m/s]
thrust_setting = 0.9 #THIS IS INPUT
C_D_jet_curr = 0.05 #current C_D value THIS IS INPUT
cV_jet = 0.20 #Climb gradient

################################# Tbp #########################################
C_fe_tbp = 0.003
S_tbp = 60                               # Adjust per concept
S_wet_tbp = 4 * S_tbp                  # Adjust per concept
A_tbp = 10                          # Adjust per concept
e_tbp = 0.8                        # Adjust per concept
eff_cruise_tbp = 0.85       # [-]
eff_loiter_tbp = 0.77       # [-]
cp_cruise_tbp = 75/1e9         #  [g/J]
cp_loiter_tbp = 90/1e9         #  [g/J]
Mach_cruise_tbp = 0.6
W_S_tbp = 2500

C_L_max_tbp_clean_min = 1.5
C_L_max_tbp_clean_max = 1.9
C_L_max_tbp_take_min = 1.7
C_L_max_tbp_take_max = 2.1
C_L_max_tbp_land_min = 1.9
C_L_max_tbp_land_max = 3.3

#take off parameter and propulsion
TOP_aquila_tbp = 139 #find from statistics THIS IS INPUT
power_setting = 0.9 #usually at 0.9 THIS IS INPUT
V_cruise_tbp = Mach_cruise_tbp * speed_of_sound #[m/s]
C_D_tbp_curr = 0.065 #current CD value THIS IS INPUT
eff_prop = 0.85 #THIS IS INPUT
cV_tbp = 0.083 #from CS23.65 climb gradient

# Weight estimation and wing loading----------------------------------------
for iter in range(1):
    if jet == True:
        ##################### CLASS I ############################
        # Weight estimations
        MTOW_jet, OEW_jet, W_fuel_jet, C_D_0_jet, f_cruise_start_jet, f_cruise_end_jet, LD_cruise_jet = Weights_Class_I_jet(W_empty_jet, W_payload, W_crew, C_fe_jet, S_jet, S_wet_jet, A_jet, e_jet, cj_loiter_jet, cj_cruise_jet, f_trapped_fuel)
        W_landing_jet = 0.98 * MTOW_jet
        #print("jet:", W_fuel_jet/MTOW_jet, MTOW_jet/g, OEW_jet/g, W_fuel_jet/g, C_D_0_jet)
        # Wing loading
        Wing_loading_jet(MTOW_jet, W_landing_jet, S_jet,  \
        C_L_max_jet_take_min, C_L_max_jet_take_max, C_L_max_jet_land_min, \
        C_L_max_jet_land_max, wing_loading_x,  \
        V_landing, V_cruise_jet, rho0, thrust_setting, weight_fraction,rho, TOP_aquila_jet, \
        C_D_0_jet, C_D_jet_curr, A_jet, e_jet, c, cV_jet, n_max_man)

        # Fuselage sizing
        length_nose, length_nosecone, length_cabin, length_tail, length_tailcone, length_fuselage, diameter_fuselage_outside, width_fuselage_outside = fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles)

        # Wing sizing
        # Required inputs
        A = A_jet
        S = MTOW_jet/W_S_jet                   # Depends on loading diagrams!
        Mach_cruise = Mach_cruise_jet
        C_L_cruise = MTOW_jet / (0.5 * rho * V_cruise_jet**2 * S)
        print("Cruise:", C_L_cruise, S, MTOW_jet)
        # Wing sweep
        sweep_chord_0_25 = det_quarter_chord_sweep(Mach_cruise, supercritical = False, delta_mach = 0.03)
        # Determine required surface area of the wing
        L_fuselage = W_S_jet*C_L_fuselage*(length_cabin*width_fuselage_outside)
        S_wing = (MTOW_jet - L_fuselage)/W_S_jet
        # Wing planform
        b, taper, root_chord, tip_chord, t_c_ratio = det_planform(S_wing, A, Mach_cruise, C_L_cruise, sweep_chord_0_25, supercritical = False, delta_mach = 0.03)
        # Wing dihedral - it requires input on wing position!
        dihedral = det_dihedral_angle(sweep_chord_0_25, low = True)
        MAC = S_wing/b
        print("Wing:", S_wing, b, taper, dihedral, t_c_ratio)
        # Engine sizing
        T_TO_jet = 0.296 * MTOW_jet          # Depends on loading diagrams!
        length_nacelle, length_f, diameter_highlight, diameter_exit_fan, diameter_gas_generator = enginedimensions_jet(rho0, n_engines, T_TO_jet, jettypeC=True)

        # Empennage sizing
        V_h = 0.7
        V_v = 0.04
        l_h = length_fuselage/2
        l_v = length_fuselage/2
        AR_h, AR_v, S_h, span_h, root_chord_h, tip_chord_h, sweepqc_h, sweepLE_h, S_v, span_v, root_chord_v, tip_chord_v, sweepLE_v = empennage(V_h, V_v, l_h, l_v, S_wing, b, MAC)
        print("Horizontal tail:", S_h, span_h, sweepqc_h)
        print("Vertical tail:", S_v, span_v)

        # Undercarriage sizing
        nose_landing_pos = 2
        main_landing_pos = 10
        wheel_height, lateral_position = undercarriage(main_landing_pos, nose_landing_pos, length_fuselage, length_tail, width_fuselage_outside)
        print("Wheels:", wheel_height, lateral_position)
        LCN = 25
        tire_pressure, P_mw, P_nw = tiresizing(MTOW_jet, LCN)
        print("Tires:", tire_pressure, P_mw, P_nw)

        ##################### CLASS II ############################

        n_max_man = ult_load_factor(MTOW_jet)

        # Determine design C_L and c_l
        q = 0.5 * rho * V_cruise_jet**2
        W_S_cruise_start = MTOW_jet * f_cruise_start_jet / S_wing
        W_S_cruise_end = MTOW_jet * f_cruise_end_jet / S_wing
        C_L_des = C_L_des(q, W_S_cruise_start, W_S_cruise_end)
        print("Design C_L:",C_L_des)

        C_l_des = C_l_des(C_L_des, sweep_chord_0_25)
        print("Design c_l", C_l_des)

        # Determine weight of the wing
        W_dg = kg_to_pounds(MTOW_jet/g)
        N_z = 1.5*n_max_man
        S_w = metersquared_to_feetsquared(S_wing)
        S_csw = 0.05*S_w
        wing_weight = g*pounds_to_kg(det_wing_weight(W_dg, N_z, S_w, A_jet, t_c_ratio, taper, sweep_chord_0_25, S_csw))
        print("Wing weight:", wing_weight)

        # Determine weight of the fuselage
        B_w = meter_to_feet(b)
        S_f = metersquared_to_feetsquared(length_fuselage*(diameter_fuselage_outside+2.4)*0.9)
        fuselage_weight = g*pounds_to_kg(det_fuselage_weight(W_dg, N_z, length_fuselage, S_f, taper, B_w, sweep_chord_0_25, LD_cruise_jet, cargo_doors = 1, fuselage_mounted_lg = False))
        print("Fuselage weight:", fuselage_weight)

        # Determine weight of the horizontal tail
        F_w = meter_to_feet(width_fuselage_outside)
        span_ht = meter_to_feet(span_h)
        S_ht = metersquared_to_feetsquared(S_h)
        L_t = meter_to_feet(10)
        quarter_chord_sweep_ht = radians(sweepqc_h)
        S_e = 0.33*S_ht
        hor_tail_weight = g*pounds_to_kg(det_hor_tail_weight(F_w, span_ht, W_dg, N_z, S_ht, L_t, quarter_chord_sweep_ht, AR_h, S_e, all_moving_unit = False, K_y = 0.3))
        print("Horizontal weight:", hor_tail_weight)

        # Determine weight of the vertical tail
        H_ht = meter_to_feet(span_v)
        H_vt = meter_to_feet(span_v)
        S_vt = metersquared_to_feetsquared(S_v)
        quarter_chord_sweep_vt = radians(sweepLE_v) - 0.1
        t_c_root_vt = 0.14
        vert_tail_weight = g*pounds_to_kg(det_vert_tail_weight(H_ht, H_vt, W_dg, N_z, L_t, S_vt, quarter_chord_sweep_vt, AR_v, t_c_root_vt, K_z = 1))
        print("Adjust line 222 for quarter chord sweep vertical tail")
        print("Vertical tail weight:", vert_tail_weight)

        # Determine weight of the landing gear
        W_l = kg_to_pounds(W_landing_jet)
        N_l = 4.5
        # Main landing gear
        L_m = meter_to_inch(wheel_height)
        N_mw = 4
        N_mss = 2
        V_stall = ms_to_knots(45)
        main_lg_weight = g*pounds_to_kg(det_main_lg_weight(W_l, N_l, L_m, N_mw, N_mss, V_stall, kneeling_main_lg = False))
        print("Adjust stall speed in line 234")
        print("Main landing gear weight:", main_lg_weight, "   Wait... This is weird...")

        # Nose landing gear
        L_n = meter_to_inch(wheel_height + 0.2)
        N_nw = 2
        nose_lg_weight = g*pounds_to_kg(det_nose_lg_weight(W_l, N_l, L_n, N_nw, kneeling_nose_lg = False))
        print("Nose landing gear weight:", nose_lg_weight)

    if tbp == True:
        # Weight estimation
        MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp, f_cruise_start_tbp, f_cruise_end_tbp, LD_cruise_tbp = Weights_Class_I_tbp(W_empty_tbp, W_payload, W_crew, C_fe_tbp, S_tbp, S_wet_tbp, A_tbp, e_tbp,  eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel)
        W_landing_tbp = 0.98 * MTOW_tbp
        print("turboprop:", W_fuel_tbp/MTOW_tbp, MTOW_tbp/g, OEW_tbp/g, W_fuel_tbp/g, C_D_0_tbp)
        # Wing loading
        Wing_loading_tbp(MTOW_tbp, W_landing_tbp, S_tbp,  \
        C_L_max_tbp_take_min, C_L_max_tbp_take_max, C_L_max_tbp_land_min, \
        C_L_max_tbp_land_max, wing_loading_x,  \
        V_landing, V_cruise_tbp, rho0, power_setting,weight_fraction,rho, TOP_aquila_tbp, \
        C_D_0_tbp, C_D_tbp_curr, A_tbp, eff_prop, e_tbp, c, cV_tbp, n_max_man)

        # Fuselage sizing
        length_nose, length_nosecone, length_cabin, length_tail, length_tailcone, length_fuselage, diameter_fuselage_outside, width_fuselage_outside = fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles)

        # Wing sizing
        A = A_tbp
        S = MTOW_tbp/W_S_tbp
        Mach_cruise = Mach_cruise_tbp
        C_L_cruise = MTOW_tbp / (0.5 * rho * V_cruise_tbp**2 * S)
        print(C_L_cruise)
        # Wing sweep
        sweep_chord_0_25 = det_quarter_chord_sweep(Mach_cruise, supercritical = False, delta_mach = 0.03)
        # Wing planform
        L_fuselage = W_S_tbp*C_L_fuselage*(length_cabin*width_fuselage_outside)
        S_wing = (MTOW_tbp - L_fuselage)/W_S_tbp
        # Wing planform
        b, taper, root_chord, tip_chord, t_c_ratio = det_planform(S_wing, A, Mach_cruise, C_L_cruise, sweep_chord_0_25, supercritical = False, delta_mach = 0.03)
        # Wing dihedral - it requires input on wing position!
        dihedral = det_dihedral_angle(sweep_chord_0_25, low = True)
        MAC = S_wing/b
        print("Wing:", length_cabin, length_fuselage, S_wing, b, taper, dihedral)

        # Engine sizing
        P_TO_tbp = MTOW_tbp/0.046           # adjust
        diameter_engine, length_engine, diameter_propeller = enginedimensions_tbp(rho0, n_engines, P_TO_tbp)

        # Empennage sizing
        V_h = 0.7
        V_v = 0.04
        l_h = length_fuselage/2
        l_v = length_fuselage/2
        AR_h, AR_v, S_h, span_h, root_chord_h, tip_chord_h, sweepqc_h, sweepLE_h, S_v, span_v, root_chord_v, tip_chord_v, sweepLE_v = empennage(V_h, V_v, l_h, l_v, S_wing, b, MAC)
        print("Horizontal tail:", S_h, span_h)
        print("Vertical tail:", S_v, span_v)

        # Undercarriage sizing
        nose_landing_pos = 2
        main_landing_pos = 10
        wheel_height, lateral_position = undercarriage(main_landing_pos, nose_landing_pos, length_fuselage, length_tail, width_fuselage_outside)
        print("Wheels:", wheel_height, lateral_position)
        LCN = 25
        tire_pressure, P_mw, P_nw = tiresizing(MTOW_tbp, LCN)
        print("Tires:", tire_pressure, P_mw, P_nw)
