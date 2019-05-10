# Class 2 weight estimation
# Equations by Raymer 1989

def ult_load_factor(MTOW):
    if MTOW <= 4100:
        n_max = 3.8
    if 4100 < MTOW <= 50000:
        n_max = 2.1 + (24000/(MTOW + 10000))
    if MTOW > 50000:
        n_max = 2.5
    return n_max

def det_wing_weight(W_dg, N_z, S_w, AR, t_c_root, taper, quarter_chord_sweep, S_csw):
    """
    Inputs:
    W_dg = design gross weigth in lb
    N_z = ultimate load factor = 1.5 * limit load factor
    S_w = total wing area in ft^2
    AR = aspect ratio
    t_c_root = thickness to chord ratio at the root
    taper = taper ratio
    quarter_chord_sweep = sweep angle at the quarter chord line in radians
    S_csw = surface area of wing mounted control surfaces

    output:
    Wing weight in lb
    """

    wing_weight = 0.0051 (W_dg * N_z)**0.557 * S_w**0.649 * AR**0.5 * t_c_root^-0.4 * (1 + taper)^0.1 * np.cos(quarter_chord_sweep)**-1 * S_csw**0.1
    return wing_weight

def det_hor_tail_weight(F_w, B_h, W_dg, N_z, S_ht, L_t, quarter_chord_sweep_ht, AR_ht, S_e, all_moving_unit = False, K_y = 0.3):
    """
    Inputs:
    F_w = fuselage width at horizontal tail intersection in ft
    B_h = horizontal tail span
    W_dg = design gross weigth in lb
    N_z = ultimate load factor = 1.5 * limit load factor
    S_ht = horizontal tail area in ft^2
    L_t = tail length (wing quarter MAC to tail quarter MAC) in feet
    quarter_chord_sweep_ht = sweep angle of the horizontal tail at the quarter chord line in radians
    AR_ht = aspect ratio of the horizontal tail
    S_e = elevator area in ft^2

    conditional inputs:
    all_moving_unit = is the horizontal tail moving as a unit (True) or only part of it (False)
    all_moving_unit = True -> 1.143
    all_moving_unit = False -> 1
    K_y = Aircraft pitching radius of gyration in ft (usually 0.3 * L_t)

    output:
    Horizontal tail weight in lb
    """
    K_uht = 1
    if all_moving_unit:
        K_uht = 1.143
    if K_y == 0.3:
        K_y = 0.3 * L_t
    hor_tail_weight = 0.0379 * K_uht * (1 + F_w/B_h)**-0.25 * W_dg**0.639 * N_z**0.10 * S_ht**0.75 * L_t^-1 * K_y**0.704 * np.cos(quarter_chord_sweep_ht)**-1 * A_h**0.166 * (1 + S_e/S_ht)**0.1
    return hor_tail_weight

def det_vert_tail_weight(H_ht, H_vt, W_dg, N_z, L_t, S_vt, quarter_chord_sweep_vt, AR_vt, t_c_root_vt, K_z = 1):
    """
    Inputs:
    H_ht = horizontal tail height above fuselage in ft
    H_vt = vertical tail height above fuselage in ft
    H_ht/H_vt is the ratio of the heigts of the empenage (0 for conventional and 1 for T-tail)
    W_dg = design gross weigth in lb
    N_z = ultimate load factor = 1.5 * limit load factor
    L_t = tail length (wing quarter MAC to tail quarter MAC) in feet
    S_vt = vertical tail area in ft^2
    quarter_chord_sweep_vt = sweep angle of the vertical tail at the quarter chord line in radians
    AR_vt = aspect ratio of the vertical tail
    t_c_root_vt = thickness to chord ratio at the root of the vertical tail

    conditional inputs:
    K_z = Aircraft yawing radius of gyration in ft (usually L_t)

    outputs:
    vertical tail weight in lb
    """
    if K_y == 1:
        K_y = L_t
    vert_tail_weight = 0.0026 * (1 + H_ht/H_vt)**0.225 * W_dg**0.556 * N_z**0.536 * L_t**-0.5 * S_vt**0.5 * K_z**0.875 * np.cos(quarter_chord_sweep_vt)**-1 * AR_vt**0.35 * t_c_root_vt**-0.5
    return vert_tail_weight

def det_fuselage_weight(W_dg, N_z, L, S_f, taper, B_w, quarter_chord_sweep, L_D_ratio, cargo_doors = 1, fuselage_mounted_lg = False):
    """
    Inputs:
    W_dg = design gross weigth in lb
    N_z = ultimate load factor = 1.5 * limit load factor
    L = fuselage structural length in feet (no radome/nosecone and tail cap)
    S_f = fuselage wetted area in ft^2
    taper = taper ratio
    B_w = wing span in feet
    quarter_chord_sweep = sweep angle at the quarter chord line in radians
    L_D_ratio = the lift over drag ratio of the aircraft

    conditional inputs:
    K_door = factor for the cargo doors (1 for no cargo door, 1.06 for one side cargo door, 1.12 for 2 side cargo doors)
    K_lg = factor for landing gear (1.12 for fuselage mounted landing gear, 1 otherwise)

    outputs:
    fuselage weight in lb
    """
    K_door = 1 + cargo_doors * 0.06
    K_lg = 1
    if fuselage_mounted_lg:
        K_lg = 1.12
    K_ws = 0.75 * (1 + 2 * taper)/(1 + taper) * (B_w * np.tan(quarter_chord_sweep)/L)
    fuselage_weight = 0.3280 * K_door * K_lg * (W_dg * N_z)**0.5 * L**0.25 * S_f**0.302 * (1 + K_ws)**0.04 * L_D_ratio**0.10
    return fuselage_weight

def det_main_lg_weight(W_l, N_l, L_m, N_mw, N_mss, V_stall, kneeling_main_lg = False):
    """
    Inputs:
    W_l = landing design weight in lb
    N_l = ultimate landing load factor = N_gear * 1.5
    L_m = length of the main landing gear in inches
    N_mw = number of main wheels
    N_mss = number of main gear shock struts
    V_stall = stall speed in knots

    conditional inputs:
    K_mp = factor for landing gear (1.126 for kneeling landing gear, 1 otherwise)

    outputs:
    main landing gear weight in lb
    """
    K_mp = 1
    if kneeling_main_lg:
        K_mp = 1.126
    main_lg_weight = 0.0106 * K_mp * W_l**0.888 * N_l**0.25 * L_m**0.4 * N_mw**0.321 * N_mss**-0.5 * V_stall**0.1
    return main_lg_weight

def det_nose_lg_weight(W_l, N_l, L_n, N_mw, kneeling_nose_lg = False):
    """
    Inputs:
    W_l = landing design weight in lb
    N_l = ultimate landing load factor = N_gear * 1.5
    L_n = length of the nose landing gear in inches
    N_mw = number of nose wheels

    conditional inputs:
    K_np = factor for landing gear (1.15 for kneeling gear, 1 otherwise)

    outputs:
    nose landing gear weight in lb
    """
    K_np = 1
    if kneeling_nose_lg:
        K_np = 1.15
    nose_lg_weight = 0.032 * K_np * W_l**0.646 * N_l**0.2 * L_n**0.5 * N_nw**0.45
    return nose_lg_weight

def det_nacelle_group_weight(N_lt, N_w, N_z, N_en, S_n, pylon_mounted = False, W_ec = 0, W_engine = 0, propeller = False, thrust_reverser = False):
    """
    Inputs:
    N_lt = nacelle length in ft
    N_w = nacelle width in ft
    N_z = ultimate load factor = 1.5 * limit load factor
    N_en = number of engines
    S_n = nacelle wetted area in ft^2

    conditional inputs:
    K_ng = factor for nacelles (1.017 for pylon mounted nacelle, 1.0 otherwise)
    W_ec = weight of the engine and contents in lb
    (approximated by 2.331 W_engine^0.901 * K_p * K_tr)
    W_engine = weight of each engine in lb
    K_p = 1.4 for propeller, 1.0 otherwise
    K_tr = 1.18 for jet with thrust reverser, 1.0 otherwise

    outputs:
    nacelle group weight in lb
    """
    K_ng = 1
    if pylon_mounted:
        K_ng = 1.017
    if W_ec == 0 and W_engine == 0:
        raise NameError("no propulsion??")
    K_p = 1
    K_tr = 1
    if propeller:
        K_p = 1.4
    if thrust_reverser:
        K_tr = 1.18
    if W_ec == 0:
        W_ec = W_engine**0.901 * K_p * K_tr
    nacelle_group_weight = 0.6724 * K_ng * N_lt**0.10 * N_w**0.294 * N_z**0.119 * W_ec**0.611 * N_en**0.984 * S_n**0.224
    return nacelle_group_weight

def det_engine_controls_weight(N_en, L_ec):
    """
    Inputs:
    N_en = number of engines
    L_ec = length from engine front to cockpit, total if multiengine, in ft

    outputs:
    engine controls weight in lb
    """
    engine_controls_weight = 5 * N_en + 0.8 * L_ec
    return engine_controls_weight

def det_starter_weight(N_en, W_engine):
    """
    for a pneumatic starter
    inputs:
    N_en = number of engines
    W_engine = weight of each engine in lb

    outputs:
    nacelle group weight in lb
    """
    starter_weight = 49.19 * (N_en * W_engine / 1000)**0.541
    return starter_weight

def det_fuel_system_weight(V_t, V_i, V_p, N_t):
    """
    Inputs:
    V_t = total fuel volume in gallons
    V_i = total integral tank volume in gallons
    V_p = total volume of self-sealing 'protected' tanks in gallons
    N_t = number of fuel tanks

    outputs:
    fuel system weight in lb
    """
    fuel_system_weight = 2.405 * V_t**0.606 * (1 + V_i/V_t)**-1 * (1 + V_p/V_t)*N_t**0.5
    return fuel_system_weight

def det_flight_controls_weight(S_cs, I_y, N_f = 6, N_m = 1):
    """
    Inputs:
    S_cs = total area of control surfaces
    I_y = yawing moment of inertia in lb * ft^2

    conditional inputs:
    N_f = number of functions performed by controls, typically 4 - 7
    N_m = number of mechanical functions, typically 0 - 2

    outputs:
    flight controls weight in lb
    """
    flight_controls_weight = 149.9* N_f**0.554 * (1 + N_m / N_f)^-1 * S_cs**0.2 *(I_y * 10**-6)**0.07
    return flight_controls_weight

def det_APU_weight(W_APU_uninstalled):
    """
    inputs:
    W_APU_uninstalled = uninstalled weight of APU

    outputs:
    the weight of the installed APU in lb
    """
    APU_weight = 2.2 * W_APU_uninstalled
    return APU_weight

def det_instruments_weight(N_c, N_en, L_f, B_w, reciprocating = False, turboprop = False):
    """
    inputs:
    N_c = number of crew
    N_en = number of engines
    L_f = total fuselage length in feet
    B_w = wing span in feet

    conditional inputs:
    K_r = 1.133 for reciprocating engine, 1.0 otherwise
    K_tbp = 0.793 for turboprop, 1.0 otherwise

    outputs:
    the instruments weight in lb
    """
    K_r = 1
    if reciprocating:
        K_r = 1.133
    K_tbp - 1
    if turboprop:
        K_tbp = 0.793
    instruments_weight = 4.509 * K_r * K_tbp * N_c**0.541 * N_en * (L_f + B_w)**0.5
    return instruments_weight

def det_hydraulics_weight(L_f, B_w, N_f = 6):
    """
    inputs:
    L_f = total fuselage length in feet
    B_w = wing span in feet

    conditional inputs:
    N_f = number of functions performed by controls, typically 4 - 7

    outputs:
    the hydraulics weight in lb
    """
    hydraulics_weight = 0.2673 * N_f * (L_f + B_w)**0.937
    return hydraulics_weight

def det_electrical_weight(R_kva = 50, N_gen = 0, N_en = 0):
    """
    inputs:
    L_a = electrical routing distance, generators to avionics to cockpit, in ft

    conditional inputs:
    R_kva = system electrical rating typically 40-60 in kV * A
    N_gen = number of generators, typically number of engines
    N_en = number of engines

    outputs:
    the total weight of the electrical system in lb
    """
    if N_gen == 0:
        N_gen = N_en
    if N_gen == 0:
        raise NameError("No generators?")
    electrical_weight = 7.291 * R_kva**0.782 * L_a**0.346 * N_gen**0.10
    return electrical

def det_avionics_weight(W_uav = 1100):
    """
    conditional inputs:
    W_uav = uninstalled avionics weight, typically 800-1400 lb

    outputs:
    the total weight of the avionics system in lb
    """
    avionics_weight = 1.73 * W_uav**0.983
    return avionics_weight

def det_furnishings_weight():
    """
    inputs:
    N_c = number of crew
    W_c = maximum cargo weight in lb
    S_f = fuselage wetted area in feet^2

    outputs:
    the total weight of the furnishings in lb
    """
    furnishings_weight = 0.0577 * N_c**0.1 * W_c**0.393 * S_f**0.75
    return furnishings_weight

def det_aircond_weight(N_p, V_pr, W_uav = 1100):
    """
    inputs:
    N_p = number of personnel on board (crew and passengers)
    V_pr = volume of pressurised section in ft^3

    conditional inputs:
    W_uav = uninstalled avionics weight, typically 800-1400 lb

    outputs:
    the total weight of the furnishings in lb
    """
    aircond_weight = 62.36 * N_p**0.25 * (V_pr/1000)**0.604 * W_uav**0.10
    return aircond_weight

def det_anti_ice_weight(W_dg):
    """
    inputs:
    W_dg = design gross weight in lb

    outputs:
    the total weight of the anti icing system in lb
    """
    anti_ice_weight = 0.002 * W_dg
    return anti_ice_weight

def det_handling_gear_weight(W_dg):
    """
    inputs:
    W_dg = design gross weight in lb

    outputs:
    the total weight of the furnishings in lb
    """
    handling_gear_weight = 3.0 * 10**-4 * W_dg
    return handling_gear_weight
