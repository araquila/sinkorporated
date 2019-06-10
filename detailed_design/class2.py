import conversion_formulas as cf
import numpy as np

def lim_load_factor(MTOM):
    """
    
    Inputs:
    MTOM = Maximum Take-Off Mass in kg
    
    Output:
    Limit Load Factor (Dimensionless)
    """
    
    MTOM = cf.kg_to_pounds(MTOM)
    
    if MTOM <= 4100:
        n_lim = 3.8
    if 4100 < MTOM <= 50000:
        n_lim = 2.1 + (24000/(MTOM + 10000))
    if MTOM > 50000:
        n_lim = 2.5
    return n_lim

def det_hor_tail_weight(F_w, B_h, W_dg, S_ht, L_t, quarter_chord_sweep_ht, AR_ht, S_e, all_moving_unit = False, K_y = 0.3):
    """
    Inputs:
    F_w = fuselage width at horizontal tail intersection in m
    B_h = horizontal tail span in m
    W_dg = design gross weigth in kg
    S_ht = horizontal tail area in m^2
    L_t = tail length (wing quarter MAC to tail quarter MAC) in m
    quarter_chord_sweep_ht = sweep angle of the horizontal tail at the quarter chord line in degrees
    AR_ht = aspect ratio of the horizontal tail
    S_e = elevator area in m^2

    conditional inputs:
    all_moving_unit = is the horizontal tail moving as a unit (True) or only part of it (False)
    all_moving_unit = True -> 1.143
    all_moving_unit = False -> 1
    K_y = Aircraft pitching radius of gyration in ft (usually 0.3 * L_t)

    output:
    Horizontal tail weight in kg
    """
    
    F_w = cf.meter_to_feet(F_w)
    B_h = cf.meter_to_feet(B_h)
    N_z = 1.5* lim_load_factor(W_dg)
    W_dg = cf.kg_to_pounds(W_dg)
    S_ht = cf.metersquared_to_feetsquared(S_ht)
    L_t = cf.meter_to_feet(L_t)
    quarter_chord_sweep_ht = quarter_chord_sweep_ht * (np.pi/180)
    S_e = cf.metersquared_to_feetsquared(S_e)
    
    K_uht = 1
    if all_moving_unit:
        K_uht = 1.143
    if K_y == 0.3:
        K_y = 0.3 * L_t
    hor_tail_weight = 0.0379 * K_uht * (1 + F_w/B_h)**-0.25 * W_dg**0.639 * N_z**0.10 * S_ht**0.75 * L_t**-1 * K_y**0.704 * np.cos(quarter_chord_sweep_ht)**-1 * AR_ht**0.166 * (1 + S_e/S_ht)**0.1
    
    hor_tail_weight = cf.pounds_to_kg(hor_tail_weight)
    return hor_tail_weight

def det_ver_tail_weight(W_dg, L_t, S_vt, quarter_chord_sweep_vt, AR_vt, t_c_root_vt, K_z = 1, TTail = True):
    """
    Inputs:
    W_dg = design gross weigth in kg
    L_t = tail length (wing quarter MAC to tail quarter MAC) in m
    S_vt = vertical tail area in m^2
    quarter_chord_sweep_vt = sweep angle of the vertical tail at the quarter chord line in degrees
    AR_vt = aspect ratio of the vertical tail
    t_c_root_vt = thickness to chord ratio at the root of the vertical tail

    conditional inputs:
    K_z = Aircraft yawing radius of gyration in ft (usually L_t)

    outputs:
    vertical tail weight in kg
    """
    
    if TTail:
        H_ht = 1
        H_vt = 1
    else:
        H_ht = 0
        H_vt = 1
    
    N_z = 1.5* lim_load_factor(W_dg)
    W_dg = cf.kg_to_pounds(W_dg)
    L_t = cf.meter_to_feet(L_t)
    S_vt = cf.metersquared_to_feetsquared(S_vt)
    quarter_chord_sweep_vt = quarter_chord_sweep_vt * (np.pi/180)

    if K_z == 1:
        K_z = L_t
        
    ver_tail_weight = 0.0026 * (1 + H_ht/H_vt)**0.225 * W_dg**0.556 * N_z**0.536 * L_t**-0.5 * S_vt**0.5 * K_z**0.875 * np.cos(quarter_chord_sweep_vt)**-1 * AR_vt**0.35 * t_c_root_vt**-0.5
    
    ver_tail_weight = cf.pounds_to_kg(ver_tail_weight)
    return ver_tail_weight

def det_main_lg_weight(W_l, L_m, N_mw, N_mss, V_stall, N_gear = 3, kneeling_main_lg = False):
    """
    Inputs:
    W_l = landing design weight in kg
    L_m = length of the main landing gear in m
    N_mw = number of main wheels
    N_mss = number of main gear shock struts
    V_stall = stall speed in m/s

    conditional inputs:
    K_mp = factor for landing gear (1.126 for kneeling landing gear, 1 otherwise)

    outputs:
    main landing gear weight in lb
    """
    
    W_l = cf.kg_to_pounds(W_l)
    N_l = N_gear * 1.5
    L_m = cf.meter_to_inch(L_m)
    V_stall = V_stall * 1.94384449
    
    K_mp = 1
    if kneeling_main_lg:
        K_mp = 1.126
        
    main_lg_weight = 0.0106 * K_mp * W_l**0.888 * N_l**0.25 * L_m**0.4 * N_mw**0.321 * N_mss**-0.5 * V_stall**0.1
    
    main_lg_weight = cf.pounds_to_kg(main_lg_weight)
    return main_lg_weight

def det_nose_lg_weight(W_l, L_n, N_nw, N_gear = 3, kneeling_nose_lg = False):
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
    
    W_l = cf.kg_to_pounds(W_l)
    N_l = N_gear * 1.5
    L_n = cf.meter_to_inch(L_n)
    
    K_np = 1
    if kneeling_nose_lg:
        K_np = 1.15
    nose_lg_weight = 0.032 * K_np * W_l**0.646 * N_l**0.2 * L_n**0.5 * N_nw**0.45
    
    nose_lg_weight = cf.pounds_to_kg(nose_lg_weight)
    return nose_lg_weight

def det_fuselage_weight(W_dg, L, S_f, taper, B_w, quarter_chord_sweep, L_D_ratio, cargo_doors = 1, fuselage_mounted_lg = True):
    """
    Inputs:
    W_dg = design gross weigth in kg
    L = fuselage structural length in m (no radome/nosecone and tail cap)
    S_f = fuselage wetted area in m^2
    taper = taper ratio
    B_w = wing span in m
    quarter_chord_sweep = sweep angle at the quarter chord line in degrees
    L_D_ratio = the lift over drag ratio of the aircraft

    conditional inputs:
    K_door = factor for the cargo doors (1 for no cargo door, 1.06 for one side cargo door, 1.12 for 2 side cargo doors)
    K_lg = factor for landing gear (1.12 for fuselage mounted landing gear, 1 otherwise)

    outputs:
    fuselage weight in kg
    """
    
    N_z = 1.5* lim_load_factor(W_dg)
    W_dg = cf.kg_to_pounds(W_dg)
    L = cf.meter_to_feet(L)
    S_f = cf.metersquared_to_feetsquared(S_f)
    B_w = cf.meter_to_feet(B_w)
    quarter_chord_sweep = quarter_chord_sweep * (np.pi/180)
    
    K_door = 1 + cargo_doors * 0.06
    K_lg = 1
    if fuselage_mounted_lg:
        K_lg = 1.12
        
    K_ws = 0.75 * (1 + 2 * taper)/(1 + taper) * (B_w * np.tan(quarter_chord_sweep)/L)
    fuselage_weight = 0.3280 * K_door * K_lg * (W_dg * N_z)**0.5 * L**0.25 * S_f**0.302 * (1 + K_ws)**0.04 * L_D_ratio**0.10
    
    fuselage_weight = cf.pounds_to_kg(fuselage_weight)
    return fuselage_weight

def det_nacelle_group_weight(N_lt, N_w, W_dg, N_en, S_n, W_engine, pylon_mounted = True, W_ec = 0, propeller = True, thrust_reverser = False):
    """
    Inputs:
    N_lt = nacelle length in m
    N_w = nacelle width in m
    W_dg = design gross weight in kg
    N_en = number of engines
    S_n = nacelle wetted area in m^2
    W_engine = weight of each engine in kg

    conditional inputs:
    K_ng = factor for nacelles (1.017 for pylon mounted nacelle, 1.0 otherwise)
    W_ec = weight of the engine and contents in lb
    (approximated by 2.331 W_engine^0.901 * K_p * K_tr)
    K_p = 1.4 for propeller, 1.0 otherwise
    K_tr = 1.18 for jet with thrust reverser, 1.0 otherwise

    outputs:
    nacelle group weight in kg
    """
    
    N_lt = cf.meter_to_feet(N_lt)
    N_w = cf.meter_to_feet(N_w)
    N_z = 1.5* lim_load_factor(W_dg)
    S_n = cf.metersquared_to_feetsquared(S_n)
    W_engine = cf.kg_to_pounds(W_engine)
    
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
        W_ec = 2.331 * W_engine**0.901 * K_p * K_tr
        
    nacelle_group_weight = 0.6724 * K_ng * N_lt**0.10 * N_w**0.294 * N_z**0.119 * W_ec**0.611 * N_en**0.984 * S_n**0.224
    
    nacelle_group_weight = cf.pounds_to_kg(nacelle_group_weight)
    return nacelle_group_weight

def det_instruments_weight(N_c, N_en, L_f, B_w, reciprocating = False, turboprop = True):
    """
    inputs:
    N_c = number of pilots
    N_en = number of engines
    L_f = total fuselage length in m
    B_w = wing span in m

    conditional inputs:
    K_r = 1.133 for reciprocating engine, 1.0 otherwise
    K_tbp = 0.793 for turboprop, 1.0 otherwise

    outputs:
    the instruments weight in kg
    """
    
    L_f = cf.meter_to_feet(L_f)
    B_w = cf.meter_to_feet(B_w)
    
    K_r = 1
    if reciprocating:
        K_r = 1.133
    K_tbp = 1
    if turboprop:
        K_tbp = 0.793
        
    instruments_weight = 4.509 * K_r * K_tbp * N_c**0.541 * N_en * (L_f + B_w)**0.5
    
    instruments_weight = cf.pounds_to_kg(instruments_weight)
    return instruments_weight

def det_hydraulics_weight(L_f, B_w, N_f = 6):
    """
    inputs:
    L_f = total fuselage length in m
    B_w = wing span in m

    conditional inputs:
    N_f = number of functions performed by controls, typically 4 - 7

    outputs:
    the hydraulics weight in kg
    """
    L_f = cf.meter_to_feet(L_f)
    B_w = cf.meter_to_feet(B_w)
    
    hydraulics_weight = 0.2673 * N_f * (L_f + B_w)**0.937
    
    hydraulics_weight = cf.pounds_to_kg(hydraulics_weight)
    return hydraulics_weight

def det_electrical_weight(L_a, N_en, R_kva = 50, N_gen = 0):
    """
    inputs:
    L_a = electrical routing distance, generators to avionics to cockpit, in m

    conditional inputs:
    R_kva = system electrical rating typically 40-60 in kV * A
    N_gen = number of generators, typically number of engines
    N_en = number of engines

    outputs:
    the total weight of the electrical system in kg
    """
    
    L_a = cf.meter_to_feet(L_a)
    
    if N_gen == 0:
        N_gen = N_en
    if N_gen == 0:
        raise NameError("No generators?")
    electrical_weight = 7.291 * R_kva**0.782 * L_a**0.346 * N_gen**0.10
    
    electrical_weight = cf.pounds_to_kg(electrical_weight)
    return electrical_weight

def det_avionics_weight(W_uav = 1100):
    """
    conditional inputs:
    W_uav = uninstalled avionics weight, typically 800-1400 lb

    outputs:
    the total weight of the avionics system in kg
    """
    avionics_weight = 1.73 * W_uav**0.983
    
    avionics_weight = cf.pounds_to_kg(avionics_weight)
    return avionics_weight

def det_furnishings_weight(N_c, W_c, S_f):
    """
    inputs:
    N_c = number of pilots
    W_c = maximum cargo weight in kg
    S_f = fuselage wetted area in m^2

    outputs:
    the total weight of the furnishings in kg
    """
    
    W_c = cf.kg_to_pounds(W_c)
    S_f = cf.metersquared_to_feetsquared(S_f)
    
    furnishings_weight = 0.0577 * N_c**0.1 * W_c**0.393 * S_f**0.75
    
    furnishings_weight = cf.pounds_to_kg(furnishings_weight)
    return furnishings_weight

def det_anti_ice_weight(W_dg):
    """
    inputs:
    W_dg = design gross weight in kg

    outputs:
    the total weight of the anti icing system in kg
    """
    
    W_dg = cf.kg_to_pounds(W_dg)
    
    anti_ice_weight = 0.002 * W_dg
    
    anti_ice_weight = cf.pounds_to_kg(anti_ice_weight)
    return anti_ice_weight

def det_aircond_weight(N_p, V_pr, W_uav = 1100):
    """
    inputs:
    N_p = number of personnel on board (crew and passengers)
    V_pr = volume of pressurised section in m^3

    conditional inputs:
    W_uav = uninstalled avionics weight, typically 800-1400 lb

    outputs:
    the total weight of the furnishings in kg
    """
    
    V_pr = cf.metercubed_to_feetcubed(V_pr)
    
    aircond_weight = 62.36 * N_p**0.25 * (V_pr/1000)**0.604 * W_uav**0.10
    
    aircond_weight = cf.pounds_to_kg(aircond_weight)
    return aircond_weight

def det_handling_gear_weight(W_dg):
    """
    inputs:
    W_dg = design gross weight in kg

    outputs:
    the total weight of the furnishings in kg
    """
    
    W_dg = cf.kg_to_pounds(W_dg)
    
    handling_gear_weight = 3.0 * 10**-4 * W_dg
    
    handling_gear_weight = cf.pounds_to_kg(handling_gear_weight)
    return handling_gear_weight