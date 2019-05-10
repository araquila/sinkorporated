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
