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

#print(det_hor_tail_weight(0.4, 7.35, 18000, 13.51, 12, 15, 4, 1.78))
#print(det_ver_tail_weight(18000, 11, 9.35, 30, 1.2, 0.12))

