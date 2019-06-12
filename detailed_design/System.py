import parameters as p
import numpy as np



#def det_hor_tail_weight(F_w, B_h, W_dg, N_z, S_ht, L_t, quarter_chord_sweep_ht, AR_ht, S_e, all_moving_unit = False, K_y = 0.3):
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
#if all_moving_unit:
K_uht = 1.143
B_h = p.b_h * 3.28
W_dg = p.mtom * 2.20
S_ht = p.S_h * 3.28 * 3.28
L_t = p.l_h * 3.28
AR_ht = p.A_h 
S_e = p.S_e
quarter_chord_sweep_ht = p.sweep_h * (2 * np.pi / 180)
F_w = p.d_fuselage_outside * 3.28
N_z = 1.5 * 2

K_y = 0.3 * L_t
hor_tail_weight = 0.0379 * K_uht * (1 + F_w/B_h)**-0.25 * W_dg**0.639 * N_z**0.10 * S_ht**0.75 * L_t**-1 * K_y**0.704 * np.cos(quarter_chord_sweep_ht)**-1 * AR_ht**0.166 * (1 + S_e/S_ht)**0.1
#return hor_tail_weight

#def det_vert_tail_weight(H_vt, W_dg, N_z, L_t, S_vt, quarter_chord_sweep_vt, AR_vt, t_c_root_vt, K_z = 1, T_tail = False):
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

H_vt = p.b_v * 3.28

#if T_tail == True:
H_ht = H_vt
S_vt = p.S_v
K_z = L_t
quarter_chord_sweep_vt = p.sweep_v * (2 * np.pi / 180)
AR_vt = p.A_v
t_c_root_vt = 0.11

vert_tail_weight = 0.0026 * (1 + H_ht/H_vt)**0.225 * W_dg**0.556 * N_z**0.536 * L_t**-0.5 * S_vt**0.5 * K_z**0.875 * np.cos(quarter_chord_sweep_vt)**-1 * AR_vt**0.35 * t_c_root_vt**-0.5
#return vert_tail_weight

#def det_instruments_weight(N_c, N_en, L_f, B_w, reciprocating = False, turboprop = False):
#"""
#inputs:
#N_c = number of pilots
#N_en = number of engines
#L_f = total fuselage length in feet
#B_w = wing span in feet
#
#conditional inputs:
#K_r = 1.133 for reciprocating engine, 1.0 otherwise
#K_tbp = 0.793 for turboprop, 1.0 otherwise
#
#outputs:
#the instruments weight in lb
#"""
N_c = p.n_pilots
L_f = p.l_fuselage * 3.28
N_en = p.n_engines
B_w = p.b * 3.28
K_r = 1
#if reciprocating:
#K_r = 1.133
K_tbp = 1
#if turboprop:
K_tbp = 0.793
instruments_weight = 4.509 * K_r * K_tbp * N_c**0.541 * N_en * (L_f + B_w)**0.5
#return instruments_weight

#def det_hydraulics_weight(L_f, B_w, N_f = 6):
"""
inputs:
L_f = total fuselage length in feet
B_w = wing span in feet

conditional inputs:
N_f = number of functions performed by controls, typically 4 - 7

outputs:
the hydraulics weight in lb
"""
#hydraulics_weight = 0.2673 * N_f * (L_f + B_w)**0.937
#return hydraulics_weight

#def det_electrical_weight(L_a, R_kva = 50, N_gen = 0, N_en = 0):
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
#if N_gen == 0:
#    N_gen = N_en
#if N_gen == 0:
#    raise NameError("No generators?")
#electrical_weight = 7.291 * R_kva**0.782 * L_a**0.346 * N_gen**0.10
#return electrical_weight
#
#def det_avionics_weight(W_uav = 1100):
#"""
#conditional inputs:
#W_uav = uninstalled avionics weight, typically 800-1400 lb
#
#outputs:
#the total weight of the avionics system in lb
#"""
#avionics_weight = 1.73 * W_uav**0.983
#return avionics_weight
#
#def det_furnishings_weight(N_c, W_c, S_f):
#"""
#inputs:
#N_c = number of pilots
#W_c = maximum cargo weight in lb
#S_f = fuselage wetted area in feet^2
#
#outputs:
#the total weight of the furnishings in lb
#"""
#furnishings_weight = 0.0577 * N_c**0.1 * W_c**0.393 * S_f**0.75
#return furnishings_weight
#
#def det_aircond_weight(N_p, V_pr, W_uav = 1100):
#"""
#inputs:
#N_p = number of personnel on board (crew and passengers)
#V_pr = volume of pressurised section in ft^3
#
#conditional inputs:
#W_uav = uninstalled avionics weight, typically 800-1400 lb
#
#outputs:
#the total weight of the furnishings in lb
#"""
#aircond_weight = 62.36 * N_p**0.25 * (V_pr/1000)**0.604 * W_uav**0.10
#return aircond_weight
#
#def det_anti_ice_weight(W_dg):
#"""
#inputs:
#W_dg = design gross weight in lb
#
#outputs:
#the total weight of the anti icing system in lb
#"""
#anti_ice_weight = 0.002 * W_dg
#return anti_ice_weight
#
#def det_handling_gear_weight(W_dg):
#"""
#inputs:
#W_dg = design gross weight in lb
#
#outputs:
#the total weight of the furnishings in lb
#"""
#handling_gear_weight = 3.0 * 10**-4 * W_dg
#return handling_gear_weight


