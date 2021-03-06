# This is where the variables come from
import parameters as p
import empennage.empennage_sizing as es
import undercarriage.undercarriage_sizing as uc

# Import class II mass estimation equations
import class2 as c

##### WEIGHT ESTIMATIONS #####

# Empennage
W_h_tail = c.det_hor_tail_weight(0.12*es.c_tip_v, es.b_h, p.mtom, es.S_h, p.l_h, es.sweep_h, es.A_h, 1.78)
W_v_tail = c.det_ver_tail_weight(p.mtom, p.l_v, es.S_v, es.sweep_v, es.A_v, 0.12)

# Landing Gear
W_nose_landing = c.det_nose_lg_weight(p.mtom, uc.wheel_height, 2)
W_main_landing = c.det_main_lg_weight(p.mtom, uc.wheel_height, 4, 2, p.V_stall)

# Fuselage
W_fuselage = c.det_fuselage_weight(p.mtom, p.l_cabin, p.S_wet_fuselage, p.taper, p.b, p.sweep_qc, p.LD_ratio)

# Engines
W_engine = 924
W_nacelle = c.det_nacelle_group_weight(p.l_engine, p.w_engine, p.mtom, p.n_engines, 4*p.l_engine*p.w_engine, W_engine/p.n_engines)

# Miscellaneous
W_instruments = c.det_instruments_weight(p.n_pilots, p.n_engines, p.l_fuselage, p.b)
W_hydraulics = c.det_hydraulics_weight(p.l_fuselage, p.b)
W_electrics = c.det_electrical_weight(p.xLEMAC+p.y_engine, p.n_engines)
W_avionics = c.det_avionics_weight()
W_furnishings = c.det_furnishings_weight(p.n_pilots, p.M_total_cargo, p.S_wet_fuselage)
W_airco = c.det_aircond_weight(p.n_crew+p.n_pilots, p.volume_fuselage)
W_anti_icing = c.det_anti_ice_weight(p.mtom)
W_handling_gear = c.det_handling_gear_weight(p.mtom)
