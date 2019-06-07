# This is where the variables come from
import parameters as p
import empennage.empennage_sizing as es

# Import class II mass estimation equations
import class2 as c

# Empennage weight estimation
W_h_tail = c.det_hor_tail_weight(0.12*es.c_tip_v, es.b_h, p.mtom, es.S_h, p.l_h, es.sweep_h, es.A_h, 1.78)
W_v_tail = c.det_ver_tail_weight(p.mtom, p.l_v, es.S_v, es.sweep_v, es.A_v, 0.12)