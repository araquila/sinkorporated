import numpy as np
import parameters as p
import empennage.horizontal_tail_sizing as h
V_R = p.V_stall
mu = 0.04
rotation_acc_TO = 7 * (2 * np.pi / 180) #[deg/s^2]
R_y = 0.34
I_yy = 1.34 * (p.l_fuselage**2 * MTOW * R_y**2) / (4 * p.g)
#LOCATIONS
x_cg = 11.06
x_mg = 12
x_ac_h = p.l_fuselage - x_cg - 0.75 * 2.41
z_drag = 0
z_cg = 0
z_mg = -1
z_thrust = 0.5

#FORCES
L_wf = p.q_TO * p.C_L_max_TO * p.S
N = p.MTOW - L_wf
F_f = mu * N
D_TO = p.q_TO * p.C_D_TO * p.S
M_ac = p.q_TO * h.C_m_ac * p.S * p.MAC

acc_TO = (p.T_TO - D_TO - 0.04 * N) / p.mtom
#MOMENTS
M_weight = p.MTOW * (x_mg - x_cg)
M_drag = D_TO * (z_drag - z_mg)
M_thrust = p.T_TO * (z_thrust - z_cg)
M_lift = L_wf * (x_mg - h.x_ac)
M_acc = p.mtom * acc_TO * (z_cg - z_mg)

L_h = (M_lift + M_ac + M_acc - M_weight + M_drag - M_thrust - I_yy * rotation_acc_TO) / (x_ac_h - x_mg)