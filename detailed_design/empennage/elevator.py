import numpy as np
import parameters as p
import empennage.horizontal_tail_sizing as h
import matplotlib.pyplot as plt
import empennage.empennage_sizing as emp

V_R = 1.3 * p.V_stall
mu = 0.04
rotation_acc_TO = 7 * (np.pi / 180) #[rad/s^2]
R_y = 0.34
I_yy = 1.34 * ((p.l_fuselage*3.28084)**2 * (p.MTOW*0.225) * R_y**2) / (4 * 32.17)
alpha_w = - np.pi 
elevator_rotation_max = -25 * np.pi / 180
span_ratio = 0.8

#LOCATIONS
x_cg = 11.06
x_mg = 12
x_ac_h = p.l_fuselage - x_cg - 0.75 * emp.c_root_v
z_drag = 0
z_cg = 0
z_mg = -1.5
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

C_L_h = L_h / (p.q_TO * p.S * span_ratio)

downwash_angle_0 = (2 * p.C_L_max_TO) / (np.pi * p.A)

de_daplha = (2 * h.C_L_alpha_w) / (np.pi * p.A)

downwash = downwash_angle_0 + (de_daplha * alpha_w)

alpha_h = alpha_w - downwash

elevator_effectiveness = (alpha_h * np.pi / 180 + (C_L_h / h.C_L_alpha_h)) / elevator_rotation_max

print(elevator_effectiveness)

x = np.arange(0,0.5,0.01)
y = -6.624 * x**4 + 12.07 * x**3 - 8.292 * x**2 + 3.295 *x + 0.004942 - elevator_effectiveness


x_0 = [0,0.5]
y_0 = [0,0]
plt.scatter(x,y)
plt.plot(x_0,y_0)
plt.show()