import numpy as np
import matplotlib.pyplot as plt
import parameters as p
import empennage.horizontal_tail_sizing as h
bank_angle_1 = 30 * (2 * np.pi / 180)
time_1 = 3.5
rho = 0.525168
aileron_angle_max = 0.28
R_x = 0.22
I_xx = (3.28*(p.b+p.l_fuselage)/2)**2 *((p.MTOW * R_x**2)/(4*32.2))

L_p = p.Clp * p.V_cruise**2 * 0.5 * p.S *rho
K_1 = -rho * p.S * p.V_cruise * p.b**2 * L_p / (2 * I_xx)



L_aileron_angle_min = p.b * L_p * bank_angle_1 / (p.V_cruise * (time_1 + (1 / K_1) * np.e**(-K_1*time_1)-1)) / aileron_angle_max
c_ratio = np.arange(0.001,0.5,0.001)
midpoint = -2 *L_aileron_angle_min / ((c_ratio**0.47+0.16/p.A) * h.C_L_alpha_w)

plt.plot(midpoint, c_ratio)
plt.show()