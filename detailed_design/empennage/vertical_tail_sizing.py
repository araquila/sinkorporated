import matplotlib.pyplot as plt
import parameters as p
import numpy as np
T_TO = p.T_TO
k_OEI = 1.3
rudder_angle_max = 35 * np.pi / 180         #[rad]
moment_OEI = p.y_engine * T_TO
q_TO = 0.5 * p.rho0 * (1.4 * p.V_stall) ** 2
c_ratio = 0.3


C_n_e = k_OEI * (p.y_engine * T_TO) / (q_TO * p.S * p.b)

S_ratio = np.arange(0.1,0.3,0.005)

N_v = 2.4 * (p.l_v / p.b) * (S_ratio)

#Rudder agle trim after OEI
N_min1 = 1.35 * C_n_e / rudder_angle_max

#Heading hold in crosswind
N_min2 = 0.375 * N_v /rudder_angle_max

#Change of heading after OEI
N_min3 = (C_n_e + 0.262 * N_v) / rudder_angle_max

#inabilty of rudder to cause fin stall
N_max = 0.5 * N_v / rudder_angle_max

N_rudder = 3 *((c_ratio) ** 0.47 + 0.08) * (p.l_v / p.b) * S_ratio

plt.plot(S_ratio, N_min2)
plt.plot(S_ratio, N_max)
plt.plot(S_ratio, N_rudder)
plt.plot(S_ratio, N_min3)
plt.show()