import numpy as np
import parameters as p
import empennage.empennage_sizing as emp

S_df = 0.164 * emp.S_v
S_df_total = 0.19 * emp.S_v
delta_sweep = 70.49 - p.sweep_v
sweep_0_dorsal = 2.244 * p.sweep_v
c_root_dorsal = 2.699 * emp.c_root_v
l_dorsal = 1.156 * c_root_dorsal
h_dorsal = 3.664 * emp.b_v

print("Dorsal fin surface area =", S_df)
print("Dorsal fin total surface area =", S_df_total)
print("Dorsal fin leading edge sweep =", sweep_0_dorsal)
print("Dorsal fin root chord =", c_root_dorsal)
print("Dorsal fin length =", l_dorsal)
print("Dorsal fin height =", h_dorsal)