import numpy as np
import parameters as p
import empennage.empennage_sizing as emp

S_df = 0.164 * emp.S_v
phi_0_v = 90 - p.sweep_v
phi_0_dorsal = 70.49 + 0.141 * phi_0_v
sweep_0_dorsal = 90 - phi_0_dorsal

phi_0_v = phi_0_v * (np.pi / 180)
phi_0_dorsal = phi_0_dorsal * (np.pi / 180)

h_dorsal = np.sqrt((2 * S_df) / (np.tan(phi_0_dorsal) - np.tan(phi_0_v)))

print("Dorsal fin surface area =", S_df)
print("Dorsal fin leading edge sweep =", sweep_0_dorsal)
print("Dorsal fin height =", h_dorsal)