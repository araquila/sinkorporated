import numpy as np
import matplotlib.pyplot as plt

import parameters as p

R = p.range_cruise

MTOM = p.mtom
MP = p.M_payload
OEM = p.OEW / p.g


# Point A


RA = 1.
PA = 1.

# Point B


RB = 1.
PB = 1.

# Point C


RC = 1.
PC = 1.

# Point D


RD = 1.
PD = 1.


rangelist = np.hstack([RA, RB, RC, RD]) / 1000.
payloadlist = np.hstack([PA, PB, PC, PD])

plt.plot(rangelist, payloadlist)

plt.xlabel("Range [km]")
plt.ylabel("Payload [kg]")
plt.xlim([0, 4500])
plt.ylim([0, 7000])

plt.show()