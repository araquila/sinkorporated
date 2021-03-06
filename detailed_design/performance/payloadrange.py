import numpy as np
import matplotlib.pyplot as plt

import parameters as p

MTOM = p.mtom
MP = p.M_payload + 300
OEM = p.OEW / p.g
MF = p.W_fuel / p.g

# Fuel fractions
f1 = 0.990      
f2 = 0.995      
f3 = 0.995      
f4 = 0.985     
f5 = 0.985      
f6 = 0.995      


# Point A
RA = 0
PA = MP

# Point B
MPB = MP
MFB = MTOM - MPB - OEM
FFB = 1 - (MFB / MTOM)

W_cruiseB = 1 / (FFB/(f1*f2*f3*f4*f5*f6))
R_cruiseB = (p.eff_cruise/(p.g*p.cp_cruise)) * p.LD_ratio * np.log(W_cruiseB)

RB = R_cruiseB
PB = MPB

# Point C
MFC = MF
MPC = MTOM - OEM - MFC
FFC = 1 - (MFC / MTOM)

W_cruiseC = 1 / (FFC/(f1*f2*f3*f4*f5*f6))
R_cruiseC = (p.eff_cruise/(p.g*p.cp_cruise)) * p.LD_ratio * np.log(W_cruiseC)

RC = R_cruiseC
PC = MPC

# Point D
MPD = 0.
MTOMD = MTOM - MP
MFD = MF
FFD = 1 - (MFD / MTOMD)

W_cruiseD = 1 / (FFD/(f1*f2*f3*f4*f5*f6))
R_cruiseD = (p.eff_cruise/(p.g*p.cp_cruise)) * p.LD_ratio * np.log(W_cruiseD)

RD = R_cruiseD
PD = MPD

rangelist = np.hstack([RA, RB, RC, RD]) / 1000.
payloadlist = np.hstack([PA, PB, PC, PD])

plt.plot(rangelist, payloadlist)

plt.plot(np.linspace(0,RC,100)/1000, np.ones(100)*PC, linestyle=':', color='black')
plt.plot(np.ones(100)*RC/1000, np.linspace(0,PC,100), linestyle=':', color='black')

plt.xlabel("Range [km]", size="large")
plt.ylabel("Payload [kg]", size="large")
plt.xlim([0, 4500])
plt.ylim([0, 7000])

plt.show()