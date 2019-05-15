from math import *
from matplotlib import pyplot as plt
from main_boxwing import MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet, W_used_fuel_jet

# Gravitional constant
g = 9.8065


# Passengers and crew
n_passenger = 60
M_passenger = 105           # (Including luggage)
n_crew = 4
M_crew_member = 100

#M_OEM = 17133.6
#M_MTOM = 25584.5

# Initial mass and fractions
M_payload = n_passenger * M_passenger
M_crew = n_crew * M_crew_member
f_trapped_fuel = 0.003      # Range 0.001-0.005
#M_empty_tbp = M_OEM-M_crew-f_trapped_fuel*M_MTOM

#M_empty_jet = M_OEM-M_crew-f_trapped_fuel*M_MTOM

# Convert to weights
W_payload = M_payload * g
#W_crew = M_crew * g
#W_empty_jet = M_empty_jet * g

## Initial jet and tbp aircraft parameters
C_fe = 0.003
S = 1
S_wet = 5 * S

# Jet
A_jet = 12
e_jet = 1.2
cj_loiter_jet = 19e-6       # (0.4-0.6) [lbs/lbs/hr]
cj_cruise_jet = 19e-6      # (0.5-0.9) [lbs/lbs/hr]
V_cruise_jet = 229
range_cruise_jet = 1850000
endurance_loiter_jet = 2700

C_D_0 = C_fe * S_wet/S
C_L_loiter_jet = sqrt(C_D_0 * pi * A_jet * e_jet)
C_D_loiter_jet = 2 * C_D_0
LD_loiter_jet = C_L_loiter_jet / C_D_loiter_jet

Wpayload = W_payload #+ 2000*g

MTOW = MTOW_jet
Wf = MTOW - OEW_jet - Wpayload

print(W_used_fuel_jet)

RA = 0
PA = Wpayload
W_used_fuel_jet
Mff = 1 - (Wf/MTOW)

f1_jet = 0.990      # W_1 / W_TO (Engine start, warm-up)
f2_jet = 0.990      # W_2 / W_1 (Taxi)
f3_jet = 0.995      # W_3 / W_2 (Take-off)
f4_jet = 0.980      # W_4 / W_3 (Climb)
f5_jet = None       # W_5 / W_4 (Cruise)
f6_jet = None       # W_6 / W_5 (Loiter)
f7_jet = 0.990      # W_7 / W_6 (Descent)
f8_jet = 0.992      # W_8 / W_7 (Landing, taxi, shutdown)

# Calculation of cruise fuel fraction
f5_jet = 1/exp(range_cruise_jet/(((V_cruise_jet/(g*cj_cruise_jet))*LD_cruise_jet)))
print(f5_jet)

f6_jet = 1/exp(endurance_loiter_jet/((1/(g*cj_loiter_jet))*LD_loiter_jet))
print(f6_jet)


W54 = Mff/(f1_jet*f2_jet*f3_jet*f4_jet*f6_jet*f7_jet*f8_jet)
W45 = 1/W54
print(W54)

R45 = (V_cruise_jet/(g*cj_cruise_jet))*LD_cruise_jet*log(W45)
print(R45)
pointB = [R45, W_payload]
RB = R45
PB = Wpayload

MTOW = MTOW_jet
W_fuel_jet = W_fuel_jet + 10000
Wpay = MTOW - OEW_jet - W_fuel_jet


Mff = 1 - (W_fuel_jet/MTOW)

f1_jet = 0.990      # W_1 / W_TO (Engine start, warm-up)
f2_jet = 0.990      # W_2 / W_1 (Taxi)
f3_jet = 0.995      # W_3 / W_2 (Take-off)
f4_jet = 0.980      # W_4 / W_3 (Climb)
f5_jet = None       # W_5 / W_4 (Cruise)
f6_jet = None       # W_6 / W_5 (Loiter)
f7_jet = 0.990      # W_7 / W_6 (Descent)
f8_jet = 0.992      # W_8 / W_7 (Landing, taxi, shutdown)

# Calculation of cruise fuel fraction
f5_jet = 1/exp(range_cruise_jet/(((V_cruise_jet/(g*cj_cruise_jet))*LD_cruise_jet)))

f6_jet = 1/exp(endurance_loiter_jet/((1/(g*cj_loiter_jet))*LD_loiter_jet))

W54 = Mff/(f1_jet*f2_jet*f3_jet*f4_jet*f6_jet*f7_jet*f8_jet)
W45 = 1/W54

R45 = (V_cruise_jet/(g*cj_cruise_jet))*LD_cruise_jet*log(W45)
pointC = [R45, Wpay]
print(pointC)
RC = R45
PC = Wpay

MTOW = W_fuel_jet + OEW_jet

Wf = W_fuel_jet

Mff = 1 - (Wf/MTOW)

f1_jet = 0.990      # W_1 / W_TO (Engine start, warm-up)
f2_jet = 0.990      # W_2 / W_1 (Taxi)
f3_jet = 0.995      # W_3 / W_2 (Take-off)
f4_jet = 0.980      # W_4 / W_3 (Climb)
f5_jet = None       # W_5 / W_4 (Cruise)
f6_jet = None       # W_6 / W_5 (Loiter)
f7_jet = 0.990      # W_7 / W_6 (Descent)
f8_jet = 0.992      # W_8 / W_7 (Landing, taxi, shutdown)

# Calculation of cruise fuel fraction
f5_jet = 1/exp(range_cruise_jet/(((V_cruise_jet/(g*cj_cruise_jet))*LD_cruise_jet)))

f6_jet = 1/exp(endurance_loiter_jet/((1/(g*cj_loiter_jet))*LD_loiter_jet))

W54 = Mff/(f1_jet*f2_jet*f3_jet*f4_jet*f6_jet*f7_jet*f8_jet)
W45 = 1/W54

R45 = (V_cruise_jet/(g*cj_cruise_jet))*LD_cruise_jet*log(W45)
pointD = [R45, 0]
RD = R45
PD = 0
print(R45)

R = [RA, RB, RC, RD]
P = [PA/g, PB/g, PC/g, PD/g]

plt.plot(R, P)
plt.show()
