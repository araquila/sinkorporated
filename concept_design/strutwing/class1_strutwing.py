W_payload = 60000. # [N]
W_crew = 0. # [N]
W_trapped_fuel_jet = 0. # [N]
W_trapped_fuel_tbp = 0. # [N]

V_cruise_jet = 0. # [m/s]
V_cruise_tbp = 0. # [m/s]

eff_prop = 0. # [-]
cp_cruise_tbp = 0. # [g/J]
cj_cruise_jet = 0. # [g/Ns]

g = 9.80665 # [m/s^2]

A = 19.5 # [-]
e = 0.75 # [-]
C_D_0 = 0.015 # [-]
C_fe = 0.0030 # [-]

# BASED ON STATISTICS: CALCULATION OF THE OEW AND MTOW_jet
payload_ratio = 0.22 # = PL/MTOW
emptyweight_ratio = 0.636 # = OEW/MTOW

MTOW_tbp = W_payload / payload_ratio
OEW_tbp = MTOW_tbp * emptyweight_ratio

print(MTOW_tbp)
print(OEW_tbp)
