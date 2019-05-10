# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from class1_conventional import Weights_Class_I
from power_wingloading_conventional import wingloading_jet, wingloading_tbp
from wingloadingfunctions import T_W_calc
import numpy as np

# Gravitional constant
g = 9.8065

# Passengers and crew
n_passenger = 60
M_passenger = 102                   #(including luggage)
n_crew = 4
M_crew_member = 100

# Initial mass and fractions
M_payload = n_passenger * M_passenger
M_crew = n_crew * M_crew_member
f_trapped_fuel = 0.003              # Range 0.001-0.005
M_empty_tbp = 14400                 
M_empty_jet = 16300                 

# Convert to weights
W_payload = M_payload * g
W_crew = M_crew * g
W_empty_tbp = M_empty_tbp * g
W_empty_jet = M_empty_jet * g

# Initial jet and tbp aircraft parameters
C_fe = 0.003
S = 1                               
S_wet = 5 * S                       

# Jet
A_jet = 10                          
e_jet = 0.8                         # Adjust per concept
cj_loiter_jet = 19e-6               # (0.4-0.6) [g/j] Propfan: 0.441
cj_cruise_jet = 19e-6               # (0.5-0.9) [g/j] Propfan: 0.441
V_cruise_jet =  200                 # [m/s]
S_jet = 61
TOP_jet = 6698

# Tbp
A_tbp = 12                          
e_tbp = 0.85                        # Adjust per concept
eff_cruise_tbp = 0.85               # [-]
eff_loiter_tbp = 0.77               # [-]
cp_cruise_tbp = 90e-9               # (0.4-0.6) [kg/ns]
cp_loiter_tbp = 90e-9               # (0.5-0.7) [kg/ns]
V_cruise_tbp = 150                  # [m/s]
S_tbp = 76
TOP_tbp = 139

# Iterator
for iter in range(1):
    ## CLASS I
    MTOW_jet, OEW_jet, W_fuel_jet, C_D_0_jet = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, jet = True)
    MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp = Weights_Class_I(W_empty_jet, W_empty_tbp, W_payload, W_crew, C_fe, S, S_wet, A_jet, A_tbp, e_jet, e_tbp, cj_loiter_jet, cj_cruise_jet, eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel, tbp = True)
    
    ## WING LOADING AND POWER LOADING
    W_S_landing_jet = wingloading_jet(MTOW_jet,OEW_jet,V_cruise_jet,e_jet,C_D_0_jet,A_jet,S_jet)
    W_S_landing_tbp = wingloading_tbp(MTOW_tbp, OEW_tbp, S_tbp, A_tbp, V_cruise_tbp, e_tbp, eff_cruise_tbp, C_D_0_tbp)
    
    ## T_W CALCULATION
    T_W_jet_range = np.zeros(len(W_S_landing_jet))
    for i in range(len(W_S_landing_jet)):
        T_W_jet_range[i] = T_W_calc(W_S_landing_jet[i],TOP_jet,1.32)
        
    S_jet = MTOW_jet/W_S_landing_jet    
    
    
MTOM_jet = MTOW_jet / g
MTOM_tbp = MTOW_tbp / g

print('Jet: ' + str(MTOW_jet))
print('Tbp: ' + str(MTOW_tbp))


