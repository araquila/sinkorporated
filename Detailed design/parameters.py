### IMPORTS
from atmosphere import atmosphere_calc

### AIRCRAFT PARAMETERS ###

## -------- CONSTANTS -------- ##
# Atmosphere
altitude = 8000
temperature0 = 288.15
temperature_gradient = -0.0065
gamma = 1.4
rho0 = 1.225
g = 9.80665
R = 287
temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, temperature0, temperature_gradient, g, R, gamma)

# Passengers and Crew
n_passenger = 60
M_passenger = 105                   
n_crew = 4
n_pilots = 2
M_crew_member = 100
M_payload = n_passenger * M_passenger
M_crew = n_crew * M_crew_member

# Convert to weights
W_payload = M_payload * g
W_crew = M_crew * g

# Mission
range_cruise = 1850000
V_loiter = 80
endurance_loiter = 2700 

# Initial aircraft parameters
C_fe = 0.003
S_ratio = 1/5.5
n_engines = 2
n_seats_abreast = 4
n_aisles = 1

# Fuel
energy_density_LNG = 53.6                                                   
energy_density_kerosene = 43                                                
chosen_fuel_energy_density = energy_density_LNG
fuel_efficiency_factor = energy_density_kerosene/chosen_fuel_energy_density

## -------- WEIGHTS AND MASSES -------- ##
# General
OEW = None
MTOW = None
MLW = None
EW = None
W_fuel = None

# Propulsion
M_engine = 450                    
W_engine = M_engine * g
W_nacelle = None
W_engine_controls = None
W_starter = None
W_APU = None   
W_fuel_system = None

# Wing
W_wing = None
W_flight_controls = None
W_anti_ice = None           

# Fuselage
W_fuselage = None
W_furnishings = None

# Empennage
W_hor_emp = None
W_ver_emp = None

# Undercarriage
W_main_landing = None
W_nose_landing = None

# Other systems
W_avionics = None
W_airco = None
W_instruments =  None
W_hydraulics = None
W_electrical = None
W_handling_gear = None

## -------- DIMENSIONS -------- ##

# Fuselage
l_fuselage = 21.118
d_fuselage_outside = None
d_fuselage_inside = None

# Wing
A = 18
S = None
b = None
sweep_qc = None
dihedral = None
taper = None
root_chord = None
tc_ratio = None

# Empennage
# Vertical
l_v = 11                       
V_v = 0.07                        
b_v = None
sweep_qc_v = None
sweep_LE_v = None
tip_chord_v = None
root_chord_v = None
taper_v = tip_chord_v/root_chord_v
A_v = None

# Horizontal
l_h = 11                           
V_h = 1.57 
b_h = None
sweep_qc_h = None
sweep_LE_h = None
tip_chord_h = None
root_chord_h = None
taper_h = tip_chord_h/root_chord_h
A_h = None

# Undercarriage
main_landing_pos = 11             
nose_landing_pos = 3   
h_wheel = None
lateral_pos = None
size_tire = None

# Propulsion
pos_engine = 10                    
n_fueltanks = 2                     
n_blades = 6 
d_engine = None
l_engine = None
d_prop = None

## -------- PERFORMANCE -------- ##
# Aerodynamic
e = 0.85 
M_cruise = 0.6 
C_L_max_land = 2.4              
C_L_max_TO = 1.4 
C_L_cruise = 0.8
V_cruise = M_cruise*speed_of_sound
V_stall = 46.3                  
C_L_max_land = 2.4              
C_L_max_TO = 1.4

# Propulsion
eff_cruise = 0.85
eff_loiter = 0.77
cp_cruise = 0.8 * fuel_efficiency_factor * 74e-9
cp_loiter = 0.8*fuel_efficiency_factor * 74e-9
P_TO = 5.8e6