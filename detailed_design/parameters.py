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

# Materials
ult_stress_carbon = 600e6

# Passengers and Crew
n_passenger = 60
M_passenger = 85
M_cargo = 20
M_total_cargo = M_cargo * n_passenger
n_crew = 4
n_pilots = 2
M_crew_member = 100
M_payload = n_passenger * M_passenger + n_passenger * M_cargo
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
seat_pitch = 31 * 0.0254

# Fuel
energy_density_LNG = 53.6
energy_density_kerosene = 43
chosen_fuel_energy_density = energy_density_LNG
fuel_efficiency_factor = energy_density_kerosene/chosen_fuel_energy_density

## -------- WEIGHTS AND MASSES -------- ##
# General
OEW = None
MTOW = 173185.74
MLW = None
EW = None
W_fuel = 7736.30
W_pod = 157.55 * g

# Propulsion
M_engine = 450
W_engine = M_engine * g
W_nacelle = None
W_engine_controls = None
W_starter = None
W_APU = None
W_fuel_system = None

# Wing
W_wing = 1288 * g
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

x_ac = 10.2

# Fuselage
l_fuselage = 21.118
d_fuselage_outside = 2.84
d_fuselage_inside = None

# Wing
A = 18
S = 52
b = 29.76
sweep_qc = 0
dihedral = None
taper = 0.4
root_chord = 1.8
tc_ratio = 0.15
strut_pos_perc = 0.5                    # % of span
MAC = (2/3) * root_chord * ((1 + taper + taper**2)/(1 + taper))
xLEMAC = 10

# Empennage
# Vertical Tail
l_v = 11
A_v = 1.5
taper_v = 0.6
sweep_v = 30
V_v = 0.05

# Horizontal Tail
l_h = 11
A_h = 4
taper_h = 0.5
sweep_h = 15
V_h = 1.57

# Undercarriage
main_landing_pos = 11
nose_landing_pos = 3
h_wheel = None
lateral_pos = None
size_tire = None

# Propulsion
n_fueltanks = 2
n_blades = 6
d_engine = None
l_engine = None
d_prop = None
engine_pos_perc = 0.27                  # % of span
pod_pos_perc = 0.5

x_engine = engine_pos_perc*b/2
x_pod = pod_pos_perc*b/2

#Strut
strut_pos_perc = 0.5
x_strut = strut_pos_perc*b/2


## -------- PERFORMANCE -------- ##
# Aerodynamic
e = 0.85
M_cruise = 0.6
C_L_max_land = 2.4
C_L_max_TO = 1.4
C_L_cruise = 0.5
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
