### IMPORTS
from atmosphere import atmosphere_calc
import numpy as np
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
#material aluminium 2014-T6 http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2014T6

# Ultimate
ult_stress_carbon = 600e6
ultimate_bending_stress_al2024 = 483e6
ultimate_shear_stress_al2024 = 290e6
ultimate_stress_al2024 = 469e6
tensile_yield_strength_al2014 = 414e6

# Yield stress
yield_stress_carbon = None
yield_stress_al2024 = 324e6

# Fatigue
fatigue_strength_al2014 = 124e6

# Poisson
poisson_ratio_al2014 = 0.33

# E-modulus
E_modulus = 414e6
E_al2014 = 72.4e9 #E-modulus

# Shear modulus
G_al2014 = 28e9 #shear modulus

# Density
density_aluminum = 2800

# Passengers and Crew
n_passenger = 60
M_passenger = 95
M_cargo = 10
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

#
pressure_inside = 100000 #N/m2
pressure_outside = 35000 #N/m2

# Forces
F_strut = 10000
R_y = 10000
R_x = 10000
M = 10000

## -------- WEIGHTS AND MASSES -------- ##
# General
OEW = None
MTOW = 173185.74
MLW = None
EW = None
W_fuel = 7736.30
W_pod = 157.55 * g
mtom = MTOW / g

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
W_fuselage = 2750 * g
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

# Safetyfactors
safetyfactor_wingloading = 2.5
safetyfactor_fuselage = 2
safetyfactor_wingbox = 1.5


## -------- DIMENSIONS -------- ##
# Fuselage
l_fuselage = 21.118
l_cabin = 13.7414
d_fuselage_outside = 2.84
d_fuselage_inside = None
l_nose = 2.8373002246584007
l_lavatory = 36 * 0.0254

# Wing
A = 20
S = 49.209
b = 31.372
sweep_qc = 0
dihedral = 1.
taper = 0.4
root_chord = 2.241
tip_chord = root_chord * taper
tc_ratio_root = 0.15
tc_ratio_tip = 0.12
strut_pos_perc = 0.5                    # % of span
MAC = (2/3) * root_chord * ((1 + taper + taper**2)/(1 + taper))
xLEMAC = 9.0816
x_ac_w = xLEMAC + 0.25*MAC

# Wingbox
# Width
w_root_wingbox = 1.008 #m
w_tip_wingbox = 0.4487 #m

# Height
h_max_root_wingbox = 0.35156
h_max_tip_wingbox = 0.08518

# Stringers
n_upper_skin_wingbox = 5
n_lower_skin_wingbox = 2

#thickness
t_sheet = 0.003 #m

#amount of ribs, excluding root and tip caps
n_ribs = 8
rib_spacing = (b/2)/(n_ribs+1)

# Strutbox
# Width
w_root_strutbox = 1 #m
w_tip_strutbox = 1 #m

# Height strutbox
h_max_root_strutbox = 0.4
h_max_tip_strutbox = 0.4

# Length strutbox
l_strutbox = d_fuselage_outside - 0.5

# Stringers
n_upper_skin_strutbox = 0
n_lower_skin_strutbox = 0

# Stringer geometry
A_stiffener = 0.001
h_stiffener = 0.03

# Strut
d_strut = 5

# Empennage
l_tail = 4.539680359453441

# Vertical Tail
l_v = 11
A_v = 1.2
taper_v = 0.6
sweep_v = 30
V_v = 0.05

# Horizontal Tail
l_h = 12
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
d_engine = 3.66
l_engine = None
d_prop = None
engine_pos_perc = 0.27                  # % of span
pod_pos_perc = 0.5
y_engine = 4.74

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
V_TO = np.sqrt((2 * MTOW) / (rho0 * S * C_L_max_TO))
q_TO = 0.5 * rho0 * V_TO ** 2
C_D_TO = 0.023 #NOT FINAL
Cl_delta_aileron = 0.11217
Clp = -25.8276
Cl_alpha = 2 * np.pi
tau = 0.57
Cd0 = 0.008


# Propulsion
eff_cruise = 0.85
eff_loiter = 0.77
cp_cruise = 0.8 * fuel_efficiency_factor * 74e-9
cp_loiter = 0.8*fuel_efficiency_factor * 74e-9
P_TO = 5.8e6
T_TO = 31000
tot_thrust = 78e3
