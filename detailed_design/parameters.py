### IMPORTS
import os
import sys
sys.path.append(os.getcwd())
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
p0 = 101325
temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, temperature0, temperature_gradient, g, R, gamma)
rho = rho0 * rho
pressure = pressure * p0

# Materials
#material aluminium 2014-T6 http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2014T6

# Ultimate
ult_stress_carbon = 600e6
ultimate_bending_stress_al2024 = 483e6
ultimate_shear_stress_al2024 = 290e6
tensile_yield_strength_al2014 = 414e6

# Yield stress
yield_stress_carbon = None
yield_stress_al2024 = 414e6

# Fatigue
fatigue_strength_al2014 = 124e6

# Poisson
poisson_ratio_al2014 = 0.33

# E-modulus
E_al2014 = 72.4e9 #E-modulus

## Trade-off aluminum fuselage skin
yield_al_2198_T8 = 407e6
yield_al_2199_T8E74 = 345e6
yield_al_2060_T8E30 = 345e6


# Shear modulus
G_al2014 = 28e9 #shear modulus

# Density
density_al_2198_T8 = 2.7 * 1000
density_carbon = 1.55 * 1000

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
pressure_inside = 81200 #N/m2
pressure_outside = 35000 #N/m2

# Forces
F_strut = 10000
R_y = 10000
R_x = 10000
M = 10000

## -------- WEIGHTS AND MASSES -------- ##
# General
OEW = 9074.716123 * g
MTOW = 17431.71612 * g
EW = None
W_fuel = 7736.30
W_pod = 300 * g
mtom = MTOW / g
W_empty = 7271.44 * g
W_wing = 13.85/100*W_empty/2
MLW = MTOW - W_fuel

# Propulsion
M_engine = 924.00 / 2
W_nacelle = (257.190498 / 2) * g
W_engine = M_engine * g  + W_nacelle
W_engine_controls = 0
W_starter = 0
W_APU = 61.23 * g
W_fuel_system = 358.7959
W_fuel_tanks = 436 * g

# Wing
W_wing = 1288 * g
W_flight_controls = 330 * g
W_anti_ice = 35.32006139 * g

# Fuselage
W_fuselage = 2750 * g
W_furnishings = 360.3467386

# Empennage
W_hor_emp = 159.1798947 * g
W_ver_emp = 83.85475628 * g

# Undercarriage
W_main_landing = 661.1444655 * g
W_nose_landing = 141.4487833 * g

# Other systems
W_avionics = 767.91 * g
W_airco = 412.35 * g
W_instruments =  62.05734245 * g
W_hydraulics = 60.49532434 * g
W_electrical = 341.580681 * g
W_handling_gear = 5.30 * g

# Safetyfactors
safetyfactor_wingloading = 2.5
safetyfactor_fuselage = 2
safetyfactor_wingbox = 1.5

m_fuselage = (W_fuselage + W_furnishings + W_avionics + W_airco + W_instruments + W_hydraulics + W_electrical + W_handling_gear + W_APU) / g
m_tail = (W_hor_emp + W_ver_emp) / g
m_wing = (W_wing + W_flight_controls + W_anti_ice) / g
m_engine = (2 * W_engine + W_fuel_system) / g
m_fuel_tanks = W_fuel_tanks / g

## -------- DIMENSIONS -------- ##
# Fuselage
l_fuselage = 21.118
l_cabin = 13.7414
d_fuselage_outside = 2.84
d_fuselage_inside = None
l_nose = 2.8373002246584007
l_lavatory = 36 * 0.0254
S_wet_fuselage = np.pi * d_fuselage_outside * l_fuselage
volume_fuselage = 0.25 * np.pi * d_fuselage_outside**2 * l_fuselage

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

n_upper_skin_wingbox = 13
n_lower_skin_wingbox = 13

#al 2099-t83 http://morita1950.info/akio/data/Al-li%20Alloy.pdf
density_stiffeners = 2780

t_hat = 0.0018
t_z = 0.0018

#2024 nrs for tradeoff
#ultimate_compressive_strength_2099 = 324*10**6
#ultimate_yield_strength_2099 = 324*10**6
#
#E_compressive_2099 = 73.1*10**9

ultimate_compressive_strength_2099 = 476*10**6
ultimate_yield_strength_2099 = 490*10**6

E_compressive_2099 = 82.1*10**9

#al2195-t84 https://www.constellium.com/sites/default/files/markets/airware_2195_t84_plate.pdf
#thickness
t_sheet = 0.0027 #m
t_sheet_strutbox = 0.003


#2024 nrs for tradeoff
#E_sheet = 73.1*10**9
#density_sheet = 2780
#ult_tensile_strength_2195 = 469*10**6
#ult_yield_strength_2195 = 324*10**6
#ult_shear_strength_2195 = 283*10**6
#fracture_toughness_2195 = 26*10**6

#E_sheet = 78*10**9
#density_sheet = 2700
#ult_tensile_strength_2195 = 595*10**6
#ult_yield_strength_2195 = 500*10**6
#ult_shear_strength_2195 = 350*10**6
#fracture_toughness_2195 = 35*10**6

#al 2055-t84 https://www.arconic.com/adip/catalog/AFE2055-factsheet.pdf

E_sheet = 76.5*10**9
density_sheet = 2710
ult_tensile_strength_2195 = 565*10**6
ult_yield_strength_2195 = 538*10**6
ult_shear_strength_2195 = 0.55*ult_yield_strength_2195
fracture_toughness_2195 = 35*10**6


#amount of ribs, excluding root and tip caps
n_ribs = 11
rib_spacing = (b/2)/(n_ribs+1)
t_rib = 0.003

#al2050-t84 https://www.constellium.com/sites/default/files/markets/airware_2050_t84_plate.pdf
E_rib = 76.5*10**9
density_rib = 2700



safety_factor_compression = 1.0
safety_factor_tension = 1.0

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
n_upper_skin_strutbox = 2
n_lower_skin_strutbox = 3

# Stringer geometry
#A_stiffener = 0.001
#h_stiffener = 0.03

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
S_v = 0.28 * S
b_v = np.sqrt(S_v * A_v)

# Horizontal Tail
l_h = 12
A_h = 4
taper_h = 0.5
sweep_h = 15
V_h = 1.57
S_h = 0.29 * S
b_h = np.sqrt(S_h * A_h)
S_e = 0.4

# Undercarriage
main_landing_pos = 11
nose_landing_pos = 3
h_wheel = None
lateral_pos = None
size_tire = None

# Propulsion
n_fueltanks = 2
n_blades = 6
w_engine = 0.66
l_engine = 2.134
h_engine = 0.838
d_prop = 3.66

strut_pos_perc = 0.5                    # % of span
engine_pos_perc = 0.27                  # % of span
pod_pos_perc = 0.55
y_engine = 4.74

x_engine = engine_pos_perc*b/2
x_pod = pod_pos_perc*b/2

#Strut
x_strut = strut_pos_perc*b/2




## -------- PERFORMANCE -------- ##
# Aerodynamic
e = 0.90
M_cruise = 0.6
C_L_max_land = 2.6
C_L_max_TO = 1.8
C_L_cruise = 0.5
V_cruise = M_cruise*speed_of_sound
V_stall = 46.3
C_L_max_land = 2.6
C_L_max_TO = 1.8
V_TO = np.sqrt((2 * MTOW) / (rho0 * S * C_L_max_TO))
q_TO = 0.5 * rho0 * V_TO ** 2
C_D_TO = 0.023 #NOT FINAL
Cl_delta_aileron = 0.11217
Clp = -25.8276
Cl_alpha = 2 * np.pi
tau = 0.57
Cd0 = 0.02 #NOT FINAL
LD_ratio = 28 #NOT FINAL


# Propulsion
eff_cruise = 0.85
eff_loiter = 0.77
cp_cruise = 0.8 * fuel_efficiency_factor * 74e-9
cp_loiter = 0.8*fuel_efficiency_factor * 74e-9
P_TO = 2.38e6
T_TO = 49.8e3
tot_thrust = 78e3
