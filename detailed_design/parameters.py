### IMPORTS
import os
import sys
sys.path.append(os.getcwd())
from detailed_design.atmosphere import atmosphere_calc
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
yield_al_2060_T8E30 = 488e6


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
OEW = 10270.31004 * g

MTOW = 18536.09649 * g
EW = 9824.301753 * g
W_fuel = 1965.786443 * g
W_pod = (471/2) * g
mtom = MTOW / g
W_empty = 9824.301753 * g
W_wing = 1375.6
MLW = MTOW - W_fuel

# Propulsion
M_engine = 836.38 / 2
W_nacelle = (357.38 / 2) * g
W_engine = M_engine * g  + W_nacelle
W_engine_controls = 0
W_starter = 0
W_APU = 61.23 * g
W_fuel_system = 358.7959 * g
W_fuel_tanks = 471 * g

# Wing
W_wing = 1508.77 * g
W_flight_controls = 328.073582 * g
W_anti_ice = 37.07219297 * g
W_strut = 79.55 * g
W_aileronflaps = 84.02 * g

# Fuselage
W_fuselage = 2516.18163 * g
W_furnishings = 161.1688986 * g

# Empennage
W_hor_emp = 159.1798947 * g
W_ver_emp = 83.85475628 * g

# Undercarriage
W_main_landing = 661.1444655 * g
W_nose_landing = 141.4487833 * g

# Other systems
W_avionics = 766.30 * g
W_airco = 324.87 * g
W_instruments =  61.93544887 * g
W_hydraulics = 90.57517011 * g
W_electrical = 340.8971918 * g
W_handling_gear = 5.56 * g

# Safetyfactors
safetyfactor_wingloading = 2.5
safetyfactor_fuselage = 2
safetyfactor_wingbox = 1.5

m_fuselage = (W_fuselage + W_furnishings + W_avionics + W_airco + W_instruments + W_hydraulics + W_electrical + W_handling_gear + W_APU) / g
m_tail = (W_hor_emp + W_ver_emp) / g
m_wing = (W_wing + W_strut + W_aileronflaps + W_flight_controls + W_anti_ice) / g
m_engine = (2 * W_engine + W_fuel_system) / g
m_fuel_tanks = W_fuel_tanks / g


W_A = OEW - M_engine*2*g # airframe weight in N
m_A = W_A/g # airframe mass in kg

## -------- DIMENSIONS -------- ##
# Fuselage
l_fuselage = 21.118
l_cabin = 13.7414
d_fuselage_outside = 2.84
d_fuselage_inside = None
l_nose = 2.8373002246584007
l_lavatory = 36 * 0.0254
l_galley = 30 * 0.0254
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

n_upper_skin_wingbox = 12
n_lower_skin_wingbox = 12



#al 2099-t83 http://morita1950.info/akio/data/Al-li%20Alloy.pdf
#https://www.smithmetal.com/2099-lithium.htm
#https://www.smithshp.com/assets/pdf/2099-aluminium-lithium.pdf


t_hat = 0.0028
t_z = 0.0028

#2024 nrs for tradeoff
#ultimate_compressive_strength_2099 = 324*10**6
#ultimate_yield_strength_2099 = 324*10**6
#
#E_compressive_2099 = 73.1*10**9
#
#density_stiffeners = 2780

#al2099-t83
#
ultimate_compressive_strength_2099 = 476*10**6
ultimate_yield_strength_2099 = 490*10**6
ultimate_tensile_strength_2099 = 560*10**6

E_compressive_2099 = 82.1*10**9

#
density_stiffeners = 2630

#al7075
#E_compressive_2099 = 71.7*10**9
#density_stiffeners = 2810
#ultimate_tensile_strength_2099 = 572*10**6
#ultimate_yield_strength_2099 = 503*10**6
#ultimate_shear_strength_2099 = 0.55*ultimate_yield_strength_2099
#fracture_toughness_2195 = 35*10**6


#CFRP
#ultimate_compressive_strength_2099 = (581-2*35.8)*10**6
#ultimate_yield_strength_2099 = (586.2-2*37.7)*10**6
#ultimate_tensile_strength_2099 = 596*10**6
#
#E_compressive_2099 = 115.9*10**9
#
#density_stiffeners = 1600
#

#al2195-t84 https://www.constellium.com/sites/default/files/markets/airware_2195_t84_plate.pdf
#thickness
t_sheet = 0.0037 #m
t_sheet_strutbox = 0.003


#2024 nrs for tradeoff
#E_sheet = 73.1*10**9
#density_sheet = 2780
#ult_tensile_strength_2195 = 469*10**6
#ult_yield_strength_2195 = 324*10**6
#ult_shear_strength_2195 = 283*10**6
#fracture_toughness_2195 = 26*10**6
#ult_compressive_strength_2195 = ult_yield_strength_2195*1.02

#al2195
#E_sheet = 78*10**9
#density_sheet = 2700
#ult_tensile_strength_2195 = 595*10**6
#ult_yield_strength_2195 = 500*10**6
#ult_shear_strength_2195 = 350*10**6
#fracture_toughness_2195 = 35*10**6
#ult_compressive_strength_2195 = 595*10**6
##
#al 2055-t84 https://www.arconic.com/adip/catalog/AFE2055-factsheet.pdf  
#https://www.smithshp.com/assets/pdf/2055-t84-extrusions.pdf
#
#
G_sheet = 27*10**9
E_sheet = 76.5*10**9
density_sheet = 2710
ult_tensile_strength_2195 = 565*10**6
ult_yield_strength_2195 = 538*10**6
ult_shear_strength_2195 = 0.55*ult_yield_strength_2195
fracture_toughness_2195 = 35*10**6
ult_compressive_strength_2195 = ult_yield_strength_2195


#al7075

#E_sheet = 71.7*10**9
#density_sheet = 2810
#ult_tensile_strength_2195 = 572*10**6
#ult_yield_strength_2195 = 503*10**6
#ult_shear_strength_2195 = 0.55*ult_yield_strength_2195
#fracture_toughness_2195 = 35*10**6
#ult_compressive_strength_2195 = ult_tensile_strength_2195
#
##CFRP

#E_sheet = 123.5*10**9
#density_sheet = 1600
#ult_tensile_strength_2195 = 600*10**6
#ult_yield_strength_2195 = (586.2-2*37.7)*10**6
#ult_shear_strength_2195 = 0.55*ult_yield_strength_2195
#fracture_toughness_2195 = 35*10**6
#ult_compressive_strength_2195 = (581-2*35.8)*10**6

#amount of ribs, excluding root and tip caps
n_ribs = 10

rib_spacing = (b/2)/(n_ribs+1)
t_rib = 0.0025



#al 2099-t83 http://morita1950.info/akio/data/Al-li%20Alloy.pdf
#https://www.smithmetal.com/2099-lithium.htm

E_rib = 78.6*10**9
density_rib = 2630
ult_shear_strength_2099 = 262*10**6

#CFRP
#E_rib = 70*10**9
#density_rib = 1600
#ult_shear_strength_2099 = 45.5*10**6


safety_factor_compression = 1.5
safety_factor_tension = 1.5


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
#d_strut = 5

# Empennage
l_tail = 4.539680359453441

# Vertical Tail
l_v = 11
A_v = 1.2
taper_v = 0.6
sweep_v = 30
V_v = 0.05
S_v = 0.27 * S
b_v = np.sqrt(S_v * A_v)
tip_chord_v = 2.5
root_chord_v = 4.16

# Horizontal Tail
l_h = 12
A_h = 4
taper_h = 0.5
sweep_h = 15
V_h = 1.57
S_h = 0.29 * S
b_h = np.sqrt(S_h * A_h)
S_e = 0.4
tip_chord_h = 1.26
root_chord_h = 2.53

# Undercarriage
main_landing_pos = 10.94
nose_landing_pos = 3
h_wheel = 0.79
lateral_pos = 1.94
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
M_cruise = 0.55
C_L_max_land = 2.6
C_L_max_TO = 1.8
C_L_cruise = 0.5
V_cruise = M_cruise*speed_of_sound
V_stall = 46.3
V_climb = 96 # Climb speed in m/s
V_descend = 122.7 # Descend speed in m/s
V_c = 169.3874 #Maximum Vc in m/s
t_cl = 0.2635 # Time to climb in hours
t_de = 0.357 # Time to descend in hours
C_L_max_land = 2.6
C_L_max_TO = 1.8
C_L_TO = C_L_max_TO / 1.21
V_TO = np.sqrt((2 * MTOW) / (rho0 * S * C_L_max_TO))
q_TO = 0.5 * rho0 * V_TO ** 2
C_D_TO = 0.023 #NOT FINAL
Cl_delta_aileron = 0.11217
Clp = -25.8276
Cl_alpha = 2 * np.pi
tau = 0.57
Cd0 = 0.027 #NOT FINAL
LD_ratio = 23.6 #NOT FINAL

f_cruise_start = 0.96542287875
f_cruise_end = 0.912778553222379

C_L_loiter_tbp = np.sqrt(3 * Cd0 * np.pi * A * e)
C_D_loiter_tbp = 4 * Cd0
LD_loiter = C_L_loiter_tbp / C_D_loiter_tbp



# Propulsion
eff_cruise = 0.85
eff_loiter = 0.77
P_TO = 2*1775e3
T_TO = 50.24e3
cp_cruise = 6.2e-8
cp_loiter = 6.5e-8
tot_thrust = 78e3

## -------- DOEKOES -------- ##
PP = 150000 #Cost per propeller  in 2019 dollar
EP = 1000000 # Engine cost per engine
#FP = 0.099 * 3# Fuel price per lbs in dollar
FP = 0.16 * 2
ASP = 5500000 #cost of avionics in 2019 dollar