### IMPORTS
from atmosphere import atmosphere_calc

### AIRCRAFT PARAMETERS ###

## CONSTANTS
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

## WEIGHTS AND MASSES
# Propulsion
mass_engine = 450                    
                  
# Fuselage

           

## DIMENSIONS
# Wing
A = 18
S = 66

# Empennage
V_h = 1.57                         
V_v = 0.07                        
l_h = 11                           
l_v = 11 

# Undercarriage
main_landing_pos = 11             
nose_landing_pos = 3   

# Propulsion
pos_engine = 10                    
n_fueltanks = 2                     
n_blades = 6 

## PERFORMANCE
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