# Gravitional constant
g = 9.80665
R = 287

# Atmospherical parameters
temperature0 = 288.15
temperature_gradient = -0.0065
gamma = 1.4
rho0 = 1.225                        # [kg/m3]

# Passengers and crew
n_passenger = 60
M_passenger = 105                   #(including luggage)
n_crew = 4
M_crew_member = 100

# Cabin layout
n_seats_abreast = 4
n_aisles = 1

# Initial mass and fractions
M_payload = n_passenger * M_passenger
M_crew = n_crew * M_crew_member
f_trapped_fuel = 0.003              # Range 0.001-0.005
M_empty_tbp = 13311
M_empty_jet = 16300

# Convert to weights
W_payload = M_payload * g
W_crew = M_crew * g
W_empty_tbp = M_empty_tbp * g
W_empty_jet = M_empty_jet * g

n_engines_jet = 2
n_engines_tbp = 2

# Fuels
energy_density_kerosene = 43 #[MJ/kg]
energy_density_LNG = 53.6 #[MJ/kg]
energy_density_HHV = 141.6 #[MJ/kg]
energy_density_LHV = 119.93 #[MJ/kg]
