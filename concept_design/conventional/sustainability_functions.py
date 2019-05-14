from constant_variables import *
def CO2_calc(fuel_per_passenger,chosen_fuel_energy_density):
    if chosen_fuel_energy_density == energy_density_kerosene:
        CO2_per_passenger = 3.00*fuel_per_passenger
    if chosen_fuel_energy_density == energy_density_LNG:
        CO2_per_passenger = 2.75*fuel_per_passenger
    return CO2_per_passenger
