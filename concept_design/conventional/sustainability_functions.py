from constant_variables import *
def CO2_calc(fuel_per_passenger,chosen_fuel_energy_density):
    if chosen_fuel_energy_density == energy_density_kerosene:
        CO2_per_passenger = 3.00*fuel_per_passenger
    if chosen_fuel_energy_density == energy_density_LNG:
        CO2_per_passenger = 2.75*fuel_per_passenger
    if chosen_fuel_energy_density == energy_density_HHV:
        CO2_per_passenger = 0
    if chosen_fuel_energy_density == energy_density_LHV:
        CO2_per_passenger = 0
    return CO2_per_passenger

def NOX_calc(fuel_per_passenger,chosen_fuel_energy_density):
    if chosen_fuel_energy_density == energy_density_kerosene:
        NOX_per_passenger = 3.00*10e-3*fuel_per_passenger
    if chosen_fuel_energy_density == energy_density_LNG:
        NOX_per_passenger =10e-3*fuel_per_passenger
    if chosen_fuel_energy_density == energy_density_HHV:
        NOX_per_passenger = 0
    if chosen_fuel_energy_density == energy_density_LHV:
        NOX_per_passenger = 0
    return NOX_per_passenger
