from constant_variables import *
import numpy as np
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

def prop_noise(D,B,n_p,P_br,N_p,c):
    # D: propeller diameter(m)
    # B: number of blades per propeller
    # n_p: propeller rotational speed (rpm)
    # P_br engine power (kW)
    # N_p number of propellers
    M_t = np.pi*D*n_p/(c*60)
    SPL_max = 83.4 + 15.3*np.log10(P_br)-20*np.log10(D)+38.5*M_t-3*(B-2)+10*np.log10(N_p)
    return SPL_max

def airframe_noise(V_cruise,W):
    SPL_airframe = 10*np.log10(V_cruise**5)+10*np.log10(W)-74
    return SPL_airframe
