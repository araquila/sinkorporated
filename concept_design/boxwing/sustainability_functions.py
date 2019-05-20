from constant_variables import *
import numpy as np
def CO2_calc(fuel_per_passenger,chosen_fuel_energy_density):
    if chosen_fuel_energy_density == energy_density_kerosene:
        CO2_per_passenger = 3.00*fuel_per_passenger
    if chosen_fuel_energy_density == energy_density_LNG:
        CO2_per_passenger = 2.75*fuel_per_passenger
    else:
        CO2_per_passenger = 0
    return CO2_per_passenger

def prop_noise(D,B,n_p,P_br,N_p,c):
    # D: propeller diameter(m)
    # B: number of blades per propeller
    # n_p: propeller rotational speed (rpm)
    # P_br engine power (kW)
    # N_p number of propellers
    M_t = np.pi*D*n_p/(c*60)
    SPL_max = 83.4 + 15.3*np.log10(P_br)-20*np.log10(D)+38.5*M_t-3*(B-2)+10*np.log10(N_p)
    return SPL_max

def turbofan_noise():
    SPL_turbofan = 100 #taken from average data
    return SPL_turbofan

def airframe_noise(V_cruise,W):
    SPL_airframe = 10*np.log10(V_cruise**5)+10*np.log10(W)-74
    return SPL_airframe

def total_noise(engine_noise,airframe_noise):
    total_noise = 10*np.log10(10**(engine_noise/10)+10**(airframe_noise/10))
    return total_noise

def noise_distance(noise_level,r1,r2):
    noise_at_distance = noise_level+20*np.log10(r1/r2)
    return noise_at_distance
