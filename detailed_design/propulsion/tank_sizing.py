import sys
import os
from sympy.solvers import solve
from sympy import Symbol
import numpy as np
from p_t_data import *

# hydrogen aircraft technologies p. 30 (allowance) and p. 154 (tank volume)

sys.path.append(os.getcwd())

number_tanks    = 2      # number of tanks
weight_CH4      = 1600      # total weight of the methane
min_amb_pres    = None
allowance       = 0.07      # allowance of the tank in fraction of 1 (usually 7%)
rho_CH4_cryo    = 425      # density of methane at cryogenic temperatures
Lambda          = 511e3          # heat of vaporation in J/kg
insulation      = 0.05 # insulation in meters
kappa_eff       = 0.003588 # effecctive insulation parameter

def LNG_volume(LNG_mass):
    """
    simple converter from mass to m^3
    input: LNG mass in kg
    output: LNG volume in m^3
    """
    return LNG_mass/rho_CH4_cryo

def det_tank_volume(LNG_mass, LNG_density, allowance):
    """
    tank volume from mass to m^3
    taken from hydrogen aircraft technologies p. 154
    input: LNG mass in kg
    output: tank volume in m^3
    """
    tank_volume = LNG_mass * (LNG_density/(1 + allowance))**-1
    return tank_volume

print(det_tank_volume(weight_CH4/number_tanks, rho_CH4_cryo, allowance))
internal_tank_volume = det_tank_volume(weight_CH4/number_tanks, rho_CH4_cryo, allowance)

def det_BOR(Q_leak, LNG_density, V_LNG, Lambda):
    """
    the boil off rate in a fraction
    input: heat input in W
    output:
    """
    BOR = (Q_leak * 60)/(LNG_density * V_LNG * Lambda)
    return BOR


# print(det_BOR(250, rho_CH4_cryo, 2, Lambda))

def det_internal_dimensions(tank_volume, shell_length_ratio =  0.5, tank_circular_ratio = 1, tank_head_ratio = 3):
    a = Symbol("a")
    r = solve((4/3 * np.pi * a * tank_circular_ratio * a * tank_head_ratio * a + np.pi * a * tank_circular_ratio * a * (2 * a * tank_head_ratio / (1 - shell_length_ratio) - 2 * a * tank_head_ratio)) - tank_volume, a)
    radius = r[0]
    length = (radius * tank_head_ratio * 2)/(1 - shell_length_ratio)
    return radius, length


def det_total_tank_volume(radius, length, thickness_insulation, shell_length_ratio =  0.5, tank_circular_ratio = 1, tank_head_ratio = 3):
    ans = length * shell_length_ratio * radius**2 *np.pi
    # ans += 4/3 * np.pi * (radius * tank_head_ratio) * radius**2
    ans += 4/3 * np.pi * (length * shell_length_ratio)/2 * radius**2
    return ans

def det_external_wetted_area(radius, length, thickness_insulation, shell_length_ratio =  0.5, tank_circular_ratio = 1, tank_head_ratio = 3):
    radius2 = radius + thickness_insulation
    ans = length * shell_length_ratio * radius2 * 2 *np.pi
    if tank_head_ratio == 1:
        ans += 1000
    if tank_head_ratio > 1:
        # prolate spheroid
        ecc = np.sqrt(1 -  radius2**2 / (tank_head_ratio * radius2)**2)
        ans += 2 * np.pi * radius2**2 * (1 + tank_head_ratio / ecc * np.arcsin(ecc))
    if tank_head_ratio < 1:
        # oblate spheroid
        ecc = np.sqrt(1 - (tank_head_ratio * radius2)**2 / radius2**2)
        ans += 2 * np.pi * radius2**2 * (1 + (1 - ecc**2)/ecc * np.arctanh(ecc))
    return ans

def det_heat_flux(kappa, delta_t, thickness_insulation, area):
    flux = kappa * delta_t / thickness_insulation
    return area * flux

def benedict_webb_rubin(T, density):
    R_u = 8314.462 # universal gas constant
    molarmass = 16.043 # molar mass of Methane
    criticaldensity = 162.7 # kg/m^3
    v_bar = molarmass / density
    B_0 = 0.04260
    A_0 = 187.91
    C_0 = 2.286e6
    b = 0.00380
    a = 5.00
    alpha = 1.244e-4
    c = 2.578e5
    gamma = 0.0060
    P = R_u * T / v_bar + (B_0 * R_u * T - A_0 - C_0 / T**2)/v_bar**2 + (b * R_u * T - a)/v_bar**3 + a * alpha / v_bar**6 + c / (v_bar**3 * T**2) * (1 + gamma/v_bar**2) * np.exp(-gamma/v_bar**2)
    return P

def pressure_build_up(total_t):


    steps = int(total_t / 60)
    pressure = []
    # initial conditions
    density = 1.8
    T_vaporise = 111
    V_gas = internal_tank_volume - internal_tank_volume/ (1 + allowance)
    liquid_mass = weight_CH4/number_tanks
    gas_mass = V_gas * density
    temperature = 111
    for t in range(steps):
        if temperature >= T_vaporise:
            pressure.append(benedict_webb_rubin(temperature, density))
            boil_off_mass = liquid_mass * det_BOR(det_heat_flux(kappa_eff, 224, insulation, det_external_wetted_area(0.4, 4.8, insulation)), rho_CH4_cryo, liquid_mass/rho_CH4_cryo, Lambda)
            liquid_mass -= boil_off_mass
            gas_mass += boil_off_mass
            V_gas += boil_off_mass / rho_CH4_cryo
            density = gas_mass / V_gas
            print(pressure)
            print("vaporise, timestep " +str(t) + "density "+ str(density))
        T_vaporise = det_T_vaporise(pressure[-1]/1000)
        temperature += det_heat_flux(kappa_eff, 224, insulation, det_external_wetted_area(0.4, 4.8, insulation)) * 60 / (2 * 1000 * liquid_mass)
    return pressure

print(pressure_build_up(432000))
# print(benedict_webb_rubin(111, 1.8))
#
# print("External Wetted Area of One Tank: " + str(det_external_wetted_area(0.4, 4.8, insulation)))
# print("Total Heat Flux of One Tank: " + str(det_heat_flux(kappa_eff, 224, insulation, det_external_wetted_area(0.4, 4.8, insulation))))
# print("Boil-Off Rate: " + str(det_BOR(det_heat_flux(kappa_eff, 224, insulation, det_external_wetted_area(0.4, 4.8, insulation)), rho_CH4_cryo, 2, Lambda)))
