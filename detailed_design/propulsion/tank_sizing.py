import sys
import os
from sympy.solvers import solve
from sympy import Symbol
import numpy as np
from p_t_data import *

# hydrogen aircraft technologies p. 30 (allowance) and p. 154 (tank volume)

sys.path.append(os.getcwd())

number_tanks    = 2      # number of tanks
weight_CH4      = 2120      # total weight of the methane
min_amb_pres    = None
allowance       = 0.07      # allowance of the tank in fraction of 1 (usually 7%)
rho_CH4_cryo    = 425      # density of methane at cryogenic temperatures
Lambda          = 511e3          # heat of vaporation in J/kg
insulation      = 0.04 # insulation in meters
kappa_eff       = 0.0032#0.011 # effecctive insulation parameter
kappa_2         = 0.032
thickness_2     = 0.012
timestep        = 30
T_outside       = 298

#weight Calculation
sigma_yield_steel     = 700e6 #[Pa]
sigma_yield_alu = 500e6 #[Pa] 7075-T651
sigma_ultimate_fibre = 600e6 #[Pa]
safety_factor   = 1.1 #[-]
p_a = 0.987*101325 #[Pa]
k = 7
n = 4
E_alu = 72e9 #[Pa] E-modulus of aluminium
E_steel = 200e9 #[Pa] E-modulus of steel
E_fibre = 70e9 #[Pa] E-modulus of carbon fibre
radius_outer = 0.5
rho_alu = 2800 #[kg/m3]
rho_aerogel = 130 #[kg/m3]
rho_last_layer = 25 #[kg/m3]

#Functions
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
    output: tank volume in m^3rho_steel = 8050 #[kg/m3]
rho_fibre = 1550 #[kg/m3]

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
    BOR = (Q_leak * timestep)/(LNG_density * V_LNG * Lambda)
    return BOR


# print(det_BOR(250, rho_CH4_cryo, 2, Lambda))

def det_internal_dimensions(tank_volume, shell_length_ratio =  0.5, tank_circular_ratio = 1, tank_head_ratio = 3):
    a = Symbol("a")
    r = solve((4/3 * np.pi * a * tank_circular_ratio * a * tank_head_ratio * a + np.pi * a * tank_circular_ratio * a * (2 * a * tank_head_ratio / (1 - shell_length_ratio) - 2 * a * tank_head_ratio)) - tank_volume, a)
    radius = r[0]
    length = (radius * tank_head_ratio * 2)/(1 - shell_length_ratio)
    return radius, length

print(det_internal_dimensions(2.67))

def det_total_tank_volume(radius, length, thickness_insulation, shell_length_ratio =  0.5, tank_circular_ratio = 1, tank_head_ratio = 3):
    ans = length * shell_length_ratio * (radius + thickness_insulation)**2 *np.pi
    ans += 4/3 * np.pi * (length * shell_length_ratio)/2 * (radius + thickness_insulation)**2
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
    """
    Compute the total heat flux in to one tank in Watt
    """
#    flux = kappa * delta_t / thickness_insulation
    flux = delta_t/(thickness_insulation/kappa+thickness_2/kappa_2)
    return area * flux

def benedict_webb_rubin(T, density):
    """
    Calculates the pressure according to benedict_webb_rubin equation of state
    valid for 2.5 times the critical density
    critical density = 162.7 kg/m^3
    """
    R_u = 8314.462 # universal gas constant
    molarmass = 16.043 # molar mass of Methane
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

    timedata = []
    steps = int(total_t / timestep)
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
            boil_off_mass = liquid_mass * det_BOR(det_heat_flux(kappa_eff, (T_outside - temperature), insulation, det_external_wetted_area(0.4, 4.8, insulation)), rho_CH4_cryo, liquid_mass/rho_CH4_cryo, Lambda)
            liquid_mass -= boil_off_mass
            gas_mass += boil_off_mass
            V_gas += boil_off_mass / rho_CH4_cryo
            density = gas_mass / V_gas
        T_vaporise = det_T_vaporise(pressure[-1]/1000)
        temperature += det_heat_flux(kappa_eff, (T_outside - temperature), insulation, det_external_wetted_area(0.4, 4.8, insulation)) * timestep / (2 * 1000 * liquid_mass)
        timedata.append([t,pressure[-1],temperature])
#        print("Timestep (" + str(timestep) + " sec): " + str(t) + " pressure " + str(pressure[-1]) + " temperature " + str(temperature))
    return pressure, timedata

#thickness functions
def t_hoop(radius,sigma_yield,p):
    t_hoop = (safety_factor*p*radius)/(2*sigma_yield)
    return t_hoop

def t_long(radius,sigma_yield,p):
    t_long = (safety_factor*p*radius)/sigma_yield
    return t_long

def t_dewar(radius_outer,p_a,n,E,k):
    t_dewar = 2*radius_outer*((p_a*n)/(E*k))**(1./3)
    return t_dewar

pressure, timedata = pressure_build_up(432000)
print("max pressure 5 days", max(pressure))
LNG_mass = weight_CH4/2
tank_volume = det_tank_volume(LNG_mass,rho_CH4_cryo,allowance)
radius, length = det_internal_dimensions(tank_volume)

#t_hoop_steel = t_hoop(0.4,sigma_yield_steel,max(pressure))
# t_long_steel = t_long(0.4, sigma_yield_steel, max(pressure))
# #t_hoop_alu = t_hoop(0.4,sigma_yield_alu,max(pressure))
# t_long_alu = t_long(0.4,sigma_yield_alu, max(pressure))
# #t_hoop_fibre = t_hoop(0.4,sigma_ultimate_fibre,max(pressure))
# t_long_fibre = t_long(0.4,sigma_ultimate_fibre, max(pressure))

t_dewar_alu = t_dewar(radius_outer,p_a,n,E_alu,k)
t_dewar_steel = t_dewar(radius_outer,p_a,n,E_steel,k)
t_dewar_fibre = t_dewar(radius_outer,p_a,n,E_fibre,k)
# print(t_long_fibre)
# print(t_long_steel)
# print(t_long_alu)
# print("Thickness dewar steel" +str(t_dewar_steel))
# print("internal skin mass of a carbon composite tank is " + str(det_external_wetted_area(0.4, 4.8, t_long_fibre) * t_long_fibre * rho_fibre))
# print("the weight of the aluminium carbon tank is " + str((rho_fibre+0.25*rho_alu)*(det_external_wetted_area(0.4, 4.8, t_long_fibre) * t_long_fibre)) + " kg")
# print("internal skin mass of a aluminium tank is " + str(det_external_wetted_area(0.4, 4.8, t_long_alu) * t_long_alu * rho_alu))
# print("external skin mass of an aluminium tank is " + str(det_external_wetted_area(0.4, 4.8, insulation) * t_dewar_alu * rho_alu))
# print("external skin mass of a steel tank is "+ str(det_external_wetted_area(0.4, 4.8, insulation) * t_dewar_steel * rho_steel))
# print("external skin mass of a fibre tank is "+ str(det_external_wetted_area(0.4, 4.8, insulation) * t_dewar_fibre * rho_fibre))

# inner_volume = det_total_tank_volume(radius,length,1.25*t_long_fibre)
# outer_volume = det_total_tank_volume(radius,length,1.25*t_long_fibre+insulation)
# outer_outer_volume = det_total_tank_volume(radius,length,1.25*t_long_fibre+insulation+thickness_2)
# outer_alu_volume = det_total_tank_volume(radius,length,1.25*t_long_fibre+insulation+thickness_2+0.0001)
# mass_last_layer = outer_outer_volume*rho_last_layer
# aerogel_mass = (outer_volume-inner_volume)*rho_aerogel
# outer_alu_layer_mass = (outer_alu_volume-outer_outer_volume)*rho_alu
# inner_tank_mass = (rho_fibre+0.25*rho_alu)*(det_external_wetted_area(0.4, 4.8, t_long_fibre) * t_long_fibre)
# total_mass = inner_tank_mass + aerogel_mass + mass_last_layer + outer_alu_layer_mass
# total_diameter = 2*(radius+1.25*t_long_fibre+insulation+thickness_2+0.0001)
# print(total_mass,total_diameter)
print(timedata[960])

#dewar calculations
# inner_volume2 = det_total_tank_volume(radius,length,1.25*t_long_fibre)
# outer_volume2 = det_total_tank_volume(radius,length,1.25*t_long_fibre+insulation)
# total_diameter2 = 2*(radius+1.25*t_long_fibre+insulation+t_dewar_fibre)
# total_mass2 = (rho_fibre+0.25*rho_alu)*(det_external_wetted_area(0.4, 4.8, t_long_fibre) * t_long_fibre)+det_external_wetted_area(0.4, 4.8, insulation) * t_dewar_fibre * rho_fibre
# print(total_mass2,total_diameter2)
# print(benedict_webb_rubin(111, 1.8))
#
# print("External Wetted Area of One Tank: " + str(det_external_wetted_area(0.4, 4.8, insulation)))
# print("Total Heat Flux of One Tank: " + str(det_heat_flux(kappa_eff, 224, insulation, det_external_wetted_area(0.4, 4.8, insulation))))
# print("Boil-Off Rate: " + str(det_BOR(det_heat_flux(kappa_eff, 224, insulation, det_external_wetted_area(0.4, 4.8, insulation)), rho_CH4_cryo, 2, Lambda)))
