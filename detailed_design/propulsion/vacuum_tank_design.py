from sympy.solvers import solve
from sympy import Symbol
import numpy as np
from tank_sizing import *

E_steel = 200e9
v_steel = 0.3
rho_steel = 8050
sigma_yield_steel = 700e6

E_alu = 73.1e9
v_alu = 0.33
rho_alu = 2780
sigma_yield_alu = 500e6

E_carbon = 70e9
v_carbon = 0.10
rho_carbon = 1550

E_aramid = 95e9
v_aramid = 0.2
rho_aramid = 1440

rho_vapor_barier = 1

r = 0.44
l = 5.28
insulation = 0.04 # meters

P_cr_buckling = 1.5e5 # Pascals
P_cr_internal = 17.5e5 # Pascals



l_b = 0.833456 # buckling length


def windenburg_trilling(E, v, P_cr):
    """
    Calculates the thickness acoording to windenburg and trilling equation for buckling
    """
    r_out = r + insulation
    t = Symbol("t")
    thickness = solve((2.42 * E * (0.5 *t /r_out)**(5/2)) / ((1-v**2)**0.75 * (0.5*l/r_out - 0.45 * (0.5 * t/r_out)**0.5)) - P_cr)
    return thickness[0]

def t_long(radius, sigma_yield, p):
    t_long = (safety_factor * p * radius) / sigma_yield
    return t_long


def rahim_buckling(E, thickness, i, l_b, v):
    P_cr = (E * thickness / r) / (1 + 0.5 * (np.pi * r / (i * l_b))**2) * ((1 / (i**2 * (1 + (i * l_b / np.pi / r)**2)) + i**2 * thickness**2 / (12 * r**2 * (1 - v**2)) * (1 + (np.pi * r /i / l_b)**2)**2))
    return P_cr


# trade-off for different materials
# Steel
# thickness_outer_steel = windenburg_trilling(E_steel, v_steel, P_cr_buckling)
# thickness_inner_steel = t_long(r, sigma_yield_steel, P_cr_internal)
# mass_inner_steel = (det_total_tank_volume(r, l, thickness_inner_steel) - det_total_tank_volume(r, l, 0)) * rho_steel
# mass_outer_steel = (det_total_tank_volume(r, l, thickness_inner_steel + insulation + thickness_outer_steel) - det_total_tank_volume(r, l, thickness_inner_steel + insulation)) * rho_steel
# mass_microspheres_steel = (det_total_tank_volume(r, l, thickness_inner_steel + insulation) - det_total_tank_volume(r, l, thickness_inner_steel)) * 69
# print("internal thickness steel ", thickness_inner_steel)
# print("external thickness steel ", thickness_outer_steel)
# print("internal mass ", mass_inner_steel)
# print("external mass ", mass_outer_steel)
# print("mass microspheres", mass_microspheres_steel)
# print("total mass ", mass_outer_steel + mass_inner_steel + mass_microspheres_steel)

# aluminium
# thickness_outer_alu = windenburg_trilling(E_alu, v_alu, P_cr_buckling)
# thickness_outer_carbon = windenburg_trilling(E_carbon, v_carbon, P_cr_buckling)
# thickness_inner_alu = t_long(r, sigma_yield_alu, P_cr_internal)
# mass_inner_alu = (det_total_tank_volume(r, l, thickness_inner_alu) - det_total_tank_volume(r, l, 0)) * rho_alu
# # mass_outer_alu = (det_total_tank_volume(r, l, thickness_inner_alu + insulation + thickness_outer_alu) - det_total_tank_volume(r, l, thickness_inner_alu + insulation)) * rho_alu
# mass_outer_carbon = (det_total_tank_volume(r, l, thickness_inner_alu + insulation + thickness_outer_carbon) - det_total_tank_volume(r, l, thickness_inner_alu + insulation)) * rho_carbon
# mass_vapor_barier = det_external_wetted_area(0.4, 4.8, thickness_inner_alu + insulation) * rho_vapor_barier
# mass_microspheres_alu = (det_total_tank_volume(r, l, thickness_inner_alu + insulation) - det_total_tank_volume(r, l, thickness_inner_alu)) * 69
# print("internal thickness aluminium ", thickness_inner_alu)
# print("external thickness carbon ", thickness_outer_carbon)
# print("internal mass ", mass_inner_alu)
# print("external mass ", mass_outer_carbon)
# print("mass microspheres", mass_microspheres_alu)
# print("mass vapor barier", mass_vapor_barier)
# print("total mass aluminium-carbon", mass_outer_carbon + mass_inner_alu + mass_microspheres_alu + mass_vapor_barier)


# aluminium internal shell and aramid external shell
thickness_outer_aramid = windenburg_trilling(E_aramid, v_aramid, P_cr_buckling)
thickness_inner_alu = t_long(r, sigma_yield_alu, P_cr_internal)
mass_inner_alu = (det_total_tank_volume(r, l, thickness_inner_alu) - det_total_tank_volume(r, l, 0)) * rho_alu
mass_outer_aramid = (det_total_tank_volume(r, l, thickness_inner_alu + insulation + thickness_outer_aramid) - det_total_tank_volume(r, l, thickness_inner_alu + insulation)) * rho_aramid
mass_vapor_barier = det_external_wetted_area(r + insulation, l + insulation, thickness_inner_alu + insulation) * rho_vapor_barier
mass_microspheres_alu = (det_total_tank_volume(r, l, thickness_inner_alu + insulation) - det_total_tank_volume(r, l, thickness_inner_alu)) * 69
print("internal thickness aluminium ", thickness_inner_alu)
print("external thickness aramid ", thickness_outer_aramid)
print("internal mass ", mass_inner_alu)
print("external mass ", mass_outer_aramid)
print("mass microspheres", mass_microspheres_alu)
print("mass vapor barier", mass_vapor_barier)
print("total mass aluminium-aramid", mass_outer_aramid + mass_inner_alu + mass_microspheres_alu + mass_vapor_barier)
