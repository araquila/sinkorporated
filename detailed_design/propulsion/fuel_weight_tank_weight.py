import numpy as np
from sympy.solvers import solve
from sympy import Symbol

rho_CH4_cryo    = 425
allowance       = 0.07
safety_factor   = 1.1 #[-]


# aluminium
E_alu = 73.1e9
v_alu = 0.33
rho_alu = 2780
sigma_yield_alu = 500e6

# aramid
E_aramid = 95e9
v_aramid = 0.2
rho_aramid = 1440

# insulation
insulation = 0.04 # meters
rho_vapor_barier = 1
rho_microspheres = 69 # kg/m^3

# critical loads
P_cr_buckling = 1.5e5 # Pascals
P_cr_internal = 17.5e5 # Pascals

fuel_masses = [1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2050, 2100, 2150, 2200, 2250, 2300, 2350, 2400]

def det_tank_volume(LNG_mass, LNG_density, allowance):
    """
    tank volume from mass to m^3
    taken from hydrogen aircraft technologies p. 154
    input: LNG mass in kg
    output: tank volume in m^3
    """
    tank_volume = LNG_mass * (LNG_density/(1 + allowance))**-1
    return tank_volume

def t_long(radius,sigma_yield,p):
    t_long = (safety_factor*p*radius)/sigma_yield
    return t_long

def windenburg_trilling(E, v, P_cr):
    """
    Calculates the thickness acoording to windenburg and trilling equation for buckling
    """
    r_out = r + insulation
    t = Symbol("t")
    thickness = solve((2.42 * E * (0.5 *t /r_out)**(5/2)) / ((1-v**2)**0.75 * (0.5*l/r_out - 0.45 * (0.5 * t/r_out)**0.5)) - P_cr)
    return thickness[0]

def det_internal_dimensions(tank_volume, shell_length_ratio =  0.5, tank_circular_ratio = 1, tank_head_ratio = 3):
    a = Symbol("a")
    r = solve((4/3 * np.pi * a * tank_circular_ratio * a * tank_head_ratio * a + np.pi * a * tank_circular_ratio * a * (2 * a * tank_head_ratio / (1 - shell_length_ratio) - 2 * a * tank_head_ratio)) - tank_volume, a)
    radius = r[0]
    length = (radius * tank_head_ratio * 2)/(1 - shell_length_ratio)
    return radius, length

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

total_tank_mass_list = []

for item in fuel_masses:
    internal_volume = det_tank_volume(item/2, rho_CH4_cryo, allowance)
    r, l = det_internal_dimensions(internal_volume)
    thickness_outer_aramid = windenburg_trilling(E_aramid, v_aramid, P_cr_buckling)
    thickness_inner_alu = t_long(r, sigma_yield_alu, P_cr_internal)
    mass_inner_alu = (det_total_tank_volume(r, l, thickness_inner_alu) - det_total_tank_volume(r, l, 0)) * rho_alu
    mass_outer_aramid = (det_total_tank_volume(r, l, thickness_inner_alu + insulation + thickness_outer_aramid) - det_total_tank_volume(r, l, thickness_inner_alu + insulation)) * rho_aramid
    mass_vapor_barier = det_external_wetted_area(r + insulation, l + insulation, thickness_inner_alu + insulation) * rho_vapor_barier
    mass_microspheres_alu = (det_total_tank_volume(r, l, thickness_inner_alu + insulation) - det_total_tank_volume(r, l, thickness_inner_alu)) * rho_aramid
    total_mass = mass_inner_alu + mass_outer_aramid + mass_vapor_barier + mass_microspheres_alu
    total_tank_mass_list.append(total_mass)

print(fuel_masses)
print(total_tank_mass_list)
