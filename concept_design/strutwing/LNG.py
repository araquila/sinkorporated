from sympy.solvers import solve
from sympy import Symbol
import numpy as np

# https://www.sciencedirect.com/topics/engineering/cryogenic-tank
# https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7564027
# https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/7815389/references#references



def LNG_volume(LNG_mass):
    """
    simple converter from mass to m^3
    input LNG mass in kg
    output LNG volume in m^3
    """
    LCH4_density = 425 # kg/m^3
    return LNG_mass/LCH4_density


def LNG_system(n_tanks, required_volume, podded = True, shell_length_ratio =  0.5, tank_circular_ratio = 1, tank_head_ratio = 3):
    """
    inputs:

    conditional inputs:

    outputs:

    sources:
    Modelling and Designing Cryogenic Hydrogen Tanks
    for Future Aircraft Applications
    """
    LCH4_density = 425 # kg/m^3
    LCH4_energy_density = 55.6 * 10**6 # Joule per kg

    max_liquid_volume = 0.97 # 3% is needed for venting
    max_liquid_volume = 0.92 # more is needed for venting at higher venting_pressure

    # volume = n_tanks * (4/3 * np.pi * a * tank_circular_ratio * a * tank_head_ratio * a + np.pi * a * tank_circular_ratio * a * (2 * a * tank_head_ratio / (1 - shell_length_ratio) - 2 * a * tank_head_ratio))
    a = Symbol("a")
    r = solve(n_tanks * (4/3 * np.pi * a * tank_circular_ratio * a * tank_head_ratio * a + np.pi * a * tank_circular_ratio * a * (2 * a * tank_head_ratio / (1 - shell_length_ratio) - 2 * a * tank_head_ratio)) - required_volume/max_liquid_volume, a)
    length_LNG_tank = (r[0] * tank_head_ratio * 2)/(1 - shell_length_ratio)
    mass_tank = 1.17 * required_volume/1.5 * 100
    mass_cryo_cooler = 150
    mass_LNG_system = mass_tank + mass_cryo_cooler

    return mass_LNG_system


# print(LNG_volume(1740))
# print(LNG_system(2, LNG_volume(1740)))
