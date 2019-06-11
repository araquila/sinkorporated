from sympy.solvers import solve
from sympy import Symbol
import numpy as np

def det_internal_dimensions(tank_volume, shell_length_ratio =  0.5, tank_circular_ratio = 1, tank_head_ratio = 3):
    a = Symbol("a")
    r = solve((4/3 * np.pi * a * tank_circular_ratio * a * tank_head_ratio * a + np.pi * a * tank_circular_ratio * a * (2 * a * tank_head_ratio / (1 - shell_length_ratio) - 2 * a * tank_head_ratio)) - tank_volume, a)
    radius = r[0]
    length = (radius * tank_head_ratio * 2)/(1 - shell_length_ratio)
    return radius, length

E_steel = 69*10**9
r = 0.4
l = 4.8
v_steel = 0.3
P_cr = 1.5e5


def windenburg_trilling():
    t = Symbol("t")
    thickness = solve((2.42 * E_steel * (0.5 *t /r)**(5/2)) / ((1-v_steel**2)**0.75 * (0.5*l/r - 0.45 * (0.5 * t/r)**0.5)) - P_cr)
    print(thickness)
    return thickness[0]

thickness = 0.003
l_b = 0.812333

for i in range(15):
    if i > 1:
        P_cr = (E_steel * thickness / r) / (1 + 0.5 * (np.pi * r / (i * l_b))**2) * ((1 / (i**2 * (1 + (i * l_b / np.pi / r)**2)) + i**2 * thickness**2 / (12 * r**2 * (1 - v_steel**2)) * (1 + (np.pi * r /i / l_b)**2)**2))
        print(str(i) + " : " + str(P_cr))
