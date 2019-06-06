import parameters as p
import numpy as np
from matplotlib import pyplot as plt
from wing.shear_and_moment_strutbox import shear_and_moment
import wing.section_properties_strutbox as scs
#from wing_deflection2 import 

# Obtain strut force, reaction forces and reaction moment
F_strut = 5
F_strut_x = 
F_



def normal_stress(M_list,x_pos, F_strut_x):
    normal_stress_list = []
    for i in len(x_pos):
        sigma = -(M_list[i] * -p.h_max_root_strutbox) / scs.I_zz_strutbox(x_pos[i]) + 2 * 
        normal_stress_list.append(sigma)    
    
def shear_stress(V_list,x_pos):
    