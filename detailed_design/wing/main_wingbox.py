# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:22:07 2019

@author: robert
"""

import parameters as p
import numpy as np
import matplotlib.pyplot as plt
import wing.section_properties as sp

def normal_stress(x):
    """Returns bending moment and normal force stress""" 
    moment_z = 10000
    moment_y = 10000
    
    normal_force = 1000 #positive along right wing span
        
    I_zz = sp.I_zz_wingbox(x)
    I_yy = sp.I_yy_wingbox(x)
    
    y = sp.y_max(x)
    z = sp.width_wingbox(x)/2
    
    #-my/I according to the formula, makes sense because for a positive Mz the top skin will be in compression
    moment_z_upperskin = -moment_z*(sp.height_wingbox(x)-y)/I_zz
    moment_z_lowerskin = -moment_z*y/I_zz
    
    moment_y_rightflange = -moment_y*z/I_yy
    moment_y_leftflange = -moment_y*-z/I_yy
    
    area = sp.cross_sectional_area(x)
    
    normal_force_stress = normal_force/area
    
    normal_ru = moment_z_upperskin + moment_y_rightflange +     normal_force_stress 
    normal_lu = moment_z_upperskin + moment_y_leftflange +     normal_force_stress 
    normal_rl = moment_z_lowerskin + moment_y_rightflange +     normal_force_stress 
    normal_ll = moment_z_lowerskin + moment_y_leftflange +  normal_force_stress 
    
    if max(normal_ru,normal_lu,normal_rl,normal_ll) == normal_ru:
        print("Max tension at right upper corner")
    elif max(normal_ru,normal_lu,normal_rl,normal_ll) == normal_lu:
        print("Max tension at left upper corner")
    elif max(normal_ru,normal_lu,normal_rl,normal_ll) == normal_rl:
        print("Max tension at right lower corner")
    elif max(normal_ru,normal_lu,normal_rl,normal_ll) == normal_ll:
        print("Max tension at left lower corner")
        
    
    if min(normal_ru,normal_lu,normal_rl,normal_ll) == normal_ru:
        print("Max compression at right upper corner")
    elif min(normal_ru,normal_lu,normal_rl,normal_ll) == normal_lu:
        print("Max compression at left upper corner")
    elif min(normal_ru,normal_lu,normal_rl,normal_ll) == normal_rl:
        print("Max compression at right lower corner")
    elif min(normal_ru,normal_lu,normal_rl,normal_ll) == normal_ll:
        print("Max compression at left lower corner")
        
        
    print("x: ",x)
    print("Neutral axis: y =",sp.centroid_y(x),"(",sp.centroid_y(x)/sp.height_wingbox(x)*100,"%)",)
    print("Maximum tension: ",max(normal_ru,normal_lu,normal_rl,normal_ll)/10**6,"MPa")
    print("Maximum compression: ", min(normal_ru,normal_lu,normal_rl,normal_ll)/10**6,"MPa")
    print("")
    
    return min(normal_ru,normal_lu,normal_rl,normal_ll)



def skin_buckling_stress(x):
    """Returns compressive stress at which buckling will occur"""
    
    k = None # to be determined based on the "final" stiffener and rib spacing
    b = sp.top_spacing
    poisson_ratio = p.poisson_ratio
    t = p.t_sheet
    
    return k*np.pi**2*p.E*(t/b)/(12*(1-poisson_ratio)**2)


def critical_colunn_buckling(x):
    """Critical column buckling force on the stiffener""" 
    
    l_eff = 0.5*p.rib_spacing
    
    return (np.pi**2*p.E*sp.I_zz_wingbox(x))/(l_eff**2)


#def critical_crippling_stiffener(x):
#    """Critical crippling top stiffener"""
#    
    
    

x = np.linspace(0,p.b/2,100)
normal_stress_values = np.zeros(len(x))
for i in range(len(x)):
    normal_stress_values[i] = normal_stress(x[i])
    
    
plt.plot(x,normal_stress_values)
plt.show()


    
    