import parameters as p
import numpy as np
from matplotlib import pyplot as plt
from wing.shear_and_moment_strutbox import shear_and_moment
import wing.section_properties_strutbox as scs
#from wing_deflection2 import 

# Obtain strut force, reaction forces and reaction moment




def normal_stress(x,y,moment_z,moment_y,normal_force,I_zz,I_yy):
    """Returns bending moment and normal force stress""" 
    moment_z = 10000
    moment_y = 10000
    
    normal_force = 1000 #positive along right wing span
        
    I_zz = sp.I_zz_wingbox(x)
    I_yy = sp.I_yy_wingbox(x)
    
    y = sp.y_max(x)
    z = sp.width_wingbox(x)/2
    
    #-my/I according to the formula, makes sense because for a positive Mz the top skin will be in compression
    moment_z_upperskin = -moment_z*(sp.height_wingbox(x)-sp.y_max(x))/I_zz
    moment_z_lowerskin = -moment_z*sp.y_max(x)/I_zz
    
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

x = np.linspace(0,p.b/2,100)
normal_stress_values = np.zeros(len(x))
for i in range(len(x)):
    normal_stress_values[i] = normal_stress(x[i])
    
plt.plot(x,normal_stress_values)
plt.show()
    