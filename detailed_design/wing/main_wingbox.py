# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:22:07 2019

@author: robert
"""

import parameters as p
import numpy as np
import matplotlib.pyplot as plt
import wing.section_properties as sp
import wing.wing_deflection2 as wd2 


### DISCRETIZATION OF THE WINGBOX ###
n = 51
x_pos = np.linspace(0,p.b/2,n)


### OBTAIN CROSSECTIONAL PROPERTIES ###
Izz_list = []
Iyy_list = []
first_moment_of_area_list = []
area_list = []
y_max_list = []

for x in x_pos:
    Izz_list.append(sp.I_zz_wingbox(x))
    Iyy_list.append(sp.I_yy_wingbox(x))
    first_moment_of_area_list.append(sp.first_moment_of_area(x))
    area_list.append(sp.cross_sectional_area(x))
    y_max_list.append(sp.y_max(x))

### OBTAIN STRUT FORCE, REACTION FORCES AND REACTION MOMENT ###
lengthdata = 50
Lift, Chord, Yle, Drag = wd2.read_aero_data("wing/aquiladata1.txt", lengthdata, p.V_cruise, p.rho)
Frx, Fry, Fs, Mrz, Frz, Fsz, Mry, momentyi, momentzi, shearyi, shearzi, vyi, vny, vzi, vnz, xi, theta = wd2.CallForces(Lift, Yle, Drag, p.tot_thrust, Iyy_list, Izz_list,70*10**9, p.engine_pos_perc, p.strut_pos_perc, p.pod_pos_perc)
alpha = np.arctan((p.strut_pos_perc * p.b/2)/p.d_fuselage_outside)        # Angle of the strut with fuselage
F_strut_y = Fs * np.cos(alpha)
F_strut_x = Fs * np.sin(alpha)

### NORMAL STRESS CALCULATOR ###
def normal_stress(x,y,moment_z,moment_y,normal_force,I_zz,I_yy,area):
    """Returns bending moment and normal force stress""" 
    
    z = sp.width_wingbox(x)/2
    
    #-my/I according to the formula, makes sense because for a positive Mz the top skin will be in compression
    moment_z_upperskin = -moment_z*(sp.height_wingbox(x)-abs(y))/I_zz
    moment_z_lowerskin = -moment_z*y/I_zz
    
    
    moment_y_rightflange = -moment_y*z/I_yy
    moment_y_leftflange = -moment_y*-z/I_yy
    
    print(moment_y_leftflange,moment_y_rightflange,moment_z_lowerskin,moment_z_upperskin)
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
    print("Normal force stress: ",normal_force_stress,"MPa")
    print("")
    
    return normal_ru,normal_lu,normal_rl,normal_ll



def skin_buckling_stress(x):
    """Returns compressive stress at which buckling will occur"""
    
    #k = None # to be determined based on the "final" stiffener and rib spacing
    b = sp.top_spacing
    t = p.t_sheet
    
    return k*np.pi**2*p.E_al2014*(t/b)/(12*(1-p.poisson_ratio_al2014)**2)

#print("Skin buckling limit: ", skin_buckling_stress(0))

def critical_column_buckling(x):
    """Critical column buckling force on the stiffener""" 
    
    l_eff = 0.5*p.rib_spacing
    
    return (np.pi**2*p.E_al2014*sp.I_zz_wingbox(x))/(l_eff**2)/10**6


def critical_crippling_stiffener(x):
    """Critical crippling stress for aluminium top stiffener in MPa"""
    
    alpha = 0.8
    n = 0.6
    yield_stress = p.tensile_yield_strength_al2014
    
    
    def stress_cc(K,b):
        return K*(np.pi**2*p.E_al2014*(sp.t_hat/b)**2/(12*(1-p.poisson_ratio_al2014**2)))
    
    t_sheet = sp.t_hat
    
    #areas
    area_a = sp.a*t_sheet
    area_b = (sp.b-t_sheet)*t_sheet
    area_c = (sp.c-t_sheet)*t_sheet
    
    #critical crippling stress
    cc_a = stress_cc(4,sp.a)
    cc_b = stress_cc(0.425,sp.b)
    cc_c = stress_cc(4,sp.c)
    
    ratio_a = alpha*(cc_a/yield_stress)**(1-n)
    ratio_b = alpha*(cc_b/yield_stress)**(1-n)
    ratio_c = alpha*(cc_c/yield_stress)**(1-n)

    total_crippling = yield_stress*(2*(ratio_a*area_a+ratio_b*area_b)+ratio_c*area_c) / (2*area_a + 2*area_b + area_c)
    
    return total_crippling/10**6


###  SHEAR STRESS CALCULATOR ###



### CALCULATE SHEAR AND NORMAL STRESS ###
normal_ru_list= []
normal_lu_list= []
normal_rl_list= []
normal_ll_list= []

for i in range(len(x_pos)):
    if x_pos[i] < p.x_strut:
        F_normal = -Frx
    else:
        F_normal = 0
    
    normal_ru, normal_lu, normal_rl, normal_ll = normal_stress(x_pos[i],y_max_list[i],momentzi[i],momentyi[i],F_normal,Izz_list[i],Iyy_list[i],area_list[i])
    normal_ru_list.append(normal_ru/10**6)
    normal_lu_list.append(normal_lu/10**6)
    normal_rl_list.append(normal_rl/10**6)
    normal_ll_list.append(normal_ll/10**6)
    
   
### MOMENT AND SHEAR DIAGRAM ###
plt.figure(2,figsize = (8,6))
plt.xlabel('Location along the length of the wingbox [m]',fontsize=13)
plt.ylabel('Moment in z-axis [Nm]',fontsize=13)
plt.plot(x_pos, momentzi, 'b') 

plt.figure(3,figsize = (8,6))
plt.xlabel('Location along the length of the wingbox [m]',fontsize=13)
plt.ylabel('Shear force in y direction[N]',fontsize=13)
plt.plot(x_pos, shearyi,'r')   

plt.figure(4,figsize = (8,6))
plt.xlabel('Location along the length of the wingbox [m]',fontsize=13)
plt.ylabel('Moment in around y-axis [Nm]',fontsize=13)
plt.plot(x_pos, momentyi, 'b') 

plt.figure(5,figsize = (8,6))
plt.xlabel('Location along the length of the wingbox [m]',fontsize=13)
plt.ylabel('Shear force in z direction [N]',fontsize=13)
plt.plot(x_pos, shearzi,'r') 


### PLOT NORMAL STRESS AT THE FOUR CORNERS
plt.figure(1,figsize = (8,6))
plt.xlabel('Location along the length of the strutbox [m]',fontsize=13)
plt.ylabel('Normal stress [MPa]',fontsize=13)
plt.plot(x_pos, normal_ru_list, 'r', label='Right upper corner')
plt.plot(x_pos, normal_lu_list, 'g', label='Left upper corner')
plt.plot(x_pos, normal_rl_list, 'b', label='Right lower corner')
plt.plot(x_pos, normal_ll_list, 'y', label='Left lower corner')
plt.legend(loc = 'upper right')      

plt.show()

print("Total weight: ",sp.total_weight,"kg")

print("Column buckling force: ",critical_column_buckling(0),"MPa")

print("Critical crippling stress of the hat stiffener: ",critical_crippling_stiffener(p.b/2/2),"MPa")








