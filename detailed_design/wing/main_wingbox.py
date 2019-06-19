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
from wing.shear_stress import max_shear_stress

### DISCRETIZATION OF THE WINGBOX ###
n = 101
x_pos = np.linspace(0,p.b/2,n)


### OBTAIN CROSSECTIONAL PROPERTIES ###
Izz_list = []
Iyy_list = []
first_moment_of_area_list = []
area_list = []
y_max_list = []
hi = []
bi = []

for x in x_pos:
    Izz_list.append(sp.I_zz_wingbox(x))
    Iyy_list.append(sp.I_yy_wingbox(x))
    first_moment_of_area_list.append(sp.first_moment_of_area_y(x))
    area_list.append(sp.cross_sectional_area(x))
    y_max_list.append(sp.y_max(x))
    hi.append(sp.height_wingbox(x))
    bi.append(sp.width_wingbox(x))

### OBTAIN STRUT FORCE, REACTION FORCES AND REACTION MOMENT ###
lengthdata = 100
Lift, Chord, Yle, Drag, AeroMoment = wd2.read_aero_data("wing/datastrut5.txt", lengthdata, p.V_cruise, p.rho)

nullen = np.zeros(len(Lift))
#Frx, Fry, Fs, Mrz, Frz, Fsz, Mry, momentyi, momentzi, shearyi, shearzi, vyi, vny, vzi, vnz, xi, theta = wd2.CallForces(nullen, Yle, nullen, 0, Iyy_list, Izz_list,70*10**9, p.engine_pos_perc, p.strut_pos_perc, p.pod_pos_perc)

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
    
    area = sp.cross_sectional_area(x)
    
    normal_force_stress = normal_force/area
    
    normal_ru = moment_z_upperskin + moment_y_rightflange +     normal_force_stress 
    normal_lu = moment_z_upperskin + moment_y_leftflange +     normal_force_stress 
    normal_rl = moment_z_lowerskin + moment_y_rightflange +     normal_force_stress 
    normal_ll = moment_z_lowerskin + moment_y_leftflange +  normal_force_stress 
    
#    if max(normal_ru,normal_lu,normal_rl,normal_ll) == normal_ru:
#        print("Max tension at right upper corner")
#    elif max(normal_ru,normal_lu,normal_rl,normal_ll) == normal_lu:
#        print("Max tension at left upper corner")
#    elif max(normal_ru,normal_lu,normal_rl,normal_ll) == normal_rl:
#        print("Max tension at right lower corner")
#    elif max(normal_ru,normal_lu,normal_rl,normal_ll) == normal_ll:
#        print("Max tension at left lower corner")
#        
#    
#    if min(normal_ru,normal_lu,normal_rl,normal_ll) == normal_ru:
#        print("Max compression at right upper corner")
#    elif min(normal_ru,normal_lu,normal_rl,normal_ll) == normal_lu:
#        print("Max compression at left upper corner")
#    elif min(normal_ru,normal_lu,normal_rl,normal_ll) == normal_rl:
#        print("Max compression at right lower corner")
#    elif min(normal_ru,normal_lu,normal_rl,normal_ll) == normal_ll:
#        print("Max compression at left lower corner")
#        
#        
#    print("x: ",x)
#    print("Neutral axis: y =",sp.centroid_y(x),"(",sp.centroid_y(x)/sp.height_wingbox(x)*100,"%)",)
#    print("Maximum tension: ",max(normal_ru,normal_lu,normal_rl,normal_ll)/10**6,"MPa")
#    print("Maximum compression: ", min(normal_ru,normal_lu,normal_rl,normal_ll)/10**6,"MPa")
#    print("Normal force stress: ",normal_force_stress,"MPa")
#    print("")
    
    return normal_ru,normal_lu,normal_rl,normal_ll

 
def von_mises(x,y,moment_z,moment_y,normal_force,I_zz,I_yy,area, tau_max):
    """Returns maximum Von Mises stress"""
    
    z = sp.width_wingbox(x)/2
    
    #-my/I according to the formula, makes sense because for a positive Mz the top skin will be in compression
    moment_z_upperskin = -moment_z*(sp.height_wingbox(x)-abs(y))/I_zz
    moment_z_lowerskin = -moment_z*y/I_zz
    
    moment_y_rightflange = -moment_y*z/I_yy
    moment_y_leftflange = -moment_y*-z/I_yy
    
    area = sp.cross_sectional_area(x)
    
    normal_force_stress = normal_force/area
    
    stress_x_lower = moment_z_lowerskin + normal_force_stress
    stress_x_upper = moment_z_upperskin + normal_force_stress
    
    stress_y_right = moment_y_rightflange
    stress_y_left = moment_y_leftflange
    
    
    vm_ll_1 = (stress_x_lower + stress_y_left)/2 + np.sqrt(((stress_x_lower-stress_y_left)/2)**2 + tau_max**2)
    vm_ll_2 = (stress_x_lower + stress_y_left)/2 - np.sqrt(((stress_x_lower-stress_y_left)/2)**2 + tau_max**2)
    vm_ll = np.sqrt(vm_ll_1**2 + vm_ll_2**2 -vm_ll_1*vm_ll_2 + 3*tau_max**2)
    
    vm_lr_1 = (stress_x_lower + stress_y_right)/2 + np.sqrt(((stress_x_lower-stress_y_right)/2)**2 + tau_max**2)
    vm_lr_2 = (stress_x_lower + stress_y_right)/2 - np.sqrt(((stress_x_lower-stress_y_right)/2)**2 + tau_max**2)
    vm_lr = np.sqrt(vm_lr_1**2 + vm_lr_2**2 -vm_lr_1*vm_lr_2 + 3*tau_max**2)
    
    vm_ur_1 = (stress_x_upper + stress_y_right)/2 + np.sqrt(((stress_x_upper-stress_y_right)/2)**2 + tau_max**2)
    vm_ur_2 = (stress_x_upper + stress_y_right)/2 - np.sqrt(((stress_x_upper-stress_y_right)/2)**2 + tau_max**2)
    vm_ur = np.sqrt(vm_ur_1**2 + vm_ur_2**2 -vm_ur_1*vm_ur_2 + 3*tau_max**2)
    
    vm_ul_1 = (stress_x_upper + stress_y_left)/2 + np.sqrt(((stress_x_upper-stress_y_left)/2)**2 + tau_max**2)
    vm_ul_2 = (stress_x_upper + stress_y_left)/2 - np.sqrt(((stress_x_upper-stress_y_left)/2)**2 + tau_max**2)
    vm_ul = np.sqrt(vm_ul_1**2 + vm_ul_2**2 -vm_ul_1*vm_ul_2 + 3*tau_max**2)
  
    
    return vm_ll,vm_lr,vm_ur,vm_ul
 


def skin_buckling_stress(x):
    """Returns compressive stress at which buckling will occur"""
    
    k = 4 # to be determined based on the "final" stiffener and rib spacing
    b = sp.top_spacing
    t = p.t_sheet
    
    return k*np.pi**2*p.E_sheet/(12*(1-p.poisson_ratio_al2014**2)) * (t/b)**2


def shear_skin_buckling_stress(x):
    """Returns compressive stress at which buckling will occur"""
    
    k = 5.35 # to be determined based on the "final" stiffener and rib spacing
    b = sp.top_spacing
    t = p.t_sheet
    
    return k*np.pi**2*p.E_sheet/(12*(1-p.poisson_ratio_al2014**2)) * (t/b)**2



#print("Skin buckling limit: ", skin_buckling_stress(0))

def critical_column_buckling(x):
    """Critical column buckling stress on the panel""" 
    
    l_eff = p.rib_spacing
    
    return (np.pi**2*p.E_sheet*sp.I_yy_wingbox(x))/(l_eff**2*sp.cross_sectional_area(x))/10**6


def column_buckling_stiffener(x):
    """Critical column buckling stress on the stiffener"""
    
    l_eff = p.rib_spacing
    
    return (np.pi**2*p.E_compressive_2099*sp.I_yy_hat)/(l_eff**2*sp.A_hat)/10**6


def critical_crippling_stiffener(x):
    """Critical crippling stress for aluminium top stiffener in MPa"""
    
    alpha = 0.8
    n = 0.6
    yield_stress = p.ultimate_yield_strength_2099
    
    
    def stress_cc(K,b):
        return K*(np.pi**2*p.E_compressive_2099*(sp.t_hat/b)**2/(12*(1-p.poisson_ratio_al2014**2)))
    
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
    
    return total_crippling


def critical_panel_buckling(x):
    """Returns critical panel buckling stress"""
    C = 6.98
    v = p.poisson_ratio_al2014
    
    crippling = critical_crippling_stiffener(x)
    
    we = p.t_sheet * np.sqrt(C*np.pi**2/(12*(1-v**2))) * np.sqrt(p.E_sheet/crippling)
    
    
    return we


def crack_length_sheet(stress):
    """Fast fracture crack length for given stress"""
    
    Kic  = p.fracture_toughness_2195
    a = Kic/((stress*10**6)**2 * np.pi)
    
    return a
    
### CALCULATE MAXIMUM SHEAR STRESS PER SECTION ###
tau_max = max_shear_stress(Lift, Drag, AeroMoment, Chord, shearyi, shearzi, hi, bi, Izz_list, Iyy_list)
 
### CALCULATE SHEAR AND NORMAL STRESS ###
normal_ru_list= []
normal_lu_list= []
normal_rl_list= []
normal_ll_list= []

vm_ll_list = []
vm_lr_list = []
vm_ur_list = []
vm_ul_list = []


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
    
    vm_ll,vm_lr,vm_ur,vm_ul = von_mises(x_pos[i],y_max_list[i],momentzi[i],momentyi[i],F_normal,Izz_list[i],Iyy_list[i],area_list[i],tau_max[i])
    vm_ll_list.append(vm_ll/10**6)
    vm_lr_list.append(vm_lr/10**6)
    vm_ul_list.append(vm_ul/10**6)
    vm_ur_list.append(vm_ur/10**6)
    
    
    
    

# Error correction
error_lower = normal_ll_list[-1]
error_upper = normal_lu_list[-1]
    
for i in range(len(x_pos)):
    normal_ru_list[i] -= error_upper
    normal_lu_list[i] -= error_upper
    normal_rl_list[i] -= error_lower
    normal_ll_list[i] -= error_lower
    
 
#plt.plot(x_pos,vm_ll_list,marker='x')
#plt.plot(x_pos,vm_lr_list)
#plt.plot(x_pos,vm_ul_list)
#plt.plot(x_pos,vm_ur_list)
#plt.show()


plt.figure()
plt.xlim([0,p.b/2])

#### MOMENT AND SHEAR DIAGRAM ###
#plt.figure(2,figsize = (8,6))
#plt.xlabel('Location along the length of the wingbox [m]',size='large')
#plt.ylabel('Moment in z-axis [kNm]',size='large')
#momentzi = np.array(momentzi)
#plt.plot(x_pos, momentzi/1000, 'b') 

##plt.figure(3,figsize = (8,6))
#plt.xlabel('Location along the length of the wingbox [m]',size='large')
#shearyi = np.array(shearyi)
#plt.ylabel('Shear force in y-direction [kN]',size='large')
#plt.plot(x_pos, shearyi/1000,'r')   

#plt.figure(4,figsize = (8,6))
#plt.xlabel('Location along the length of the wingbox [m]',size='large')
#momentyi = np.array(momentyi)
#plt.ylabel('Moment in around y-axis [kNm]',size='large')
#plt.plot(x_pos, momentyi/1000,'b')

#plt.figure(5,figsize = (8,6))
#plt.xlabel('Location along the length of the wingbox [m]',size='large')
#shearzi = np.array(shearzi)
#plt.ylabel('Shear force in z-direction [kN]',size='large')
#plt.plot(x_pos, shearzi/1000,'r') 

#### PLOT SHEAR STRESS ###
#plt.figure(6,figsize = (8,6))
#plt.xlabel('Location along the length of the wingbox [m]',size='large')
#plt.ylabel('Shear stress [MPa]',size='large')
#plt.plot(x_pos, tau_max,'r')

### PLOT NORMAL STRESS AT THE FOUR CORNERS
plt.figure(7,figsize = (8,6))
plt.xlabel('Location along the length of the wingbox [m]',size='large')
plt.ylabel('Normal stress [MPa]',size='large')
plt.plot(x_pos, normal_ru_list, 'r', label='Right-up')
plt.plot(x_pos, normal_lu_list, 'g', label='Left-up')
plt.plot(x_pos, normal_rl_list, 'b', label='Right-bottom')
plt.plot(x_pos, normal_ll_list, 'y', label='Left-bottom')
plt.legend(loc="best", fontsize="large")
plt.xlim([0,p.b/2])


plt.show()

max_compressive_stress = min(normal_ru_list) #*p.safety_factor_compression
max_tensile_stress = max(normal_ll_list) #*p.safety_factor_tension


print("Top stiffeners: ",p.n_upper_skin_wingbox)
print("Bottom stiffeners: ",p.n_lower_skin_wingbox)
print("Skin thickness: ",p.t_sheet*1000, "mm")
print("Total weight: ",sp.total_weight,"kg")
print("")


print("Max shear: ",max(tau_max),"MPa")
print("Max compressive: ",max_compressive_stress,"MPa")
print("Max tensile: ",max_tensile_stress,"MPa")
print("Max Von Mises: ", max(max(vm_ll_list),max(vm_lr_list),max(vm_ul_list),max(vm_ur_list)))
print("")

print("Skin buckling limit: ",skin_buckling_stress(p.strut_pos_perc*p.b/2)/10**6,"MPa")

print("Panel column buckling limit: ",critical_column_buckling(p.strut_pos_perc*p.b/2),"MPa")

print("Stiffener column buckling limit: ",column_buckling_stiffener(p.strut_pos_perc*p.b/2),"MPa")

print("Critical crippling stress of the hat stiffener: ",critical_crippling_stiffener(p.strut_pos_perc*p.b/2)/10**6,"MPa")

print("Shear stress buckling limit: ",shear_skin_buckling_stress(p.strut_pos_perc*p.b/2)/10**6,"MPa (still to be determined where the maximum shear stress occurs spanwise)")

print("")
print("Tests: ")
#print("Crack length: ",crack_length_sheet(max(max_tensile_stress,max_compressive_stress)),"m")

if abs(max_compressive_stress) < abs(skin_buckling_stress(p.strut_pos_perc*p.b/2)/10**6):
    print("Skin buckling passed")
else:
    print("Failure on skin buckling")
    
if abs(max_compressive_stress) < abs(critical_column_buckling(p.strut_pos_perc*p.b/2)):
    print("Panel column buckling passed")
else:
    print("Failure on panel column buckling")
    
if abs(max_compressive_stress) < abs(column_buckling_stiffener(p.strut_pos_perc*p.b/2)):
    print("Stiffener column buckling passed")
else:
    print("Failure on stiffener column buckling")
    
if abs(max_compressive_stress)*10**6 < abs(critical_crippling_stiffener(p.strut_pos_perc*p.b/2)):
    print("Cripple limit of the stiffener passed")
else:
    print("Failure on stiffener crippling")
    
if shear_skin_buckling_stress(p.strut_pos_perc*p.b/2)/4 > max(tau_max):
    print("Shear buckling passed")
else:
    print("Failure on shear buckling")

print("")

if (max_tensile_stress > p.ult_yield_strength_2195/10**6) :
    print("Failure on tensile yielding at the strut",max_tensile_stress/(p.ult_yield_strength_2195/10**6)*100)
else:
    print("Max tensile stress",max_tensile_stress/(p.ult_yield_strength_2195/10**6)*100,"% of the yield strength")

    
if abs(max_compressive_stress) > (p.ult_compressive_strength_2195/10**6):
    print("Failure on compressive yielding at the strut",-max_compressive_stress/(p.ult_compressive_strength_2195/10**6)*100)
else:
    print("Max compressive stress",-max_compressive_stress/(p.ult_compressive_strength_2195/10**6)*100,"% of the yield strength")

if p.ultimate_shear_stress_al2024 > max(tau_max)*10**6:
    print("Shear limit of the material passed: ",max(tau_max)*10**6/p.ult_shear_strength_2195 *100,"% of the max")
else: 
    print("Shear failure of the material")


