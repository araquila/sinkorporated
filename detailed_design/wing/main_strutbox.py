import parameters as p
import numpy as np
from matplotlib import pyplot as plt
from wing.shear_and_moment_strutbox import shear_and_moment
import wing.section_properties_strutbox as scs
import wing.wing_deflection2 as wd2 


### DISCRETIZATION OF THE STRUTBOX ###
n = 51
x_pos = np.linspace(0,p.l_strutbox,n)

### OBTAIN CROSSECTIONAL PROPERTIES ###
Izz_list = []
Iyy_list = []
first_moment_of_area_list_y = []
area_list = []
y_max_list = []

for x in x_pos:
    Izz_list.append(scs.I_zz_strutbox(x))
    Iyy_list.append(scs.I_yy_strutbox(x))
    first_moment_of_area_list_y.append(scs.first_moment_of_area(x))
    area_list.append(scs.cross_sectional_area(x))
    y_max_list.append(scs.y_max(x))
    
### OBTAIN STRUT FORCE, REACTION FORCES AND REACTION MOMENT ###
lengthdata = 50
Lift, Chord, Yle, Drag = wd2.read_aero_data("wing/aquiladata1.txt", lengthdata, p.V_cruise, p.rho)
Frx, Fry, Fs, Mrz, Frz, Fsz, Mry, momentyi, momentzi, shearyi, shearzi, vyi, vny, vzi, vnz, xi, theta = wd2.CallForces(Lift, Yle, Drag, p.tot_thrust, Iyy_list, Izz_list,70*10**9, p.engine_pos_perc, p.strut_pos_perc, p.pod_pos_perc)

alpha = np.arctan((p.strut_pos_perc * p.b/2)/p.d_fuselage_outside)        # Angle of the strut with fuselage
F_strut_x = Fs * np.sin(alpha)

### OBTAIN SHEAR AND MOMENT DIAGRAM ###
V_list_y, V_list_z, M_list_z, M_list_y = shear_and_moment(Fs,Fsz,n)

### NORMAL STRESS CALCULATOR ###
def normal_stress(x,y,moment_z,moment_y,normal_force,I_zz,I_yy,area):
    """Returns bending moment and normal force stress""" 
    
    z = scs.width_strutbox(x)/2
    
    #-my/I according to the formula, makes sense because for a positive Mz the top skin will be in compression
    moment_z_upperskin = -moment_z*(scs.height_strutbox(x)-y)/I_zz
    moment_z_lowerskin = -moment_z*y/I_zz
    
    moment_y_rightflange = -moment_y*z/I_yy
    moment_y_leftflange = -moment_y*-z/I_yy
    
    normal_force_stress = 2*normal_force/area
    
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
    print("Neutral axis: y =",scs.centroid_y(x),"(",scs.centroid_y(x)/scs.height_strutbox(x)*100,"%)",)
    print("Maximum tension: ",max(normal_ru,normal_lu,normal_rl,normal_ll)/10**6,"MPa")
    print("Maximum compression: ", min(normal_ru,normal_lu,normal_rl,normal_ll)/10**6,"MPa")
    print("")
    
    return normal_ru,normal_lu,normal_rl,normal_ll

###  SHEAR STRESS CALCULATOR ###


def skin_buckling_stress(x):
    """Returns compressive stress at which buckling will occur"""
    
    k = 4 # to be determined based on the "final" stiffener and rib spacing
    b = scs.top_spacing
    poisson_ratio = p.poisson_ratio_al2014
    t = p.t_sheet
    
    return k*np.pi**2*p.E_al2014*(t/b)/(12*(1-poisson_ratio)**2)


def critical_column_buckling(x):
    """Critical column buckling force on the stiffener""" 
    
    l_eff = 0.5*p.rib_spacing
    
    return (np.pi**2*p.E_al2014*scs.I_yy_strutbox(x))/(l_eff**2)


def critical_crippling_stiffener(x):
    """Critical crippling stress for aluminium top stiffener in MPa"""
    
    alpha = 0.8
    n = 0.6
    yield_stress = p.tensile_yield_strength_al2014
    
    
    def stress_cc(K,b):
        return K*(np.pi**2*p.E_al2014*(scs.t_hat/b)**2/(12*(1-p.poisson_ratio_al2014**2)))
    
    t_sheet = scs.t_hat
    
    #areas
    area_a = scs.a*t_sheet
    area_b = (scs.b-t_sheet)*t_sheet
    area_c = (scs.c-t_sheet)*t_sheet
    
    #critical crippling stress
    cc_a = stress_cc(4,scs.a)
    cc_b = stress_cc(0.425,scs.b)
    cc_c = stress_cc(4,scs.c)
    
    ratio_a = alpha*(cc_a/yield_stress)**(1-n)
    ratio_b = alpha*(cc_b/yield_stress)**(1-n)
    ratio_c = alpha*(cc_c/yield_stress)**(1-n)

    total_crippling = yield_stress*(2*(ratio_a*area_a+ratio_b*area_b)+ratio_c*area_c) / (2*area_a + 2*area_b + area_c)
    
    return total_crippling/10**6

### CALCULATE REQUIRED THICKNESS ###

def required_shear_thickness(V_list_y,V_list_z,Iyy_list,Izz_list,first_moment_of_area_list_y,moment_of_area_list_z):
    t = 0
    for i in range(len(x_pos)):
        t_temp = 1
        if t_temp > t:
            t = t_temp
        return t


### CALCULATE NORMAL STRESS ###
normal_ru_list= []
normal_lu_list= []
normal_rl_list= []
normal_ll_list= []

for i in range(len(x_pos)):
    normal_ru, normal_lu, normal_rl, normal_ll = normal_stress(x_pos[i],y_max_list[i],M_list_z[i],M_list_y[i],F_strut_x,Izz_list[i],Iyy_list[i],area_list[i])
    normal_ru_list.append(normal_ru/10**6)
    normal_lu_list.append(normal_lu/10**6)
    normal_rl_list.append(normal_rl/10**6)
    normal_ll_list.append(normal_ll/10**6)
    
### PLOT NORMAL STRESS AT THE FOUR CORNERS
plt.figure(3,figsize = (8,6))
plt.xlabel('Location along the length of the strutbox [m]',fontsize=13)
plt.ylabel('Normal stress [MPa]',fontsize=13)
plt.plot(x_pos, normal_ru_list, 'r', label='Right upper corner')
plt.plot(x_pos, normal_lu_list, 'g', label='Left upper corner')
plt.plot(x_pos, normal_rl_list, 'b', label='Right lower corner')
plt.plot(x_pos, normal_ll_list, 'y', label='Left lower corner')
plt.legend(loc = 'upper right')

max_compressive_stress = min(normal_ru_list)*p.safety_factor_compression
max_tensile_stress = max(normal_ll_list)*p.safety_factor_tension

print("Max compressive: ",max_compressive_stress,"MPa")
print("Max tensile: ",max_tensile_stress,"MPa")

print("Total weight: ",scs.total_weight,"kg")

print("Skin buckling limit: ",skin_buckling_stress(p.b/2/2)/10**6,"MPa")

print("Column buckling limit: ",critical_column_buckling(p.b/2/2),"MPa")

print("Critical crippling stress of the hat stiffener: ",critical_crippling_stiffener(p.b/2/2),"MPa")

print("")


if abs(max_compressive_stress) < abs(skin_buckling_stress(8)/10**6):
    print("Skin buckling passed")
else:
    print("Failure on skin buckling")
    
if abs(max_compressive_stress) < abs(critical_column_buckling(8)):
    print("Column buckling passed")
else:
    print("Failure on column buckling")
    
if abs(max_compressive_stress) < abs(critical_crippling_stiffener(p.b/2/2)):
    print("Cripple limit of the stiffener passed")
else:
    print("Failure on stiffener crippling")
    
if max_tensile_stress > p.tensile_yield_strength_al2014/10**6:
    print("Yielded at the strut")
else:
    print("Max stress",max_tensile_stress/(p.tensile_yield_strength_al2014/10**6)*100,"% of the yield strength")
    






