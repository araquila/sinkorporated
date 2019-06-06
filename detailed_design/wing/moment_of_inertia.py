import numpy as np
import parameters as p
import matplotlib.pyplot as plt

#wingbox discretization
discretizations = 26
stepsize = p.b/2/discretizations
safety_factor = p.safetyfactor_wingbox


#spanwise locations for each of the discretizations
steps = np.zeros(discretizations)
for i in range(len(steps)):
    steps[i] = (i+.5)*stepsize


#height wing box
h_root = 0.929*p.h_max_root_wingbox
h_tip = 0.825*p.h_max_tip_wingbox #arbitrary

def width_wingbox(x):
    return p.w_root_wingbox - (2/p.b) * x * (p.w_root_wingbox - p.w_tip_wingbox) 

def height_wingbox(x):
    return h_root - (2/p.b) * x * (h_root - h_tip)

def y_neutral_line(x):
    return (p.n_upper_skin*(height_wingbox(x)-p.h_stiffener/2)*p.A_stiffener + p.n_lower_skin*(p.h_stiffener/2)*p.A_stiffener) / ((p.n_upper_skin+p.n_lower_skin)*p.A_stiffener)


#arbitrary loading diagrams
def M(x):
    return 100000 - 100000 * x /(p.b/2)
    
def T(x):
    return 300000 - 300000 * x /(p.b/2)


def enclosed_area(x):
    return width_wingbox(x)*height_wingbox(x)

def effective_airfoil_depth(x):
    bending_eff_factor = (1./3.)*(1+(height_wingbox(x+stepsize)/height_wingbox(x-stepsize))**2+(height_wingbox(x-stepsize)/height_wingbox(x-stepsize))**2)
    max_height = height_wingbox(x+stepsize)
    return bending_eff_factor*max_height
    
def Ix(x):
    moment_contribution = M(x)/(width_wingbox(x)*effective_airfoil_depth(x)*p.ultimate_bending_stress_al2024)
    torsion_contribution = T(x)/(2*enclosed_area(x)*p.ultimate_shear_stress_al2024)
    return safety_factor*width_wingbox(x)*effective_airfoil_depth(x)**2 / 2 * (moment_contribution+torsion_contribution)

def equivalent_sheet_thickness(x):
    moment_contribution = safety_factor * M(x) /(width_wingbox(x)*effective_airfoil_depth(x)*p.ultimate_bending_stress_al2024)
    torsion_contribution = safety_factor * T(x) / (2*enclosed_area(x)*p.ultimate_shear_stress_al2024)
    return moment_contribution + torsion_contribution


print(Ix(equivalent_sheet_thickness(steps[-1])))
print(equivalent_sheet_thickness(steps[-1]))

#import numpy as np
#import parameters as p
#import matplotlib.pyplot as plt
#
##wingbox discretization
#discretizations = 26
#stepsize = p.b/2/discretizations
#safety_factor = p.safetyfactor_wingbox
#
#
##spanwise locations for each of the discretizations
#steps = np.zeros(discretizations)
#for i in range(len(steps)):
#    steps[i] = (i+.5)*stepsize
#
#
##width wingbox
#w_root = 1.27758 #m
#w_tip = 0.51075 #m
#
##height wing box
#h_max_root = 0.35156
#h_root = 0.929*h_max_root 
#
#h_max_tip = 0.09704
#h_tip = 0.825*h_max_tip #arbitrary
#
#def width_wingbox(x):
#    return w_root - (2/p.b) * x * (w_root - w_tip) 
#
#def height_wingbox(x):
#    return h_root - (2/p.b) * x * (h_root - h_tip)
#
##stiffener geometry
#A_stiffener = 0.001 #m^2
#h_stiffener = 0.03#m
#
##number of stiffeners
#n_upper_skin = 5
#n_lower_skin = 2
#
#def y_neutral_line(x):
#    return (n_upper_skin*(height_wingbox(x)-h_stiffener/2)*A_stiffener + n_lower_skin*(h_stiffener/2)*A_stiffener) / ((n_upper_skin+n_lower_skin)*A_stiffener)
#
#
#
##Aluminium 2014-T6
#E_modulus = 414*10**6
#ultimate_bending_stress = 483*10**6
#ultimate_shear_stress = 290*10**6
#
#
##arbitrary loading diagrams
#def M(x):
#    return 100000 - 100000 * x /(p.b/2)
#    
#def T(x):
#    return 300000 - 300000 * x /(p.b/2)
#
#
#
#def enclosed_area(x):
#    return width_wingbox(x)*height_wingbox(x)
#
#def effective_airfoil_depth(x):
#    bending_eff_factor = (1./3.)*(1+(height_wingbox(x+stepsize)/height_wingbox(x-stepsize))**2+(height_wingbox(x-stepsize)/height_wingbox(x-stepsize))**2)
#    max_height = height_wingbox(x+stepsize)
#    return bending_eff_factor*max_height
#    
#def Ix(x):
#    moment_contribution = M(x)/(width_wingbox(x)*effective_airfoil_depth(x)*ultimate_bending_stress)
#    torsion_contribution = T(x)/(2*enclosed_area(x)*ultimate_shear_stress)
#    return safety_factor*width_wingbox(x)*effective_airfoil_depth(x)**2 / 2 * (moment_contribution+torsion_contribution)
#
#def equivalent_sheet_thickness(x):
#    moment_contribution = safety_factor * M(x) /(width_wingbox(x)*effective_airfoil_depth(x)*ultimate_bending_stress)
#    torsion_contribution = safety_factor * T(x) / (2*enclosed_area(x)*ultimate_shear_stress)
#    return moment_contribution + torsion_contribution
#
#
#print(Ix(equivalent_sheet_thickness(steps[-1])))
#print(equivalent_sheet_thickness(steps[-1]))