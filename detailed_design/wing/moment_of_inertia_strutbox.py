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
h_root = 0.929*p.h_max_root_strutbox 
h_tip = h_root

def width_wingbox(x):
    return p.w_root_strutbox - (2/p.b) * x * (p.w_root_strutbox - p.w_tip_strutbox) 

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