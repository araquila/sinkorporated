# Aerodynamic design for the wing

# Imports
import numpy as np


# --------------------------------------------------------


# ------------------- Aerodynamic formulas -------------------

def induced_angle_of_attack(A, e, C_L_cruise):
    # Inputs:
#    A = aspect ratio 
#    e = Owswald efficiency factor
#    C_L_cruise = lift coefficient during cruise
    
    # Determine the induced angle of attack
    induced_alpha = C_L_cruise / (np.pi * A * e)
    
    # Return outputs:
    # induced_alpha is the induced angle of attack
    return induced_alpha

print(induced_angle_of_attack(18,0.85,0.41))

# ------------------- Planform formula -------------------
def planform(S, A, M_cruise, C_L_cruise, quarter_chord_sweep, t_c, delta_mach = 0.03, t_c_forced = False):
    # Inputs: 
#    S = surface area in m^2
#    A = aspect ratio
#    M_cruise = Mach number during cruise
#    C_L_cruise = lift coefficient during cruise
#    quarter_chord_sweep = sweep at quarter chord in radians
#    t_c = a forced thickness over chord ratio if t_c_forced is set to True
#    Delta_mach can range from 0 to 0.05 but is given as 0.03 in ADSEE Slides

    # Determine span
    b = np.sqrt(A * S)
    
    # Determine taper ratio
    taper = 0.2 * (2 - quarter_chord_sweep)
    
    # Determine root chord
    root_chord = (2 * S)/((1 + taper) * b)
    
    # Determine tip chord
    tip_chord = taper * root_chord
    
    # Determine sweep at half chord point
    half_chord_sweep = np.arctan(np.tan(quarter_chord_sweep) - (4 / A) * (0.25 * (1 - taper)/(1 + taper)))
    
    # Determine thickness over chord ratio: forced or unforced
    if t_c_forced == True:
        t_c_ratio = t_c
    else:
        if 0 < M_cruise < 0.8:
                t_c_ratio = min(0.18, (np.cos(half_chord_sweep)**3 * (1 - (M_cruise + delta_mach) * np.cos(half_chord_sweep)) - 0.115 * C_L_cruise**1.5) / np.cos(half_chord_sweep)**2)
    
    # Check speed
    if M_cruise <= 0:
        raise NameError("Plane flying backwards??")
    if M_cruise >= 0.8:
        raise NameError("Going supercritical now, are we?")
        
    # Return outputs:
#    b = span in m
#    taper = taper ratio
#    root_chord = root chord in m
#    tip_chord = tip chord in m
#    t_c_ratio = thickness over chord ratio
    return b, taper, root_chord, tip_chord, t_c_ratio


b, taper, root, tip, tc = planform(49.21,18,0.6,0.5,0,0,delta_mach=0.03,t_c_forced=False)
print(b, taper, root, tip, tc)


# ------------------- MAC formula -------------------
def MAC(root_chord,t_c_ratio):
    # Inputs:
#    root_chord = root chord in m
#    t_c_ratio = thickness over chord ratio
    
    # Determine MAC
    MAC = root_chord*(2/3)*((1+t_c_ratio+t_c_ratio**2)/(1+t_c_ratio))
    
    # Return outputs:
    # MAC = mean aerodynamic chord in m
    return MAC


# ------------------- Dihedral angle -------------------
def det_dihedral_angle(quarter_chord_sweep, high = False, mid = False, low = False):
    # Inputs:
#    quarter_chord_sweep = sweep angle at quarter chord line in radians
#    high, mid, low = wing position on fuselage
    
    # Determine effect of sweep
    dihedral = quarter_chord_sweep * 18 / np.pi
    
    # Determine dihedral dependent on wing position on fuselage
    # Return outputs:
    # angle is the dihedral angle in degrees
    if high:
        dihedral_angle = 1 - dihedral
        return dihedral_angle
    if mid:
        dihedral_angle = 3 - dihedral
        return dihedral_angle
    if low:
        dihedral_angle = 5 - dihedral
        return dihedral_angle
    else:
        raise NameError("Where is the wing?")


