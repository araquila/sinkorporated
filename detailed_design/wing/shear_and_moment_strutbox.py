import parameters as p
from matplotlib import pyplot as plt
import numpy as np

def shear_and_moment(F_strut,n):
    alpha = np.arctan((p.strut_pos_perc * p.b/2)/p.d_fuselage_outside)        # Angle of the strut with fuselage
    F_strut_y = F_strut * np.cos(alpha)
    F_support = F_strut_y
    
    x_support = 0.25
    
    x_pos = np.linspace(0,p.l_strutbox,n)
    V_list = []
    M_list = []
    
    for x in x_pos:
        if x < x_support:
            V = F_strut_y
            M = F_strut_y * x
        elif x_support <= x < (p.l_strutbox - x_support):
            V = F_strut_y - F_support
            M = F_strut_y * x - (x - x_support) * F_support
        else:
            V = F_strut_y - 2*F_support
            M = F_strut_y * x - ((x - x_support) + (x - (p.l_strutbox - x_support))) * F_support
        V_list.append(V)
        M_list.append(M)
    
    plt.figure(2,figsize = (8,6))
    plt.xlabel('Location along the length of the strutbox [m]',fontsize=13)
    plt.ylabel('Moment [Nm]',fontsize=13)
    plt.plot(x_pos, M_list, 'b')
    
    x_pos_shear = np.linspace(0,p.l_strutbox,n+2)
    V_list_0 = V_list
    V_list_0.insert(0,0)
    V_list_0.append(0)  
    
    plt.figure(1,figsize = (8,6))
    plt.xlabel('Location along the length of the strutbox [m]',fontsize=13)
    plt.ylabel('Shear force [N]',fontsize=13)
    plt.plot(x_pos_shear, V_list_0, 'r')
    return V_list, M_list