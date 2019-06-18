import parameters as p
from matplotlib import pyplot as plt
import numpy as np

def shear_and_moment(F_strut,F_strut_z,n):
    alpha = np.arctan((p.strut_pos_perc * p.b/2)/p.d_fuselage_outside)        # Angle of the strut with fuselage
    F_strut_y = F_strut * np.cos(alpha)
    F_support_y = F_strut_y
    F_support_z = F_strut_z
    
    
    x_support = 0.25
    
    x_pos = np.linspace(0,p.l_strutbox,n)
    V_list_y = []
    V_list_z = []
    M_list_z = []
    M_list_y = []
    
    
    for x in x_pos:
        if x < x_support:
            V_y = F_strut_y
            V_z = F_strut_z
            M_z = F_strut_y * x
            M_y = F_strut_z * x
        elif x_support <= x < (p.l_strutbox - x_support):
            V_y = F_strut_y - F_support_y
            V_z = F_strut_z - F_support_z
            M_z = F_strut_y * x - (x - x_support) * F_support_y
            M_y = F_strut_z * x - (x - x_support) * F_support_z
        else:
            V_y = F_strut_y - 2*F_support_y
            V_z = F_strut_z - 2*F_support_z
            M_z = F_strut_y * x - ((x - x_support) + (x - (p.l_strutbox - x_support))) * F_support_y
            M_y = F_strut_z * x - ((x - x_support) + (x - (p.l_strutbox - x_support))) * F_support_z
        
        V_list_y.append(-V_y)
        V_list_z.append(V_z)
        M_list_z.append(M_z)
        M_list_y.append(-M_y)
        
    
    plt.figure(2,figsize = (8,6))
    plt.xlabel('Location along the length of the strutbox [m]',fontsize=13)
    plt.ylabel('Moment around the z-axis [Nm]',fontsize=13)
    plt.plot(x_pos, M_list_z, 'b')
    
    plt.figure(4,figsize = (8,6))
    plt.xlabel('Location along the length of the strutbox [m]',fontsize=13)
    plt.ylabel('Moment around the y-axis [Nm]',fontsize=13)
    plt.plot(x_pos, M_list_y, 'b')
    
    x_pos_shear = np.linspace(0,p.l_strutbox,n+2)
    V_list_y_0 = V_list_y
    V_list_z_0 = V_list_z
#    V_list_z_0.insert(0,0)
#    V_list_z_0.append(0)
#    V_list_y_0.insert(0,0)
#    V_list_y_0.append(0)  
    
    plt.figure(1,figsize = (8,6))
    plt.xlabel('Location along the length of the strutbox [m]',fontsize=13)
    plt.ylabel('Shear force in the y-direction [N]',fontsize=13)
    plt.plot(x_pos, V_list_y_0, 'r')
    
    plt.figure(5,figsize = (8,6))
    plt.xlabel('Location along the length of the strutbox [m]',fontsize=13)
    plt.ylabel('Shear force in the z-direction [N]',fontsize=13)
    plt.plot(x_pos, V_list_z_0, 'r')
    
    return V_list_y, V_list_z, M_list_z, M_list_y