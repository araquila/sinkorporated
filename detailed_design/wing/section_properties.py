# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:48:26 2019

@author: robert
"""

import parameters as p
import numpy as np
import matplotlib.pyplot as plt


#actual root enclosed area 0.316
#actual tip enclosed area 0.045


#height at root of the wing box
h_root = 0.929*p.h_max_root_wingbox
h_tip = 0.825*p.h_max_tip_wingbox #arbitrary


def width_wingbox(x):
    return p.w_root_wingbox - (2/p.b) * x * (p.w_root_wingbox - p.w_tip_wingbox) 

def height_wingbox(x):
    return h_root - (2/p.b) * x * (h_root - h_tip)



t_sheet = 0.003 #m
t_hat = 0.002 #m
t_z = 0.002 #m


#hat geometry
a = 0.06
b = 0.04
c = 0.05

width_hat = 2*b+c-t_hat

A_hat = t_hat*(2*b+2*a+c)

y_centroid_hat = a + 2*t_hat - ( 2*b*t_hat*(a+3*t_hat/2) + 2*a*t_hat*(a/2+t_hat)+c*t_hat*t_hat/2 ) / (t_hat*(2*b+2*a+c))
z_centroid_hat = b - t_hat + c/2
 
I_zz_hat = 2*(b*t_hat*(a+3*t_hat/2)**2 + a**3*t_hat/12 + a*t_hat*(a/2+t_hat)**2)
I_yy_hat = 2*b**3*t_hat/12 + c**3*t_hat/12 + 2*b*t_hat*(b-t_hat+c/2) + c*t_hat*(c/2-t_hat/2)


#Z-stiffener geometry
d = 0.04
e = 0.06 
f = 0.045

width_z = d+f-t_z

A_z = t_z*(d+e+f)

y_centroid_z = (d*t_z**2/2 + e*t_z*(t_z+e/2)+f*t_z*(e+3*t_z/2))/(t_z*(d+e+f))
z_centroid_z = ((d/2)*d*t_z + (d-t_z/2)*e*t_z + (d-t_z + f/2)*f*t_z) / (t_z*(d+e+f) )

I_zz_z = e**3*t_z*(e+t_z)**2 +f*t_z*(e+3*t_z/2)**2
I_yy_z = d**3*t_z/12 + (d/2-z_centroid_z)**2*d*t_z + e*t_z*(d-t_z/2-z_centroid_z)**2 + f**3*t_z/12 + f*t_z*(d-t_z+f/2-z_centroid_z)**2

#section properties
n_top = 0 # at the upper skin
n_bottom = 0 # at the lower skin


#select top and bottom stiffener
width_top = width_hat
width_bottom = width_z

y_top = y_centroid_hat
y_bottom = y_centroid_z

z_top = z_centroid_hat
z_bottom = z_centroid_z

A_top = A_hat
A_bottom = A_z

I_zz_top = I_zz_hat
I_zz_bottom = I_zz_z

I_yy_top = I_yy_hat
I_yy_bottom = I_yy_z



#stiffener spacing
top_spacing = (width_wingbox(0) - n_top * (width_top))/(n_top+1)
lower_spacing = (width_wingbox(0) - n_bottom * (width_bottom))/(n_bottom+1)

#centroid y
def centroid_y(x):
    return (2*(height_wingbox(x)/2*(height_wingbox(x)-2*t_sheet)*t_sheet) + width_wingbox(x)*t_sheet*((height_wingbox(x)-t_sheet/2)+t_sheet/2) + n_top*(height_wingbox(x)-t_sheet-y_top)*A_top + n_bottom*(t_sheet+y_bottom)*A_bottom) / (2*((height_wingbox(x)-2*t_sheet)*t_sheet + width_wingbox(x)*t_sheet) + n_top*A_top + n_bottom*A_bottom)    



#at the root
centroid_y_root = (2*((t_sheet+height_wingbox(0)/2)*(height_wingbox(0)-2*t_sheet)*t_sheet) + (-t_sheet/2+height_wingbox(0))*width_wingbox(0)*t_sheet + t_sheet**2/2*width_wingbox(0) + n_top*(3*t_sheet/2+height_wingbox(0)-y_top)*A_top + n_bottom*(t_sheet+y_bottom)*A_bottom) / (2*(height_wingbox(0)-2*t_sheet)*t_sheet + 2*width_wingbox(0)*t_sheet + n_top * A_top + n_bottom * A_bottom)


#moment of inertia about z-axis
dy_vertical_flange = np.abs(height_wingbox(0)/2 - centroid_y_root)
dy_upper_flange = np.abs((height_wingbox(0)-t_sheet/2)-centroid_y_root)
dy_bottom_flange = np.abs(centroid_y_root - t_sheet/2)

dy_top_stiffener = np.abs(centroid_y_root - (height_wingbox(0)-t_sheet-y_top))
dy_bottom_stiffener = np.abs(centroid_y_root - (t_sheet + y_bottom))

Izz = t_sheet*(height_wingbox(0)-2*t_sheet)**3/12 + n_top*I_zz_top + n_bottom*I_zz_bottom + 2*(height_wingbox(0)-2*t_sheet)*t_sheet*dy_vertical_flange**2 + width_wingbox(0)*t_sheet*dy_upper_flange**2 + width_wingbox(0)*t_sheet*dy_bottom_flange**2 + n_top*A_top*dy_top_stiffener**2 + n_bottom*A_bottom*dy_bottom_stiffener**2

#moment of inertia about y-axis without stiffeners
I_yy_nostiffeners = 2*((height_wingbox(0)-2*t_sheet)*t_sheet*((width_wingbox(0)-t_sheet)/2)**2 + width_wingbox(0)**3*t_sheet/12)

    

def I_zz_wingbox(x):
    """Returns moment of inertia around the z-axis as a function of the spanwise position"""
    #neutral axis
    centroid_y = (2*(height_wingbox(x)/2*(height_wingbox(x)-2*t_sheet)*t_sheet) + width_wingbox(x)*t_sheet*((height_wingbox(x)-t_sheet/2)+t_sheet/2) + n_top*(height_wingbox(x)-t_sheet-y_top)*A_top + n_bottom*(t_sheet+y_bottom)*A_bottom) / (2*((height_wingbox(x)-2*t_sheet)*t_sheet + width_wingbox(x)*t_sheet) + n_top*A_top + n_bottom*A_bottom)    
    
    #moment of inertia
    dy_vertical_flange = np.abs(height_wingbox(x)/2 - centroid_y)
    dy_upper_flange = np.abs((height_wingbox(x)-t_sheet/2)-centroid_y)
    dy_bottom_flange = np.abs(centroid_y - t_sheet/2)
    
    dy_top_stiffener = np.abs(centroid_y - (height_wingbox(x)-t_sheet-y_top))
    dy_bottom_stiffener = np.abs(centroid_y - (t_sheet + y_bottom))
    
    Izz = t_sheet*(height_wingbox(x)-2*t_sheet)**3/12 + n_top*I_zz_top + n_bottom*I_zz_bottom + 2*(height_wingbox(x)-2*t_sheet)*t_sheet*dy_vertical_flange**2 + width_wingbox(x)*t_sheet*dy_upper_flange**2 + width_wingbox(x)*t_sheet*dy_bottom_flange**2 + n_top*A_top*dy_top_stiffener**2 + n_bottom*A_bottom*dy_bottom_stiffener**2
    
    return Izz



def I_yy_wingbox(x):
    """Returns moment of inertia around the y-axis as a function of the spanwise position"""
    top_spacing = (width_wingbox(0) - n_top * (width_top))/(n_top+1)
    bottom_spacing = (width_wingbox(0) - n_bottom * (width_bottom))/(n_bottom+1)

    #moment of inertia
    dz_vertical_flange = np.abs(width_wingbox(x)-t_sheet)/2
    dz_horizontal_flange = 0
    
    #top stiffeners
    if n_top%2==0:
        Adz_top = 0
        for i in range(1,int(n_top/2)+1):
            Adz_top += A_top*(i*((top_spacing + z_top)/2))**2
    else:
        Adz_top = 0
        for i in range(1,int(n_top/2)+1):
            Adz_top += A_top*((i-1)*(top_spacing+z_top))**2
     
    #bottom stiffeners
    if n_bottom%2==0:
        Adz_bottom = 0
        for i in range(1,int(n_bottom/2)+1):
            Adz_bottom += A_bottom*(i*((bottom_spacing + z_bottom)/2))**2
    else:
        Adz_bottom = 0
        for i in range(1,int(n_bottom/2)+1):
            Adz_bottom += A_bottom*((i-1)*(bottom_spacing+z_bottom))**2
     
    
    I_yy_nostiffeners = 2*((height_wingbox(x)-2*t_sheet)*t_sheet*dz_vertical_flange**2+width_wingbox(x)**3*t_sheet/12)
    I_yy = I_yy_nostiffeners + 2*Adz_top + 2*Adz_bottom + n_top*I_yy_top + n_bottom*I_yy_bottom
    
    return I_yy
    

def first_moment_of_area(x):
    """Returns Q along spanwise direction of the wingbox"""

    A_top = A_hat
    A_bottom = A_z
       
    #neutral axis
    centroid_y = (2*(height_wingbox(x)/2*(height_wingbox(x)-2*t_sheet)*t_sheet) + width_wingbox(x)*t_sheet*((height_wingbox(x)-t_sheet/2)+t_sheet/2) + n_top*(height_wingbox(x)-t_sheet-y_top)*A_top + n_bottom*(t_sheet+y_bottom)*A_bottom) / (2*((height_wingbox(x)-2*t_sheet)*t_sheet + width_wingbox(x)*t_sheet) + n_top*A_top + n_bottom*A_bottom)    
   
    
    #areas above and below neutral axis
    A_toppart = width_wingbox(x)*t_sheet + (height_wingbox(x) - centroid_y - t_sheet)*t_sheet*2 + n_top*A_top
    A_bottompart = width_wingbox(x)*t_sheet + (centroid_y - t_sheet)*t_sheet*2 + n_bottom*A_bottom
        
    
    #Q's for top and bottom, one as verification of the other, as Q is max at n.a. so Qbottom should be equal to Qtop
    sumproduct_yA_toppart = width_wingbox(x)*t_sheet*(height_wingbox(x)-centroid_y-t_sheet/2) + 0.5*(height_wingbox(x) - centroid_y - t_sheet)**2*t_sheet*2 + n_top*A_top*(height_wingbox(x)-centroid_y-t_sheet-y_top)
    sumproduct_yA_bottompart = width_wingbox(x)*t_sheet*(centroid_y-t_sheet/2) + 0.5*(centroid_y - t_sheet)**2*t_sheet*2 + n_bottom*A_bottom*(centroid_y-t_sheet-y_bottom)
    

    Q_toppart = sumproduct_yA_toppart/A_toppart
    Q_bottompart = sumproduct_yA_bottompart/A_bottompart    
    
    return Q_toppart

def cross_sectional_area(x):
    """Returns cross sectional area at spanwise position x"""
    centroid_y = (2*(height_wingbox(x)/2*(height_wingbox(x)-2*t_sheet)*t_sheet) + width_wingbox(x)*t_sheet*((height_wingbox(x)-t_sheet/2)+t_sheet/2) + n_top*(height_wingbox(x)-t_sheet-y_top)*A_top + n_bottom*(t_sheet+y_bottom)*A_bottom) / (2*((height_wingbox(x)-2*t_sheet)*t_sheet + width_wingbox(x)*t_sheet) + n_top*A_top + n_bottom*A_bottom)    

    A_toppart = width_wingbox(x)*t_sheet + (height_wingbox(x) - centroid_y - t_sheet)*t_sheet*2 + n_top*A_top
    A_bottompart = width_wingbox(x)*t_sheet + (centroid_y - t_sheet)*t_sheet*2 + n_bottom*A_bottom
    
    A_total = A_bottompart + A_toppart
    return A_total

def y_max(x):
    """Returns maximum y-distance from the neutral axis for a given cross section"""
    centroid_y = (2*(height_wingbox(x)/2*(height_wingbox(x)-2*t_sheet)*t_sheet) + width_wingbox(x)*t_sheet*((height_wingbox(x)-t_sheet/2)+t_sheet/2) + n_top*(height_wingbox(x)-t_sheet-y_top)*A_top + n_bottom*(t_sheet+y_bottom)*A_bottom) / (2*((height_wingbox(x)-2*t_sheet)*t_sheet + width_wingbox(x)*t_sheet) + n_top*A_top + n_bottom*A_bottom)    
   
    if height_wingbox(x)/2 > centroid_y:
#        print('y max is at the top, i.e. has a positive value: ') 
        return height_wingbox(x) - centroid_y
    else:
#        print('y max is at the bottom, i.e. has a negative value: ') 
        return -centroid_y
   
    
#print(first_moment_of_area(0))

#x = np.linspace(0,16,50)
#plt.plot(x,first_moment_of_area(x))
#plt.show()
        
    
    
    
    


