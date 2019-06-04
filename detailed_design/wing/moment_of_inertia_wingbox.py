# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 14:48:26 2019

@author: robert
"""

import parameters as p
import numpy as np
import matplotlib.pyplot as plt


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

I_zz_hat = 2*(b*t_hat*(a+3*t_hat/2)**2 + a**3*t_hat/12 + a*t_hat*(a/2+t_hat)**2)


#Z-stiffener geometry
d = 0.04
e = 0.06 
f = 0.045

width_z = d+f-t_z

A_z = t_z*(d+e+f)

y_centroid_z = (d*t_z**2/2 + e*t_z*(t_z+e/2)+f*t_z*(e+3*t_z/2))/(t_z*(d+e+f))

I_zz_z = e**3*t_z*(e+t_z)**2 +f*t_z*(e+3*t_z/2)**2

#section properties
n_top = 0 # at the upper skin
n_bottom = 0 # at the lower skin

width_top = width_hat
width_bottom = width_z

y_top = y_centroid_hat
y_bottom = y_centroid_z

A_top = A_hat
A_bottom = A_z

I_zz_top = I_zz_hat
I_zz_bottom = I_zz_z

#height at root of the wing box
h_root = 0.929*p.h_max_root_wingbox
h_tip = 0.825*p.h_max_tip_wingbox #arbitrary

def width_wingbox(x):
    return p.w_root_wingbox - (2/p.b) * x * (p.w_root_wingbox - p.w_tip_wingbox) 

def height_wingbox(x):
    return h_root - (2/p.b) * x * (h_root - h_tip)

def y_neutral_line(x):
    return (p.n_upper_skin*(height_wingbox(x)-p.h_stiffener/2)*p.A_stiffener + p.n_lower_skin*(p.h_stiffener/2)*p.A_stiffener) / ((p.n_upper_skin+p.n_lower_skin)*p.A_stiffener)

#stiffener spacing
top_spacing = (width_wingbox(0) - n_top * (width_top))/(n_top+1)

lower_spacing = (width_wingbox(0) - n_bottom * (width_bottom))/(n_bottom+1)

centroid_root = (2*((t_sheet+height_wingbox(0)/2)*(height_wingbox(0)-2*t_sheet)*t_sheet) + (-t_sheet/2+height_wingbox(0))*width_wingbox(0)*t_sheet + t_sheet**2/2*width_wingbox(0) + n_top*(3*t_sheet/2+height_wingbox(0)-y_top)*A_top + n_bottom*(t_sheet+y_bottom)*A_bottom) / (2*(height_wingbox(0)-2*t_sheet)*t_sheet + 2*width_wingbox(0)*t_sheet + n_top * A_top + n_bottom * A_bottom)

#moment of inertia
dy_vertical_flange = np.abs(height_wingbox(0)/2 - centroid_root)
dy_upper_flange = np.abs((height_wingbox(0)-t_sheet/2)-centroid_root)
dy_bottom_flange = np.abs(centroid_root - t_sheet/2)

dy_top_stiffener = np.abs(centroid_root - (height_wingbox(0)-t_sheet-y_top))
dy_bottom_stiffener = np.abs(centroid_root - (t_sheet + y_bottom))

Izz = t_sheet*(height_wingbox(0)-2*t_sheet)**3/12 + n_top*I_zz_top + n_bottom*I_zz_bottom + 2*(height_wingbox(0)-2*t_sheet)*t_sheet*dy_vertical_flange**2 + width_wingbox(0)*t_sheet*dy_upper_flange**2 + width_wingbox(0)*t_sheet*dy_bottom_flange**2 + n_top*A_top*dy_top_stiffener**2 + n_bottom*A_bottom*dy_bottom_stiffener**2




def I_zz_wingbox(x,n_top,n_bottom):

    y_top = y_centroid_hat
    y_bottom = y_centroid_z
    
    A_top = A_hat
    A_bottom = A_z
    
    I_zz_top = I_zz_hat
    I_zz_bottom = I_zz_z
    
    centroid_root = (2*((t_sheet+height_wingbox(x)/2)*(height_wingbox(x)-2*t_sheet)*t_sheet) + (-t_sheet/2+height_wingbox(x))*width_wingbox(x)*t_sheet + t_sheet**2/2*width_wingbox(x) + n_top*(3*t_sheet/2+height_wingbox(x)-y_top)*A_top + n_bottom*(t_sheet+y_bottom)*A_bottom) / (2*(height_wingbox(x)-2*t_sheet)*t_sheet + 2*width_wingbox(x)*t_sheet + n_top * A_top + n_bottom * A_bottom)
    
    #moment of inertia
    dy_vertical_flange = np.abs(height_wingbox(x)/2 - centroid_root)
    dy_upper_flange = np.abs((height_wingbox(x)-t_sheet/2)-centroid_root)
    dy_bottom_flange = np.abs(centroid_root - t_sheet/2)
    
    dy_top_stiffener = np.abs(centroid_root - (height_wingbox(x)-t_sheet-y_top))
    dy_bottom_stiffener = np.abs(centroid_root - (t_sheet + y_bottom))
    
    Izz = t_sheet*(height_wingbox(x)-2*t_sheet)**3/12 + n_top*I_zz_top + n_bottom*I_zz_bottom + 2*(height_wingbox(x)-2*t_sheet)*t_sheet*dy_vertical_flange**2 + width_wingbox(x)*t_sheet*dy_upper_flange**2 + width_wingbox(x)*t_sheet*dy_bottom_flange**2 + n_top*A_top*dy_top_stiffener**2 + n_bottom*A_bottom*dy_bottom_stiffener**2
    
    return Izz




