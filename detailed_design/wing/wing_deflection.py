# -*- coding: utf-8 -*-
"""
Created on Wed May 29 09:14:33 2019

@author: robert
"""

import matplotlib.pyplot as plt
import numpy as np
import parameters as p
from scipy.integrate import quad


def func(x):
    return np.cos((np.pi*x) / (2*L))

L = p.b/2

F_sy = 100000
safetyfactor = 2.5

lift = p.MTOW*safetyfactor/2 #per wing

lift_dist = lift / quad(func,0,L)[0]


def defl(x):
    """Deflection function without 1/EI term"""
    defl_W_eng = p.W_engine*p.x_engine**3/(3) + p.W_engine*p.x_engine**2*(x-p.x_engine)/2

    defl_W_fuel_pod = (p.W_fuel/2+p.W_pod)*p.x_pod**3/3 + (p.W_fuel/2+p.W_pod)*p.x_pod**2*(x-p.x_pod)/2

    defl_F_strut = F_sy*p.x_strut**3/3 + F_sy*p.x_strut**2*(x-p.x_strut)/2

    defl_W_wing = p.taper*p.W_wing/L*x**2*(6*L**2-4*L*x+x**2)/24 + (1-p.taper)*p.W_wing/L*x**2*(20*L**3-10*L**2*x+5*L*x**2-x**3)/(120*L)
    
    defl_lift = lift_dist*L/(3*np.pi**4)*(48*L**3*np.cos(np.pi*x/(2*L))-48*L**3+3*np.pi**3*L*x**2-np.pi**3*x**3)
    
    return -(defl_W_eng + defl_W_fuel_pod + defl_F_strut + defl_W_wing - defl_lift)
 #   return "engine",defl_W_eng,"pod",defl_W_pod,"fuel", defl_W_fuel ,"strut", defl_F_strut,"wing", defl_W_wing,"lift",defl_lift



    
    

print(defl(0))
print(defl(p.b/2))
x = np.linspace(0,p.b/2,100)
plt.figure()
plt.plot(x,defl(x))

plt.show()