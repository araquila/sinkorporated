# -*- coding: utf-8 -*-
"""
Created on Wed May 29 09:14:33 2019

@author: rober
"""

import matplotlib.pyplot as plt
import numpy as np
import parameters as p



L = p.b/2
F_sy = 10000
safetyfactor = 1.1
lift = p.MTOW*safetyfactor/2 #per wing

def defl(x):
    """Deflection function without 1/EI term"""
    defl_W_eng = p.W_engine*p.x_engine**3/(3) + p.W_engine*p.x_engine**2*(x-p.x_engine)/2

    defl_W_pod = p.W_pod*p.x_pod**3/3 + p.W_pod*p.x_pod**2*(x-p.x_pod)/2

    defl_W_fuel = (p.W_fuel+p.W_pod)*p.x_pod**3/3 + (p.W_fuel+p.W_pod)*p.x_pod**2*(x-p.x_pod)/2

    defl_F_strut = F_sy*p.x_strut**3/3 + F_sy*p.x_strut**2*(x-p.x_strut)/2

    defl_W_wing = p.taper*p.W_wing*x**2*(6*L**2-4*L*x+x**2)/24 + (1-p.taper)*p.W_wing*x**2*(20*L**3-10*L**2*x+5*L*x**2-x**3)/(120*L)
    
    defl_lift = lift*L/(3*np.pi**4)*(48*L**3*np.cos(np.pi*x/(2*L))-48*L**3+3*np.pi**3*L*x**2-np.pi**3*x**3)
    
    return -(defl_W_eng + defl_W_pod + defl_W_fuel + defl_F_strut + defl_W_wing - defl_lift)

print(defl(p.b/2))
x = np.linspace(0,p.b/2,100)
plt.figure()
plt.plot(x,defl(x))

plt.show()