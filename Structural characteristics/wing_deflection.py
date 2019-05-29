# -*- coding: utf-8 -*-
"""
Created on Wed May 29 09:14:33 2019

@author: rober
"""

import matplotlib.pyplot as plt
import numpy as np
import parameters as p

#p0 = 20
#L = 14.88
#base = 0
#
#
#def f(x):
#    return p0*np.cos(np.pi*x/(2*L))
#
#
#x = np.linspace(0,L,100)
#plt.plot(x,f(x)+base)
#plt.xlim(right=20)
#plt.ylim(top=100)
#plt.ylim(bottom=0)
#plt.show()


def defl(x):
    """Deflection function without 1/EI term"""
    defl_W_eng = W_eng*x_eng**3/(3) + W_eng*x_eng**2*(x-x_eng)/2
    
    defl_W_pod = W_pod*x_pod**3/3 + W_pod*x_pod**2*(x-x_pod)/2
    
    defl_W_fuel = W_fuel*x_fuel**3/3 + W_fuel*x_fuel**2*(x-x_fuel)/2
    
    defl_F_strut = F_sy*x_strut**3/3 + F_sy*x_strut**2*(x-x_strut)/2
    
    defl_W_wing = taper*W_wing*x**2*(6*L**2-4*L*x+x**2)/24 + (1-taper)*W_wing*x**2*(20*L**3-10*L**2*x+5*L*x**2-x**3)/(120*L)
    
    defl_lift = lift*L/(3*np.pi**4)*(48*L**3*np.cos(np.pi*x/(2*L))-48*L**3+3*np.pi**3*L*x**2-np.pi**3*x**3)
    
    return defl_W_eng + defl_W_pod + defl_W_fuel + defl_F_strut + defl_W_wing + defl_lift
