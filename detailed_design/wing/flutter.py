# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 15:47:29 2019

@author: robert
"""

import parameters as p
import numpy as np
import control.matlab as control
import wing.section_properties as sp
import matplotlib.pyplot as plt


m = p.W_wing
Clalpha = 5.65 #per radian
I_theta = 0.2 #kgm^2
x_theta = 0.15
chord = (p.taper-1)*p.root_chord * 0.75 + p.root_chord
b = chord/2
S = 2*b
L = 0.75*p.b/2
m = p.W_wing / p.g

S_theta = m*x_theta*b

E = p.E_sheet
I_zz = sp.I_zz_wingbox(p.b/2*0.75)

G = p.G_sheet
J = sp.I_zz_wingbox(L) + sp.I_yy_wingbox(L)

K_h = (3*E*I_zz)/((L)**3)
K_theta = G*J/L

e = 0.125*chord

q = 0.5*p.rho*p.V_cruise**2

length = 20000
qrange = np.linspace(0,2*q,length)

real = np.zeros(4*length)
imag = np.zeros(4*length)

def roots_flutter(q):
    a0 = K_h*(K_theta-2*e*b*q*S*Clalpha)
    a2 = m*K_theta + I_theta*K_h - (2*m*e*b+S_theta)*q*S*Clalpha
    a4 = m*I_theta - S_theta**2
    
    coeff = [a4, 0, a2, 0, a0]
    
    roots = np.roots(coeff)
    
    return roots
    
for i in range(length):
    q = qrange[i]
    roots = roots_flutter(q)

    real[i] = np.real(roots[0])
    imag[i] = np.imag(roots[0])
    
    real[i*2] = np.real(roots[1])
    imag[i*2] = np.imag(roots[1])
    
    real[i*3] = np.real(roots[2])
    imag[i*3] = np.imag(roots[2])
    
    real[i*4] = np.real(roots[3])
    imag[i*4] = np.imag(roots[3])
    
plt.plot(imag)
plt.show()
