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

n_ribs = p.n_ribs
t_wing_skin = 0.001
t_sheet = p.t_sheet
t_rib = p.t_rib
chord_avg = (p.taper-1)*p.root_chord * 0.5 + p.root_chord


I_xx_rib = t_rib*chord_avg**3/12 + t_rib*chord_avg*(0.12*chord_avg)**12

I_xx_wing = t_wing_skin * p.b/2 * (0.38*chord_avg)**2 + 2 * (t_sheet*p.b/2 * (0.23*chord_avg)**2) + n_ribs * (t_rib*chord_avg**3/12 + t_rib*chord_avg*(0.12*chord_avg)**2) + t_wing_skin*p.b/2 * (0.62*chord_avg)**2  
area_xx_wing = 2 * t_wing_skin * p.b/2 + 2 * t_sheet *p.b/2 + n_ribs * t_rib * chord_avg


m = p.W_wing
Clalpha = 5.65 #per radian

x_theta = 0.015
chord = (p.taper-1)*p.root_chord * 0.75 + p.root_chord
chord_avg = (p.taper-1)*p.root_chord * 0.5 + p.root_chord
b = chord/2
S = 2*b
L = 0.75*p.b/2
m = p.W_wing / p.g

r_theta = np.sqrt(I_xx_rib/area_xx_wing)
I_theta = (r_theta*b)**2*m
S_theta = m*x_theta*b

E = p.E_sheet
I_zz = sp.I_zz_wingbox(p.b/2*0.75)

G = p.G_sheet
J = sp.I_zz_wingbox(L) + sp.I_yy_wingbox(L)

K_h = (8*E*I_zz)/(L**4)
K_theta = G*J/L


e = 0.125

q = 0.5*p.rho*p.V_cruise**2
V_cruise = np.sqrt(2*q/p.rho)
print("V at cruise; ",V_cruise)
print("")

length = 1000
qrange = np.linspace(0,90*q,length)

q_divergence_torsion = K_theta/(Clalpha*e*chord*S)
print("Torsional divergence at q:",q_divergence_torsion)
print("Corresponding V:", np.sqrt(2*q_divergence_torsion/p.rho))
print("")

real_1 = np.zeros(length)
real_2 = np.zeros(length)
real_3 = np.zeros(length)
real_4 = np.zeros(length)

imag_1 = np.zeros(length)
imag_2 = np.zeros(length)
imag_3 = np.zeros(length)
imag_4 = np.zeros(length)

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
    
    real_1[i] = np.real(roots[0])
    imag_1[i] = np.imag(roots[0])
    
    if i > 2:
        if real_1[i] - real_1[i-1] > 0.1:
            print("Flutter at q: ",q)
            q_f = q
            break
    
    real_2[i] = np.real(roots[1])
    imag_2[i] = np.imag(roots[1])
    
    real_3[i] = np.real(roots[2])
    imag_3[i] = np.imag(roots[2])
    
    real_4[i] = np.real(roots[3])
    imag_4[i] = np.imag(roots[3])
    
#plt.scatter(qrange,real_1)
##plt.scatter(qrange,real_2)
##plt.scatter(qrange,real_3)
##plt.scatter(qrange,real_4)
#
##plt.scatter(qrange,imag_1)
##plt.scatter(qrange,imag_2)
##plt.scatter(qrange,imag_3)
##plt.scatter(qrange,imag_4)
#
#plt.show()

V_f_cruise = np.sqrt(2*q/p.rho)
print("Corresponding V: ",V_f_cruise)
print("")

def uncoupled_bending_frequency(x):
    K_h = (8*E*sp.I_zz_wingbox(x))/((L)**4)
    return np.sqrt(K_h/m)

def uncoupled_torsional_frequency(x):
    K_theta = K_theta = G*(sp.I_yy_wingbox(x)+sp.I_zz_wingbox(x))/L
    return np.sqrt(K_theta/I_theta)

#print(uncoupled_bending_frequency(0.75*L))
#print(uncoupled_torsional_frequency(0.75*L))


condition_1 = x_theta*(x_theta + 2*e - 2*e*(uncoupled_bending_frequency(L)/uncoupled_torsional_frequency(L))**2*(1+2*e*x_theta/(r_theta**2)))
condition_2 = x_theta + 2*e + (uncoupled_bending_frequency(L)/uncoupled_torsional_frequency(L))**2*(x_theta-2*e+4*e*(x_theta/r_theta)**2)

print("Pine's flutter condition based on wing geometry, if both positive no flutter:")
print(condition_1)
print(condition_2)
print("")




length = 1000
rho_list = np.linspace(0.1,1.225,length)
V_f_list = np.zeros(length)


for i in range(length):
    V_f_list[i] = np.sqrt(2*q_f/(rho_list[i]))
    
rho_list = rho_list[::-1]
plt.plot(V_f_list, rho_list)
plt.axhline(p.rho)
plt.ylim(1.225, 0.1)
plt.xlim(0,1000)
plt.show()
    

















