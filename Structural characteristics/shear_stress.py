# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 10:01:19 2019

@author: Stijn
"""
from matplotlib import pyplot as plt

t = 0.002
xi = [0, 1]
h = [0.2, 0.2]
b = [2, 2]
Izz = [10**-5, 10**-5]
Iyy = [10**-5, 10**-5]

sheary_root = 87096.88196835331
sheary_1 = -11901.700560542904
shearz_root = 0
shearz_1 = 0
sheary = [sheary_root, sheary_1]
shearz = [shearz_root, shearz_1]
T = [0, 0]
taumax = []

for i in range(len(xi)):
    step = 0.001
    s1 = np.arange(0, b[i]/2 + step , step)
    s2 = np.arange(0, h[i] + step, step)
    s3 = np.arange(0, b[i]/2 + step, step)
    dz1 = (1/8*t*b[i]**3 + 1/4*t*b[i]*h[i]**2 + 3/8*t*b[i]**2*h[i])
    dzeta = 1/Iyy[i]*(1/4*t*b[i]*h[i]**2 + 1/8*t*b[i]**2*h[i] - h[i]/(h[i]+b[i])*dz1 + 1/4*b[i]**3*t - 1/8*t*b[i]**3 + 1/2*t*b[i]**2*h[i] + 1/8*t*b[i]**3 - b[i]/(b[i]+h[i])*dz1)
    eta = 0
    ShC = [eta, dzeta]
    
    qs0 = -sheary[i]/Izz[i]*(1/8*b[i]*t + 1/4*h[i]*t + 1/12*h[i]**2*t/b[i]) - shearz[i]/Iyy[i]*(1/8*b[i]**2/h[i]*t + 3/8*b[i]*t + 1/4*h[i]*t) + shearz[i]*dzeta/(2*h[i]*b[i])
    qT = T[i]/(2*(h[i]*b[i]))
    
    qb12 = sheary[i]/Izz[i]*(1/2*h[i]*t*s1) + shearz[i]/Iyy[i]*(1/2*t*(s1)**2)
    qb23 = -sheary[i]/Izz[i]*(1/2*t*(s2)**2 - 1/2*t*h[i]*s2 - 1/4*h[i]*b[i]*t) - shearz[i]/Iyy[i]*(-1/2*t*b[i]*s2 - 1/8*t*b[i]**2)
    qb34 = -sheary[i]/Izz[i]*(h[i]/2*t*s3 - 1/4*h[i]*b[i]*t) - shearz[i]/Iyy[i]*(-b[i]/2*t*s3 + 1/2*t*(s3)**2 - 1/2*t*b[i]*h[i] - 1/8*t*b[i]**2)
    
    q12 = qb12 + qs0 + qT
    q23 = qb23 + qs0 + qT
    q34 = qb34 + qs0 + qT
    
    tau12 = q12/t
    tau23 = q23/t
    tau34 = q34/t
    
    maxtau12 = max(tau12), min(tau12)
    maxtau23 = max(tau23), min(tau23)
    maxtau34 = max(tau34), min(tau34)
    
    maxtau = max(maxtau12[0], maxtau23[0], maxtau34[0]), min(maxtau12[1], maxtau23[1], maxtau34[1])
    taumax.append(max(abs(maxtau[0]), abs(maxtau[1])))

#plt.figure()
#plt.plot(s1, q12)
#plt.plot(s2, q23)
#plt.plot(s3, q34)
#plt.show()