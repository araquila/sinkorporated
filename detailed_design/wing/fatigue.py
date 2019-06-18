# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 12:39:14 2019

@author: robert
"""

import parameters as p
import numpy as np

strain_f =
T_s = 
E = 

stress_amplitude = 



b_4p = -.1785*np.log10(2.78*(1+strain_f))
M_1 = T_s/E
M_2 = np.log10(2.5*(M_1*(1+strain_f))
C_E_4p = 0.5*10**(0.301*b_4p+M_2) * 2**b_4p
M_3 = 10**(4.602*b_4p+M_2)
M_4 = np.log10(0.25*strain_f**0.75)

c_4p = 0.333*(np.log10(0.00691-0.52356*M_3)-M_4)
C_P_4p = 0.5*10**(-1.301*c_4p+M_4) * 2**(c_4p)

N_f = 72,000



strain_amplitude = C_E_4p*N_F**b_4p + C_P_4p*N_f**c_4p

