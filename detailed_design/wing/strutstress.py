# -*- coding: utf-8 -*-
"""
Created on Wed May 29 16:44:05 2019

@author: robert
"""

import numpy as np
import matplotlib.pyplot as plt
import parameters as p


M_root = 1000
Ry = 20000

q_w = p.W_wing/(p.b/2)

#root-strut
V_rs = Ry - p.W_engine - q_w*p.x_strut
M_rs = Ry*p.x_strut - q_w*p.x_strut**2/2 - p.W_engine*(p.x_strut - p.x_engine) - M_root

#strut-tip
V_st = (p.W_fuel+p.W_pod) + q_w*(p.b/2 - p.x_strut)
M_st = (p.W_fuel+p.W_pod)*(p.x_pod - p.x_strut) + q_w*(L-p.x_strut)**2/2


print(M_rs)
print(M_st)