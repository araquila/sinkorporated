# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 17:17:50 2019

@author: Stijn
"""
from math import *
import numpy as np
#import parameters as p
from matplotlib import pyplot as plt
from wing_deflection2 import read_aero_data, CallForces
#from section_proporties import 

V_cruise = 184.84
rho_cruise = 0.5258
lengthdata = 50
Lift, Chord, Yle, Drag = read_aero_data("aquiladata1.txt", lengthdata, V_cruise, rho_cruise)
Frx, Fry, Fs, Mr, momentyi, momentzi, shearyi, shearzi, vii = CallForces(Lift, Chord, Yle, Drag, np.ones(len(Lift)+1)*10**(-4), 70*10**9, 0.2, 0.4, 0.4)
