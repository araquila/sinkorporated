# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 14:34:18 2019

@author: robert
"""

import parameters as p


batch_size = 300

#cost of aircraft components

C_wing = 1730 * p.W_wing**0.766 * batch_size ** -.218
C_tail = 1820 * p.W_wing**0.766 * batch_size ** -.218
C_fuselage = 2060 * p.W_fuselage**0.766 * batch_size ** -.218
C_lg = 1180 * 
