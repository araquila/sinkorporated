# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 15:01:52 2019

@author: Stijn
"""

import numpy as np
import parameters as p



b = p.b

noseangle = np.arange(30, 71, 5)

NLGpos = 2
MLGpos = 10.49
latW = 1.49


distTurnpoint = (MLGpos-NLGpos)*np.tan((90-noseangle)/180*np.pi)

radiusLwheel = distTurnpoint - latW
radiusRwheel = distTurnpoint + latW
radiusNwheel = np.sqrt(distTurnpoint**2 + (MLGpos - NLGpos)**2)
radiusLwing = b/2 - distTurnpoint
radiusRwing = b/2 + distTurnpoint

