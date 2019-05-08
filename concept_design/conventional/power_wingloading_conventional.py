#pieter responsible
from math import *
from wingloadingfunctions import *
import numpy as np
#load data
n_max_flap = 2
n_max_clean = 2.5
n_min = -1
s_landing = 1400 #[m]
rho0 = 1.225 #[kg/m3]

#graph data
W_S_x = np.linspace(0,4000,200)
#turboprop data
W = 200000 #[N]
S = 55 #[m2]

#data props
CLmax_turboprop_clean_min = 1.5
CLmax_turboprop_clean_max = 1.9
CLmax_turboprop_take_min = 1.7
CLmax_turboprop_take_max = 2.1
CLmax_turboprop_land_min = 1.9
CLmax_turboprop_land_max = 3.3

#take off parameter turboprop
TOP_aquila_turboprop = 580
#data jets
CLmax_jet_clean_min = 1.2
CLmax_jet_clean_max = 1.8
CLmax_jet_take_min = 1.6
CLmax_jet_take_max = 2.2
CLmax_jet_land_min = 1.8
CLmax_jet_land_max = 2.8

#take off parameter jet
TOP_aquila_jet_single = 240
TOP_aquile_jet_double = 100

#props, for jets scroll DOWN################
#calculate stall speeds and the wing loading
V_stall = V_stall_calc(W,rho0,CLmax_turboprop_take_max,S)
W_S_stall = W_S_calc(rho0,V_stall,CLmax_turboprop_take_max)

##########take-off################
k = TOP_aquila_turboprop
CL_TO_min = CL_TO_calc(CLmax_turboprop_take_min)
CL_TO_max = CL_TO_calc(CLmax_turboprop_take_max)
CL_TO_range = np.linspace(CL_TO_min,CL_TO_max,5)
TOP_takeoff = np.zeros(shape=(len(CL_TO_range),len(W_S_x)))
for i in range(len(TOP_takeoff)):
    for j in range(len(W_S_x)):
        TOP_takeoff[i,j] = W_P_calc(W_S_x[j],k,CL_TO_range[i])
print(TOP_takeoff)
#jets
