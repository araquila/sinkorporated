#pieter responsible
from math import *
from wingloadingfunction.py import *
#load data
n_max_flap = 2
n_max_clean = 2.5
n_min = -1
s_landing = 1400 #[m]
rho0 = 1.225 #[kg/m3]

#turboprop data
W = 10000 #[N]
S = 10 #[m2]

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
V_stall = V_stall_calc(W,rho0,CLmax_turboprop_take_max,S)
W_S_stall = W_S_calc(rho0,V_stall,CLmax_turboprop_take_max)
print W_S_stall
