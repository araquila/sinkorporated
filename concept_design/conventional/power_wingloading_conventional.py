#pieter responsible
from math import *
#load data
n_max_flap = 2
n_max_clean = 2.5
n_min = -1
s_landing = 1400 #[m]

#data props
Clmax_turboprop_clean_min = 1.5
Clmax_turboprop_clean_max = 1.9
Clmax_turboprop_take_min = 1.7
Clmax_turboprop_take_max = 2.1
Clmax_turboprop_land_min = 1.9
Clmax_turboprop_land_max = 3.3

#take off parameter turboprop
TOP_aquila_turboprop = 580
#data jets
Clmax_jet_clean_min = 1.2
Clmax_jet_clean_max = 1.8
Clmax_jet_take_min = 1.6
Clmax_jet_take_max = 2.2
Clmax_jet_land_min = 1.8
Clmax_jet_land_max = 2.8

#take off parameter jet
TOP_aquila_jet_single = 240
TOP_aquile_jet_double = 100

#props, for jets scroll DOWN
V_landing = sqrt(s_landing/0.5847)
print V_landing
#jets
