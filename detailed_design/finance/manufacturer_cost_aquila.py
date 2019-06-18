# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 14:46:30 2019

@author: victo
"""
import numpy as np


#--------- INPUTS --------------#
R_eng = 120 #cost per hour for engineer in 2019 dollar
R_m = 70 #cost per hour for manufacturing in 2019 dollar
R_t = 90 #cost per hour for tooling in 2019 dollar
MTOW_kg = 18025.5 #Maximum Take-Off Weight in kg
MTOW_lbs = MTOW_kg*2.2046226 #Maximum Take-Off Weight in pounds
V_max_ms = 175 #Maximum Vc in m/s
V_max_kts = 1.94384449*V_max_ms #Maximum Vc in knots
N_rdte = 6 #Test airplanes in rdte
N_st = 3 #Number of static test aircraft
F_diff = 1.8 #Complexity factor
F_cad = 0.8 #Skill of CAD developers
CEF = 6.66 #Cost escalation factor
CEF90 = 3.37 #Cost escalation factor FOR 1990
C_e = 1000000 #Cost per engine in 2019 dollar
N_e = 2 #Number of engines in 2019 dollar
C_p = 150000 #Cost per propeller  in 2019 dollar
N_p = 2 #Number of propellers
C_avionics = 7800000#cost of avionics in 2019 dollar
F_mat = 2.2 #Correction factor for type of materials
N_rr = 0.33 #RDTE production rate in units per month
F_obs = 1.0 #Factor for stealth (1.0 for commercial aircraft)
F_tsf = 0.20 #Test and facilities estimation factor
F_pror = 0.10 #Profit margin
F_fin = 0.10 #Cost to finance RDTE


N_m = 630 # Number of airplanes built to standard
N_program = N_rdte +  N_m #Expected airplanes build in the whole program
F_int = 1000 #interior cost factor
N_pax = 60
N_rm = 5 #Manufacturing production rate in units per month


#---------------------------------------------------#
# RESEARCH, DEVELOPMENT, TEST AND EVALUATION COST   #
#---------------------------------------------------#

# AIRFRAME ENGINEERING AND DESIGN COST (RESEARCH)
W_ampr = 10**(0.1936+0.8645*np.log10(MTOW_lbs)) #Aeronautical Manufacturers Planning Report Weight (what the airframe manufacturer builds)
MHR_aedr = 0.0396 * (W_ampr**0.791) * (V_max_kts**1.526) * (N_rdte**0.183) * F_diff * F_cad #Engineering manhours required for phase 1-3
C_aedr = MHR_aedr * R_eng # Engineering cost for phase 1-3


# DEVELOPMENT SUPPORT AND TESTING COST
C_dst = 0.008325 * (W_ampr**0.873) * (V_max_kts**1.890) * (N_rdte**0.346) * F_diff * CEF # Cost for development support and testing cost


# FLIGHT TEST AIRPLANES COST
C_ear = (C_e * N_e + C_p * N_p + C_avionics)*(N_rdte - N_st) # Cost engines and avionics

MHR_manr = 28.984*(W_ampr**0.740)*(V_max_kts**0.543)*(N_rdte**0.524)*F_diff #Manufacturing manhours required for phase 1-3
C_manr = MHR_manr * R_m #Manufacturing cost of flight test airplanes

C_matr = 37.632*F_mat*(W_ampr**0.689)*(V_max_kts**0.624)*(N_rdte**0.792)*CEF #Material cost of flight test airplanes

MHR_toolr = 4.0127*(W_ampr**0.764)*(V_max_kts**0.899)*(N_rdte**0.178)*(N_rr*0.066)*F_diff #Tooling manhours required for phase 1-3
C_toolr = MHR_toolr * R_t #Tooling cost of flight test airplanes

C_qcr = 0.13*C_manr # #Quality control cost of flight test airplanes

C_fta = C_ear + C_manr + C_matr + C_toolr + C_qcr #Flight Test Airplanes cost


# FLIGHT TEST OPERATIONS COST
C_ftor = 0.001244 * (W_ampr**1.160) * (V_max_kts**1.371) * ((N_rdte-N_st)**1.281) * CEF * F_diff * F_obs # Flight test operations cost


C_rdte = (C_aedr + C_dst + C_fta + C_ftor)/(1-F_tsf-F_pror-F_fin)

# TEST AND SIMULATION FACILITIES COST
C_tsfr = F_tsf * C_rdte


# RDTE PROFIT
C_pror = F_pror * C_rdte


# COST TO FINANCE RDTE
C_finr = F_fin * C_rdte

assert round((C_aedr + C_dst + C_fta + C_ftor+C_tsfr+C_pror+C_finr),2) == round(C_rdte,2)

"""
print('Cost in million USD')
print('Airframe engineering and design cost: ',round(C_aedr/1000000,2))
print('Development support and testing cost: ',round(C_dst/1000000,2))
print('Flight test airplanes cost:           ',round(C_fta/1000000,2))
print('Flight test operations cost:          ',round(C_ftor/1000000,2))
print('Test and simulation facilities cost:  ',round(C_tsfr/1000000,2))
print('RDTE profit:                          ',round(C_pror/1000000,2))
print('Cost to finance RDTE phase:           ',round(C_finr/1000000,2))
print('--------------------------------------')
print('Total RDTE cost:                      ',round(C_rdte/1000000,2))
"""

#-----------------------------------------------------#
# MANUFACTURING AND ACQUISITION                       #
#-----------------------------------------------------#



# AIRFRAME ENGINEERING AND DESIGN COST (RESEARCH)
MHR_aedprogram = 0.0396*(W_ampr**0.791)*(V_max_kts**1.526)*(N_program**0.183)*F_diff*F_cad # Engineering manhours required for the entire program
C_aedm = MHR_aedprogram *R_eng - C_aedr # Engineering cost for phase 4


# AIRPLANE PROGRAM PRODUCTION COST
C_eam = (C_e * N_e + C_p * N_p + C_avionics)*N_m # Cost engines and avionics

C_intm = F_int * N_pax * N_m * CEF/CEF90 # Cost interior

MHR_manprogram = 28.984 * (W_ampr**0.74)*(V_max_kts**0.543)*(N_program**0.524)*F_diff #Total number of manhours for the production of N_program airplanes
C_manm = MHR_manprogram * R_m - C_manr # Labor cost for manufacturing N_m airplanes to standard

C_matprogram = 37.632*F_mat*(W_ampr**0.689)*(V_max_kts**0.624)*(N_program**0.792)*CEF # Total material cost of the program
C_matm = C_matprogram - C_matr # Cost of material for N_m airplanes 

MHR_toolprogram = 4.0127*(W_ampr**0.764)*(V_max_kts**0.899)*(N_program**0.178)*(N_rm**0.066)*F_diff
C_toolm = MHR_toolprogram * R_t - C_toolr

C_qcm = 0.13*C_manm # Quality control

C_apcm = C_eam + C_intm + C_manm + C_matm + C_toolm + C_qcm


# MANUFACTURING AND ACQUISITION COST
C_MAN = (C_aedm + C_apcm) / (1 - F_fin)
C_finm = C_MAN * F_fin
assert round((C_aedm + C_apcm+C_finm),2) == round(C_MAN,2)
C_PRO = 0.1 * C_MAN
C_ACQ = C_MAN + C_PRO
AEP = (C_MAN + C_PRO + C_rdte)/N_m

"""
print('Manufacturing Cost per aircraft:',round(C_MAN/N_m/1000000,2))
print('Acquisition Cost per aircraft  :',round((C_ACQ)/N_m/1000000,2))
print('Unit Cost                      :',round(AEP/1000000,2))
"""
