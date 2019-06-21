# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 10:55:29 2019

@author: victo
"""

#import matplotlib as mpl
#mpl.rcParams['font.size'] = 20.0
import matplotlib.pyplot as plt
import numpy as np
import detailed_design.parameters as p
import detailed_design.finance.manufacturer_cost_aquila as manu_cost
mtomlist=np.arange(0.8,1.3,0.05)


doclist=[]
for i in range(len(mtomlist)):
    #-----------------------------------------------------#
    # INPUTS                                              #
    #-----------------------------------------------------#
    
    #--------- AIRCRAFT INPUTS --------------#
    MTOW_kg = mtomlist[i] #Maximum Take-Off Weight in kg
    V_cr_ms = p.V_cruise #Cruise speed in m/s
    V_cl_ms = p.V_climb # Climb speed in m/s
    V_de_ms = p.V_descend # Descend speed in m/s
    t_cl = p.t_cl # Time to climb in hours
    t_de = p.t_de # Time to descend in hours
    n_c1 = 1 # Number of captains
    n_c2 = 1 # Number of co-pilots
    W_Fbl_kg = 1496 # Block fuel in kg for 1000 km trip
    
    EP = p.EP # Engine cost per engine
    FP = p.FP # Fuel price per lbs in dollar
    N_e = p.n_engines #Number of engines
    PP = p.PP #Cost per propeller  in 2019 dollar
    N_p = p.n_engines #Number of propellers
    AEP = 24500000 # Unit cost airplane
    AFP = AEP - EP* N_e # Airframe cost
    ASP = p.ASP #cost of avionics in 2019 dollar
    P_TO_1ENG = p.P_TO/2 # Shaft Power per engine in kW
    W_A_kg = p.m_A # airframe weight in kg
    F_disp = 0.02
    
    
    ATF = 1.0 # Airplane type factor (1.0 for MTOW > 10,000 lbs)
    
    #--------- UNIT CONVERSIONS --------------#
    MTOW_lbs = MTOW_kg*2.2046226 #Maximum Take-Off Weight in pounds
    W_Fbl = W_Fbl_kg*2.2046226 # Block fuel in lbs
    W_A = W_A_kg*2.2046226 # airframe weight in lbs
    V_cr_kts = 1.94384449*V_cr_ms #Cruise speed in kts
    V_cl_kts = 1.94384449*V_cl_ms # Climb speed in kts
    V_de_kts = 1.94384449*V_de_ms # Descend speed in kts
    SHP_TO = P_TO_1ENG*0.0013410220888 # Shaft Horse Power per engine
    #--------- GENERAL AIRCRAFT INPUTS --------------#
    N_yr = 20 # operating years aircraft
    R_bl = 540 # block distance in NM
    
    #--------- STATISTICAL INPUTS --------------#
    # COST ESCALATION FACTORS
    CEF = 6.66 #Cost escalation factor
    CEF89 = 3.20 #Cost escalation factor FOR 1989
    CEF90 = 3.37 #Cost escalation factor FOR 1990
    
    # CREW FINANCIALS
    K_j = 0.26 # Extra factor for vacation and shit
    AH_j = 900 # Flight hours a year for crew member
    TEF_j = 7 * CEF/CEF90 # Travel Expense Factor
    SAL_1 = 86500# Yearly salary captain
    SAL_2 = 47500# Yearly salary co-pilot
    
    # OIL PRICES
    OLP = 53.05 # Lubricating oil price per gallon in dollar
    OD = 7.74 # Lubricating oil density in lbs/gallon
    
    f_inshull = 0.0175 # Insurance rate
    H_em = 4000 # Time in hours between engine overhauls
    R_l = 22.44 # Mechanic labor rate dollar per hour
    
    ESPPF = 1.3 # Engine spare parts price factor
    f_amblab = 1.2 #Overhead distribution factors labor
    f_ambmat = 0.55 #Overhead distribution factors materials 
    
    F_dap = 0.85 # Depreciation factor airframe
    DP_ap = 10 #Suggested depreciation period airframe
    F_deng = 0.85 # Depreciation factor engines
    DP_eng = 7 #Suggested depreciation period engines
    F_dprp = 0.85 # Depreciation factor propellers
    DP_prp = 7 #Suggested depreciation period propellers
    F_dav = 1.00 # Depreciation factor avionics
    DP_av = 5 #Suggested depreciation period avionics
    F_dapsp = 0.85 # Depreciation factor airframe spares
    F_apsp = 0.10 # Airframe spares factor
    DP_apsp = 10 #Suggested depreciation period airplane spares
    F_dengsp = 0.85 # Depreciation factor engine spares
    F_engsp = 0.30 # Engine spares factor
    DP_engsp = 7 #Suggested depreciation period engine spares
    
    #-----------------------------------------------------#
    # TOTAL ANNUAL BLOCK MILES FLOWN                      #
    #-----------------------------------------------------#
    
    # BLOCK TIME
    t_gm = 0.51 * (10**-6) * MTOW_lbs + 0.125 # Time for ground maneuvers in hours
    
    R_cl = V_cl_kts*t_cl # Horizontal climb distance in NM
    R_de = V_de_kts*t_de # Horizontal descend distance in NM
    V_man = V_cr_kts # Maneuvering speed for ATC commands
    t_man = 0.25 * (10**-6) * MTOW_lbs + 0.0625 # Time spend for ATC commands
    R_man = V_man * t_man # Horizontal distance for ATC commands
    t_cr = (1.06 * R_bl - R_cl - R_de + R_man)/V_cr_kts # Time for cruise
    t_flt = t_cl + t_cr + t_de # Time in flight
    t_bl = t_gm + t_cl + t_cr + t_de # block time in hours
    
    
    # BLOCK SPEED
    V_bl = R_bl/t_bl # Block speed in kts
    
    
    # ANNUAL UTILIZATION
    #U_annbl = (10**3) * (3.4546 * t_bl + 2.994 - (12.289 * t_bl**2 - 5.6626 * t_bl + 8.964)**0.5) # Annual utilization in block hours
    U_annbl = 3240
    
    # TOTAL ANNUAL BLOCK MILES FLOWN   
    R_blann = V_bl * U_annbl # Total annual block miles
    
    
    #-----------------------------------------------------#
    # DIRECT OPERATING COST PER NAUTICAL MILE             #
    #-----------------------------------------------------#
    
    
    # COST OF FLYING
    C_crew = ((n_c1*((1 + K_j)/V_bl)*(SAL_1/AH_j)) + (TEF_j/V_bl)) + ((n_c2*((1 + K_j)/V_bl)*(SAL_2/AH_j)) + (TEF_j/V_bl)) # Cost of pilot crew
    
    W_olbl = 0.70 * N_e * t_bl # Lubricants consumption in lbs
    C_olbl = (W_olbl/R_bl)*(OLP/OD) 
    C_fuel = (W_Fbl/R_bl)*(FP) #cost of fuel per NM
    C_pol = C_olbl + C_fuel # Fuel and lubricants cost in dollar per NM
    
    C_ins = (f_inshull*AEP)/(U_annbl*V_bl) # Insurance cost in dollar per NM
    
    DOC_flt = C_crew + C_pol + C_ins #cost of flying in dollar per NM
    
    # COST OF MAINTENANCE
    MHR_mapbl = 3.0 + 0.067*W_A/1000
    C_labap = 1.03 * MHR_mapbl * R_l / V_bl # Cost of airframe maintenance per NM
    
    MHR_mengbl = (0.4956+0.0532*(SHP_TO/N_e)/1000)*(1100/H_em)+0.10 # Turbine engine maintenance manhours
    C_labeng = 1.03 * 1.3 * N_e * MHR_mengbl * R_l / V_bl # Cost of engine maintenance per NM 
    
    C_matapblhr = 30.0 * CEF/CEF89 * ATF + 0.79 * (10**-5) * AFP # Material cost airframe maintenance per block hour
    C_matap = 1.03 * C_matapblhr/V_bl # Material cost airframe maintenance per NM
    
    K_Hem = 0.021 * (H_em/100) + 0.769 # Attained period between engine overhauls factor
    C_matengblhr = ((5.43*(10**-5)*EP)*ESPPF-0.47)/K_Hem # Material cost engine maintenance per block hour
    C_mateng = 1.03*1.3*N_e*C_matengblhr/V_bl # Material cost engine maintenance per NM
    
    C_amb = 1.03 * (f_amblab * (MHR_mapbl*R_l + N_e * MHR_mengbl * R_l) + f_ambmat * (C_matapblhr + N_e * C_matengblhr)) / V_bl # cost of applied maintenance burder
    
    DOC_maint = C_labap + C_labeng + C_matap + C_mateng + C_amb # Cost of maintenance
    
    # COST OF DEPRECIATION
    C_dap = (F_dap * (AEP-(N_e*EP)-N_p*PP) - ASP)/(DP_ap*U_annbl*V_bl) #Depreciation cost airframe
    
    C_deng = (F_deng*N_e*EP)/(DP_eng*U_annbl*V_bl) #Depreciation cost engines
    
    C_dprp = (F_dprp*N_p*PP)/(DP_prp*U_annbl*V_bl) #Depreciation cost propellers
     
    C_dav = (F_dav*ASP)/(DP_av*U_annbl*V_bl)#Depreciation cost avionics
      
    C_dapsp = (F_dapsp*F_apsp*AEP-N_e*EP)/(DP_apsp*U_annbl*V_bl) #Depreciation cost airframe spares
       
    C_dengsp = (F_dengsp*F_engsp*N_e*EP*ESPPF)/(DP_engsp*U_annbl*V_bl) #Depreciation cost engine spares
    
    DOC_depr = C_dap + C_deng + C_dprp + C_dav + C_dapsp + C_dengsp #Total epreciation cost
    
    
    # COST OF LANDING FEES, NAVIGATION FEES
    C_lf = 0.005 * MTOW_kg / (V_bl/t_bl) # Landing per NM in USD (Taken from airbus) 
    C_nf = 0.006 * (MTOW_kg**0.5) # Navigation fees per NM in USD (Taken from airbus) 
    
    DOC_lnr = C_lf + C_nf
    
    
    # COST OF FINANCING
    TI = (AFP*(1+F_apsp)+N_e*EP*(1+F_engsp))
    DOC_fin = ((TI * t_bl)/(14 * U_annbl) + (0.03 * TI * t_bl)/U_annbl)/R_bl
    
    
    # TOTAL DIRECT OPERATIONAL COST
    DOC = DOC_flt + DOC_maint + DOC_depr + DOC_lnr + DOC_fin
    LCC = (AEP + (DOC-DOC_depr) * R_blann * N_yr)/(1-F_disp)
    doclist.append(DOC)
plt.clf()
plt.scatter(mtomlist,doclist)
plt.ylim([13.386,13.390])
plt.show()
#plt.clf()
#labels = ['Insurance', 'Capital cost', 'Navigation cost', 'Landing fees','Flight crew','Maintenance','Fuel']
#sizes = [round(C_ins*R_bl,2), round(DOC_fin*R_bl,2),round(C_nf*R_bl,2),round(C_lf*R_bl,2),round(C_crew*R_bl,2),round(DOC_maint*R_bl,2), round(C_pol*R_bl,2)]
#patches, texts, pcts = plt.pie(sizes,labels=sizes, startangle=90,autopct='%1.1f%%',pctdistance=0.80, labeldistance=1.03)
#plt.legend(patches, labels, loc='center left')
#plt.axis('equal')
#plt.tight_layout()
#plt.title('Operational costs 1000 km Aquila, Total cost: '+str(round((DOC-DOC_depr)*R_bl,2))+' USD')
#plt.show()

#labels = ['Flight cost', 'Maintenance cost', 'Depreciation cost', 'Landing fees and navigation fees','Financing costs']
#sizes = [DOC_flt, DOC_maint, DOC_depr, DOC_lnr,DOC_fin]
#patches, texts, pcts = plt.pie(sizes, shadow=True, startangle=90,autopct='%1.1f%%')
#plt.legend(patches, labels, loc="best")
#plt.axis('equal')
#plt.tight_layout()
#plt.show()

#labels = ['Crew', 'Fuel and oil', 'Insurance cost']
#sizes = [C_crew, C_pol, C_ins]
#patches, texts, pcts = plt.pie(sizes, shadow=True, startangle=90,autopct='%1.1f%%')
#plt.legend(patches, labels, loc="best")
#plt.axis('equal')
#plt.tight_layout()
#plt.show()

#labels = ['Oil', 'Fuel']
#sizes = [C_olbl, C_fuel]
#patches, texts, pcts = plt.pie(sizes, shadow=True, startangle=90,autopct='%1.1f%%')
#plt.legend(patches, labels, loc="best")
#plt.axis('equal')
#plt.tight_layout()
#plt.show()