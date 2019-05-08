import math
tbp=False
jet =True
P_TO_tbp=4000000
T_TO_jet=50000
n_engines=2
rho_0=1.225
                    #TURBOPROP
if tbp==True:
    #turboshaft dimensions
    diameter_engine=0.2*(P_TO_tbp/(1000*n_engines))**0.18
    length_engine=0.1*(P_TO_tbp/(1000*n_engines))**0.4
    #propeller dimensions
    diameter_propeller=0.55*(P_TO_tbp/(1000*n_engines))**0.25
    #engine envelope dimensions
    height_engine_envelope=1.5*diameter_engine
    width_engine_envelope=1.1*diameter_engine
    length_engine_envelope=length_engine

if jet==True:
    bypass_ratio=8
    a_0=340.3 #[m/s]
    e_nozzle=0.97
    e_tf=0.75
    T_t4=1500 #1350-1650 [K]
    G=(T_t4/600)-1.25

    massflow=(T_TO_jet*(1+bypass_ratio))/(n_engines*a_0*math.sqrt(5*e_nozzle*G*(1+e_tf*bypass_ratio)))
    #Intake dimensions
    inlet_spinner_ratio=0.05*(1+(rho_0*a_0)/massflow+3*bypass_ratio/(1+bypass_ratio))
    diameter_inlet=1.65*math.sqrt((massflow/(rho_0*a_0)+0.0050)/(1-(inlet_spinner_ratio)**2))
    diameter_highlight=diameter_inlet
    #Fan cowling dimensions
    if jettypeB==True:
        lenth_nacelle=c_l_nacelle*(math.sqrt((massflow*(1+0.2*bypass_ratio))/(rho_0*a_0*(1+bypass_ratio)))+delta_l_nacelle)
