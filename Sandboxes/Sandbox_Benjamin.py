import math
import matplotlib.pyplot as plt
import numpy as np

tbp=False
jet =True
jettypeC=True
jettypeB=False
P_TO_tbp=4000000
T_TO_jet=50000
n_engines=2
rho_0=1.225
b=15
                    #TURBOPROP
if tbp:
    #turboshaft dimensions
    diameter_engine=0.2*(P_TO_tbp/(1000*n_engines))**0.18
    length_engine=0.1*(P_TO_tbp/(1000*n_engines))**0.4
    #propeller dimensions
    diameter_propeller=0.55*(P_TO_tbp/(1000*n_engines))**0.25
    #engine envelope dimensions
    height_engine_envelope=1.5*diameter_engine
    width_engine_envelope=1.1*diameter_engine
    length_engine_envelope=length_engine

if jet:
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
    phi_fan_cowling=0.6 #0.5-0.75
    if jettypeB:
        c_l_nacelle=9.8
        delta_l_nacelle=0.05
        beta_nacelle=0.35
    if jettypeC:
        c_l_nacelle=7.8
        delta_l_nacelle=0.10
        beta_nacelle=0.21+(0.12/math.sqrt(phi_fan_cowling-0.3))
    #nacelle length
    length_nacelle=c_l_nacelle*(math.sqrt((massflow*(1+0.2*bypass_ratio))/(rho_0*a_0*(1+bypass_ratio)))+delta_l_nacelle)
    #fan cowl length
    length_fan_cowl=phi_fan_cowling*length_nacelle
    #maximum nacelle diameter
    diameter_nacelle=diameter_inlet+0.06*phi_fan_cowling*length_nacelle+0.03
    #exit diameter fan
    diameter_exit_fan=diameter_nacelle*(1-(phi_fan_cowling**2)/3)
                #gas generator cowling dimensions
    #exposed length of gas generator
    length_generator=(1-phi_fan_cowling)*length_nacelle
    #gas generator cowling at fan exit diameter
    gasgenerator_coefficient=(massflow*bypass_ratio)/(rho_0*a_0)
    diameter_gas_generato_fan=diameter_exit_fan*((0.089+4.5)/(0.067+5.8))**2
    #gas generator cowling at gas generator exit diameter
    diameter_gas_generator=0.55*diameter_gas_generato_fan
                #integration
    #lateral disposition
    y_engine=0.35*(b/2)
    #longitudinal position
    def f(x_f_c):
        return 0.07+0.03*np.cos(15*(x_f_c+0.03))
    x_f_c=np.arange(-0.2,0.18,0.02)
    plt.plot(x_f_c,f(x_f_c))
    plt.show()
