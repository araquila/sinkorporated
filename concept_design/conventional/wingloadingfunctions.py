#a list of the functions that are used to calculate the wing loading
from math import *
import numpy as np
def V_stall_calc(W,rho,CL,S):
    V_stall = sqrt(2*W/(rho*CL*S))
    return V_stall

def CL_TO_calc(CL_max_TO):
    CL_TO = CL_max_TO/(1.1*1.1)
    return CL_TO

def V_TO_calc(V_stall):
    V_TO = 1.1*V_stall
    return V_stall

def TOP_jet_calc(W_S,T_W,CL_TO):
    TOP_jet = W_S/T_W*1/CL_TO
    return TOP_jet

def TOP_prop_calc(W_S,BHP_W,CL_TO):
    TOP_prop = W_S/BHP_W*1/CL_TO

def W_S_calc(rho,V_stall,CL_max):
    W_S = 0.5*rho*V_stall**2*CL_max
    return W_S

def W_P_calc(x,k,CL_TO):
    y = k/x*CL_TO
    return y

def T_W_calc(x,k,CL_TO):
    y = x/k/CL_TO
    return y

def W_S_landing_calc(CL_max,rho,V_landing,f):
    W_S_landing = CL_max*rho*V_landing**2/(2*f)
    return W_S_landing

def W_P_cruise_tbp_calc(power_setting,weight_fraction,eff_prop,rho,rho0,C_D_0_tbp,x,A_tbp,V_cruise_tbp,e_tbp):
    W_P_cruise_tbp = (power_setting/weight_fraction*eff_prop*(rho/rho0)**0.75)*((C_D_0_tbp*0.5*rho*V_cruise_tbp**3)/(weight_fraction*x)+(weight_fraction*x)/(np.pi*e_tbp*A_tbp*0.5*rho*V_cruise_tbp))**-1
    return W_P_cruise_tbp

def T_W_cruise_jet_calc(thrust_setting,weight_fraction,rho,rho0,C_D_0_jet,x,A_jet,V_cruise_jet,e_jet):
    T_W_cruise_jet = weight_fraction/thrust_setting*(rho0/rho)**0.75*((C_D_0_jet*0.5*rho*V_cruise_jet**2)/(weight_fraction*x)+(weight_fraction*x)/(pi*A_jet*e_jet*0.5*rho*V_cruise_jet**2))
    return T_W_cruise_jet

def W_P_climb_calc(eff_prop,c,x,rho,A_tbp,e_tbp,C_D_0_tbp):
    W_P_climb = eff_prop/(c+sqrt(x)*sqrt(2/rho))/((1.345*((A_tbp*e_tbp)**0.75))/(C_D_0_tbp**0.25))
    return W_P_climb
