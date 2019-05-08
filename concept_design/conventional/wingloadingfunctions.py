#a list of the functions that are used to calculate the wing loading
from math import *
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
