#a list of the functions that are used to calculate the wing loading
from math import *
def stall_calc(W,rho,CL,S):
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


print TOP_jet_calc(1000,0.2,1.2)
