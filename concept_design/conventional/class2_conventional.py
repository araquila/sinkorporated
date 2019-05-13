def C_L_des(q,W_S_cruise_start,W_S_cruise_end):
    C_L_des = 1.1*1/q*(0.5*(W_S_cruise_start+W_S_cruise_end))
    return C_L_des

def C_l_des(C_L_des,sweep):
    C_l_des = C_L_des/cos(sweep)**2
    return C_ldes
