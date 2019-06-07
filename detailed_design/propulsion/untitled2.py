#radius = 0.4 #[m]
#sigma_yield = 700e6 #[Pa]
#p = 12e6 #[Pa]
#safety_factor = 1.5 #[-]
#p_a = 101325 #[Pa]
#k = 12
#n = 4
#E = 72e9 #[Pa] E-modulus of aluminium
#radius_outer = 0.5

def t_hoop(radius,sigma_yield,p):
    t_hoop = (safety_factor*p*radius)/(2*sigma_yield)
    return t_hoop

def t_long(radius,sigma_yield,p):
    t_long = (safety_factor*p*radius)/sigma_yield
    return t_long

def t_dewar(radius_outer,p_a,n,E,k):
    t_dewar = 2*radius_outer*((p_a*n)/(E*k))**(1./3)
    return t_dewar