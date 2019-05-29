radius = 1.42 #m
pressure_inside = 100000 #N/m2
pressure_outside = 35000 #N/m2
j = 2 #safety factor
sigma_yield = 138e6
density = 2700 #[kg/m3]

def t_circ(sigma_yield,p_in,p_out,R):
    t = j*(p_in-p_out)*R/sigma_yield
    return t

def t_long(sigma_yield,p_in,p_out,R):
    t = j*(p_in-p_out)*R/(2*sigma_yield)
    return t

t1 = t_circ(sigma_yield,pressure_inside,pressure_outside,radius)
t2 = t_long(sigma_yield,pressure_inside,pressure_outside,radius)

print(1000*t1,1000*t2)
weight_per_m = radius*t2*density
print(weight_per_m)
