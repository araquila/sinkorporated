import numpy as np
import parameters as p
import matplotlib.pyplot as plt

#### FUSELAGE GROUP ###
cgx_fuselage = 0.4 * p.l_fuselage
cgx_tail = 0.9 * p.l_fuselage
cgx_nl = 2 / p.l_fuselage

m_fusgroup = p.m_fuselage + p.m_tail + p.W_nose_landing / p.g
cgx_fusgroup = (cgx_fuselage * p.m_fuselage + cgx_tail * p.m_tail + cgx_nl * (p.W_nose_landing / p.g)) / m_fusgroup


### WING GROUP ###
cgx_wing = 0.35
cgx_fueltank = 0.5
cgx_engine = -0.1
cgx_ml = (10.94 - 0.4283 * p.l_fuselage) / p.MAC

m_winggroup = (p.m_wing + p.m_engine + p.m_fuel_tanks + p.W_main_landing / p.g)
cgx_winggroup = (cgx_wing * p.m_wing + cgx_fueltank * p.m_fuel_tanks + cgx_ml * (p.W_main_landing / p.g)) / m_winggroup


n_rows = 15
begin_cabin = p.l_nose + p.l_lavatory
end_cabin = begin_cabin + n_rows * p.seat_pitch

x_cargo = end_cabin + p.l_galley + 0.5
x_seats = np.arange(begin_cabin, end_cabin, p.seat_pitch)
### OEW CG ###

x_lemac = np.arange(8,15,2)
for i in range(len(x_lemac)):
    cgx_oew = (m_fusgroup * cgx_fusgroup + (x_lemac[i] + p.MAC * cgx_winggroup) * m_winggroup) / (m_fusgroup + m_winggroup)
    mass_oew = m_fusgroup + m_winggroup

    cgx.append((x_cargo * p.M_total_cargo + cgx[-1] * mass[-1]) / (mass[-1] + p.M_total_cargo))
    mass.append(p.M_total_cargo + mass[-1])
    for i in range(len(x_seats)):
        mass.append(mass[-1] + 4 * p.M_passenger)
        cgx.append((4 * p.M_passenger * x_seats[-(i + 1)] + cgx[-1] * mass[-2]) / mass[-1])
    
    plt.scatter(cgx, mass)
    plt.show()
### cargo ###
#cg_cargo = 



### loading pax aft ###



#

