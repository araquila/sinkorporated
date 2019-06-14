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

cg_margin = 0.02

m_fuel = p.W_fuel / p.g

x_lemac = np.arange(9,9.5,0.1)

cgx_min = []
cgx_max = []
xlemac = []
for i in range(len(x_lemac)):
    ### OEW CG ###
    cgx_oew = (m_fusgroup * cgx_fusgroup + (x_lemac[i] + p.MAC * cgx_winggroup) * m_winggroup) / (m_fusgroup + m_winggroup)
    mass_oew = m_fusgroup + m_winggroup
#    plt.plot(cgx_oew, mass_oew)

    ### CARGO LOADING ###
    cgx_cargo = [cgx_oew]
    mass_cargo = [mass_oew]
    cgx_cargo.append((x_cargo * p.M_total_cargo + cgx_cargo[-1] * mass_cargo[-1]) / (mass_cargo[-1] + p.M_total_cargo))
    mass_cargo.append(p.M_total_cargo + mass_cargo[-1])
#    plt.plot(cgx_cargo, mass_cargo)
    
    ### WINDOW SEATS ###
    cgx_windowback = [cgx_cargo[-1]]
    cgx_windowfront = [cgx_cargo[-1]]
    
    mass_window = [mass_cargo[-1]]
    for wseat in range(len(x_seats)):
        mass_window.append(mass_window[-1] + 2 * p.M_passenger)
        cgx_windowback.append((2 * p.M_passenger * x_seats[-(wseat + 1)] + cgx_windowback[-1] * mass_window[-2]) / mass_window[-1])
        cgx_windowfront.append((2 * p.M_passenger * x_seats[wseat] + cgx_windowfront[-1] * mass_window[-2]) / mass_window[-1])
#    plt.scatter(cgx_windowback, mass_window)
#    plt.scatter(cgx_windowfront, mass_window)
    
    ### AISLE SEATS ###
    cgx_aisleback = [cgx_windowback[-1]]
    cgx_aislefront = [cgx_windowfront[-1]]

    mass_aisle = [mass_window[-1]]
    for aseat in range(len(x_seats)):
        mass_aisle.append(mass_aisle[-1] + 2 * p.M_passenger)
        cgx_aisleback.append((2 * p.M_passenger * x_seats[-(aseat + 1)] + cgx_aisleback[-1] * mass_aisle[-2]) / mass_aisle[-1])
        cgx_aislefront.append((2 * p.M_passenger * x_seats[aseat] + cgx_aislefront[-1] * mass_aisle[-2]) / mass_aisle[-1])    
#    plt.scatter(cgx_aisleback, mass_aisle)
#    plt.scatter(cgx_aislefront, mass_aisle)
    
    ### FUEL ###
    x_fuel = x_lemac[i]
    cgx_fuel = [cgx_aisleback[-1]]
    
    mass_fuel = [mass_aisle[-1]]
    for pod in range(p.n_fueltanks):
        mass_fuel.append(mass_fuel[-1] + m_fuel / 2)
        cgx_fuel.append((x_fuel * m_fuel / 2 + mass_fuel[-2] * cgx_fuel[-1]) / mass_fuel[-1])

    cgx_min.append((1 - cg_margin) * min(cgx_aislefront))
    cgx_max.append((1 + cg_margin) * max(cgx_windowback))
    xlemac.append(x_lemac[i])

    if x_lemac[i] < 9.3 and x_lemac[i] > 9.2:
        plt.plot(cgx_oew, mass_oew)
        plt.plot(cgx_cargo, mass_cargo)
        plt.scatter(cgx_windowback, mass_window)
        plt.scatter(cgx_windowfront, mass_window)
        plt.scatter(cgx_aisleback, mass_aisle)
        plt.scatter(cgx_aislefront, mass_aisle)
        plt.scatter(cgx_fuel, mass_fuel)
        
CGmacmin = (cgx_min-x_lemac) / p.MAC
CGmacmax = (cgx_max-x_lemac) / p.MAC
xlemac = (x_lemac / p.l_fuselage)
plt.ylim([10000, 20000])
plt.show()