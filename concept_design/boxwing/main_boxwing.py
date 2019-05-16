# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from class1_boxwing import Weights_Class_I
from wingsizing import iterempty
from sustainability_functions import CO2_calc
from constant_variables import *
from class1sizing import *
from math import *
from class2_boxwing import *
from atmosphere import atmosphere_calc
from cost_equations import *
import numpy as np

chosen_fuel_energy_density = energy_density_kerosene
fuel_efficiency_factor = energy_density_kerosene/chosen_fuel_energy_density

def class1box(M_empty_jet):

    # Gravitional constant
    g = 9.8065

    # Passengers and crew
    n_passenger = 60
    M_passenger = 105           # (Including luggage)
    n_crew = 4
    M_crew_member = 100

    #M_OEM = 17133.6
    #M_MTOM = 25584.5


    # Initial mass and fractions
    M_payload = n_passenger * M_passenger
    M_crew = n_crew * M_crew_member
    f_trapped_fuel = 0.003      # Range 0.001-0.005
    #M_empty_tbp = M_OEM-M_crew-f_trapped_fuel*M_MTOM


    #M_empty_jet = M_OEM-M_crew-f_trapped_fuel*M_MTOM


    # Convert to weights
    W_payload = M_payload * g
    W_crew = M_crew * g

    W_empty_jet = M_empty_jet * g

    ## Initial jet and tbp aircraft parameters
    C_fe = 0.003
    S = 1
    S_wet = 5 * S

    # Jet
    A_jet = 12
    e_jet = 1.2
    cj_loiter_jet = 19e-6/1.16       # (0.4-0.6) [lbs/lbs/hr]
    cj_cruise_jet = 19e-6/1.16    # (0.5-0.9) [lbs/lbs/hr]

    V_cruise_jet = 229

    range_cruise_jet = 1850000

    endurance_loiter_jet = 2700

        # CLASS I ESTIMATION

    MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet, W_used_fuel_jet, f_cruise_start_jet, f_cruise_end_jet =  Weights_Class_I(W_empty_jet, W_payload, W_crew, C_fe, S, S_wet, A_jet, e_jet, cj_loiter_jet, cj_cruise_jet, f_trapped_fuel, V_cruise_jet, range_cruise_jet, endurance_loiter_jet, jet = True)



    W_empty_jet = (OEW_jet-W_crew)-f_trapped_fuel*MTOW_jet
    MTOW_jet_1000, OEW_jet_1000, W_fuel_jet_1000, LD_cruise_jet, W_used_fuel_jet, f_cruise_start_jet, f_cruise_end_jet = Weights_Class_I(W_empty_jet, W_payload, W_crew, C_fe, S, S_wet, A_jet, e_jet, cj_loiter_jet, cj_cruise_jet, f_trapped_fuel, V_cruise_jet, 1000*1000, 0, jet = True, tbp = False)

    fuel_per_passenger_jet_1000 = (W_fuel_jet_1000/n_passenger)/g
    CO2_jet = fuel_per_passenger_jet_1000*3.0
    NOX_jet = 3e-3*fuel_per_passenger_jet_1000
    #print('MTOM jet: ' + str(round(MTOW_jet/g,2)))
    fuelinit = 30.52339
    co2init = 86.54437377
    noxinit = 0.0865443737
    print()
    #print('Fuel per passenger per 1000 km turbofan: ',fuel_per_passenger_jet_1000)
    #print((fuel_per_passenger_jet_1000-fuelinit)/fuelinit)
    print('NOX per passanger per 1000 km turbofan: ',NOX_jet)
    print((NOX_jet-noxinit)/noxinit)
    print('CO2 per passanger per 1000 km turbofan: ',CO2_jet)
    print((CO2_jet-co2init)/co2init)

    return MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet, W_used_fuel_jet, f_cruise_start_jet, f_cruise_end_jet


M_empty_jet = 13874.75
g = 9.80665


for i in range(50):
    MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet, W_used_fuel_jet, f_cruise_start_jet, f_cruise_end_jet = class1box(M_empty_jet)
    M_empty_jet, geometrylistfuselage, geometrylistwings, geometrylistvtail, mainlg_cg, wing_weight1_box, wing_weight2_box, x_cg, noselg_cg, S = iterempty(MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet)

    #print(M_empty_jet, mainlg_cg, MTOW_jet/9.81)
#print(geometrylistfuselage)
#print(wing_weight1_box,wing_weight2_box, x_cg)

Fn = (mainlg_cg-x_cg)/(mainlg_cg-noselg_cg)
#print(W_fuel_jet/g)
#print(W_used_fuel_jet/g/60, S)


def C_L_des(q,W_S_cruise_start,W_S_cruise_end):
    C_L_des = 1.1*1/q*(0.5*(W_S_cruise_start+W_S_cruise_end))
    return C_L_des

def C_l_des(C_L_des,sweep):
    C_l_des = C_L_des/cos(sweep)**2
    return C_l_des
altitude = 8000
gamma = 1.4
temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, temperature0, temperature_gradient, g, R, gamma)


S_jet = S
V_cruise_jet = 229
g =  9.80665
altitude = 10000

TW = 0.25
#input parameters for wing area
WS = 3080
S = MTOW_jet/WS
S = S *1.2

s1frac = 0.5
s2frac = 1 - s1frac
length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside, length_nosecone, length_tailcone = fuselage(60,4,4,1)
#print(fuselage(60,4,4,1))
S1 = S * s1frac
S2 = S * s2frac
AR1 = 12
AR2 = 12
Fus_len = length_fuselage
frac_qtrchord_fus = 0.35

#undercarriage initial Values
wheel_height, lateral_position_undercarridge = undercarriage(0.58*Fus_len, 0.085*Fus_len, Fus_len, length_tail, diameter_fuselage_outside)
tire_pressure, P_mw, P_nw = tiresizing(MTOW_jet/g, 25)
#print(wheel_height, lateral_position_undercarridge)



taper1 = 0.3
sweep1 = (2-taper1/0.2)*180/pi+3
b = sqrt(S1*AR1)
cr1 = 2*S1/(1+taper1)/b
ct1 = taper1 * cr1
#print(b,cr1,ct1)

ct2 = ct1
cr2 = (2*S2-ct2*b)/b
taper2 = ct2/cr2
AR_v = 0.8
sweepv = 30
taperv = 0.7
ctv = cr2
crv = ctv/taperv
bv = AR_v*(ctv+crv)/2
Sv = bv**2 / AR_v
vtail_xoffset = -crv + sin(sweepv*pi/180)*bv + ctv
sweep2 = -atan((Fus_len + vtail_xoffset - Fus_len*frac_qtrchord_fus - b/2*tan(sweep1/180*pi) - 0.75*cr2 )/(b/2))*180/pi

kgtolbs = 2.20462262
mtoft = 3.2808399
MTOWlbs = MTOW_jet/g*kgtolbs
n_lim = ult_load_factor(MTOWlbs)
n_ult = n_lim*1.5
V_cruise = 229
V_diving = 1.4*V_cruise
alt_cruise = 35000/mtoft
a_cruise = 296
M_cruise = V_cruise/a_cruise
M_cross = 0.935
M_dd = M_cruise + 0.03
qhat = 0.5*1.4*23842*M_cruise**2
Clhat1 = MTOW_jet/(qhat*S1)
Clhat2 = MTOW_jet/(qhat*S2)
Clhatv = 0

sweephalfchord1 = atan((sin(sweep1/180*pi)*b/2+ct1/4)/(b/2))
sweephalfchord2 = atan((sin(sweep2/180*pi)*b/2+ct2/4)/(b/2))
sweephalfchordv = atan((sin(sweepv/180*pi)*b/2+ctv/4)/(bv/2))
t_c1 = min((0.18, ((cos(sweephalfchord1)**3*(M_cross-M_dd*cos(sweephalfchord1))-0.115*Clhat1**1.5)/(cos(sweephalfchord1)**2))))
t_c2 = min((0.18, ((cos(sweephalfchord2)**3*(M_cross-M_dd*cos(sweephalfchord2))-0.115*Clhat2**1.5)/(cos(sweephalfchord2)**2))))
t_cv = min((0.18, ((cos(sweephalfchordv)**3*(M_cross-M_dd*cos(sweephalfchordv))-0.115*Clhatv**1.5)/(cos(sweephalfchordv)**2))))

#desing cruise speed, diving speed, t_c
S_control_ail = 0.05*S
S_control_elev = 0.03*S



C_L_des1 = C_L_des(rho*(V_cruise_jet**2)/2,f_cruise_start_jet*MTOW_jet/2/S1,f_cruise_end_jet*MTOW_jet/2/S1)
C_L_des2 = C_L_des(rho*(V_cruise_jet**2)/2,f_cruise_start_jet*MTOW_jet/2/S2,f_cruise_end_jet*MTOW_jet/2/S2)

C_l_des1 = C_l_des(C_L_des1,sweep1*pi/180)
C_l_des2 = C_l_des(C_L_des2,sweep2*pi/180)
#print(S1,'S1')
#print(S2,'S2')
#print(b,'b')
#print(sweep1,'sweep1')
#print(sweep2,'sweep2')
#print(taper1,'taper1')
#print(taper2,'taper2')
#print(cr1,'cr1')
#print(cr2,'cr2')
#print(t_c1,'tc1')
#print(t_c2,'tc2')
#print(AR1,'AR1')
#print(AR2,'AR2')
#print()
#print(bv,'bv')
#print(sweepv,'sweepv')
#print(taperv,'taperv')
#print(crv,'crv')
#print(t_cv,'t_cv')
#print(AR_v,'arv')
mtow_init = 199349.285
wfuel_init = 27338.952

print('mtow is',MTOW_jet)
print((MTOW_jet-mtow_init)/mtow_init)
#print('fuel weight is', W_fuel_jet/g)
#print((W_fuel_jet-wfuel_init)/wfuel_init)
print()
print()
