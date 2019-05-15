# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from class1_boxwing import Weights_Class_I
from wingsizing import iterempty
from sustainability_functions import CO2_calc
from constant_variables import *

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
    cj_loiter_jet = 19e-6       # (0.4-0.6) [lbs/lbs/hr]
    cj_cruise_jet = 19e-6      # (0.5-0.9) [lbs/lbs/hr]

    V_cruise_jet = 229

    range_cruise_jet = 1850000

    endurance_loiter_jet = 2700

        # CLASS I ESTIMATION

    MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet, W_used_fuel_jet =  Weights_Class_I(W_empty_jet, W_payload, W_crew, C_fe, S, S_wet, A_jet, e_jet, cj_loiter_jet, cj_cruise_jet, f_trapped_fuel, V_cruise_jet, range_cruise_jet, endurance_loiter_jet, jet = True)



    W_empty_jet = (OEW_jet-W_crew)-f_trapped_fuel*MTOW_jet
    MTOW_jet_1000, OEW_jet_1000, W_fuel_jet_1000, LD_cruise_jet, W_used_fuel_jet  = Weights_Class_I(W_empty_jet, W_payload, W_crew, C_fe, S, S_wet, A_jet, e_jet, cj_loiter_jet, cj_cruise_jet, f_trapped_fuel, V_cruise_jet, 1000*1000, 0, jet = True, tbp = False)

    fuel_per_passenger_jet_1000 = (W_fuel_jet_1000/n_passenger)/g
    CO2_jet = fuel_per_passenger_jet_1000*3.0
    print('MTOM jet: ' + str(round(MTOW_jet/g,2)))
    print('Fuel per passenger per 1000 km propfan: ' + str(round(fuel_per_passenger_jet_1000,2)))
    print('CO2 per passanger per 1000 km propfan: ' + str(round(CO2_jet,2)))

    return MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet, W_used_fuel_jet


M_empty_jet = 13874.75
g = 9.80665
for i in range(50):
    MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet, W_used_fuel_jet = class1box(M_empty_jet)
    M_empty_jet, geometrylistfuselage, geometrylistwings, geometrylistvtail, mainlg_cg, wing_weight1_box, wing_weight2_box, x_cg, noselg_cg, S= iterempty(MTOW_jet, OEW_jet, W_fuel_jet, LD_cruise_jet)
    print(M_empty_jet, mainlg_cg, MTOW_jet/9.81)
print(geometrylistfuselage)
print(wing_weight1_box,wing_weight2_box, x_cg)

Fn = (mainlg_cg-x_cg)/(mainlg_cg-noselg_cg)
print(W_fuel_jet/g)
print(W_used_fuel_jet/g/60, S)
