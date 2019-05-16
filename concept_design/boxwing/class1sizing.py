import math
import numpy as np
#FUSELAGE SIZING
def fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles):

    mass_perpassenger=195*0.45359237
    mass_peroverheadluggage=30*0.45359237
    masspl=n_passenger*(mass_perpassenger+mass_peroverheadluggage)

    #general fuselage shape = cylindrical

                #CROSS SECTION
    #cabin layout
    n_rows=n_passenger/n_seats_abreast
    width_seats=17*0.0254 #16-17
    pitch_seats=31*0.0254 #30-32
    head_room=1.7
    height_shoulder=1.3
    width_aisles=0.50 #min 38cm
    n_trolleys=6
    n_lavatories=1
    n_galleys=1
    width_armrest=0.05
    s_clearance=0.02
    #cabin width up to shoulder
    width_cabin=(n_seats_abreast*width_seats)+(n_seats_abreast)*width_armrest+n_aisles*width_aisles+2*s_clearance
    #cabin width at head of passengers
    width_headroom=width_cabin-2*(width_armrest+s_clearance)-width_seats
    #overhead storage
    k_os=1
    n_os=2
    area_overhead=0.1
    length_seats=n_rows*pitch_seats+(10*0.0254)
    volume_overhead=(n_os*area_overhead)*length_seats*k_os

    #cargo compartment
    density_luggage=170
    density_cargo=160
    mass_cargo=0
    mass_luggage=n_passenger*mass_peroverheadluggage
    volume_luggage=mass_luggage/density_luggage
    volume_cargo=mass_cargo/density_cargo
    volume_cargocompartment=volume_cargo+(volume_luggage-volume_overhead)
    #length_cargo=volume_cargocompartment/(width_cabin*height_shoulder)
    length_cargocompartment=2


    #structural dimensions
    thickness_floor=0.150 #0.100-0.300
    thickness_fuselage_skin_frame=0.150

    diameter_fuselage_inside=2*math.sqrt(((height_shoulder/2)**2)+((width_cabin/2)**2))

    diameter_fuselage_outside=1.045*diameter_fuselage_inside+0.084



    height_cargo=(diameter_fuselage_inside/2)-((height_shoulder/2)+thickness_floor)
    height_aisle=diameter_fuselage_inside-(thickness_floor+height_cargo)

                    #FUSELAGE TOP VIEW
    width_lavatory=36*0.0254
    length_lavatory=36*0.0254

    width_galley=36*0.0254
    length_galley=30*0.0254

    length_cabin=length_seats+length_lavatory+length_galley+length_cargocompartment

    nose_fineness=1 #can be altered using aerodynamic data
    nosecone_fineness=2 #from data
    length_nose=nose_fineness*diameter_fuselage_outside
    length_nosecone=nosecone_fineness*diameter_fuselage_outside


    tail_fineness=1.6
    tailcone_fineness=3 #2-4
    length_tail=tail_fineness*diameter_fuselage_outside
    length_tailcone=tailcone_fineness*diameter_fuselage_outside

                    #SIDE VIEW
    overnose_angle=11 #11-20 #dependent on approach situation
    overside_angle=35 #from table
    length_flightdeck=2.5

    length_fuselage=length_nose+length_cabin+length_tail

    return length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside, length_nosecone, length_tailcone

def enginedimensions(n_engines, T_TO_jet):


    bypass_ratio=12
    a_0=340.3 #[m/s]
    e_nozzle=0.97
    e_tf=0.75
    T_t4=1500 #1350-1650 [K]
    G=(T_t4/600)-1.25
    rho_0 = 1.225

    massflow=(T_TO_jet*(1+bypass_ratio))/(n_engines*a_0*math.sqrt(5*e_nozzle*G*(1+e_tf*bypass_ratio)))
                #Intake dimensions
    inlet_spinner_ratio=0.05*(1+(rho_0*a_0)/massflow+3*bypass_ratio/(1+bypass_ratio))
    diameter_inlet=1.65*math.sqrt((massflow/(rho_0*a_0)+0.0050)/(1-(inlet_spinner_ratio)**2))
    diameter_highlight=diameter_inlet
                #Fan cowling dimensions
    phi_fan_cowling=0.6 #0.5-0.75
    jettypeB = True
    jettypeC = False
    if jettypeB:
        phi_fan_cowling = 1.
        c_l_nacelle=9.8
        delta_l_nacelle=0.05
        beta_nacelle=0.35
    if jettypeC:
        c_l_nacelle=7.8
        delta_l_nacelle=0.10
        beta_nacelle=0.21+(0.12/math.sqrt(phi_fan_cowling-0.3))
    #nacelle length
    length_nacelle=c_l_nacelle*(math.sqrt((massflow*(1+0.2*bypass_ratio))/(rho_0*a_0*(1+bypass_ratio)))+delta_l_nacelle)
    #fan cowl length
    length_fan_cowl=phi_fan_cowling*length_nacelle
    #maximum nacelle diameter
    diameter_nacelle=diameter_inlet+0.06*phi_fan_cowling*length_nacelle+0.03
    #exit diameter fan
    diameter_exit_fan=diameter_nacelle*(1-(phi_fan_cowling**2)/3)
                #gas generator cowling dimensions
    #exposed length of gas generator
    length_generator=(1-phi_fan_cowling)*length_nacelle
    #gas generator cowling at fan exit diameter
    gasgenerator_coefficient=(massflow*bypass_ratio)/(rho_0*a_0)
    diameter_gas_generato_fan=diameter_exit_fan*((0.089+4.5)/(0.067+5.8))**2
    #gas generator cowling at gas generator exit diameter
    diameter_gas_generator=0.55*diameter_gas_generato_fan


    #print(length_nacelle,'nacelle length')
    #print(diameter_nacelle, 'nacelle diameter')
    #print(length_fan_cowl, 'fan cowling length')
    #print(diameter_exit_fan, 'fan exit diameter')
    #print(diameter_gas_generator,'gas generator diameter')
    return length_fan_cowl, diameter_nacelle, length_generator, length_nacelle, diameter_highlight, diameter_exit_fan, diameter_gas_generator, diameter_inlet

def empennage(V_v, l_v, S, b):
    S_v = (V_v * S * b) / l_v

    #Statistical Values
    AR_v = 1.5
    taper_v = 0.7
    sweepLE_v = 35

    span_v = np.sqrt(AR_v * S_v)
    av_chord_v = span_v / AR_v
    root_chord_v = av_chord_v / (0.5*(1+taper_v))
    tip_chord_v = root_chord_v * taper_v

def undercarriage(main_landing_pos, nose_landing_pos, length_fuselage, length_tail, diameter_fuselage_outside):
    dist_to_tail = length_fuselage - main_landing_pos - length_tail
    scrap_angle = np.radians(15)
    wheel_height = dist_to_tail * np.tan(scrap_angle)

    lateral_position = (main_landing_pos + nose_landing_pos) / np.sqrt(((nose_landing_pos**2 * np.tan(np.radians(55))**2)/(wheel_height+0.3*diameter_fuselage_outside))-1)
    return wheel_height, lateral_position

def tiresizing(MTOW, LCN):
    # LCN is in the range of 25 when looking at reference aircraft
    if 0 < LCN <= 100:
        tire_pressure = 430*np.log(LCN) - 630
    else:
        print("That is not possible")
    # This is jsut fixed: Our aircraft is not big enough for more wheels and not small enough for less wheels.
    N_mw = 4
    N_nw = 2

    # If MTOW is entered in Newton: This will bring it back to kg
    if MTOW > 100000:
        MTOW = MTOW / 9.81

    # Load on each wheel
    P_mw = (0.92 * MTOW) / N_mw
    P_nw = (0.08 * MTOW) / N_nw

    # LOOK IN THE SLIDES FOR THE WHEEL DIMENSIONS MATCHING THE TIRE LOAD
    # Probably the best dimensions will be:
    # Wheel = outer dimension x width - inner dimension
    # Main wheel = 0.84 x 0.25 - 0.41
    # Nose wheel = 0.46 x 0.11 - 0.25

    return tire_pressure, P_mw, P_nw
