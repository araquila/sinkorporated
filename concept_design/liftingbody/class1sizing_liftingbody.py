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
    height_shoulder=1.30
    width_aisles=0.38 #min 38cm
    n_trolleys=6
    n_lavatories=1
    n_galleys=1
    width_armrest=0.05
    s_clearance=0.02

    #overhead storage
    k_os=1
    n_os=2
    area_overhead=0.3       # verified with CATIA model
    length_seats=n_rows*pitch_seats+(10*0.0254)
    volume_overhead=(n_os*area_overhead)*length_seats*k_os

    # Data measured from CATIA
    diameter_fuselage_inside = 2.9646
    diameter_fuselage_outside = 3.182
    diff_outside_inside = diameter_fuselage_outside - diameter_fuselage_inside
    width_fuselage_inside = 3.931
    width_fuselage_outside = 4.382

    #structural dimensions
    thickness_floor=0.150 #0.100-0.300
    thickness_fuselage_skin_frame=0.150

                    #FUSELAGE TOP VIEW
    width_lavatory=36*0.0254
    length_lavatory=36*0.0254

    width_galley=36*0.0254
    length_galley=30*0.0254

    length_cabin=length_seats+length_lavatory+length_galley

    #cargo compartment
    density_luggage=170
    density_cargo=160
    mass_cargo=n_passenger*23
    mass_luggage=n_passenger*mass_peroverheadluggage
    volume_luggage=mass_luggage/density_luggage
    volume_cargo=mass_cargo/density_cargo
    if volume_luggage > volume_overhead:
        volume_cargocompartment=volume_cargo+(volume_luggage-volume_overhead)
    else:
        volume_cargocompartment = volume_cargo
    height_cargo=(diameter_fuselage_inside/2)-((height_shoulder/2)+thickness_floor)
    area_cargocompartment = volume_cargocompartment / (0.8*length_cabin)

    height_container=64*0.0254
    width_container=61.5*0.0254
    base_container=60.4*0.0254

    nose_fineness=1 #can be altered using aerodynamic data
    nosecone_fineness=2 #from data
    length_nose=nose_fineness*width_fuselage_outside
    length_nosecone=nosecone_fineness*width_fuselage_outside


    tail_fineness=1.6
    tailcone_fineness=3 #2-4
    length_tail=tail_fineness*diameter_fuselage_outside
    length_tailcone=tailcone_fineness*diameter_fuselage_outside

                    #SIDE VIEW
    overnose_angle=11 #11-20 #dependent on approach situation
    overside_angle=35 #from table
    length_flightdeck=2.5

    length_fuselage=length_nose+length_cabin+length_tail

    return length_nose, length_nosecone, length_cabin, length_tail, length_tailcone, length_fuselage, diameter_fuselage_inside, diameter_fuselage_outside, width_fuselage_outside

                    #WING SIZING
def det_quarter_chord_sweep(M_cruise, supercritical = False, delta_mach = 0.03):
    """
    determines the quarter chord sweep in radians
    """
    # Delta_mach can range from 0 to 0.05 but is given as 0.03 in ADSEE Slides
    if 0.7 < M_cruise < 1:
        if supercritical:
            sweep = 0.75 * 0.935 / (M_cruise + delta_mach) # 0.935 from statistical data from Torenbeek
            return np.arccos(sweep)
        else:
            sweep = 0.75 / (M_cruise + delta_mach)
            return np.arccos(sweep)
    if M_cruise <= 0:
        raise NameError("Plane flying backwards??")
    if M_cruise >= 1:
        raise NameError("Going supersonic now, are we?")
    else:
        return np.arccos(1)

def det_planform(S, AR, M_cruise, C_L_cruise, sweep, supercritical = False, delta_mach = 0.03):
    # Delta_mach can range from 0 to 0.05 but is given as 0.03 in ADSEE Slides
    b = np.sqrt(AR * S)
    taper = 0.2 * (2 - sweep)
    root_chord = (2 * S)/((1 + taper) * b)
    tip_chord = taper * root_chord
    half_chord_sweep = np.arctan(np.tan(sweep) - (4 / AR) * (0.25 * (1 - taper)/(1 + taper)))
    if 0 < M_cruise < 1:
        if supercritical:
            t_c_ratio = min(0.18, (np.cos(half_chord_sweep)**3 * (0.935 - (M_cruise + delta_mach) * np.cos(half_chord_sweep)) - 0.115 * C_L_cruise**1.5) / np.cos(half_chord_sweep)**2)
        else:
            t_c_ratio = min(0.18, (np.cos(half_chord_sweep)**3 * (1 - (M_cruise + delta_mach) * np.cos(half_chord_sweep)) - 0.115 * C_L_cruise**1.5) / np.cos(half_chord_sweep)**2)
    if M_cruise <= 0:
        raise NameError("Plane flying backwards??")
    if M_cruise >= 1:
        raise NameError("Going supersonic now, are we?")
    return b, taper, root_chord, tip_chord, t_c_ratio

def det_dihedral_angle(sweep, high = False, mid = False, low = False):
    dihedral = sweep * 18 / np.pi
    if high:
        angle = 1 - dihedral
        return angle
    if mid:
        angle = 3 - dihedral
        return angle
    if low:
        angle = 5 - dihedral
        return angle
    else:
        raise NameError("Where is the wing?")

def MAC(root_chord,t_c_ratio):
    MAC = root_chord*(2/3)*((1+t_c_ratio+t_c_ratio**2)/(1+t_c_ratio))
    return MAC



                    #ENGINE DIMENSIONS
def enginedimensions_tbp(rho0, n_engines, P_TO_tbp):
    #turboshaft dimensions
    diameter_engine=0.2*(P_TO_tbp/(1000*n_engines))**0.18
    length_engine=0.1*(P_TO_tbp/(1000*n_engines))**0.4
    #propeller dimensions
    diameter_propeller=0.55*(P_TO_tbp/(1000*n_engines))**0.25
    #engine envelope dimensions
    height_engine_envelope=1.5*diameter_engine
    width_engine_envelope=1.1*diameter_engine
    length_engine_envelope=length_engine
    return diameter_engine, length_engine, diameter_propeller



def enginedimensions_jet(rho0, n_engines, T_TO_jet, jettypeB=False, jettypeC=False):
    bypass_ratio=12
    a_0=340.3 #[m/s]
    e_nozzle=0.97
    e_tf=0.75
    T_t4=1500 #1350-1650 [K]
    G=(T_t4/600)-1.25

    massflow=(T_TO_jet*(1+bypass_ratio))/(n_engines*a_0*math.sqrt(5*e_nozzle*G*(1+e_tf*bypass_ratio)))
                #Intake dimensions
    inlet_spinner_ratio=0.05*(1+(rho0*a_0)/massflow+3*bypass_ratio/(1+bypass_ratio))
    diameter_inlet=1.65*math.sqrt((massflow/(rho0*a_0)+0.0050)/(1-(inlet_spinner_ratio)**2))
    diameter_highlight=diameter_inlet
                #Fan cowling dimensions
    phi_fan_cowling=0.6 #0.5-0.75
    if jettypeB:
        c_l_nacelle=9.8
        delta_l_nacelle=0.05
        beta_nacelle=0.35
    if jettypeC:
        c_l_nacelle=7.8
        delta_l_nacelle=0.10
        beta_nacelle=0.21+(0.12/math.sqrt(phi_fan_cowling-0.3))
    #nacelle length
    length_nacelle=c_l_nacelle*(math.sqrt((massflow*(1+0.2*bypass_ratio))/(rho0*a_0*(1+bypass_ratio)))+delta_l_nacelle)
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
    gasgenerator_coefficient=(massflow*bypass_ratio)/(rho0*a_0)
    diameter_gas_generato_fan=diameter_exit_fan*((0.089+4.5)/(0.067+5.8))**2
    #gas generator cowling at gas generator exit diameter
    diameter_gas_generator=0.55*diameter_gas_generato_fan

    return length_nacelle, length_fan_cowl, diameter_highlight, diameter_exit_fan, diameter_gas_generator, diameter_nacelle

def empennage(V_h, V_v, l_h, l_v, S, b, c):
    S_h = (V_h * S * c) / l_h
    S_v = (V_v * S * b) / l_v

    #Statistical Values
    AR_v = 1.5
    taper_v = 0.45
    sweepLE_v = 40

    AR_h = 4.5
    taper_h = 0.35
    sweepqc_h = 28

    span_v = np.sqrt(AR_v * S_v)
    av_chord_v = span_v / AR_v
    root_chord_v = av_chord_v / (0.5*(1+taper_v))
    tip_chord_v = root_chord_v * taper_v

    span_h = np.sqrt(AR_h * S_h)
    av_chord_h = span_h / AR_h
    root_chord_h = av_chord_h / (0.5*(1+taper_h))
    tip_chord_h = root_chord_h * taper_h
    sweepLE_h = sweepqc_h + np.degrees(np.arctan((0.25*root_chord_h-0.25*tip_chord_h)/(span_h/2)))

    return AR_h, AR_v, S_h, span_h, root_chord_h, tip_chord_h, sweepqc_h, sweepLE_h, S_v, span_v, root_chord_v, tip_chord_v, sweepLE_v

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
