import math
#FUSELAGE SIZING
def fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles):

    mass_perpassenger=195*0.45359237
    mass_peroverheadluggage=30*0.45359237
    masspl=n_passenger*(mass_perpassenger+mass_peroverheadluggage)

    #general fuselage shape = cylindrical

                #CROSS SECTION
    #cabin layout
    n_rows=n_passenger/n_seats_abreast
    width_seats=16*0.0254 #16-17
    pitch_seats=31*0.0254 #30-32
    head_room=1.5
    height_shoulder=1.20
    width_aisles=0.38 #min 38cm
    n_trolleys=6
    n_lavatories=1
    n_galleys=1
    width_armrest=0.05
    s_clearance=0.02
    #cabin width up to shoulder
    width_cabin=(n_seats_abreast*width_seats)+(n_seats_abreast+n_aisles+1)*width_armrest+n_aisles*width_aisles+2*s_clearance
    #cabin width at head of passengers
    width_headroom=width_cabin-2*(width_armrest+s_clearance)-width_seats
    #overhead storage
    k_os=1
    n_os=2
    area_overhead=0.2
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

    height_container=64*0.0254
    width_container=61.5*0.0254
    base_container=60.4*0.0254
    #structural dimensions
    thickness_floor=0.150 #0.100-0.300
    thickness_fuselage_skin_frame=0.150

    diameter_fuselage_inside=2*math.sqrt(((height_shoulder/2)**2)+((width_cabin/2)**2))

    diameter_fuselage_outside=1.045*diameter_fuselage_inside+0.084

    height_cargo=(diameter_fuselage_inside/2)-((height_shoulder/2)+thickness_floor)

                    #FUSELAGE TOP VIEW
    width_lavatory=36*0.0254
    length_lavatory=36*0.0254

    width_galley=36*0.0254
    length_galley=30*0.0254

    length_cabin=length_seats+length_lavatory+length_galley

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

    return length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside

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
    if 0.7 < M_cruise < 1:
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

                    #ENGINE DIMENSIONS

def enginedimensions(n_engines, P_TO_tbp, tbp=True):
    pass
                    #TURBOPROP
    if tbp:
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

print(fuselage(60,4,4,1))
print(wing(0.72, 55, 19.5, 0.75))
print(enginedimensions(2,5*10**6,))
