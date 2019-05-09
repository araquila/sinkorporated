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
def wing(Mach_cruise, S, A, C_L, high=True, mid=False, low=False ):
    #AERODYNAMICS
    Mach_dd=Mach_cruise+0.03
    Mach_t=0.935
    #quarter chord sweep
    if Mach_cruise <0.7:
        sweep_chord_0_25=0
    else :
        sweep_chord_0_25=math.acos(0.75*(Mach_t/Mach_dd)) #[rad]
    #geometric parameters
    #taper ratio
    taper=0.2*(2-sweep_chord_0_25)
    #wingspan
    b=math.sqrt(S*A)
    #Chord lengths
    rootchord=(2*S)/((1+taper)*b)
    tipchord=taper*rootchord
    #half wing sweep
    sweep_chord_0_5=10*(math.pi/180)
    #thickness-to-chord ratio
    cos=math.cos(sweep_chord_0_5)
    thickness_chord_ratio=((cos**3)*(Mach_t-Mach_dd*cos)-(0.115*C_L**1.5))/(cos**2)
    if thickness_chord_ratio > 0.18:
        thickness_chord_ratio=0.18
    #dihedral angle
    dihedral=3-(sweep_chord_0_5*(180/math.pi)/10)
    if high:
        dihedral=dihedral+2 #0-1[deg]
    if low:
        dihedral=dihedral+2
    return taper, b, rootchord, tipchord, sweep_chord_0_5, sweep_chord_0_25, thickness_chord_ratio, dihedral

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
