import math
def fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles):
    #can be updated by including the mach number in order to alter the nose length
    #n_passenger=60
    #List of crew+operational items
    #n_crew=4
    mass_perpassenger=195*0.45359237
    mass_peroverheadluggage=30*0.45359237
    masspl=n_passenger*(mass_perpassenger+mass_peroverheadluggage)

    #general fuselage shape = cylindrical

    #Cross section
    #cabin layout
    #n_seats_abreast=4
    #n_aisles=1
    n_rows=n_passenger/n_seats_abreast
    width_seats=16*0.0254 #16-17
    pitch_seats=31*0.0254 #30-32
    head_room=1.5
    height_shoulder=1.20
    width_aisles=0.38 #min 38cm
    #width_aisles=0.51
    height_floor=1.9
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

    #fuselage top view
    width_lavatory=36*0.0254
    length_lavatory=36*0.0254

    width_galley=36*0.0254
    length_galley=30*0.0254

    length_cabin=length_seats+length_lavatory+length_galley

    nose_fineness=1
    nosecone_fineness=2 #from data
    length_nose=nose_fineness*diameter_fuselage_outside
    length_nosecone=nosecone_fineness*diameter_fuselage_outside


    tail_fineness=1.6
    tailcone_fineness=3 #2-4
    length_tail=tail_fineness*diameter_fuselage_outside
    length_tailcone=tailcone_fineness*diameter_fuselage_outside

    #side view
    overnose_angle=11 #11-20 #dependent on approach situation
    overside_angle=35 #from table
    length_flightdeck=2.5

    length_fuselage=length_nose+length_cabin+length_tail
    return length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside
print(fuselage(60, 4, 4, 1))
