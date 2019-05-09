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
    area_overhead=0.2
    length_seats=n_rows*pitch_seats+(10*0.0254)
    volume_overhead=(n_os*area_overhead)*length_seats*k_os

    #cargo compartment
    density_luggage=170
    density_cargo=160
    mass_cargo=n_passenger*23
    mass_luggage=n_passenger*mass_peroverheadluggage
    volume_luggage=mass_luggage/density_luggage
    volume_cargo=mass_cargo/density_cargo
    volume_cargocompartment=volume_cargo+(volume_luggage-volume_overhead)
    area_cargocompartment = volume_cargocompartment / (0.8*length_cabin)

    height_container=64*0.0254
    width_container=61.5*0.0254
    base_container=60.4*0.0254
    #structural dimensions
    thickness_floor=0.150 #0.100-0.300
    thickness_fuselage_skin_frame=0.150

    # Data measured from CATIA
    diameter_fuselage_inside = 2.9646
    diameter_fuselage_outside = 1.045*diameter_fuselage_inside + 0.084
    diff_outside_inside = diameter_fuselage_outside - diameter_fuselage_inside
    width_fuselage_inside = 3.931
    width_fuselage_outside = width_fuselage_inside + diff_outside_inside

    height_cargo=(diameter_fuselage_inside/2)-((height_shoulder/2)+thickness_floor)

                    #FUSELAGE TOP VIEW
    width_lavatory=36*0.0254
    length_lavatory=36*0.0254

    width_galley=36*0.0254
    length_galley=30*0.0254

    length_cabin=length_seats+length_lavatory+length_galley
    print(length_cabin)

    nose_fineness=1.5 #can be altered using aerodynamic data
    nosecone_fineness=2.25 #from data
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

    return length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside, width_fuselage_outside

                    #WING SIZING
def wing(Mach_cruise, S, A, C_L, high=False, mid=False, low=False):
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
        dihedral=dihedral-2 #0-1[deg]
    if low:
        dihedral=dihedral+2
    return taper, b, rootchord, tipchord, sweep_chord_0_5, sweep_chord_0_25, thickness_chord_ratio, dihedral

                    #ENGINE DIMENSIONS

def enginedimensions(n_engines, P_TO_tbp, T_TO_jet, tbp=False, jet=True, jettypeB=False, jettypeC=False):
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
    if jet:
        bypass_ratio=8
        a_0=340.3 #[m/s]
        e_nozzle=0.97
        e_tf=0.75
        T_t4=1500 #1350-1650 [K]
        G=(T_t4/600)-1.25

        massflow=(T_TO_jet*(1+bypass_ratio))/(n_engines*a_0*math.sqrt(5*e_nozzle*G*(1+e_tf*bypass_ratio)))
                    #Intake dimensions
        inlet_spinner_ratio=0.05*(1+(rho_0*a_0)/massflow+3*bypass_ratio/(1+bypass_ratio))
        diameter_inlet=1.65*math.sqrt((massflow/(rho_0*a_0)+0.0050)/(1-(inlet_spinner_ratio)**2))
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

    return length_nacelle, length_f, diameter_highlight, diameter_exit_fan, diameter_gas_generator
