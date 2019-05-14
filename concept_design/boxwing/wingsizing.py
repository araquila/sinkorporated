from math import *
from class1sizing import enginedimensions, undercarriage, fuselage, tiresizing
from class2_boxwing import *
from wingloadingfunctions import V_stall_calc

def iterempty(MTOW, OEW, WF, LD):
    g =  9.80665


    TW = 0.25
    length_fan_cowl, diameter_nacelle, length_generator, length_nacelle, diameter_highlight, diameter_exit_fan, diameter_gas_generator, diameter_inlet = enginedimensions(2, TW*MTOW)
    #input parameters for wing area
    WS = 3080
    S = MTOW/WS
    S = S * 1.36

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
    tire_pressure, P_mw, P_nw = tiresizing(MTOW/g, 25)
    #print(wheel_height, lateral_position_undercarridge)



    taper1 = 0.3
    sweep1 = (2-taper1/0.2)*180/pi
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

    #print(sweep1, sweep2)
    #print(b, cr2, ct2)

    kgtolbs = 2.20462262
    mtoft = 3.2808399
    MTOWlbs = MTOW/g*kgtolbs
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
    Clhat1 = MTOW/(qhat*S1)
    Clhat2 = MTOW/(qhat*S2)


    sweephalfchord1 = atan((sin(sweep1/180*pi)*b/2+ct1/4)/(b/2))
    sweephalfchord2 = atan((sin(sweep2/180*pi)*b/2+ct2/4)/(b/2))
    t_c1 = min((0.18, ((cos(sweephalfchord1)**3*(M_cross-M_dd*cos(sweephalfchord1))-0.115*Clhat1**1.5)/(cos(sweephalfchord1)**2))))
    t_c2 = min((0.18, ((cos(sweephalfchord2)**3*(M_cross-M_dd*cos(sweephalfchord2))-0.115*Clhat2**1.5)/(cos(sweephalfchord2)**2))))

    #desing cruise speed, diving speed, t_c
    S_control_ail = 0.05*S
    S_control_elev = 0.03*S


    #wing weights as conventional and box method
    wing_weight1_conv = det_wing_weight(MTOWlbs, n_ult, S1*mtoft**2, AR1, 0.13, taper1, sweep1/180*pi, S_control_ail*mtoft**2)/kgtolbs
    wing_weight2_conv = det_wing_weight(MTOWlbs, n_ult, S2*mtoft**2, AR2, 0.13, taper2, sweep2/180*pi, S_control_elev*mtoft**2)/kgtolbs
    wing_weight1_box = det_wing_weight_box(b, S1, taper1, sweep1, MTOW/g, n_ult, V_diving, t_c1, 0.028)
    wing_weight2_box = det_wing_weight_box(b, S2, taper2, sweep2, MTOW/g, n_ult, V_diving, t_c2, 0.028)

    #vertical tail weight
    vert_tail_weight = det_vert_tail_weight(bv*mtoft, bv*mtoft, MTOWlbs, n_ult, sin(sweep1/180*pi)*b/2+ct1/4, Sv*mtoft**2, sweepv, AR_v, 0.18, 1)/kgtolbs


    #fuselage_weight
    fuselage_weight = det_fuselage_weight(MTOWlbs, n_ult, (Fus_len - length_nosecone - length_tailcone)*mtoft, (Fus_len - length_nosecone - length_tailcone)*mtoft*(pi*(diameter_fuselage_outside*mtoft)), taper1, b*mtoft, sweep1/180*pi, 16, 1, 1)/kgtolbs


    #landing gear
    V_stall = V_stall_calc(MTOW, 1.225, 2.1, S)*1.94384
    main_lg_weight = det_main_lg_weight(MTOWlbs, 4.5, wheel_height/0.0254, 4, 2, V_stall, False)/kgtolbs
    nose_lg_weight = det_nose_lg_weight(MTOWlbs, 4.5, wheel_height/0.0254, 2, False)/kgtolbs


    #nacell group weights
    nacelle_group_weight = det_nacelle_group_weight(length_nacelle*mtoft, diameter_nacelle*mtoft, n_ult, 2, pi*diameter_nacelle*length_nacelle, True, 0, 820*kgtolbs, False, True)/kgtolbs


    #engine control weights
    engine_control_weight = det_engine_controls_weight(2, 0.75*Fus_len*2*mtoft)

    #starter weights
    starter_weight = det_starter_weight(2, 820*kgtolbs)

    #Fuel system weights
    fuel_system_weight = det_fuel_system_weight(WF/g/6.6732*kgtolbs, WF/g/6.6732*kgtolbs, 0, 4)/kgtolbs

    #fligth controls
    flight_controls_weight = det_flight_controls_weight(0.08*S + 0.2*Sv, ((length_fuselage*mtoft)**2*MTOWlbs*0.34**2)/(4*32.19), 6, 1)/kgtolbs

    #APU weights
    apu_weight = det_APU_weight(200)

    #insturments weights
    instruments_weight = det_instruments_weight(2, 2, length_fuselage*mtoft, b*mtoft, False, False)/kgtolbs

    #hydraulics
    hydraulics_weight = det_hydraulics_weight(length_fuselage*mtoft, b*mtoft, 6)/kgtolbs

    #electrical
    electrical_weight =  det_electrical_weight(length_fuselage*mtoft, 50, 2, 2)/kgtolbs

    #avionics
    avionics_weight = det_avionics_weight(1100)/kgtolbs

    #furnishings
    furnishings_weight = det_furnishings_weight(2, 13.608*60*kgtolbs, (Fus_len - length_nosecone - length_tailcone)*mtoft*(pi*(diameter_fuselage_outside*mtoft)) )/kgtolbs

    #aircondinintingiogniongionig
    aircond_weight = det_aircond_weight(64, (length_cabin*mtoft)*pi*(diameter_fuselage_outside*mtoft/2)**2, 1100)/kgtolbs

    #anti_ice_weight
    anti_ice_weight = det_anti_ice_weight(MTOWlbs)/kgtolbs

    #handling_gear_weight
    handling_gear_weight = det_handling_gear_weight(MTOWlbs)/kgtolbs

    Empty_mass_new= 800*2 + wing_weight1_box + wing_weight2_box + fuselage_weight + vert_tail_weight + main_lg_weight + nose_lg_weight + nacelle_group_weight + engine_control_weight + starter_weight + fuel_system_weight + flight_controls_weight + apu_weight + instruments_weight + hydraulics_weight + electrical_weight + avionics_weight+ furnishings_weight + aircond_weight + anti_ice_weight + handling_gear_weight

    geometrylistfuselage = (('Length fuselage nose', length_nose), ('Length cabin', length_cabin), ('Lenght fuselage tail', length_tail), ('Length fuselage', length_fuselage), ('Length nosecone', length_nosecone), ('Length tailcone', length_tailcone), ('Diameter outside fuselage', diameter_fuselage_outside))
    geometrylistwings = (('Total wing area', S), ('Wing span', b), ('Fore wing area', S1) ,('Aft wing area', S2), ('Fore root chord', cr1), ('Fore tip chord', ct1), ('Aft root chord', cr2), ('Aft tip chord', ct2), ('Fore qc sweep', sweep1), ('Aft qc sweep', sweep2), ('Fore wing position', Fus_len*frac_qtrchord_fus) )
    geometrylistvtail = (('Height vtail', bv), ('Area vtail', Sv), ('Sweep vtail', sweepv), ('Root chord vtail', crv), ('Tip chord vtail', ctv)  )



    return Empty_mass_new, geometrylistfuselage, geometrylistwings, geometrylistvtail

def cg_prandtl():
    fus_cg = 0.47*Fus_len
    wing1_cg = frac_qtrchord_fus*Fus_len + sin(sweep1*pi/180)*0.35*b/2 + (cr1-(cr1-ct1)*0.35)*0.35
    wing2_cg = Fus_len - crv + sin(sweep2*pi/180)*0.35*b/2 + (cr2-(cr2-ct2)*0.35)*0.35 + sin(sweepv*pi/180)*bv + 0.25*ctv
    tailv_cg = Fus_len - crv + sin(sweepv*pi/180)*0.55*bv/2 + (crv-(crv-ctv)*0.55)*0.42
    engine_cg = 0.75*Fus_len + 0.4*length_nacelle
