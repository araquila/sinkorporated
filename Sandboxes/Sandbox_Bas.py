det_wing_weight(kg_to_pounds(MTOW_tbp), 1.5*ult_load_factor(kg_to_pounds(MTOW_tbp)), metersquared_to_feetsquared(S_tbp), A_tbp, t_c_ratio, taper, sweepqc, metersquared_to_feetsquared(0.05*S_tbp))

det_hor_tail_weight(meter_to_feet(diameter_fuselage_outside), meter_to_feet(span_h), kg_to_pounds(MTOW_tbp),  1.5*ult_load_factor(kg_to_pounds(MTOW_tbp)), meter_to_feet(S_h), meter_to_feet(l_h), sweepqc_h, AR_h, meter_to_feet(0.3*S_h))

det_vert_tail_weight(meter_to_feet(span_v), meter_to_feet(span_v), kg_to_pounds(MTOW_tbp), 1.5*ult_load_factor(kg_to_pounds(MTOW_tbp)), l_v, S_v, sweepLE_v, AR_v, t_c_ratio)

det_fuselage_weight(kg_to_pounds(MTOW_tbp), 1.5*ult_load_factor(kg_to_pounds(MTOW_tbp)), meter_to_feet(length_cabin), metersquared_to_feetsquared(np.pi*diameter_fuselage_outside*length_fuselage), taper, b, sweepqc, LD_cruise_tbp, fuselage_mounted_lg=True)

det_main_lg_weight(kg_to_pounds(MTOW_tbp), 3, meter_to_inch(wheel_height), 4, 2, ms_to_knots(V_stall_tbp))

det_nose_lg_weight(kg_to_pounds(MTOW_tbp), 4.5, meter_to_inch(wheel_height), 2)
