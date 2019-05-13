#data is taken from the Fokker F-27-200
mf_structure_tbp = 0.284
mf_empty_tbp = 0.537
mf_wing_group_tbp = 0.118
mf_empennage_tbp = 0.024
mf_fuselage_tbp = 0.099
mf_nacelle_tbp = 0.015
mf_landgear_tbp = 0.042

#data jets, taken from DC-9-10
mf_structure_jet = 0.310
mf_empty_jet = 0.562
mf_wing_group_jet = 0.128
mf_empennage_jet = 0.023
mf_fuselage_jet = 0.093
mf_nacelle_jet = 0.016
mf_landgear_jet = 0.051

cgm_oew = 0.25 #with respect to the MAC
def x_lemac_tbp(cg_fus,c_tbp):
    x_lemac = cg_fus+c_tbp*(0.4*(0.118/0.099)-0.25*(1+0.118/0.099))
    return x_lemac

def x_lemac_jet(cg_fus,c_jet):
    x_lemac = cg_fus+c_jet*(0.4*(0.128/0.093)-0.25*(1+0.128/0.093))
    return x_lemac
