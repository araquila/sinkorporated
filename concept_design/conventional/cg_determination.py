from main_conventional.py import length_fuselage
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
cg_fus = 0.4 * length_fuselage
def x_lemac(cg_fus,c,cgm_wing,mf_wing_group,mf_fuselage,cgm_oew):
    x_lemac = cg_fus+c*(cgm_wing*(mf_wing_group/mf_fuselage)-cgm_oew*(1+mf_wing_group/mf_fuselage))
    return x_lemac
