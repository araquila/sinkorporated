# MAIN OF THE CONVENTIAL AIRCRAFT SIZING PROGRAM

# Import modules
from class1_liftingbody import *
from class1sizing_liftingbody import *
from power_wingloading_liftingbody import *

# Decide whether you like jet or turboprop:
jet = False
tbp = True

# Gravitional constant
g = 9.8065

# Passengers and crew
n_passenger = 60
M_passenger = 102                   #(including luggage)
n_crew = 4
M_crew_member = 100

# Initial mass and fractions
M_payload = n_passenger * M_passenger
M_crew = n_crew * M_crew_member
f_trapped_fuel = 0.003              # Range 0.001-0.005
M_empty_tbp = 10000                 # Adjust per concept
M_empty_jet = 15000                 # Adjust per concept

# Convert to weights
W_payload = M_payload * g
W_crew = M_crew * g
W_empty_tbp = M_empty_tbp * g
W_empty_jet = M_empty_jet * g

# More general data
n_max_flap = 2
n_max_clean = 2.5
n_min = -1
n_max_man = 4.4
s_landing = 1400 #[m]
rho0 = 1.225 #[kg/m3]
V_landing = 48.93 #[m/s] maximum landing speed that is allowed on a runway of 1400 m this is set for all aircraft
#INPUTS!!!!!!!!!!!!###
rho = 0.4 #[kg/m3] altitude for cruise flight THIS IS INPUT
c = 10 #[m/s] climb rate THIS IS INPUT
weight_fraction = 0.8 #weight fraction of aircraft in cruise wrt MTOW THIS IS INPUT

#graph data
wing_loading_x = np.linspace(0.1,6000,200)

# Initial jet and tbp aircraft parameters
C_L = 0.4                            # during cruise

################################## Jet #######################################
C_fe_jet = 0.003
S_jet = 60                               # Adjust per concept
S_wet_jet = 3.5 * S_jet                  # Adjust per concept
A_jet = 7                          # Adjust per concept
e_jet = 0.8                        # Adjust per concept
cj_loiter_jet = 0.5         # (0.4-0.6) [lbs/lbs/hr]
cj_cruise_jet = 0.6         # (0.5-0.9) [lbs/lbs/hr]
Mach_cruise_jet = 0.8

#Coefficients
C_L_max_jet_clean_min = 1.2
C_L_max_jet_clean_max = 1.8
C_L_max_jet_take_min = 1.6
C_L_max_jet_take_max = 2.2
C_L_max_jet_land_min = 1.8
C_L_max_jet_land_max = 2.8

#take off parameter jet
TOP_aquila_jet_single = 6000 #Take from statistics THIS IS INPUT
TOP_aquile_jet_double = 6000 #Take from statistics THIS IS INPUT
V_cruise_jet = 200 #[m/s] THIS IS INPUT
thrust_setting = 0.9 #THIS IS INPUT
C_D_jet_curr = 0.05 #current C_D value THIS IS INPUT
cV_jet = 0.20 #Climb gradient

################################# Tbp #########################################
C_fe_tbp = 0.003
S_tbp = 60                               # Adjust per concept
S_wet_tbp = 3.5 * S_tbp                  # Adjust per concept
A_tbp = 7                          # Adjust per concept
e_tbp = 0.8                        # Adjust per concept
eff_cruise_tbp = 0.85       # [-]
eff_loiter_tbp = 0.77       # [-]
cp_cruise_tbp = 0.5         # (0.4-0.6) [lbs/hp/hr]
cp_loiter_tbp = 0.6         # (0.5-0.7) [lbs/hp/hr]
Mach_cruise_tbp = 0.6

C_L_max_tbp_clean_min = 1.5
C_L_max_tbp_clean_max = 1.9
C_L_max_tbp_take_min = 1.7
C_L_max_tbp_take_max = 2.1
C_L_max_tbp_land_min = 1.9
C_L_max_tbp_land_max = 3.3

#take off parameter and propulsion
TOP_aquila_tbp = 500 #find from statistics THIS IS INPUT
power_setting = 0.9 #usually at 0.9 THIS IS INPUT
V_cruise_tbp = 150 #[m/s] THIS IS INPUT
C_D_tbp_curr = 0.065 #current CD value THIS IS INPUT
eff_prop = 0.85 #THIS IS INPUT
cV_tbp = 0.083 #from CS23.65 climb gradient

# Weight estimation and wing loading----------------------------------------
for iter in range(1):
    if jet == True:
        # Weight estimations
        MTOW_jet, OEW_jet, W_fuel_jet, C_D_0_jet = Weights_Class_I_jet(W_empty_jet, W_payload, W_crew, C_fe_jet, S_jet, S_wet_jet, A_jet, e_jet, cj_loiter_jet, cj_cruise_jet, f_trapped_fuel)
        W_landing_jet = 0.98 * MTOW_jet
        # Wing loading
        Wing_loading_jet(MTOW_jet, W_landing_jet, S_jet,  \
        C_L_max_jet_take_min, C_L_max_jet_take_max, C_L_max_jet_land_min, \
        C_L_max_jet_land_max, wing_loading_x,  \
        V_landing, V_cruise_jet, rho0, thrust_setting,weight_fraction,rho, TOP_aquila_jet_single, \
        C_D_0_jet, C_D_jet_curr, A_jet, e_jet, c, cV_jet, n_max_man)
        # Wing sizing
        A = A_jet
        S = 60
        Mach_cruise = Mach_cruise_jet
        taper, b, rootchord, tipchord, sweep_chord_0_5, sweep_chord_0_25, thickness_chord_ratio, dihedral = wing(Mach_cruise, S, A, C_L, low=True)
        print(taper, b, rootchord, tipchord, sweep_chord_0_5, sweep_chord_0_25, thickness_chord_ratio, dihedral)

    if tbp == True:
        MTOW_tbp, OEW_tbp, W_fuel_tbp, C_D_0_tbp = Weights_Class_I_tbp(W_empty_tbp, W_payload, W_crew, C_fe_tbp, S_tbp, S_wet_tbp, A_tbp, e_tbp,  eff_loiter_tbp, eff_cruise_tbp, cp_loiter_tbp, cp_cruise_tbp, f_trapped_fuel)
        W_landing_tbp = 0.98 * MTOW_tbp
        # Wing loading
        Wing_loading_tbp(MTOW_tbp, W_landing_tbp, S_tbp,  \
        C_L_max_tbp_take_min, C_L_max_tbp_take_max, C_L_max_tbp_land_min, \
        C_L_max_tbp_land_max, wing_loading_x,  \
        V_landing, V_cruise_tbp, rho0, power_setting,weight_fraction,rho, TOP_aquila_tbp, \
        C_D_0_tbp, C_D_tbp_curr, A_tbp, eff_prop, e_tbp, c, cV_tbp, n_max_man)
        # Wing sizing
        A = A_tbp
        S = 40
        Mach_cruise = Mach_cruise_jet
        taper, b, rootchord, tipchord, sweep_chord_0_5, sweep_chord_0_25, thickness_chord_ratio, dihedral = wing(Mach_cruise, S, A, C_L, low=True)
        print(taper, b, rootchord, tipchord, sweep_chord_0_5, sweep_chord_0_25, thickness_chord_ratio, dihedral)





# Sizing ---------------------------------------------------
# Fuselage
n_seats_abreast = 6
n_aisles = 2

length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside, width_fuselage_outside = fuselage(n_passenger, n_crew, n_seats_abreast, n_aisles)

print(length_nose, length_cabin, length_tail, length_fuselage, diameter_fuselage_outside, width_fuselage_outside)

# Wing
A = A_jet
S = S_jet
Mach_cruise = Mach_cruise_jet
taper, b, rootchord, tipchord, sweep_chord_0_5, sweep_chord_0_25, thickness_chord_ratio, dihedral = wing(Mach_cruise, S, A, C_L, low=True)

#print(taper, b, rootchord, tipchord, sweep_chord_0_5, sweep_chord_0_25, thickness_chord_ratio, dihedral)
