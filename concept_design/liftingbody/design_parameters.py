from conversion_formulas import *
from atmosphere import atmosphere_calc
from constant_variables import *

# Flight parameters
s_landing = 1400                    #[m]
altitude = 10000 #feet_to_meter(25000)     # 8000 m for tbp, 10000 m for jet
V_landing = 48.93                   #[m/s] maximum landing speed that is allowed on a runway of 1400 m this is set for all aircraft
range_cruise_jet = 1850000          # [m]
range_cruise_tbp = 1850000          # [m]

# Atmospherical parameters at cruise altitude
temperature, pressure, rho, speed_of_sound = atmosphere_calc(altitude, temperature0, temperature_gradient, g, R, gamma)

# General aircraft parameters
C_L_fuselage = 0.08
n_engines = 2

# Jet parameters
C_fe_jet = 0.003
A_jet = 14
e_jet = 0.85                         # Adjust per concept
S_jet = 61
S_wet_jet = 4.5 * S_jet
TOP_jet = 6698
endurance_loiter_jet = 2700

# Cruise speed: first with a given V, then with a given Mach number
# V_cruise_jet =  236.11               # [m/s]
# M_cruise_jet = V_cruise_jet/speed_of_sound
M_cruise_jet = 0.8                   # Use only if V_cruise_jet is deactivated
V_cruise_jet = M_cruise_jet * speed_of_sound

C_L_cruise_jet = 0.4
C_L_max_jet = 2.3
C_L_max_land_jet = 2.6
C_L_max_TO_jet = 1.9

# Empennage jet
V_h_jet = 1.07                        # [-]
V_v_jet = 0.085                          # [-]
nose_landing_pos_jet = 3                # [m]

# Turboprop parameters
C_fe_tbp = 0.003
A_tbp = 14
e_tbp = 0.9                        # Adjust per concept
V_loiter_tbp = 60                   # [m/s]
endurance_loiter_tbp = 2700

# Cruise speed: first with a given V, then with a given Mach number
# V_cruise_tbp = 180                  # [m/s]
# M_cruise_tbp = V_cruise_tbp/speed_of_sound
M_cruise_tbp = 0.6                   # Use only if V_cruise_jet is deactivated
V_cruise_tbp = M_cruise_tbp * speed_of_sound

C_L_cruise_tbp = 0.8
S_tbp = 60
S_wet_tbp = 4.5 * S_tbp
TOP_tbp = 139
C_L_max_tbp = 2.6
C_L_max_land_tbp = 2.6
C_L_max_TO_tbp = 1.9

# Empennage tbp
V_h_tbp = 1.57                          # [-]
V_v_tbp = 0.07                          # [-]
nose_landing_pos_tbp = 3                # [m]

# Engine characteristics
engine_gear_mass = pounds_to_kg(250)    # additional mass due to geared turbofan
thrust_to_weight_jet = 73.21        # [N/kg]
cj_loiter_jet = 19e-6 / 1.16              # (0.4-0.6) [g/j] Propfan: 0.441
cj_cruise_jet = 19e-6 / 1.16              # (0.5-0.9) [g/j] Propfan: 0.441

power_to_weight_tbp = 4000          # [W/kg]
eff_cruise_tbp = 0.85               # [-]
eff_loiter_tbp = 0.77               # [-]
cp_cruise_tbp = 90e-9               # (0.4-0.6) [kg/ns]
cp_loiter_tbp = 90e-9               # (0.5-0.7) [kg/ns]
