# wingloading and stuff
# Authors: Bas and Toon

rho_0 = 1.225


## input
# OEW
# MTOW
# Fuel weight
# stall speed

## wing loading function [W/S]
def det_wing_loading(C_L_max, V_stall, rho = rho_0):
    # make sure stall speed is higher than regulations

    wing_loading = 0.5 * rho_0 * V_stall**2 * C_L_max
    return wing_loading

# test
# print(det_wing_loading(1.2, 50))

def det_TOP_jet(wing_loading, thrust_loading_jet, C_L_to):

    TOP = wing_loading/thrust_loading_jet
    return TOP

def det_TOP_tbp(wing_loading, )

## power loading function [W/P]

## thrust loading function [T/W]
