# wingloading and stuff
# Authors: Bas and Toon

rho_0 = 1.225


## input
# OEW
# MTOW
# Fuel weight
# stall speed

## wing loading function [W/S]
def det_wing_loading(C_L_max, V_Stall, rho = rho_0):
    # make sure stall speed is higher than regulations

    wing_loading = 0.5 * rho_0 * V_stall^2 * C_L_max
    return wing_loading

print(det_wing_loading(1.2, 50))

## power loading function [W/P]

## thrust loading function [T/W]
