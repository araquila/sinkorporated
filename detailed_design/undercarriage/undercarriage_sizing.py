import parameters as p
import numpy as np

# Insert Scrap Angle and Tip-Over Angle Here
scrap_angle = 8.
tip_lat = 55.

# Vertical Position of the CG in the Fuselage (From the Bottom)
cg_height = 0.4

# Nose Landing Gear Parameters
load_nl = 0.08
l_nl = 2        

# Number of Wheels
N_ml = 4.
N_nl = 2.

##### OUTPUT VALUES #####

# Centre of Gravity Calculation
x_cg = 0.71 * p.MAC + 0.4283 * p.l_fuselage
y_cg = cg_height * p.d_fuselage_outside


# Calculate Main Landing Gear Parameters
load_ml = 1- load_nl
l_ml = ((load_nl * (x_cg-l_nl)) / load_ml) + x_cg
distance_to_tail = p.l_fuselage - p.l_tail - l_ml
wheel_height = distance_to_tail * np.tan(np.radians(scrap_angle))

# Check if the Aircraft will not tip over
z_dist_to_cg = l_ml - x_cg
y_dist_to_cg = y_cg + wheel_height
tip_lon = np.degrees(np.arctan(z_dist_to_cg / y_dist_to_cg))

if tip_lon < scrap_angle:
    print("WARNING, AIRCRAFT WILL TIP OVER")
    print("")
    
# Calculate Lateral Position
y_mlg = ((l_ml-x_cg) + (x_cg-l_nl)) / (np.sqrt((((x_cg-l_nl)**2 * np.tan(np.radians(tip_lat))**2) / ((y_cg+wheel_height)**2)) - 1))

print("AMOUNT OF WHEELS")
print("Number of Nose Landing Gear Wheels:", N_nl)
print("Number of Main Landing Gear Wheels:", N_ml)

print("")

print("WHEEL DIMENSIONS")
print("Wheel = Outer Diameter x Width - Inner Diameter")
print("Nose Landing Gear = 0.46 x 0.11 - 0.25")
print("Main Landing Gear = 0.84 x 0.25 - 0.41")

print("")

print("WHEEL POSITION")
print("Wheel Height:", np.round(wheel_height, decimals=2), "m")
print("Wheel Lateral Position:", np.round(y_mlg, decimals=2), "m")



