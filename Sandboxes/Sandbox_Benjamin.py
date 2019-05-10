import math
import matplotlib.pyplot as plt
import numpy as np
n_passenger=60
n_aisles=1
n_seats_abreast=4
n_crew=4

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
height_shoulder=1.3
width_aisles=0.50 #min 38cm
n_trolleys=6
n_lavatories=1
n_galleys=1
width_armrest=0.05
s_clearance=0.02
#cabin width up to shoulder
width_cabin=(n_seats_abreast*width_seats)+(n_seats_abreast)*width_armrest+n_aisles*width_aisles+2*s_clearance
#cabin width at head of passengers
width_headroom=width_cabin-2*(width_armrest+s_clearance)-width_seats
#overhead storage
k_os=1
n_os=2
area_overhead=0.1
length_seats=n_rows*pitch_seats+(10*0.0254)
volume_overhead=(n_os*area_overhead)*length_seats*k_os

#cargo compartment
density_luggage=170
density_cargo=160
mass_cargo=0
mass_luggage=n_passenger*mass_peroverheadluggage
volume_luggage=mass_luggage/density_luggage
volume_cargo=mass_cargo/density_cargo
volume_cargocompartment=volume_cargo+(volume_luggage-volume_overhead)
length_cargocompartment=volume_cargocompartment/(width_cabin*height_shoulder)
#structural dimensions
thickness_floor=0.150 #0.100-0.300
thickness_fuselage_skin_frame=0.150

diameter_fuselage_inside=2*math.sqrt(((height_shoulder/2)**2)+((width_cabin/2)**2))

diameter_fuselage_outside=1.045*diameter_fuselage_inside+0.084



height_cargo=(diameter_fuselage_inside/2)-((height_shoulder/2)+thickness_floor)
height_aisle=diameter_fuselage_inside-(thickness_floor+height_cargo)
                #FUSELAGE TOP VIEW
width_lavatory=36*0.0254
length_lavatory=36*0.0254

width_galley=36*0.0254
length_galley=30*0.0254

length_cabin=length_seats+length_lavatory+length_galley+length_cargocompartment

nose_fineness=1 #can be altered using aerodynamic data
nosecone_fineness=2 #from data
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
print(length_fuselage)
