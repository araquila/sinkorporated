import math
n_passenger=60
#List of crew+operational items
n_crew=4
mass_perpassenger=195*0.45359237
mass_peroverheadluggage=30*0.45359237
#general fuselage shape = cylindrical

#Cross section
#cabin layout
n_seats_abreast=4
n_aisles=1
n_rows=n_passenger/n_seats_abreast
width_seats=17*0.0254 #16-17
pitch_seats=31*0.0254 #30-32
head_room=1.7
width_aisles=0.51
height_floor=1.9
n_trolleys=6
n_lavatories=2
width_armrest=0.05
s_clearance=0.05
#cabin width up to shoulder
width_cabin=(n_seats_abreast*width_seats)+(n_seats_abreast+n_aisles+1)*width_armrest+n_aisles*width_aisles+2*s_clearance
#cabin width at head of passengers
width_headroom=width_cabin-2*(width_armrest+s_clearance)-width_seats
#overhead storage
k_os=0.74
n_os=2
area_overhead=0.2
length_cabin=n_rows*pitch_seats+10*0.0254
volume_overhead=(n_os*area_overhead)*length_cabin*k_os

#cargo compartment
density_luggage=170
density_cargo=160
mass_cargo=0
mass_luggage=n_passenger*mass_peroverheadluggage
volume_luggage=mass_luggage/density_luggage
volume_cargo=mass_cargo/density_cargo
volume_cargocompartment=volume_cargo+(volume_luggage-volume_overhead)
print(volume_cargocompartment)
height_cargo=64*0.0254
width_cargo=61.5*0.0254
base_cargo=60.4*0.0254
#structural dimensions
height_floor=0.150 #0.100-0.300
thickness_fuselage_skin_frame=0.150
#diameter_fuselage_inside=math.sqrt((height_houlder/2)**2)+((width_cabin/2)**))
