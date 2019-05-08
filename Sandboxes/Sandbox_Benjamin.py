import math
Mach_cruise=0.8
A=11
S=55
C_L=0.55
high=False
mid=False
low=True
#AERODYNAMICS
Mach_dd=Mach_cruise+0.03
Mach_t=0.935
if Mach_cruise <0.7:
    sweep_chord_0_25=0
else :
    sweep_chord_0_25=math.acos(0.75*(Mach_t/Mach_dd)) #[rad]
#geometric parameters
#taper ratio
taper=0.2*(2-sweep_chord_0_25)
#wingspan
b=math.sqrt(S*A)
#Chord lengths
rootchord=(2*S)/((1+taper)*b)
tipchord=taper*rootchord
#half wing sweep
sweep_chord_0_5=10*(math.pi/180)
#thickness-to-chord ratio
cos=math.cos(sweep_chord_0_5)
tip_chord_ratio=((cos**3)*(Mach_t-Mach_dd*cos)-(0.115*C_L**1.5))/(cos**2)
if tip_chord_ratio > 0.18:
    tip_chord_ratio=0.18
#dihedral angle
dihedral=3-(sweep_chord_0_5*(180/math.pi)/10)
if high:
    dihedral=dihedral+2 #0-1[deg]
if low:
    dihedral=dihedral+2
print(dihedral)
