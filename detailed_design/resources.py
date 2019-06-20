#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 12:33:22 2019

@author: erikpeeters
"""

import matplotlib.pyplot as plt



strut = 89.15
wingbox	 = 1106
htail	= 135.487902
vtail= 151.8452678
fuselage = 	2766.807516
mlg =	587.8276219
nlg= 142.3155805
nacelle=	357.38
engines	=836.38
props=	330
fueltanks=	471
fuelsystem	=358.7959
flightcontrol=	328.073582
apu=	61.23
flightinstr=	61.93544887
hydr=	90.57517011
elecsyst=	340.8971918
avionics	=766.30
furnish=	161.1688986
aircond=324.87
antiice	=37.07219297
handlinggear=	5.56
wingskin	=229.6
aileronflaps	=84.02


wing = wingbox + wingskin + antiice + aileronflaps + strut
fuselage = fuselage + furnish + handlinggear + aircond
engine = engines + props + nacelle
fuel = fueltanks + fuelsystem 
tail = htail + vtail
control = avionics + flightcontrol + flightinstr + hydr
landinggear = mlg + nlg
aux = apu + elecsyst


plt.clf()
#labels = ['Strut and Strutbox', 'Wingbox', 'Horizontal Tail', 'Vertical Tail', 'Fuselage', 'Main Landing Gear','Nose Landing Gear','Engine Nacelle Group','Engines','Propellers','Fuel Tanks','Fuel System','Flight Controls','APU','Flight Instruments','Hydraulics','Electrical Systems','Avionics','Furnishings','Air Conditioning','Anti-Icing','Handling Gear','Wing Skin','Aileron and Flaps']
#sizes = [round(strut,2), round(wingbox,2),round(htail,2),round(vtail,2),round(fuselage,2),round(mlg,2), round(nlg,2),round(nacelle,2), round(engines,2),round(props,2),round(fueltanks,2),round(fuelsystem,2),round(flightcontrol,2), round(apu,2),round(flightinstr,2), round(hydr,2),round(elecsyst,2),round(avionics,2),round(furnish,2),round(aircond,2), round(antiice,2),round(handlinggear,2), round(wingskin,2),round(aileronflaps,2)]
labels = ['Wing and strut','Fuselage','Engine','Fuel systems and fuel tanks','Tail','Control','Landing gear','Electrical system and APU']
sizes = [round(wing,2), round(fuselage,2),round(engine,2),round(fuel,2),round(tail,2),round(control,2), round(landinggear,2),round(aux,2)]


patches, texts, pcts = plt.pie(sizes, startangle=90,autopct='%1.1f%%',pctdistance=0.80, labeldistance=1.03)
plt.legend(patches, labels, loc='center left')
plt.axis('equal')
plt.tight_layout()

plt.show()








