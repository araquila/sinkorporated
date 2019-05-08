import math
tbp=True
P_TO_tbp=4000000
n_engines=2
if tbp==True:

    diameter_engine=0.2*(P_TO_tbp/(1000*n_engines))**0.18
    length_engine=0.1*(P_TO_tbp/(1000*n_engines))**0.4

    diameter_propeller=0.55*(P_TO_tbp/(1000*n_engines))**0.25
print(diameter_engine)
print(diameter_propeller)
print( length_engine)
