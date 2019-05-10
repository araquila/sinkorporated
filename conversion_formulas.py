            #SI tu US

def meter_to_inch(length_in_meter):
    length_in_inch = 39.3701 * length_meter
    return length_in_inch

def meter_to_feet(length_in_meter):
    length_in_feet = 3.28084 * length_in_meter
    return length_in_feet

def metersquared_to_feetsquared(area_in_metersquared):
    area_in_feetsquared = 10.7639 * area_in_metersquared
    return area_in_feetsquared

def metercubed_to_feetcubed(volume_in_metercubed):
    volume_in_feetcubed = 35.3147 * volume_in_metercubed
    return volume_in_feetcubed

def kg_to_pounds(mass_in_kg):
    mass_in_pounds = 2.20462 * mass_in_kg
    return mass_in_pounds

def SIdensity_to_USdensity(SI_density):
    USdensity = 0.062428 * SI_density
    return USdensity

def ms_to_knots(velocity_in_ms):
    velocity_in_knots = 1.94384 * velocity_in_ms
    return velocity_in_knots

                #US to SI

def feet_to_meter(length_in_feet):
    length_in_meter = length_in_feet / 3.28084
    return length_in_meter

def feetsquared_to_metersquared(area_in_feetsquared):
    area_in_metersquared = area_in_feetsquared / 10.7639
    return area_in_metersquared

def feetcubed_to_metercubed(volume_in_feetcubed):
    volume_in_metercubed = 35.3147 * volume_in_feetcubed
    return volume_in_metercubed

def pounds_to_kg(mass_in_pounds):
    mass_in_kg = 2.20462 * mass_in_pounds
    return mass_in_kg

def USdensity_to_ISdensity(US_density):
    ISdensity = 0.062428 * US_density
    return ISdensity

def ms_to_knots(velocity_in_knots):
    velocity_in_ms = 1.94384 * velocity_in_knots
    return velocity_in_ms
