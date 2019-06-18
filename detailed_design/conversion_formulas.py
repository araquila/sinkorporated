            #SI tu US

def meter_to_inch(length_in_meter):
    length_in_inch = 39.3701 * length_in_meter
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

def metercubed_to_USgallon(volume_in_metercubed):
    volume_in_USgallon = 264.172 * volume_in_metercubed
    return volume_in_USgallon

def kilogram_square_meter_to_pound_square_foot(moi_in_kilogram_square_meter):
    moi_in_pound_square_foot = 23.73036 * moi_in_kilogram_square_meter
    return moi_in_pound_square_foot

                #US to SI

def feet_to_meter(length_in_feet):
    length_in_meter = length_in_feet / 3.28084
    return length_in_meter

def feetsquared_to_metersquared(area_in_feetsquared):
    area_in_metersquared = area_in_feetsquared / 10.7639
    return area_in_metersquared

def feetcubed_to_metercubed(volume_in_feetcubed):
    volume_in_metercubed = volume_in_feetcubed / 35.3147
    return volume_in_metercubed

def pounds_to_kg(mass_in_pounds):
    mass_in_kg = mass_in_pounds / 2.20462
    return mass_in_kg

def USdensity_to_ISdensity(US_density):
    ISdensity = US_density / 0.062428
    return ISdensity

def ms_to_knots(velocity_in_knots):
    velocity_in_ms = velocity_in_knots / 1.94384
    return velocity_in_ms

def USgallon_to_metercubed(volume_in_USgallon):
    volume_in_metercubed = volume_in_USgallon / 264.172
    return volume_in_metercubed

def pound_square_foot_to_kilogran_square_meter(moi_in_pound_square_foot):
    moi_in_kilogram_square_meter = moi_in_pound_square_foot / 23.73036
    return moi_in_kilogram_square_meter
