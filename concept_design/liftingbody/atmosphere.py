import numpy as np

def temperature_calc(altitude, temperature0, temperature_gradient):
    temperature = temperature0 + temperature_gradient * altitude
    return temperature

def pressure_calc(temperature, temperature0, temperature_gradient, g, R):
    pressure = (temperature/temperature0)**(-g/(temperature_gradient*R))
    return pressure

def density_calc(temperature, temperature0, temperature_gradient, g, R):
    rho = (temperature/temperature0)**(-g/(temperature_gradient*R) - 1)
    return rho

def speed_of_sound_calc(rho, R, temperature):
    speed_of_sound = np.sqrt(rho * R * temperature)
    return speed_of_sound

def atmosphere_calc(altitude, temperature0, temperature_gradient, g, R):
    temperature = temperature_calc(altitude, temperature0, temperature_gradient)
    pressure = pressure_calc(temperature, temperature0, temperature_gradient, g, R)
    rho = density_calc(temperature, temperature0, temperature_gradient, g, R)
    speed_of_sound = speed_of_sound_calc(rho, R, temperature)
    return temperature, pressure, rho, speed_of_sound
