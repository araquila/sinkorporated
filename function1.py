import numpy as np

def det_quarter_chord_sweep(M_cruise, supercritical = False, delta_mach = 0.03):
    """
    determines the quarter chord sweep in radians
    """
    # Delta_mach can range from 0 to 0.05 but is given as 0.03 in ADSEE Slides
    if 0.7 < M_cruise < 1:
        if supercritical:
            sweep = 0.75 * 0.935 / (M_cruise + delta_mach) # 0.935 from statistical data from Torenbeek
            return np.arccos(sweep)
        else:
            sweep = 0.75 / (M_cruise + delta_mach)
            return np.arccos(sweep)
    if M_cruise <= 0:
        raise NameError("Plane flying backwards??")
    if M_cruise >= 1:
        raise NameError("Going supersonic now, are we?")
    else:
        return np.arccos(1)

from Sandboxes.Sandbox_Toon import det_quarter_chord_sweep
