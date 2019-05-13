# plotting and stuff

# import matplotlib.pyplot as plt
# import numpy as np
#
# def plotfiller(ax, xlim, ylim, x_data = 0, data = 0, vline = 0, direction = "right", alpha = 0.5, color = 'red'):
#     """
#     Inputs
#     ax = axis object
#     xlim/ylim = limits in x and y
#     x_data = steps in x direction
#     data = line for which you want to plotting
#     vline = the vertical line
#     """
#     if direction == "right":
#         ax.axvspan(vline, xlim, alpha = alpha, facecolor = color)
#         return
#     if direction == "left":
#         ax.axvspan(0, vline, alpha = alpha, facecolor = color)
#         return
#     if direction == "down":
#         ax.fill_between(x_data, data, alpha = alpha, facecolor = color)
#         return
#     if direction == "up":
#         topline = np.linspace(ylim, ylim, len(data))
#         ax.fill_between(x_data, topline, data, alpha = alpha, facecolor = color)
#         return
#
#
# # Example
#
# # the data
# l = np.linspace(0, 0.8, 200)
# x = np.linspace(0, 4000, 200)
# vertical = 3000
#
# # initialise the figure
# fig, ax1 = plt.subplots(1,1)
# xlim = 4000
# ylim = 1
#
# # plot lines
# ax1.plot(x, l)
# ax1.axvline(vertical)
#
# # plot filled parts of the graph
# plotfiller(ax1, xlim, ylim, x_data = x, data = l, direction = "down")
# plotfiller(ax1, xlim, ylim, vline = vertical, direction = "right")
#
# # plot cosmetics (add some legends/labels/title)
# ax1.set_ylim([0, ylim])
# ax1.set_xlim([0, xlim])
#
# plt.show()

# Wing Sizing
# import numpy as np
#
# def det_quarter_chord_sweep(M_cruise, supercritical = False, delta_mach = 0.03):
#     """
#     determines the quarter chord sweep in radians
#     """
#     # Delta_mach can range from 0 to 0.05 but is given as 0.03 in ADSEE Slides
#     if 0.7 < M_cruise < 1:
#         if supercritical:
#             sweep = 0.75 * 0.935 / (M_cruise + delta_mach) # 0.935 from statistical data from Torenbeek
#             return np.arccos(sweep)
#         else:
#             sweep = 0.75 / (M_cruise + delta_mach)
#             return np.arccos(sweep)
#     if M_cruise <= 0:
#         raise NameError("Plane flying backwards??")
#     if M_cruise >= 1:
#         raise NameError("Going supersonic now, are we?")
#     else:
#         return np.arccos(1)
#
# def det_planform(S, AR, M_cruise, C_L_cruise, sweep, supercritical = False, delta_mach = 0.03):
#     # Delta_mach can range from 0 to 0.05 but is given as 0.03 in ADSEE Slides
#     b = np.sqrt(AR * S)
#     taper = 0.2 * (2 - sweep)
#     root_chord = (2 * S)/((1 + taper) * b)
#     tip_chord = taper * root_chord
#     half_chord_sweep = np.arctan(np.tan(sweep) - (4 / AR) * (0.25 * (1 - taper)/(1 + taper)))
#     if 0.7 < M_cruise < 1:
#         if supercritical:
#             t_c_ratio = min(0.18, (np.cos(half_chord_sweep)**3 * (0.935 - (M_cruise + delta_mach) * np.cos(half_chord_sweep)) - 0.115 * C_L_cruise**1.5) / np.cos(half_chord_sweep)**2)
#         else:
#             t_c_ratio = min(0.18, (np.cos(half_chord_sweep)**3 * (1 - (M_cruise + delta_mach) * np.cos(half_chord_sweep)) - 0.115 * C_L_cruise**1.5) / np.cos(half_chord_sweep)**2)
#     if M_cruise <= 0:
#         raise NameError("Plane flying backwards??")
#     if M_cruise >= 1:
#         raise NameError("Going supersonic now, are we?")
#     return b, taper, root_chord, tip_chord, t_c_ratio
#
# def det_dihedral_angle(sweep, high = False, mid = False, low = False):
#     dihedral = sweep * 18 / np.pi
#     if high:
#         angle = 1 - dihedral
#         return angle
#     if mid:
#         angle = 3 - dihedral
#         return angle
#     if low:
#         angle = 5 - dihedral
#         return angle
#     else:
#         raise NameError("Where is the wing?")
#
#
# def wing(Mach_cruise, S, A, C_L, high=False, mid=False, low=False):
#     #AERODYNAMICS
#     Mach_dd=Mach_cruise+0.03
#     Mach_t=0.935
#     #quarter chord sweep
#     if Mach_cruise <0.7:
#         sweep_chord_0_25=0
#     else :
#         sweep_chord_0_25=math.acos(0.75*(Mach_t/Mach_dd)) #[rad]
#     #geometric parameters
#     #taper ratio
#     taper=0.2*(2-sweep_chord_0_25)
#     #wingspan
#     b=math.sqrt(S*A)
#     #Chord lengths
#     rootchord=(2*S)/((1+taper)*b)
#     tipchord=taper*rootchord
#     #thickness-to-chord ratio
#     cos=math.cos(sweep_chord_0_5)
#     thickness_chord_ratio=((cos**3)*(Mach_t-Mach_dd*cos)-(0.115*C_L**1.5))/(cos**2)
#     if thickness_chord_ratio > 0.18:
#         thickness_chord_ratio=0.18
#
#     return taper, b, rootchord, tipchord, sweep_chord_0_5, sweep_chord_0_25, thickness_chord_ratio
import os
import sys
import numpy as np

sys.path.append(os.getcwd())

from concept_design.strutwing.class2_strutwing import *
from conversion_formulas import *

MTOW_tbp = 20000

handling_gear_weight = det_handling_gear_weight(kg_to_pounds(MTOW_tbp))
anti_ice_weight = det_anti_ice_weight(kg_to_pounds(MTOW_tbp))
pres_vol = np.pi / 4 * diameter_fuselage_inside**2 * (length_nose + length_nose)
aircond_weight = aircond_weight(n_passenger, pres_vol)
