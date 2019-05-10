#pieter responsible
from math import *
from wingloadingfunctions import *
import numpy as np
import matplotlib.pyplot as plt
from main_boxwing import class1box

# plotting function
def plotfiller(ax, xlim, ylim, x_data = 0, data = 0, vline = 0, direction = "right", alpha = 0.3, color = 'red'):
    """
    Inputs
    ax = axis object
    xlim/ylim = limits in x and y
    x_data = steps in x direction
    data = line for which you want to plotting
    vline = the vertical line
    """
    if direction == "right":
        ax.axvspan(vline, xlim, alpha = alpha, facecolor = color)
        return
    if direction == "left":
        ax.axvspan(0, vline, alpha = alpha, facecolor = color)
        return
    if direction == "down":
        ax.fill_between(x_data, data, alpha = alpha, facecolor = color)
        return
    if direction == "up":
        topline = np.linspace(ylim, ylim, len(data))
        ax.fill_between(x_data, topline, data, alpha = alpha, facecolor = color)
        return

#general
n_max_flap = 2
n_max_clean = 2.5
n_min = -1
n_max_man = 4.4
s_landing = 1400 #[m]
rho0 = 1.225 #[kg/m3]
rho = 0.4 #[kg/m3]
c = 20 #[m/s]
V_landing = 48.93 #[m/s] maximum landing speed that is allowed on a runway of 1400 m
weight_fraction = 1. #weight fraction of MTOW
S = 48000 #[m2]
#graph data
wing_loading_x = np.linspace(0.1,6000,200)

MTOW_jet, OEW_jet, W_fuel_jet = class1box()
#############DATA JETS##################
#MTOW_jet = 275365.44 #[N]
#OEW_jet = 140811.41 #[N]
W_landing_jet = MTOW_jet*weight_fraction #[N]

#Coefficients
C_L_max_jet_clean_min = 1.4
C_L_max_jet_clean_max = 1.8
C_L_max_jet_take_min = 1.8
C_L_max_jet_take_max = 2.662
C_L_max_jet_land_min = 1.6
C_L_max_jet_land_max = 2.6

#take off parameter jet
TOP_aquila_jet_single = 6698
TOP_aquila_jet_double = 6698
V_cruise_jet = 230
e_jet = 1.2
C_D_0_jet = 0.0145
thrust_setting = 0.9
A_jet = 12
C_D_jet_curr = 4*C_D_0_jet #current C_D value
cV_jet = 0.2
C_L_jet_curr = sqrt(3*C_D_0_jet*pi*A_jet*e_jet)


########jets################
#calculate stall speeds and the wing loading
V_stall_jet = V_stall_calc(MTOW_jet,rho0,C_L_max_jet_take_max,S)
W_S_stall = W_S_calc(rho0,V_stall_jet,C_L_max_jet_take_max)

##########take-off################
k = TOP_aquila_jet_double
C_L_TO_min_jet = CL_TO_calc(C_L_max_jet_take_min)
C_L_TO_max_jet = CL_TO_calc(C_L_max_jet_take_max)
C_L_TO_range_jet = np.linspace(C_L_TO_min_jet,C_L_TO_max_jet,5)
TOP_takeoff_jet = np.zeros(shape=(len(C_L_TO_range_jet),len(wing_loading_x)))
for i in range(len(TOP_takeoff_jet)):
    for j in range(len(wing_loading_x)):
        TOP_takeoff_jet[i,j] = T_W_calc(wing_loading_x[j],k,C_L_TO_range_jet[i])

######landing###########
#in this parth the wing loading during landing is calculated
#this is done by using the maximum allowed landing speed (same for all aircraft)
#the maximum minimum allowed wing loading is plotted in the wing loading diagram
f = W_landing_jet/MTOW_jet
C_L_landing_range_jet = np.linspace(C_L_max_jet_land_min,C_L_max_jet_land_max,3)
W_S_landing_jet = [0,0,0]
for i in range(len(C_L_landing_range_jet)):
    W_S_landing_jet[i] = W_S_landing_calc(C_L_landing_range_jet[i],rho0,V_landing,f)

######CRUISE########
T_W_cruise_jet= np.zeros(len(wing_loading_x))
for i in range(len(wing_loading_x)):
    T_W_cruise_jet[i] = T_W_cruise_jet_calc(thrust_setting,weight_fraction,rho,rho0,C_D_0_jet,wing_loading_x[i],A_jet,V_cruise_jet,e_jet)

########Climb Rate#########
T_W_climb_jet = np.zeros(len(wing_loading_x))
for i in range(len(wing_loading_x)):
    T_W_climb_jet[i] = T_W_climb_calc(c,wing_loading_x[i],rho,C_L_jet_curr,C_D_jet_curr)

########Climb Gradient#########
T_W_climb_grad_jet = np.zeros(len(wing_loading_x))
for i in range(len(wing_loading_x)):
    T_W_climb_grad_jet[i] = T_W_climb_grad_calc(cV_jet,C_D_0_jet,A_jet,e_jet)

#####Maneauvring#########
T_W_maneuvring_jet = np.zeros(len(wing_loading_x))
for i in range(len(wing_loading_x)):
    T_W_maneuvring_jet[i] = T_W_maneuvring_jet_calc(C_D_0_jet,rho,V_cruise_jet,wing_loading_x[i],n_max_man,A_jet,e_jet)
# the data
l = np.linspace(0, 0.8, 200)
x = np.linspace(0, 4000, 200)

# initialise the figure
fig, ax1 = plt.subplots(1,1)
xlim = 4000
ylim = 0.4

# plot lines
ax1.plot(wing_loading_x,TOP_takeoff_jet[0,:])
#ax1.plot(wing_loading_x,TOP_takeoff_jet[1,:])
ax1.plot(wing_loading_x,TOP_takeoff_jet[2,:])
#ax1.plot(wing_loading_x,TOP_takeoff_jet[3,:])
ax1.plot(wing_loading_x,TOP_takeoff_jet[4,:])
ax1.axvline(W_S_landing_jet[0])
ax1.axvline(W_S_landing_jet[1])
ax1.axvline(W_S_landing_jet[2])
ax1.plot(wing_loading_x,T_W_cruise_jet)
ax1.plot(wing_loading_x,T_W_climb_jet)
ax1.plot(wing_loading_x,T_W_climb_grad_jet)

# plot filled parts of the graph
plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = TOP_takeoff_jet[4,:], direction = "down")
plotfiller(ax1, xlim, ylim, vline = W_S_landing_jet[0], direction = "right")
plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = T_W_cruise_jet, direction = "down")
plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = T_W_climb_jet, direction = "down")
plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = T_W_climb_grad_jet, direction = "down")
# plot cosmetics (add some legends/labels/title)
ax1.set_ylim([0, ylim])
ax1.set_xlim([0, xlim])
plt.xlabel('W/S')
plt.ylabel('T/W')
ax1.legend(["takeoff CL =" + str(round(C_L_TO_range_jet[0],2)), "takeoff CL =" + str(round(C_L_TO_range_jet[2],2)), "takeoff CL =" + str(round(C_L_TO_range_jet[4],2)),
"landing CL =" + str(round(C_L_landing_range_jet[0],2)),"landing CL =" + str(round(C_L_landing_range_jet[1],2)),"landing CL =" + str(round(C_L_landing_range_jet[2],2)),
"Cruise A =" + str(round(A_jet,2)), "Climb Rate A =" + str(round(A_jet,2)), "Climb Grad A =" + str(round(A_jet,2))])
plt.show()
