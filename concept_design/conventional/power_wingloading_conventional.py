#pieter responsible
from math import *
from wingloadingfunctions import *
import numpy as np
import matplotlib.pyplot as plt

# plotting function
def plotfiller(ax, xlim, ylim, x_data = 0, data = 0, vline = 0, direction = "right", alpha = 0.5, color = 'red'):
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

#load data
n_max_flap = 2
n_max_clean = 2.5
n_min = -1
s_landing = 1400 #[m]
rho0 = 1.225 #[kg/m3]
rho = 1.225 #[kg/m3]
V_landing = 48.93 #[m/s]

#graph data
W_S_x = np.linspace(0,4000,200)
#turboprop data
W_TO = 200000. #[N]
W_L = 120000.
S = 55 #[m2]

#data props
MTOW_tbp = 200000 #[N]
C_L_max_tbp_clean_min = 1.5
C_L_max_tbp_clean_max = 1.9
C_L_max_tbp_take_min = 1.7
C_L_max_tbp_take_max = 2.1
C_L_max_tbp_land_min = 1.9
C_L_max_tbp_land_max = 3.3

#take off parameter turboprop
TOP_aquila_tbp = 580
#data jets
MTOW_jet = 230000 #[N]
C_L_max_jet_clean_min = 1.2
C_L_max_jet_clean_max = 1.8
C_L_max_jet_take_min = 1.6
C_L_max_jet_take_max = 2.2
C_L_max_jet_land_min = 1.8
C_L_max_jet_land_max = 2.8

#take off parameter jet
TOP_aquila_jet_single = 6000
TOP_aquile_jet_double = 6000

#props, for jets scroll DOWN################
#calculate stall speeds and the wing loading
V_stall = V_stall_calc(W_TO,rho0,C_L_max_tbp_take_max,S)
W_S_stall = W_S_calc(rho0,V_stall,C_L_max_tbp_take_max)

##########take-off################
k = TOP_aquila_tbp
C_L_TO_min_tbp = CL_TO_calc(C_L_max_tbp_take_min)
C_L_TO_max_tbp = CL_TO_calc(C_L_max_tbp_take_max)
C_L_TO_range_tbp = np.linspace(C_L_TO_min_tbp,C_L_TO_max_tbp,5)
TOP_takeoff_tbp = np.zeros(shape=(len(C_L_TO_range_tbp),len(W_S_x)))
for i in range(len(TOP_takeoff_tbp)):
    for j in range(len(W_S_x)):
        TOP_takeoff_tbp[i,j] = W_P_calc(W_S_x[j],k,C_L_TO_range_tbp[i])

# the data
l = np.linspace(0, 0.8, 200)
x = np.linspace(0, 4000, 200)
vertical = 3000

# initialise the figure
fig, ax1 = plt.subplots(1,1)
xlim = 4000
ylim = 1

# plot lines
ax1.plot(W_S_x,TOP_takeoff_tbp[0,:])
ax1.plot(W_S_x,TOP_takeoff_tbp[1,:])
ax1.plot(W_S_x,TOP_takeoff_tbp[2,:])
ax1.plot(W_S_x,TOP_takeoff_tbp[3,:])
ax1.plot(W_S_x,TOP_takeoff_tbp[4,:])
ax1.axvline(vertical)

# plot filled parts of the graph
plotfiller(ax1, xlim, ylim, x_data = W_S_x, data = TOP_takeoff_tbp[4,:], direction = "up")
plotfiller(ax1, xlim, ylim, vline = vertical, direction = "right")

# plot cosmetics (add some legends/labels/title)
ax1.set_ylim([0, ylim])
ax1.set_xlim([0, xlim])

plt.show()


##########landing#################
f = W_L/W_TO
W_S_landing = W_S_landing_calc(C_L_max_tbp_land_max,rho,V_landing,f)
print(W_S_landing)

########jets################
#calculate stall speeds and the wing loading
V_stall = V_stall_calc(W_TO,rho0,C_L_max_jet_take_max,S)
W_S_stall = W_S_calc(rho0,V_stall,C_L_max_jet_take_max)

##########take-off################
k = TOP_aquila_jet_single
C_L_TO_min_jet = CL_TO_calc(C_L_max_jet_take_min)
C_L_TO_max_jet = CL_TO_calc(C_L_max_jet_take_max)
C_L_TO_range_jet = np.linspace(C_L_TO_min_jet,C_L_TO_max_jet,5)
TOP_takeoff_jet = np.zeros(shape=(len(C_L_TO_range_jet),len(W_S_x)))
for i in range(len(TOP_takeoff_jet)):
    for j in range(len(W_S_x)):
        TOP_takeoff_jet[i,j] = T_W_calc(W_S_x[j],k,C_L_TO_range_jet[i])
# the data
l = np.linspace(0, 0.8, 200)
x = np.linspace(0, 4000, 200)
vertical = 3000

# initialise the figure
fig, ax1 = plt.subplots(1,1)
xlim = 4000
ylim = 1

# plot lines
ax1.plot(W_S_x,TOP_takeoff_jet[0,:])
ax1.plot(W_S_x,TOP_takeoff_jet[1,:])
ax1.plot(W_S_x,TOP_takeoff_jet[2,:])
ax1.plot(W_S_x,TOP_takeoff_jet[3,:])
ax1.plot(W_S_x,TOP_takeoff_jet[4,:])
ax1.axvline(vertical)

# plot filled parts of the graph
plotfiller(ax1, xlim, ylim, x_data = W_S_x, data = TOP_takeoff_jet[4,:], direction = "down")
plotfiller(ax1, xlim, ylim, vline = vertical, direction = "right")

# plot cosmetics (add some legends/labels/title)
ax1.set_ylim([0, ylim])
ax1.set_xlim([0, xlim])

plt.show()
