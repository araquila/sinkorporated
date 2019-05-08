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

#graph data
W_S_x = np.linspace(0,4000,200)
#turboprop data
W = 200000 #[N]
S = 55 #[m2]

#data props
CLmax_turboprop_clean_min = 1.5
CLmax_turboprop_clean_max = 1.9
CLmax_turboprop_take_min = 1.7
CLmax_turboprop_take_max = 2.1
CLmax_turboprop_land_min = 1.9
CLmax_turboprop_land_max = 3.3

#take off parameter turboprop
TOP_aquila_turboprop = 580
#data jets
CLmax_jet_clean_min = 1.2
CLmax_jet_clean_max = 1.8
CLmax_jet_take_min = 1.6
CLmax_jet_take_max = 2.2
CLmax_jet_land_min = 1.8
CLmax_jet_land_max = 2.8

#take off parameter jet
TOP_aquila_jet_single = 6000
TOP_aquile_jet_double = 6000

#props, for jets scroll DOWN################
#calculate stall speeds and the wing loading
V_stall = V_stall_calc(W,rho0,CLmax_turboprop_take_max,S)
W_S_stall = W_S_calc(rho0,V_stall,CLmax_turboprop_take_max)

##########take-off################
k = TOP_aquila_turboprop
CL_TO_min_turboprop = CL_TO_calc(CLmax_turboprop_take_min)
CL_TO_max_turboprop = CL_TO_calc(CLmax_turboprop_take_max)
CL_TO_range_turboprop = np.linspace(CL_TO_min_turboprop,CL_TO_max_turboprop,5)
TOP_takeoff_turboprop = np.zeros(shape=(len(CL_TO_range_turboprop),len(W_S_x)))
for i in range(len(TOP_takeoff_turboprop)):
    for j in range(len(W_S_x)):
        TOP_takeoff_turboprop[i,j] = W_P_calc(W_S_x[j],k,CL_TO_range_turboprop[i])

# the data
l = np.linspace(0, 0.8, 200)
x = np.linspace(0, 4000, 200)
vertical = 3000

# initialise the figure
fig, ax1 = plt.subplots(1,1)
xlim = 4000
ylim = 1

# plot lines
ax1.plot(W_S_x,TOP_takeoff_turboprop[0,:])
ax1.plot(W_S_x,TOP_takeoff_turboprop[1,:])
ax1.plot(W_S_x,TOP_takeoff_turboprop[2,:])
ax1.plot(W_S_x,TOP_takeoff_turboprop[3,:])
ax1.plot(W_S_x,TOP_takeoff_turboprop[4,:])
ax1.axvline(vertical)

# plot filled parts of the graph
plotfiller(ax1, xlim, ylim, x_data = W_S_x, data = TOP_takeoff_turboprop[4,:], direction = "up")
plotfiller(ax1, xlim, ylim, vline = vertical, direction = "right")

# plot cosmetics (add some legends/labels/title)
ax1.set_ylim([0, ylim])
ax1.set_xlim([0, xlim])

plt.show()


##########landing#################

########jets################
#calculate stall speeds and the wing loading
V_stall = V_stall_calc(W,rho0,CLmax_jet_take_max,S)
W_S_stall = W_S_calc(rho0,V_stall,CLmax_jet_take_max)

##########take-off################
k = TOP_aquila_jet_single
CL_TO_min_jet = CL_TO_calc(CLmax_jet_take_min)
CL_TO_max_jet = CL_TO_calc(CLmax_jet_take_max)
CL_TO_range_jet = np.linspace(CL_TO_min_jet,CL_TO_max_jet,5)
TOP_takeoff_jet = np.zeros(shape=(len(CL_TO_range_jet),len(W_S_x)))
for i in range(len(TOP_takeoff_jet)):
    for j in range(len(W_S_x)):
        TOP_takeoff_jet[i,j] = T_W_calc(W_S_x[j],k,CL_TO_range_jet[i])
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
