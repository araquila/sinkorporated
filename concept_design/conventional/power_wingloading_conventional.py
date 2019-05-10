#pieter responsible
from math import *
from wingloadingfunctions import *
import numpy as np
import matplotlib.pyplot as plt

##############README############################################################
#READ THIS FIRST!!!!!!!!!!!!!!!!!!!!
#READ THIS FIRST!!!!!!!!!!!!!!!!!!!!
#This code needs to be used together with wingloadingfunctions.py, so make sure
#that you put them together in the same folder.
#This code can indicate the allowed wing loadings, Cl, A so that one can size
#the wing and HLD according to the requirements that have been set. By changing
#the input variables can be changed according to your initial sizing.
#running the script will give you diagrams. The white parts are the parts were
#the requirements are met, red areas must be avoided.


###############plotting function################################################
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

#graph data
wing_loading_x = np.linspace(0.1,6000,200)

#############################General Data#######################################
n_max_flap = 2
n_max_clean = 2.5
n_min = -1
n_max_man = 4.4
s_landing = 1400 #[m]
rho0 = 1.225 #[kg/m3]
V_landing = 48.93 #[m/s] maximum landing speed that is allowed on a runway of 1400 m this is set for all aircraft
g = 9.80665 
    
# GENERAL INPUTS!!!!!!!!!!!!###
rho = 0.4292 #[kg/m3] altitude for cruise flight THIS IS INPUT
c = 10 #[m/s] climb rate THIS IS INPUT

weight_fraction = 0.9 #weight fraction during cruise of MTOW THIS IS INPUT


def wingloading_tbp(MTOW_tbp, OEW_tbp, S_tbp, A_tbp, V_cruise_tbp, e_tbp, eff_prop, C_D_0_tbp):
    
    #########DATATURBOPROP##########################################################
#    S_tbp = 84 #[m2] THIS IS INPUT
#    A_tbp = 12 #THIS IS INPUT
    
#    MTOW_tbp = 27925 * g    #[N] fill in your MTOW for turboprop THIS IS INPUT
#    OEW_tbp =     #[N] fill in your OEW for turboprop THIS IS INPUT
    W_landing_tbp = MTOW_tbp #[N] fill in your landing weight for turboprop THIS IS INPUT
    
    #Coefficients
    C_L_max_tbp_clean_min = 1.5
    C_L_max_tbp_clean_max = 1.9
    C_L_max_tbp_take_min = 1.7
    C_L_max_tbp_take_max = 2.1
    C_L_max_tbp_land_min = 1.9
    C_L_max_tbp_land_max = 3.3
    
    C_D_tbp_curr = 0.065 #current CD value THIS IS INPUT
    #take off parameter and propulsion
    TOP_aquila_tbp = 500 #find from statistics THIS IS INPUT
    power_setting = 0.9 #usually at 0.9 THIS IS INPUT
#    V_cruise_tbp = 100 #[m/s] THIS IS INPUT
#    C_D_0_tbp = 0.015 #THIS IS INPUT
#    e_tbp = 0.85 #oswald efficiency factor THIS IS INPUT
#    eff_prop = 0.85 #THIS IS INPUT
    cV_tbp = 0.083 #from CS23.65 climb gradient
    
    
    ####################PROP CALCULATIONS /FOR JETS SCROLL DOWN#####################
    #calculate stall speeds and the wing loading
    V_stall_tbp = V_stall_calc(MTOW_tbp,rho0,C_L_max_tbp_take_max,S_tbp)
    W_S_stall = W_S_calc(rho0,V_stall_tbp,C_L_max_tbp_take_max)
    
    ##########take-off################
    k = TOP_aquila_tbp
    C_L_TO_min_tbp = CL_TO_calc(C_L_max_tbp_take_min)
    C_L_TO_max_tbp = CL_TO_calc(C_L_max_tbp_take_max)
    C_L_TO_range_tbp = np.linspace(C_L_TO_min_tbp,C_L_TO_max_tbp,5)
    TOP_takeoff_tbp = np.zeros(shape=(len(C_L_TO_range_tbp),len(wing_loading_x)))
    for i in range(len(TOP_takeoff_tbp)):
        for j in range(len(wing_loading_x)):
            TOP_takeoff_tbp[i,j] = W_P_calc(wing_loading_x[j],k,C_L_TO_range_tbp[i])
    
    ##########landing#################
    #in this parth the wing loading during landing is calculated
    #this is done by using the maximum allowed landing speed (same for all aircraft)
    #the maximum minimum allowed wing loading is plotted in the wing loading diagram
    f = W_landing_tbp/MTOW_tbp
    C_L_landing_range_tbp = np.linspace(C_L_max_tbp_land_min,C_L_max_tbp_land_max,3)
    W_S_landing_tbp = [0,0,0]
    for i in range(len(C_L_landing_range_tbp)):
        W_S_landing_tbp[i] = W_S_landing_calc(C_L_landing_range_tbp[i],rho0,V_landing,f)
    
    ########Cruise###########
    W_P_cruise_tbp = np.zeros(len(wing_loading_x))
    for i in range(len(wing_loading_x)):
        W_P_cruise_tbp[i]=W_P_cruise_tbp_calc(power_setting,weight_fraction,eff_prop,rho,rho0,C_D_0_tbp,wing_loading_x[i],A_tbp,V_cruise_tbp,e_tbp)
    
    ########Climb Rate#########
    W_P_climb_tbp = np.zeros(len(wing_loading_x))
    for i in range(len(wing_loading_x)):
        W_P_climb_tbp[i] = W_P_climb_calc(eff_prop,c,wing_loading_x[i],rho,A_tbp,e_tbp,C_D_0_tbp)
    
    #########Climb Gradient#####
    W_P_climb_grad_tbp = np.zeros(len(wing_loading_x))
    for i in range(len(wing_loading_x)):
        W_P_climb_grad_tbp[i] = W_P_climb_grad_calc(eff_prop,wing_loading_x[i],cV_tbp,C_D_tbp_curr,C_L_max_tbp_take_min,rho)
    
    #########Maneuvring########
    W_P_maneuvring_tbp = np.zeros(len(wing_loading_x))
    for i in range(len(wing_loading_x)):
        W_P_maneuvring_tbp[i] = W_P_maneuvring_calc(eff_prop,C_D_0_tbp,rho,V_cruise_tbp,wing_loading_x[i],n_max_man,A_tbp,e_tbp)
    
    #####plotting the data#############
    l = np.linspace(0, 0.8, 200)
    x = np.linspace(0, 6000, 200)
    
    # initialise the figure
    fig, ax1 = plt.subplots(1,1)
    xlim = 6000
    ylim = 0.4
    
    # plot lines
    ax1.plot(wing_loading_x,TOP_takeoff_tbp[0,:], label= 'Inline label')
    #ax1.plot(wing_loading_x,TOP_takeoff_tbp[1,:])
    ax1.plot(wing_loading_x,TOP_takeoff_tbp[2,:], label = 'Inline label')
    #ax1.plot(wing_loading_x,TOP_takeoff_tbp[3,:])
    ax1.plot(wing_loading_x,TOP_takeoff_tbp[4,:])
    ax1.axvline(W_S_landing_tbp[0])
    ax1.axvline(W_S_landing_tbp[1])
    ax1.axvline(W_S_landing_tbp[2])
    ax1.plot(wing_loading_x,W_P_cruise_tbp)
    ax1.plot(wing_loading_x,W_P_climb_tbp)
    ax1.plot(wing_loading_x,W_P_climb_grad_tbp)
    #ax1.plot(wing_loading_x,W_P_maneuvring_tbp)
    
    # plot filled parts of the graph
    plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = TOP_takeoff_tbp[4,:], direction = "up")
    plotfiller(ax1, xlim, ylim, vline = W_S_landing_tbp[0], direction = "right")
    plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = W_P_cruise_tbp, direction = "up")
    plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = W_P_climb_tbp, direction = "up")
    plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = W_P_climb_grad_tbp, direction = "up")
    #plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = W_P_maneuvring_tbp, direction = "up")
    
    # plot cosmetics (add some legends/labels/title)
    ax1.set_ylim([0, ylim])
    ax1.set_xlim([0, xlim])
    ax1.legend(["CL_TO =" + str(round(C_L_TO_range_tbp[0],2)), "CL_TO=" + str(round(C_L_TO_range_tbp[2],2)), "CL_TO =" + str(round(C_L_TO_range_tbp[4],2)),
    "landing CL =" + str(round(C_L_landing_range_tbp[0],2)),"landing CL =" + str(round(C_L_landing_range_tbp[1],2)),"landing CL =" + str(round(C_L_landing_range_tbp[2],2)),
    "Cruise A =" + str(round(A_tbp,2)), "Climb Rate A =" + str(round(A_tbp,2)), "Climb Grad A =" + str(round(A_tbp,2))])
    
    plt.show()
    return
    
def wingloading_jet(MTOW_jet,OEW_jet,V_cruise_jet,e_jet,C_D_0_jet,A_jet,S_jet):
    ##################################jets##########################################
    
    #############################DATA JETS##########################################
#    MTOW_jet = 230000 #[N] THIS IS INPUT
#    OEW_jet = 150000 #[N] THIS IS INPUT
    W_landing_jet = 0.9*MTOW_jet #[N] THIS IS INPUT
    
    #Coefficients
    C_L_max_jet_clean_min = 1.2
    C_L_max_jet_clean_max = 1.8
    C_L_max_jet_take_min = 1.6
    C_L_max_jet_take_max = 2.2
    C_L_max_jet_land_min = 1.8
    C_L_max_jet_land_max = 2.8
    
    #take off parameter jet
    TOP_aquila_jet_single = 6000 #Take from statistics THIS IS INPUT
    TOP_aquila_jet_double = 6000 #Take from statistics THIS IS INPUT
#    V_cruise_jet = 200 #[m/s] THIS IS INPUT
#    e_jet = 0.85 #THIS IS INPUT
#    C_D_0_jet = 0.0145 #THIS IS INPUT
    thrust_setting = 0.9 #THIS IS INPUT
#    A_jet = 10 #THIS IS INPUT
    C_D_jet_curr = 0.02 #current C_D value THIS IS INPUT
    cV_jet = 0.20 #Climb gradient
    
    
    #calculate stall speeds and the wing loading
    V_stall_jet = V_stall_calc(MTOW_jet,rho0,C_L_max_jet_take_max,S_jet)
    W_S_stall = W_S_calc(rho0,V_stall_jet,C_L_max_jet_take_max)
    
    ##########take-off################
    k = TOP_aquila_jet_single
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
        T_W_climb_jet[i] = T_W_climb_calc(c,wing_loading_x[i],rho,C_L_max_jet_take_max,C_D_jet_curr)
    
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
    ylim = 0.5
    
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
    #ax1.plot(wing_loading_x,T_W_maneuvring_jet)
    # plot filled parts of the graph
    plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = TOP_takeoff_jet[4,:], direction = "down")
    plotfiller(ax1, xlim, ylim, vline = W_S_landing_jet[0], direction = "right")
    plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = T_W_cruise_jet, direction = "down")
    plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = T_W_climb_jet, direction = "down")
    plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = T_W_climb_grad_jet, direction = "down")
    #plotfiller(ax1, xlim, ylim, x_data = wing_loading_x, data = T_W_maneuvring_jet, direction = "down")
    # plot cosmetics (add some legends/labels/title)
    ax1.set_ylim([0, ylim])
    ax1.set_xlim([0, xlim])
    ax1.legend(["CL_TO =" + str(round(C_L_TO_range_jet[0],2)), "CL_TO =" + str(round(C_L_TO_range_jet[2],2)), "CL_TO =" + str(round(C_L_TO_range_jet[4],2)),
    "landing CL =" + str(round(C_L_landing_range_jet[0],2)),"landing CL =" + str(round(C_L_landing_range_jet[1],2)),"landing CL =" + str(round(C_L_landing_range_jet[2],2)),
    "Cruise A =" + str(round(A_jet,2)), "Climb Rate A =" + str(round(A_jet,2)), "Climb Grad A =" + str(round(A_jet,2))])
    plt.show()
    return
