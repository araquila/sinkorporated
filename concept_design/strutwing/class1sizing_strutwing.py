import numpy


def struttedpropengine(cp_tbp = 0.5, eff_prop = 0.85, C_D_0 = 0.015, e = 0.85, C_L_max = 1.5, N = 2., A = 19.5, MTOW = 274680., powerloading = 0.4, wingloading = 1000.):
    P_TO_tbp = MTOW * powerloading
    S_tbp = MTOW / wingloading
    print(P_TO_tbp)

    D_engine_tbp = 0.2*(P_TO_tbp/(1000.*N))**(0.18)
    l_engine_tbp = 0.1*(P_TO_tbp/(1000.*N))**(0.40)
    h_ee_tbp = 1.5*D_engine_tbp
    w_ee_tbp = 1.1*D_engine_tbp
    l_ee_tbp = l_engine_tbp

    D_propeller = 0.55*(P_TO_tbp/(1000*N))**(0.25)

    return h_ee_tbp, w_ee_tbp, l_ee_tbp, D_propeller

print(struttedpropengine())
