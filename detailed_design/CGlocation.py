import numpy as np
import parameters as p

#### FUSELAGE GROUP ###
cgx_fuselage = 0.4 * p.l_fuselage
cgx_tail = 0.9 * p.l_fuselage
cgx_nl = 2 / p.l_fuselage

m_fusgroup = p.m_fuselage + p.m_tail + p.W_nose_landing / p.g
cgx_fusgroup = (cgx_fuselage * p.m_fuselage + cgx_tail * p.m_tail + cgx_nl * (p.W_nose_landing / p.g)) / m_fusgroup


### WING GROUP ###
cgx_wing = 0.35
cgx_fueltank = 0.5
cgx_engine = -0.1
cgx_ml = (10.94 - 0.4283 * p.l_fuselage) / p.MAC

m_winggroup = (p.m_wing + p.m_engine + p.m_fuel_tanks + p.W_main_landing / p.g)
cgx_winggroup = (cgx_wing * p.m_wing + cgx_fueltank * p.m_fuel_tanks + cgx_ml * (p.W_main_landing / p.g)) / m_winggroup

### OEW CG ###
cgx_oew = (cgx_fusgroup * cgx_fusgroup + cgx_winggroup * m_winggroup) / (m_fusgroup + m_winggroup)

