import numpy as np

def non_recurring_cost(m_wing,m_empennage,m_fuselage,m_gear,m_engines, m_systems, m_payloads):
    """Returns the total non-recurring cost in million USD(April 2019) based on aircraft mass in kg"""
    mass_vector_kg = np.array([m_wing,m_empennage,m_fuselage,m_gear,m_engines, m_systems, m_payloads])
    mass_vector_lbs = mass_vector_kg * 2.20462262
    cost_per_lbs_vector = np.array([17731,52156,32093,2499,8691,34307,10763])
    cost_vector = mass_vector_lbs*cost_per_lbs_vector

    c_total_nonrecurring_2002 = np.sum(cost_vector)
    c_total_nonrecurring_2019 = c_total_nonrecurring_2002*1.44

    return round(c_total_nonrecurring_2019/1000000,2)

def recurring_cost(n_aircraft,m_wing,m_empennage,m_fuselage,m_gear,m_engines, m_systems, m_payloads,m_assembly):
    """Returns the total recurring cost in million USD(April 2019) based on aircraft mass in kg"""
    mass_vector_kg = np.array([m_wing,m_empennage,m_fuselage,m_gear,m_engines, m_systems, m_payloads,m_assembly])
    mass_vector_lbs = mass_vector_kg * 2.20462262
    cost_per_lbs = np.array([[609,1614,679,107,248,315,405,58],[204,484,190,98,91,91,100,4],[88,233,98,16,36,46,59,3]])
    TFU_cost = np.dot(cost_per_lbs,mass_vector_lbs)
    labor  = np.array([])
    materials = []
    other = []
    for n in range(1,n_aircraft+1):
        labor = np.append(labor, n**(np.log(0.85)/np.log(2)))
        materials = np.append(materials, n**(np.log(0.95)/np.log(2)))
        other = np.append(other, n**(np.log(0.95)/np.log(2)))

    costs = (TFU_cost[0]*labor,TFU_cost[1]*materials,TFU_cost[2]*other)

    c_total_recurring_2002 = np.sum(costs)
    c_total_recurring_2019 = c_total_recurring_2002*1.44

    return round(c_total_recurring_2019/1000000,2)

def total_cost(n_aircraft,m_wing,m_empennage,m_fuselage,m_gear,m_engines, m_systems, m_payloads,m_assembly):
    "Return the total cost per aircraft in in million USD(April 2019) based on aircraft mass in kg"
    total_cost = non_recurring_cost(m_wing,m_empennage,m_fuselage,m_gear,m_engines, m_systems, m_payloads) + recurring_cost(n_aircraft,m_wing,m_empennage,m_fuselage,m_gear,m_engines, m_systems, m_payloads,m_assembly)
    PU_cost = total_cost/n_aircraft
    return round(PU_cost,2)

