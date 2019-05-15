import numpy as np

def non_recurring_cost(m_wing,m_empennage,m_fuselage,m_gear,m_engines, m_systems, m_payloads):
    """Returns the total non-recurring cost in million USD(April 2019) based on aircraft mass in kg"""
    mass_vector_kg = np.array([m_wing,m_empennage,m_fuselage,m_gear,m_engines, m_systems, m_payloads])
    mass_vector_lbs = mass_vector_kg * 2.20462262
    cost_per_lbs_vector = np.array([17731,52156,32093,2499,8691,34307,10763])
    cost_vector = mass_vector_lbs*cost_per_lbs_vector

    c_total_nonrecurring_2002 = np.sum(cost_vector)
    c_total_nonrecurring_2019 = c_total_nonrecurring_2002*1.44

    return c_total_nonrecurring_2019/1000000


print(non_recurring_cost(2045,266,2718,1198,1200,2624,6120))
