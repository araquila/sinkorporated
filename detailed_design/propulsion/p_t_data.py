import numpy as np
p_t = [[105,56.4],
        [110,88.2],
        [111.7,101.3],
        [115,132.3],
        [120,191.6],
        [125,269.0],
        [130,367.6],
        [135,490.7],
        [140,641.6],
        [145,823.7],
        [150,1040.5],
        [155,1295.6],
        [160,1592.8],
        [165,1935.9],
        [170,2329.3],
        [175,2777.6],
        [180,3286.4],
        [185,3863.2],
        [190,4520.5]]

p_t_gradient = np.zeros(len(p_t)-1)
for i in range(len(p_t_gradient)):
    p_t_gradient[i] = (p_t[i+1][1]-p_t[i][1])/(p_t[i+1][0]-p_t[i][0])
    # print(i,p_t_gradient[i])

def det_T_vaporise(pressure):
    for item in p_t:
        if pressure <= item[1]:
            return ans
        if pressure > item[1]:
            ans = item[0] + p_t_gradient[p_t.index(item)]**-1 * (pressure - item[1])
