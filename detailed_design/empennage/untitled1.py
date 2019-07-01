# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 14:16:19 2019

@author: benja
"""

import matplotlib.pyplot as plt

t = ['12-03-18','12-04-18' , '12-05-18', '12-06-18', '12-07-18']
s = [2,5,2,3,2]

plt.plot(t,s, marker='o')
for a,b in zip(t, s): 
    plt.text(a, b, str(b))

plt.title("Line Chart")
plt.xlabel("date")
plt.ylabel("count")
plt.show()