# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 16:14:53 2019

@author: benja
"""

import pandas as pd
import matplotlib.pyplot as plt
data = pd.read_excel(r'C:\Users\benja\Documents\TUD\DSE\sinkorporated\Noise_Calculations.xlsx', sheet_name='Approach') 
df = pd.DataFrame(data, columns= ['Frequency', 'SPL_wing', 'SPL_hor_tail', 'SPL_ver_tail', 'SPL_flaps', 'SPL_nose', 'SPL_main', 'SPL_strut','Total','ATR' ])

### APPROACH ###

x = df['Frequency'].values.tolist()
#x.pop(0)

y1 = df['SPL_wing'].values.tolist()
#y1.pop(0)

y2 = df['SPL_hor_tail'].values.tolist()
#y2.pop(0)

#y3 = df['SPL_ver_tail'].values.tolist()
#y3.pop(0)

y4 = df['SPL_flaps'].values.tolist()
#y4.pop(0)

y5 = df['SPL_nose'].values.tolist()
#y5.pop(0)

y6 = df['SPL_main'].values.tolist()
#y6.pop(0)

y7 = df['SPL_strut'].values.tolist()
#y7.pop(0)

y8 = df['Total'].values.tolist()
#y8.pop(0)

atrapp = df['ATR'].values.tolist()

plt.figure(1)
plt.ylim([0,100])
plt.xlim([10,10000])
plt.semilogx(x,y1, label="SPL wing")
plt.semilogx(x,y2, label="SPL hor tail")
#plt.plot(x,y3)
plt.semilogx(x,y4, label="SPL flaps")
plt.semilogx(x,y5, label="SPL nose")
plt.semilogx(x,y6, label="SPL main landing gear")
plt.semilogx(x,y7, label="SPL strut")
plt.semilogx(x,y8,marker="*", label="Total SPL")

plt.semilogx(x,atrapp, marker="o", label="Total ATR 72")


plt.legend()

plt.ylabel("SPL [dB]")
plt.xlabel("Frequency [Hz]")

#plt.show()

### FLY OVER ###
data2 = pd.read_excel(r'C:\Users\benja\Documents\TUD\DSE\sinkorporated\Noise_Calculations.xlsx', sheet_name='Fly-over') 
df2 = pd.DataFrame(data2, columns= ['Frequency', 'SPL_wing', 'SPL_hor_tail', 'SPL_ver_tail', 'SPL_flaps', 'SPL_nose', 'SPL_main', 'SPL_strut','Total', 'ATR' ])

x2 = df2['Frequency'].values.tolist()
#x2.pop(0)

y12 = df2['SPL_wing'].values.tolist()
#y12.pop(0)

y22 = df2['SPL_hor_tail'].values.tolist()
#y22.pop(0)

y32 = df2['SPL_ver_tail'].values.tolist()
#y32.pop(0)

y42 = df2['SPL_flaps'].values.tolist()
#y42.pop(0)

y52 = df2['SPL_nose'].values.tolist()
#y52.pop(0)

y62 = df2['SPL_main'].values.tolist()
#y62.pop(0)

y72 = df2['SPL_strut'].values.tolist()
#y72.pop(0)

y82 = df2['Total'].values.tolist()
#y82.pop(0)

atrfo = df2['ATR'].values.tolist()

plt.figure(2)
plt.ylim([0,80])
plt.xlim([10,10000])
plt.semilogx(x2,y12, label="SPL wing")
plt.semilogx(x2,y22, label="SPL hor tail")
#plt.plot(x,y3)
plt.semilogx(x2,y42, label="SPL flaps")
plt.semilogx(x2,y52, label="SPL nose")
plt.semilogx(x2,y62, label="SPL main landing gear")
plt.semilogx(x2,y72, label="SPL strut")
plt.semilogx(x2,y82,marker="*", label="Total SPL")

plt.semilogx(x2,atrfo, marker="o", label="Total ATR 72")

plt.legend()

plt.ylabel("SPL [dB]")
plt.xlabel("Frequency [Hz]")

### PROPELLER ###
data3 = pd.read_excel(r'C:\Users\benja\Documents\TUD\DSE\sinkorporated\Noise_Calculations.xlsx', sheet_name='Propeller') 
df3 = pd.DataFrame(data3, columns= ['Frequency', 'SPL', 'ATR'])

x3 = df3['Frequency'].values.tolist()
#x2.pop(0)

y13 = df3['SPL'].values.tolist()

y23 = df3['ATR'].values.tolist()

plt.figure(3)
#plt.ylim([0,80])
#plt.xlim([10,10000])
plt.semilogx(x3,y13, label="Aquila")
plt.semilogx(x3,y23, label = "ATR 72")


plt.legend()
#
plt.ylabel("SPL propeller [dB]")
plt.xlabel("Frequency [Hz]")