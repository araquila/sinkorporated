import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_excel(r'C:\Users\benja\Documents\TUD\DSE\market_price.xlsx') 
df = pd.DataFrame(data, columns= ['Pax', 'Price', 'Aircraft', 'AQ', 'AQPax', 'AQPrice'])

x = df['Pax'].values.tolist()
y = df['Price'].values.tolist()
n = df['Aircraft'].values.tolist()

a1 = df['AQ'].values.tolist()
a2 = df['AQPax'].values.tolist()
a3 = df['AQPrice'].values.tolist()



z = np.polyfit(x, y, 1)
p = np.poly1d(z)
fig, ax = plt.subplots()
ax.scatter(x, y)

ax.plot(x,p(x), "g")

for i, txt in enumerate(n):
    ax.annotate(txt, (x[i], y[i]))
    

plt.xlabel("Number of passengers [-]", size="large")
plt.ylabel("Unit cost [million USD]", size="large")

ax.scatter(a2, a3)
ax.annotate('AQ60', (a2[0],a3[0]))


plt.plot()