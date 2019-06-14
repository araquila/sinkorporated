import matplotlib.pyplot as plt
import numpy as np

x = np.arange(0, 11, 1)
y = np.arange(0, 11, 1)
z = np.arange(0, 22, 2)

#####################################################################
# READ ME                                                           #
#                                                                   #
# You should not make the plot in-line, but in a separate sceen.    #
#                                                                   #                                                         
# If you haven't already:                                           #
# 1. Ctrl + Alt + Shift + P                                         #
# 2. Go to: IPython console                                         #
# 3. Go to: Graphics                                                #
# 4. Put "Graphics Backend" to: Automatic                           #
# 5. Restart Spyder to make this work                               #
#                                                                   #
# Copy the lines below and put them in your own script. Do NOT      #
# full-screen the plot. Simply keep the window small and press      #
# the save-button. This will keep the font-size the same.           #
#####################################################################

# Plot lines
plt.plot(x, y, label="y-values")
plt.plot(x, z, label="z-values")

# Plot Points
plt.scatter(x, y, label="y-values")
plt.scatter(x, z, label="z-values")

# Label of the axes
plt.xlabel("x-values", size="large")
plt.ylabel("y-values and z-values", size="large")

# Limits of the axes
plt.xlim([0,10])
plt.ylim([0,20])

# Create Legend
plt.legend(loc="best", fontsize="large")

# Show plot
plt.show()