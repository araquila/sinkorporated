# plotting and stuff

import matplotlib.pyplot as plt
import numpy as np

def plotfiller(ax, xlim, ylim, x_data = 0, data = 0, vline = 0, direction = "right", alpha = 0.5, color = 'green'):
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
        topline = np.linspace(ylim, ylim, len(data), facecolor = color)
        ax.fill_between(x_data, topline, data, alpha = alpha)
        return


# Example

# the data
l = np.linspace(0, 0.8, 200)
x = np.linspace(0, 4000, 200)
vertical = 3000

# initialise the figure
fig, ax1 = plt.subplots(1,1)
xlim = 4000
ylim = 1

# plot lines
ax1.plot(x, l)
ax1.axvline(vertical)

# plot filled parts of the graph
plotfiller(ax1, xlim, ylim, x_data = x, data = l, direction = "down")
plotfiller(ax1, xlim, ylim, vline = vertical, direction = "right")

# plot cosmetics (add some legends/labels/title)
ax1.set_ylim([0, ylim])
ax1.set_xlim([0, xlim])

plt.show()
