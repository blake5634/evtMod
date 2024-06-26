import numpy as np
import matplotlib.pyplot as plt

# Generate sample data (you can replace this with your own data)
NC = 4
x = np.linspace(0, NC*2*np.pi, 1000)
y = np.sin(x)

# Compute the direction of the curve (tangent vectors)
dx = np.gradient(x)
dy = np.gradient(y)
directions = np.arctan2(dy, dx)
m = np.sqrt(dx**2+dy**2)
# Normalize the tangent vectors to get unit vectors
#dx_unit = np.cos(directions)
#dy_unit = np.sin(directions)
dx_unit = dx/m
dy_unit = dy/m

# Plot the curve
plt.plot(x, y)

# Plot arrows along the curve at specific intervals
interval = 10  # Adjust this based on how many arrows you want
for i in range(0, len(x), interval):
    scale_factor = 0.5 * np.sqrt((x[-1] - x[0])**2 + (y[-1] - y[0])**2)  # Adjust scale factor as needed
    plt.arrow(x[i], y[i], dx_unit[i], dy_unit[i], head_width=0.1, head_length=0.1, fc='red', ec='red')

plt.show()
