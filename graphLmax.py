import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import et_lib as et
import glob, os

#
#   graph the F_lim function
#

Psource = 120E3  # Pascals
Patmo   = 1.0132000E+05

r_tube = 0.0125    # m
Lmax = 0.6

Fsteadystate = 0.5 * ( Psource - Patmo) * np.pi*r_tube**2

F=[]

X = np.linspace(0, 1.0, 100)

for L in X:
    F.append( (L/Lmax)**17 * Fsteadystate )  # a very steep increase

fig,ax = plt.subplots(1,1,figsize=(6,2))

ax.plot(X,F)
ax.set_ylim((0, 10))
ax.set_ylabel('Length Limit Force (N)')
ax.set_xlabel('Length (m)')
# draw lines
x1l = [0, 0.9]
y1l = [Fsteadystate, Fsteadystate]
x2l = [Lmax, Lmax]
y2l = [0, 1.5*Fsteadystate]
ax.plot(x1l,y1l, x2l, y2l)

# Adjust layout
fig.tight_layout()
plt.show()
