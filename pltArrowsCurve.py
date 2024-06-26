import numpy as np
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle

def plot_curve_with_arrows(x, y, axis, interval=10, arrow_scale=1.0, arrow_color='red'):
    """
    Plot a curve with arrows indicating the direction of the curve.

    Parameters:
        x (array-like): x-coordinates of the data points.
        y (array-like): y-coordinates of the data points.
        axis (matplotlib axis): The axis on which to plot the curve and arrows.
        interval (int): Spacing between arrows along the curve.
        arrow_scale (float): Scale factor for the arrow lengths.
        arrow_color (str): Color of the arrows.
    """

    # scale factors for x and y
    #xsf = 0.03* arrow_scale * np.abs(np.max(x)-np.min(x))
    #ysf = 0.03* arrow_scale * np.abs(np.max(y)-np.min(y))


    # Compute the direction of the curve (tangent vectors)
    dx = np.gradient(x)
    dy = np.gradient(y)


    #axis.plot(x,dx)
    dirs = np.arctan2(dy, dx)

    # Normalize the tangent vectors to get unit vectors
    dx_unit = np.cos(dirs)
    dy_unit = np.sin(dirs)

    # plot orig curve
    #axis.plot(x,y,marker=MarkerStyle(">", "full"))
    axis.plot(x,y)

    xsf = .1
    ysf = .1

    ## Plot arrows along the curve at specific intervals
    for i in range(0, len(x), interval):
        axis.arrow(x[i], y[i], dx_unit[i] * xsf, dy_unit[i] * ysf,
                   head_width=0.1, head_length=0.1, fc=arrow_color, ec=arrow_color)


def plot_curve_with_arrows2(x, y, axis, interval, HW =0.1, arrow_scale=1.0, arrow_color='red'):
    """
    Plot a curve with arrows indicating the direction of the curve.

    Parameters:
        y (array-like): y-coordinates of the data points.
        axis (matplotlib axis): The axis on which to plot the curve and arrows.
        interval (int): Spacing between arrows along the curve.
        arrow_scale (float): Scale factor for the arrow lengths.
        arrow_color (str): Color of the arrows.
    """
    # scale factors for x and y
    #if max(ranges) == 0:
    xsf = 0.03* arrow_scale / np.abs(np.max(x)-np.min(x))
    ysf = 0.03* arrow_scale / np.abs(np.max(y)-np.min(y))
    #else:
        #xsf = 1/np.abs(ranges[1]-ranges[0])
        #ysf = 1/np.abs(ranges[3]-ranges[2])

    # Compute the direction of the curve (tangent vectors)
    dx = np.gradient(x)
    dy = np.gradient(y)
    th = np.arctan2(dy, dx)

    #  get unit vectors
    dx_unit = np.cos(th)
    dy_unit = np.sin(th)

    # Plot the curve
    axis.plot(x, y)

    xsf = .01
    ysf = .01

    xq = []
    yq = []
    uq = []
    vq = []
    ## Plot arrows along the curve at specific intervals
    for i in range(0, len(x), interval):

        #print('plotting arrow: ', x[i], y[i], dx_unit[i]*xsf, dy_unit[i]*xsf)
        #axis.arrow(x[i], y[i], dx_unit[i] * xsf, dy_unit[i] * ysf)
                #head_width=0.1, head_length=0.1, fc=arrow_color, ec=arrow_color)
        dx = dx_unit[i]
        dy = dy_unit[i]
        #start = (x[i],y[i])
        #end   = (x[i]+dx*xsf, y[i]+dy*ysf)
        #end   = (x[i]+dx_unit[i], y[i]+dy_unit[i])

        w = 1

        xsf=0.002
        ysf=0.002

        # fill the arrow quiver starts
        xq.append(x[i])
        yq.append(y[i])
        uq.append(dx * xsf)
        vq.append(dy * ysf)

        #axis.arrow(x[i], y[i], dx_unit[i] * xsf, dy_unit[i] * ysf,
                #head_width=HW, head_length=20*HW, fc=arrow_color, ec=arrow_color)
                #fc=arrow_color, ec=arrow_color)

        #arrow = mpatches.FancyArrow(x[i],y[i],dx_unit[i]*xsf,dy_unit[i]*ysf,  length_includes_head=True)

        #arrow = mpatches.FancyArrowPatch(start,end, mutation_scale=10000, transform=axis.transAxes)
#arrow(x_tail + 1, y_tail - .4, dx, dy, width=.1, length_includes_head=True, color="C2")

        #arrow = mpatches.Arrow(start[0],start[1], dx,dy , width=0.1, color='blue')

        #axis.add_patch(arrow)

    #axis.quiver(xq,yq,uq,vq, scale_units='xy', scale=1, width=0.01)
    axis.quiver(xq,yq,uq,vq, angles='xy',)

    return


# Example usage:
NC = 3
x = np.linspace(0, NC*2*np.pi, 100*NC)
y = np.sin(x)

#  a circle
th = np.linspace(0, 2*np.pi, 100)
x  = np.cos(th)
y = np.sin(th)

# an elipse:
xsc = 5
x *= xsc

Interval = int(15)

fig, ax = plt.subplots()
#plot_curve_with_arrows(x, y, ax, arrow_scale=1.0)
plot_curve_with_arrows2(x, y, ax, Interval, arrow_scale=1.0)

plt.show()
