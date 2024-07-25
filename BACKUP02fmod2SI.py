import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import et_lib as et
import glob,os

def plot_curve_with_arrows2(x, y, axis, interval, color='b',HW =0.1, arrow_scale=1.0, arrow_color='red'):
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
    axis.plot(x, y,color=color)

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
    #axis.quiver(xq,yq,uq,vq, scale_units='xy', scale=1, width=0.01)
    axis.quiver(xq,yq,uq,vq, angles='xy',color=color, scale=0.1)

    return


# unit constants
sec_per_min = 60
kPa_per_Pa = 0.001
Pa_per_kPa = 1.0/kPa_per_Pa
min_per_sec = 1/sec_per_min
Gal_per_Liter = 0.2642
Liter_per_Gal = 3.7854
Liter_per_m3  = 1000.0
Liter_per_mm3 = Liter_per_m3 / 1000**3
Gal_per_mm3 = Liter_per_mm3 *Gal_per_Liter
mm3_per_Gal = 1.0/Gal_per_mm3
MM3perLiter = 1.0 / Liter_per_mm3
# Ideal Gas Law  https://pressbooks.uiowa.edu/clonedbook/chapter/the-ideal-gas-law/
m3_per_mole = 0.02241 # m3/mol of Air
moles_per_m3 = 1.0/m3_per_mole
Pa_per_PSI  = 6894.76
atmos_Pa = 14.5 * Pa_per_PSI
Patmosphere = 101325.0    # Pascals

m3_per_Liter =  1.0 / Liter_per_m3  # m3

dataDirName = 'dataAndyMay24'

pd = {}  # let's keep track of our params and their units
pu = {}
def print_param_table(pd,pu):
    print('\n Parameter Table ')
    for k in pd.keys():
        print(f'{k:12}  {pd[k]:8.4E}  {pu[k]:15}')
    print('')

def print_param_table2(pd,pdo,pu):  # highlight changes in params due to hacks
    print('\n Parameter Table ')
    for k in pd.keys():
        f = ' '
        if pd[k] != pdo[k]:
            f = '*'
        print(f'{k:12} {f} {pd[k]:8.4E}  {pu[k]:15}')
    print('')

#
#   What to plot
#
PLOT_TYPE = 'SIMULATION' # only simulation
PLOT_TYPE = 'OVERLAY'   # includes experimental data
FPLOT = True             # make a force plot as well

pd['moles_per_m3'] = moles_per_m3
pu['moles_per_m3'] = 'moles / m3 of Air'

#  Data
Psource_SIu = Patmosphere + 3.0 * Pa_per_PSI # pascals

LLine_meas = 0.6   # L/ min / kPa  (page 5 col 2)
LLine_SIu = LLine_meas * m3_per_Liter /(sec_per_min * Pa_per_kPa)  # m3 / sec / Pa

Fopen = LLine_SIu * Psource_SIu

Rsource_SIu = Psource_SIu / Fopen

pd['Patmosphere'] = Patmosphere
pu['Patmosphere'] = 'Pascals'
pd['Psource_SIu'] = Psource_SIu
pu['Psource_SIu'] = 'Pascals'
pd['LLine_SIu']  = LLine_SIu
pu['LLine_SIu']  = 'm3/sec / Pascal'
pd['Rsource_SIu'] = Rsource_SIu
pu['Rsource_SIu'] = 'Pa/m3/sec'


# Reel Params
J = 1.0E-04 * [4.67, 5.10, 5.64][1]  # kg/m2: Reel rot. Inertia  Table I
rReel_SIu = 10.9 / 1000  # meters
#  Reel friction (Coulomb type)
Tau_coulomb = [0.0029, 0.0174, 0.0694][1] # N m :   Table II

pd['J'] = J
pu['J'] = 'kg/m2'
#Vhousing_m3 = 0.5 * m3_per_Liter  # 0.5L in m3

r_hous = 0.075 / 2 # 75mm diam.
px_p_m = (441-322)/(5*2.54/100.0)
l_hous = (328-18)/px_p_m       # m - measured off photo
print(f'Housing r: {r_hous:5.3f} Lh: {l_hous:5.3f}')
Vhousing_m3 = np.pi * r_hous**2 * l_hous

pd['Vhousing_m3'] = Vhousing_m3
pu['Vhousing_m3'] = 'm3'


#Tube
Diam_MM = 25
Radius_MM = Diam_MM/2
Radius_m = Radius_MM/1000.0
#area_mm2 = np.pi * Radius_MM**2
area_m2  = np.pi * Radius_m**2

Kdrag = 8.0 # N / m / sec2  (pure speculation!!)
K2drag =  0    # N m
pd['Kdrag'] = Kdrag           # velocity coeff
pu['Kdrag'] = 'N / m2 / sec'
pd['K2drag'] = K2drag         # inverse Length coeff
pu['K2drag'] = 'Nm'

def Fdrag(L,Ldot):
    v = 2*Ldot # eversion kinematics
    return Kdrag * v + K2drag / max(L, 0.100)


pd['area_m2'] = area_m2
pu['area_m2'] = 'm2'

#
#  FLOW load line for source
#

# Presure intercept:
Pintercept = Psource_SIu
pd['Pintercept'] = Pintercept
pu['Pintercept'] = 'Pascals'

# Velocity intercept:
FlowMax_SIu = (Psource_SIu - Patmosphere)/Rsource_SIu
Fintercept = FlowMax_SIu
pd['Fintercept'] =  FlowMax_SIu
pu['Fintercept'] = 'm3/sec'



# source VELOCITY load line
Vintercept = FlowMax_SIu / area_m2
pd['Vintercept'] =  Vintercept
pu['Vintercept'] = 'm/sec'


x1 = 0
y1 = Pintercept
x2 = Fintercept
y2 = Patmosphere
#
#def P2V_SourceLL_SIu(P):  # get evert Vel from Pressure via Source LL
    #dydxLL = (y2-y1) / (x2-x1)
    #return dydxLL

print('Pressure and velocity intercept points (SI units)')
print(Pintercept, 'Pascals')
print(Vintercept, 'meters/sec')


#
#    Eversion Thresholds
#

PBA_static = Patmosphere + 2.0 * Pa_per_PSI   #  Pascals  Static Breakaway Pressure
PHalt_dyn  = Patmosphere + 1.25 * Pa_per_PSI  #  Pascals   dynamic stopping pressure
P1 = PHalt_dyn
P2 = PBA_static


#  TESTING: taper the thresholds together depending time
dP1dL = 0.5*(P2-P1)/1.1   # pa/m
dP2dL = -dP1dL

dF1dL = 0.5*(26-18)/1.1
dF2dL = -dF1dL


pd['Threshold Taper *'] = dP1dL
pu['Threshold Taper *'] = 'Pa /m'
pd['PBA_static'] = PBA_static
pu['PBA_static'] = 'Pascals'
pd['PHalt_dyn'] = PHalt_dyn
pu['PHalt_dyn'] = 'Pascals'

# States
PRESSURE_TEST = 2
GROWING = 1
STUCK = 0

# Constant Eversion forces
f_Brake_SIu = 1.0  # Newtons (positive opposes motion)
fEverForce_SIu = 1.0 # Newtons (could be prop to L!)


# PV = NRT !  # need SI units here!
# https://en.wikipedia.org/wiki/Ideal_gas_law
#  P  Pascals
#  V  m3
#  N  moles of Air
#  T  degK
#  R = 8.314 (ideal gas constant)

# Physical Constants
R = 8.314 # wikipedia Gas Constant
T= 295.4  # 72F in deg K

print_param_table(pd,pu)  # before hacks
pd_orig = pd.copy()

##########################################################
#
#  Parameter hacks needed to match results (Fig 3)
#
#   TESTING  HACK

# these two speed up to about 2 bursts per sec.
#Rsource_SIu *= 0.65
# Andy: flow(m3/sec) = -2E-8 P (Pa) + 0.0005
#    ->    R = 1/2E-8 =
#Rsource_SIu = 1.0 / (1.0*10**(-8))
Rsource_SIu *= 2.0
pd['Rsource_SIu'] = Rsource_SIu

Vhousing_m3 *= 0.5
pd['Vhousing_m3'] = Vhousing_m3

#Kdrag *= 2
#pd['Kdrag'] = Kdrag
#pd['Vhousing_m3'] = Vhousing_m3

J *=4
pd['J'] = J
#########
#
##  Initial Conditions
#

P = Patmosphere  # 1 atmosphere in Pa
#

# Starting STATE
state = STUCK
#state = PRESSURE_TEST

# state var initial values

N = P*Vhousing_m3/(R*T)    #   N = PV/(RT), ideal gas law
L = 0       # ET length (m)
Ldot  = 0
Lddot = 0

# declare data storage
p = []      # Pressure (Pa)
vol = []    # Volume (m3)
f = []      # SourceFlow (m3/sec)
fet = []    # Flow out due to everting
l = []      # L (m)
ldot = []   # velocity (m/sec)

F_e = []   # eversion force
F_c = []   # coulomb force
F_d = []   # drag force
F_j = []   # reel inertia force

#######################################################3
#
# model the burst and stick behavior
#
##   Simulation Parameters
dt = 0.005  # sec

TMIN = 0.0
TMAX = 8.0

# Time range for plots
PltTMIN = TMIN
PltTMAX = TMAX

#PltTMIN = 0
#PltTMAX = 0.5



time = np.arange(TMIN,TMAX,dt)

for t in time:

    # Eq 1
    Voltot = Vhousing_m3 + L*area_m2

    # Eq 2
    P = N*R*T / Voltot

    # Eq 5.1
    Fcoulomb = Tau_coulomb / rReel_SIu
    F_c.append(Fcoulomb)

    # Eq 5.2
    #Fdrag = L * Kdrag * 2.0 * Ldot  # material speed is 2x Ldot
    F_d.append(Fdrag(L,Ldot))
    F_j.append(Lddot*J/rReel_SIu**2)

    # Eq 6


    F_ever = 0.5 * P*area_m2 -Fdrag(L,Ldot) -Fcoulomb
    F_e.append(F_ever)   # eversion force
    if state == GROWING:
        #  0.5 due to eversion kinematics
        Lddot =  (rReel_SIu**2/J) * max(0, F_ever )
    else: # smooth slow down
        alpha = 2000 * dt  # empirical fit
        Lddot = -1 * max(0, alpha * Ldot)
        #Ldot  = 0


    # Eq 2.5
    sourceFlowIn = (Psource_SIu - P)/Rsource_SIu

    # Eq 3
    #ETGrowthFlow = Ldot*area_m2    # vol is growing so N is growing
    Ndot =+ (sourceFlowIn) * moles_per_m3

    # Eq 3.5
    #  thresholds seem to grow closer together with eversion
    P1 = PHalt_dyn  + dP1dL * L
    P2 = PBA_static + dP2dL * L


    F1 = 18  + dF1dL * L * 0.5 # Newtons
    F2 = 26  + dF2dL * L * 0.5

    PRESSURE_BREAK = False

    # Eq 4
    if PRESSURE_BREAK:
        if P > P2:
            state = GROWING
        if P < P1:
            state = STUCK
    else:
        if F_ever > F2:
            state = GROWING
        if F_ever < F1:
            state = STUCK


    # Eq 7
    Ldot += Lddot * dt
    L    += Ldot  * dt
    N    += Ndot  * dt



    # Record data
    # tube length
    l.append(L)
    # tube eversion velocity
    ldot.append(Ldot)
    # flow from source
    f.append(sourceFlowIn)
    fet.append(Ldot * area_m2)  #flow out into the tube
    # Pressure
    p.append(P) # Pa
    # Volume
    vol.append(Voltot)


if FPLOT:
    # Create force plot
    fig = plt.figure()
    ax = fig.gca()
    plt.plot(time, F_e)
    plt.plot(time, F_c)
    plt.plot(time, F_d)
    plt.plot(time, F_j)
    plt.legend(['F_Eversion','F_Coulomb','F_Drag','F_Inertia'])



# Create simulation output figure with subplots
fig, axs = plt.subplots(3, 2, figsize=(8, 10))

# Plot 1
#axs[0].plot(vel, p)
vel1 = []
pr1 = []



if state==PRESSURE_TEST:
    axs[0].text(0.5, 10000, 'PRESSURE_TEST (no eversion)')

#############################################################################
#
#   Plot the simulation outputs
#

clrs = ['b','g','r','c','m']
#trajectory:
# "ideal" loadline
axs[0,0].plot([x1,x2], [y1,y2], color='k')  # x1,..y2 defined above
#axs[0,0].plot(fet,p) #pressure is Pa
plot_curve_with_arrows2(fet,p, axs[0,0], 50,color=clrs[0])

axs[0,0].legend(['Source Load Line',  'Trajectory'])
axs[0,0].set_xlabel('Source Flow (m3/sec)')
axs[0,0].set_ylabel('Pressure (Pa)')
axs[0,0].set_xlim(0,0.001)
#plt.sca(axs[0,0])
#plt.xticks([0.0001, 0.0002, 0.0003])
axs[0,0].set_ylim(Patmosphere, Psource_SIu)

# Plot 2   # PRESSURE
axs[1,0].plot(time, p)
axs[1,0].legend(['Pressure (Pa)' ])
axs[1,0].set_xlabel('Time (sec)')
axs[1,0].set_ylabel('Pressure (Pa)')
axs[1,0].set_xlim(PltTMIN, PltTMAX)
axs[1,0].set_ylim(Patmosphere, Psource_SIu)
axs[1,0].plot([0, PltTMAX] , [PHalt_dyn, P1] ,'g', linestyle='dashed',color=clrs[3])
axs[1,0].plot([0, PltTMAX] , [PBA_static, P2] ,'r', linestyle='dashed',color=clrs[4] )


# Plot 3
axs[2,0].plot(time, vol )
axs[2,0].legend(['Vol (m3)'])
axs[2,0].set_xlabel('Time (sec)')
axs[2,0].set_ylabel('Volume (m3)')
axs[2,0].set_xlim(PltTMIN, PltTMAX)
axs[2,0].set_ylim( 0.000, 0.0015 )

axs[0,1].plot(time, f)
axs[0,1].set_xlabel('Time (sec)')
axs[0,1].set_ylabel('Flow (m3/sec)')
axs[0,1].set_xlim(PltTMIN, PltTMAX)
axs[0,1].set_ylim(      0, 0.00020 )

axs[1,1].plot(time, l)
axs[1,1].set_xlabel('Time (sec)')
axs[1,1].set_ylabel('Tube Len (m)')
axs[1,1].set_xlim(PltTMIN, PltTMAX)
axs[1,1].set_ylim(      0, 0.600 )

axs[2,1].plot(time, ldot)
axs[2,1].set_xlabel('Time (sec)')
axs[2,1].set_ylabel('Tube Velocity (m/sec)')
axs[2,1].set_xlim(PltTMIN, PltTMAX)
axs[2,1].set_ylim(      0, 0.5 )



if PLOT_TYPE == 'OVERLAY':
    #open file

    dataDirNames = ['dataAndyMay24']


    #files = ['eversion_flow-hi-inr_hi-fric_tube-1_trial-2.csv', 'eversion_flow-hi-inr_lo-fric_tube-1_trial-1.csv', 'eversion_flow-hi-inr_lo-fric_tube-1_trial-2.csv', 'eversion_flow-hi-inr_lo-fric_tube-1_trial-3.csv', 'eversion_flow-hi-inr_hi-fric_tube-2_trial-1.csv', 'eversion_flow-hi-inr_hi-fric_tube-2_trial-2.csv', 'eversion_flow-hi-inr_hi-fric_tube-3_trial-1.csv', 'eversion_flow-hi-inr_hi-fric_tube-3_trial-2.csv', 'eversion_flow-hi-inr_hi-fric_tube-3_trial-3.csv',
    #]

    files = []
    mdfiles = []
    for datadir in dataDirNames:
            tfiles = list(glob.glob(datadir + '/' + "*"))
            tfiles.sort(key=lambda x: os.path.getmtime(x),reverse=True) # newest first
            if len(tfiles) ==0:
                cto.error('No brl_data files found')
            dirfiles = []
            dirmdfiles = []
            for f in tfiles:
                #print('found: ',f)
                if '.csv' in f :
                    dirfiles.append(f)
                elif  '_meta.json' in f:
                    dirmdfiles.append(f)

            # dict.fromkeys better than set because preserves order (thanks google)
            #filenameroots = list(dict.fromkeys(filenameroots)) # elim dupes
            files += dirfiles
            mdfiles += dirmdfiles

    print('Choose a dataset to plot with the simulation: ')
    for i,fn in enumerate(files):
        print(f'{i:5}', fn)

    n = int(input('Select file by number: '))
    if (n>=len(files) or len(files) <0):
        print('illegal file name index')
        quit()
    fn = files[n]

    fig.suptitle(fn.split('/')[-1])


    #fn = dataDirName + '/' + files[n]
    fn = files[n]

    #print('opening: ',fn)
    #x=input('       ... OK?? <cr>')

    ed = et.get_data_from_AL_csv(fn)
    ed = et.convert_units(ed)


    #ed['time'] = data[:,0]
    #ed['flow'] = data[:,idxflow]
    #ed['P'] = data[:,idxpress]
    #ed['L'] = data[:,idxL]

    # HACK
    ed['flow'] *= 10.0

    # phase plot with arrows






    HW= 0.00001
    HL= 0.0001
    Interval = 25

    x = ed['flow']
    y = ed['P']

    if False:
        ######################################3
        # TESTING
        #  a circle
        th = np.linspace(0, 2*np.pi, 100)
        x  = np.cos(th)
        y = np.sin(th)
        # an elipse:
        xsc = 0.03
        x *= xsc
        #############################################


    Plt_ranges = [0, 0.001, Patmosphere, Psource_SIu]

    #plt.figure()
    #ax = plt.gca()
    ax = axs[0,0]
    #plot_curve_with_arrows(xdata, ydata, Interval, ax, arrow_scale=1.0 )
    plot_curve_with_arrows2(x, y, ax, Interval,color=clrs[1])

    axs[1,0].plot(np.array(ed['time']), np.array(ed['P']), '--', color=clrs[1])  # Experimental Data
    #axs[2,0].plot(ed['time'], ed['vol'], '--',color=clrs[1])  # Experimental Data
    axs[0,1].plot(ed['time'], ed['flow'], '--',color=clrs[1])  # Experimental Data
    axs[1,1].plot(ed['time'], ed['L'], '--',color=clrs[1])  # Experimental Data
    axs[2,1].plot(ed['time'], ed['Ldot'], '--',color=clrs[1])  # Experimental Data


# Adjust layout
plt.tight_layout()

print_param_table2(pd,pd_orig, pu)

plt.show()

