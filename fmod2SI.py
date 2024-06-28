import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import simulation as sim
import et_lib as et
from et_lib import error


# Unit Conversions
uc = {
    "sec_per_min": 60,
    "kPa_per_Pa": 0.001,
    "Pa_per_kPa": 1.0 / 0.001,
    "min_per_sec": 1 / 60,
    "Gal_per_Liter": 0.2642,
    "Liter_per_Gal": 3.7854,
    "Liter_per_m3": 1000.0,
    "Liter_per_mm3": 1000.0 / 1000**3,
    "Gal_per_mm3": (1000.0 / 1000**3) * 0.2642,
    "mm3_per_Gal": 1.0 / ((1000.0 / 1000**3) * 0.2642),
    "MM3perLiter": 1.0 / (1000.0 / 1000**3), # Ideal Gas Law  https://pressbooks.uiowa.edu/clonedbook/chapter/the-ideal-gas-law/
    "m3_per_mole": 0.02241,  # m3/mol of Air
    "moles_per_m3": 1.0 / 0.02241,
    "Pa_per_PSI": 6894.76,
    "atmos_Pa": 14.5 * 6894.76,
    "m3_per_Liter": 1.0 / 1000.0  # m3
}


dataDirName = 'dataAndyMay24'


#
#   What to plot
#
PLOT_TYPE = 'SIMULATION' # only simulation
PLOT_TYPE = 'OVERLAY'   # includes experimental data
FPLOT = True             # make a force plot as well

#
#  Load parameters (argv)
#
pd,pu = et.load_params(['name']) # defaults for now.


print('Pressure and velocity intercept points (SI units)')
print(pd['Pintercept'], 'Pascals')
print(pd['Vintercept'], 'meters/sec')


#
#    Eversion Thresholds
#

PBA_static = pd['Patmosphere'] + 2.0 * uc['Pa_per_PSI']  # Pascals  Static Breakaway Pressure
pd['PBA_static'] = PBA_static
pu['PBA_static'] = 'Pascals'

PHalt_dyn = pd['Patmosphere'] + 1.25 * uc['Pa_per_PSI']  # Pascals   dynamic stopping pressure
pd['PHalt_dyn'] = PHalt_dyn
pu['PHalt_dyn'] = 'Pascals'

#PBA_static = Patmosphere + 2.0  * Pa_per_PSI   #  Pascals  Static Breakaway Pressure
#PHalt_dyn  = Patmosphere + 1.25 * Pa_per_PSI  #  Pascals   dynamic stopping pressure
#P1 = PHalt_dyn
#P2 = PBA_static


#  TESTING: taper the thresholds together depending time
dP1dL = 0.5*(pd['PBA_static'] - pd['PHalt_dyn'])/0.8  # pa/m  (0.8m is fully everted)
dP2dL = -dP1dL
pd['Threshold Taper'] = dP1dL
pu['Threshold Taper'] = 'Pa /m'

# force threshold taper
pd['dF1dL'] =  0.5*(26-18)/1.1
pu['dF1dL'] = 'N/m??' #??


# States
PRESSURE_TEST = 2
GROWING = 1
STUCK = 0

## Constant Eversion forces
#f_Brake_SIu = 1.0  # Newtons (positive opposes motion)
#fEverForce_SIu = 1.0 # Newtons (could be prop to L!)


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

#########################################################################################33

# print and save orig baseline params
et.print_param_table(pd,pu)  # before hacks and fn changes
pd_orig = pd.copy()


if PLOT_TYPE == 'OVERLAY':
    files, mdfiles = et.get_files()

    print('Choose a dataset to plot with the simulation: ')
    for i,fn in enumerate(files):
        print(f'{i:5}', fn)

    n = int(input('Select file by number: '))
    if (n>=len(files) or len(files) <0):
        print('illegal file name index')
        quit()
    fn = files[n]

# get inertia and friction from data file name
Ji, Fric_i = et.get_inr_fric_from_name(fn)
J = 1.0E-4*[ 4.67, 5.10, 5.64][Ji]
Tau_coulomb = [0.0029, 0.0174, 0.0694][Fric_i]

pd['J'] = J
pd['Tau_coulomb'] = Tau_coulomb
pu['Tau_coulomb'] = 'Nm'


##########################################################
#
#  Parameter hacks needed to match results (Fig 3)
#
#   TESTING  HACK


pd['Rsource_SIu'] *= 1.25

pd['Vhousing_m3'] *= 0.5

#Kdrag *= 2
#pd['Kdrag'] = Kdrag

pd['J'] *= 4

##################################################  Run Simulation
t1 = 0.0
t2 = 8.0
state = STUCK
(time,l,lc,ldot,f, fet, p, pbat, pstt, vol, F_e,F_c,F_d,F_j) = sim.simulate(pd,uc,t1,t2)

###################################################

FPLOT = False
PltTMIN = t1
PltTMAX = t2

if FPLOT:
    # Create force plot
    fig = plt.figure()
    fig.suptitle(fn.split('/')[-1])

    ax = fig.gca()
    plt.plot(time, F_e)
    plt.plot(time, F_c)
    plt.plot(time, F_d)
    plt.plot(time, F_j)
    plt.plot(time,lc)
    plt.legend(['F_Eversion','F_Coulomb','F_Drag','F_Inertia','Crumple'])
    ax = plt.gca()
    ax.set_xlim(PltTMIN, PltTMAX)
    ax.set_ylim(-100,300)



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

x1 = 0.0
y1 = pd['Pintercept']
x2 = pd['Fintercept']
y2 = 0.0


axs[0,0].plot([x1,x2], [y1,y2], color='k')  # x1,..y2 defined above
#axs[0,0].plot(fet,p) #pressure is Pa
#
# fet = flow into everting tube (not measured in reality)
# f  = flow from source (experimental data e.g.)
#plot_curve_with_arrows2(fet, p, axs[0,0], 50,color=clrs[0])
et.plot_curve_with_arrows2(f, p, axs[0,0], 500, color=clrs[0])

axs[0,0].legend(['Source Load Line',  'Trajectory'])
axs[0,0].set_xlabel('Source Flow (m3/sec)')
axs[0,0].set_ylabel('Pressure (Pa)')
axs[0,0].set_xlim(0,0.0003)
plt.sca(axs[0,0])
plt.xticks([0.0001, 0.0002, 0.0003])
axs[0,0].set_ylim(pd['Patmosphere'], pd['Psource_SIu'])

# Plot 2   # PRESSURE
axs[1,0].plot(time, p)
axs[1,0].legend(['Pressure (Pa)' ])
axs[1,0].set_xlabel('Time (sec)')
axs[1,0].set_ylabel('Pressure (Pa)')
axs[1,0].set_xlim(PltTMIN, PltTMAX)
axs[1,0].set_ylim(pd['Patmosphere'], pd['Psource_SIu'])
#  plot the eversion thresholds (function of L)
axs[1,0].plot(time, pstt, linestyle='dashed', color=clrs[3])
axs[1,0].plot(time, pbat, linestyle='dashed', color=clrs[4])

# Plot 3
axs[2,0].plot(time, vol )
axs[2,0].legend(['Vol (m3)'])
axs[2,0].set_xlabel('Time (sec)')
axs[2,0].set_ylabel('Volume (m3)')
axs[2,0].set_xlim(PltTMIN, PltTMAX)
axs[2,0].set_ylim( 0.000, 0.0015 )

axs[0,1].plot(time, f)
axs[0,1].set_xlabel('Time (sec)')
axs[0,1].set_ylabel('Source Flow (m3/sec)')
axs[0,1].set_xlim(PltTMIN, PltTMAX)
axs[0,1].set_ylim(      0, 0.00020 )

axs[1,1].plot(time, l, time, lc)
axs[1,1].set_xlabel('Time (sec)')
axs[1,1].set_ylabel('Length (m)')
axs[1,1].legend(['Tube Length', 'Crumple Length'])
axs[1,1].set_xlim(PltTMIN, PltTMAX)
axs[1,1].set_ylim(      0, 0.600 )

axs[2,1].plot(time, ldot)
axs[2,1].set_xlabel('Time (sec)')
axs[2,1].set_ylabel('Tube Velocity (m/sec)')
axs[2,1].set_xlim(PltTMIN, PltTMAX)
axs[2,1].set_ylim(      0, 0.5 )



if PLOT_TYPE == 'OVERLAY':

    # fn is chosen above prior to sim

    fig.suptitle(fn.split('/')[-1])


    #fn = dataDirName + '/' + files[n]
    fn = files[n]

    #print('opening: ',fn)
    #x=input('       ... OK?? <cr>')

    ed = et.get_data_from_AL_csv(fn)
    ed = et.convert_units(ed,uc)


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

    if False:   # testing for phase plotting with data
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

    Plt_ranges = [0, 0.001, pd['Patmosphere'], pd['Psource_SIu'] ]

    #plt.figure()
    #ax = plt.gca()
    ax = axs[0,0]
    #plot_curve_with_arrows(xdata, ydata, Interval, ax, arrow_scale=1.0 )
    et.plot_curve_with_arrows2(x, y, ax, Interval,color=clrs[1])

    axs[1,0].plot(np.array(ed['time']), np.array(ed['P']), '--', color=clrs[1])  # Experimental Data

    dtexp = ed['dtexp']
    dtsim = pd['dt']
    print('dt Sim:', dtsim, 'dt Exp:', dtexp)

    #axs[2,0].plot(ed['time'], ed['vol'], '--',color=clrs[1])  # Experimental Data
    axs[0,1].plot(ed['time'], ed['flow'], '--',color=clrs[1])  # Experimental Data
    axs[1,1].plot(ed['time'], ed['L'], '--',color=clrs[1])  # Experimental Data
    axs[2,1].plot(ed['time'], ed['Ldot'], '--',color=clrs[1])  # Experimental Data

    if False:
        print('Pressure Simulation Losses: ')
        print('Time Domain: ',  et.TD_loss(p, ed['P']))
        print('Freq Domain: ',  et.FD_loss2(p, ed['P'], dtsim, dtexp, test=True, ptitle='Pressure') )


        print('Velocity Simulation Losses: ')
        print('Time Domain: ',  et.TD_loss(ldot,  ed['Ldot']))
        print('Freq Domain: ',  et.FD_loss2(ldot, ed['Ldot'],dtsim,dtexp,test=True,ptitle='Velocity') )


    # Adjust layout
    plt.tight_layout()

    et.print_param_table2(pd,pd_orig, pu)

plt.show()

