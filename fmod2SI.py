import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import simulation as sim
import et_lib as et
import sys
import glob
import re

modeltypes = ['1Comp','2Comp','23-Jul']

args = sys.argv
cl = ' '.join(args)
paramDir = 'evtParams/2Comp/'

if len(args) == 3:  #  >fmod2SI   nn  [1Comp]   nn=fileNo  '1Comp' for non-default one comp model
    if args[2] not in modeltypes:
            et.error(f'fmod2SI: unknown command line arg: \n     >{cl}')
    if args[2] == '1Comp':     # not same as pd['Compartments']
        paramDir = 'evtParams/1Comp/'
    elif args[2] == '2Comp':   # not same as pd['Compartments']
        paramDir = 'evtParams/2Comp/'   # default is 2Compartment model
    elif args[2] == '23-Jul':
        paramDir = 'evtParams/23-Jul-FlowData/'

unitsConvfilename = 'unitConv.txt'
defaultParamName = 'InitialParams.txt'
defaultUnitsName = 'units_'+defaultParamName
defaultUnitsDir = 'evtParams/'

##set up first version of params and units files
## just run these once (then comment out)
##
#et.saveUnitConv(unitsConvfilename, uc)  # liters to cubic meters e.g.
## first param file
#pd,pu = et.setup_params()
#et.saveParams(defaultParamName,pd)
#et.savePUnits(defaultUnitsName,pu)     # same keys as parameter file

# load parameter file(s)

if len(args) == 1:
    paramDir = 'evtParams/'   # this is only generic with no data (keep??) default is 2Compartment model
    paramFileName = defaultParamName
    print('Usage: fmod2SI <parameterfile>')
    print('   loading nominal default parameters from: ' + defaultParamName)
    pd = et.loadParams(paramDir, defaultParamName)
    pu = et.loadPUnits(paramDir, 'units_'+defaultParamName)

else:
    paramFileNoStr = args[1]
    paramFileName = 'Set'+paramFileNoStr+'Params.txt'  # look up SetxxParams.txt
    if not bool(re.fullmatch(r"\d+", paramFileNoStr)):
        print('unknown param file: ',paramFileName, '  ...  quitting')
        quit()
    paramFileNo = int(paramFileNoStr)

    print('loading ',paramFileName)
    pd = et.loadParams(paramDir, paramFileName)

    try:
        pu = et.loadPUnits(paramDir, 'units_'+paramFileName)
    except:
        # units change only when a new param is added so give option for one
        #   default units file.
        print('loading default parameter units from: ' + paramDir + defaultUnitsName)
        pu = et.loadPUnits(defaultUnitsDir, defaultUnitsName)

# default the comp 1 role if not in param file
try:
    COMP1 = pd['COMP1']
except:
    # define what is compartment 1 (also determines comp2)
    #COMP1 = 'housing'       # comp2 = et vol
    COMP1 = 'supply_tubing'  # comp2 = housing+et vol
    pd['COMP1'] = COMP1  # store it.

# load unit conversions and phys constants (they never change!)
uc = et.loadUnitConv(defaultUnitsDir, unitsConvfilename)
print('loaded unit conversions')

dataDirNames = ['dataAndyMay24',
                'dataFlowTests23Jul/flowEversion/processed_data',
                'dataFlowTests23Jul/flowFixedVol/processed_data']

# generate string for plot headings etc.
try:
    if pd['Compartments'] == 1:
        compModName = 'ONE Compartment'
    else:
        compModName = 'TWO Compartment'
except:
    compModName = 'Two Compartment'

compModName += ', '+ pd['COMP1']  # what is compartment 1?

#
#   What to plot
#
PLOT_TYPE = 'SIMULATION' # only simulation
PLOT_TYPE = 'OVERLAY'   # includes experimental data
FPLOT = False             # make a force plot as well
REYNOLDSPLOT = False
PHASEPLOT = False
STATEPLOT = False

# States
PRESSURE_TEST = 2
GROWING = 1
STUCK = 0

pd['radius_modes'] = ['constant','box','constrict', 'gap', 'ramp']
pu['radius_modes'] = 'text list'

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


#   Get data file to compare
#

df_hashes = [
    "15g42423",   # these are the data file hashes if needed.
    "891a0abc",   # note: these hashes cannot look like floats or ints!!!
    "52f8bea7",   # examples: 12E404 and 87654321 are NOT allowed
    "e50137ee",
    "eb610645",
    "c4f507b1",
    "4bf58872",
    "0ddc3bfa",
    "5d0f0862",
    "fc284fa7"
]

if PLOT_TYPE == 'OVERLAY':
    #files, mdfiles = et.get_files()
    #for i,fn in enumerate(files):
        #print(f'{i:5}', fn)

    print('Simulating Dataset: ', pd['DataFile'])
    #fn = dataDirName + '/' + pd['DataFile']

    dataFileNames = []
    # DataFile parameter now is just an 8char hash code
    for ddn in dataDirNames:
        dataFileNames += glob.glob(ddn + '/' + '*' + pd['DataFile'] + '*.csv')


    if len(dataFileNames) < 1:
        et.error('Overlay plot: file not found: '+ pd['DataFile'])
    if len(dataFileNames) > 1:
        et.error('Multiple files found: ' + pd['DataFile'] + str(dataFileNames))

    print('Simulating Dataset: ', dataFileNames[0])

    fn =  dataFileNames[0]  # should be just one!!

    print('Opening data file: ', fn)


##################################################  Run Simulation
t1 = 0.0
t2 = 8.0
state = STUCK


(time, th, thdot, state_seq, l, lc, ldot, f, ft, pc1,pc2, pbat, pstt, vol1, vol2, F_e,F_c,F_d,F_j) = sim.simulate(pd,uc,t1,t2)

###################################################
#
#    Load plotting parameters
#
###################################################
fnpr = 'plottingRanges.txt'
fpp = open(fnpr,'r')
prd = {}  #plotting range dict
for line in fpp:
    par, n1, n2 = line.split(':')
    prd[par.strip()] = (float(n1), float(n2))
print(prd)

PltTMIN  = prd['Time'][0]
PltTMAX  = prd['Time'][1]
PltPrMIN = prd['Pressure'][0]
PltPrMAX = prd['Pressure'][1]
#
# TEMP HACK
#PltPrMAX = 140000  #Pa
print('Pressure lims: ', PltPrMIN, PltPrMAX)
PltFlMIN = prd['Flow'][0]
PltFlMAX = prd['Flow'][1]
#
#  HACK
#PltFlMAX = 25 * uc['m3_per_Liter']/uc['sec_per_min']
print('Flow lims    : ', PltFlMIN, PltFlMAX)
PltLeMIN = prd['Length'][0]
PltLeMAX = prd['Length'][1]
PltSpMIN = prd['Speed'][0]
PltSpMAX = prd['Speed'][1]
PltVoMIN = prd['Volume'][0]
PltVoMAX = prd['Volume'][1]



if FPLOT:
    # Create force plot
    fig = plt.figure()
    #fig.suptitle(fn.split('/')[-1] + ' ' + paramFileName)
    fig.suptitle(fn.split('/')[-1] + '\n       ' + paramFileName + ', ' + compModName)

    ax = fig.gca()
    plt.plot(time, F_e)
    plt.plot(time, F_c)
    plt.plot(time, F_d)
    plt.plot(time, F_j)
    plt.plot(time,lc)
    plt.legend(['F_Eversion','F_Coulomb','F_Drag','F_NET','Crumple'])
    ax = plt.gca()
    ax.set_xlim(PltTMIN, PltTMAX)
    ax.set_ylim(-1,7)
    ax.set_ylabel('|F|, (N)')

    print('\nFmax/min Report: ')
    names = ['Feversion','Fcoulomb','Fdrag','F_NET']
    maxs = [max(F_e), max(F_c), max(F_d), max(F_j)]
    mins = [min(F_e), min(F_c), min(F_d), min(F_j)]

    for i,n in enumerate(names):
        print(f'{n:15}   {mins[i]:10.2f}   {maxs[i]:10.2f} ')
    print('')

#REYNOLDSPLOT = False
                   # https://en.wikipedia.org/wiki/Reynolds_number
if REYNOLDSPLOT:   # https://en.wikipedia.org/wiki/Density_of_air
    Re = []    # store Reynolds number
    Rspec = 8.12446  #Joules/K/mol
    M = 0.0289652    #kg/mol of air
    for i,t in enumerate(time):
        pressAbs = pc1[i]  # Pa
        rho = pressAbs*M / (pd['RT'])
        V = ldot[i]  # m/sec
        L = pd['ET_radius'] * 2.0
        mu = 1.81E-05   #Pa sec
        Re.append(rho*V*L/mu)
    fig = plt.figure()
    fig.suptitle(fn.split('/')[-1] + '\n       ' + paramFileName + ', ' + compModName)
    plt.plot(time, Re)
    ax = plt.gca()
    ax.set_xlim(PltTMIN, PltTMAX)
    ax.set_ylim(0,3000)
    ax.set_xlabel('Time')
    ax.set_ylabel('Reynolds Number')

# Create simulation output figure with subplots
fig, axs = plt.subplots(3, 2, figsize=(8, 10))

# Plot 1
#axs[0].plot(vel, p)
vel1 = []
pr1 = []



if state==PRESSURE_TEST:
    axs[0].text(0.5, 10000, 'PRESSURE_TEST (no eversion)')

#
#   Load in the exp. data and briefly inspect/correct it
#

try:
    ed = et.get_data_from_AL_csv(fn)
except:
    et.error("I don't know how to read from data file: "+fn)

#if paramFileNo < 10:
ed = et.convert_units(ed,uc)

# HACK
if paramFileNo < 10:  # correct old flow bug
#if True:  # correct old flow bug
    ed['flow'] *= 10.0
else:
    pass

# adapt plotting limits for high pressure
pexpMax = max(ed['P'])
try:
    x = pd['Pmax']  # if there is a param file entry it takes over
    PltPrMAX = x
except:
    pass

if pexpMax > PltPrMAX:
    rnum = (pexpMax // 1000 + 1)*1000
    print('resetting Pr Max plotting to ', rnum)
    PltPrMAX = rnum


#############################################################################
#
#   Plot the simulation outputs
#

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color'] # match the defaults

clrs = [colors[0],colors[1],'g','c','m']
#trajectory:

# constrain the resistances (5-Aug-24)
r_source, r_tube = et.constrainR(pd)

# "ideal" loadline endpoints
x1 = 0.0
y1 = pd['Psource_SIu']
x2 = (pd['Psource_SIu'] - pd['Patmosphere']) / r_source
y2 = pd['Patmosphere']

#   Presure vs FLow Plot
#

#Presure- Flow phase plot

# special "avg flow" to distribute 1st compartment(??)
favg = 0.5*(np.array(f)+np.array(ft))

axs[0,0].plot([x1,x2], [y1,y2], color='k', linestyle='-.')  # x1,..y2 defined above
# simulated housing and tube pressures
#axs[0,0].plot(f, pc1, favg, pc2)
axs[0,0].plot(f, pc1,  color=clrs[0])
axs[0,0].plot(ft, pc2, color=clrs[1])
if pd['COMP1'] == 'housing':
    axs[0,0].legend(['Source Load Line',  'Housing', 'Ev. Tube'])
if pd['COMP1'] == 'supply_tubing':
    axs[0,0].legend(['Source Load Line',  'Supply Tubing', 'Housing + ET'])

axs[0,0].set_xlabel('Flow (m3/sec)')
axs[0,0].set_ylabel('Pressure (Pa)')
axs[0,0].set_xlim(PltFlMIN, PltFlMAX)
axs[0,0].set_ylim(PltPrMIN, PltPrMAX)
plt.sca(axs[0,0])
ax00 = plt.gca()
xpressticks = ticker.MaxNLocator(3)
ax00.xaxis.set_major_locator(xpressticks)
#plt.xticks([0.0001, 0.0002, 0.0003])

#plt.sca(axs[0,0])
# Plot 2   # PRESSURE
axs[1,0].plot(time, pc1, color=clrs[0])
axs[1,0].plot(time, pc2, color=clrs[1])

if pd['COMP1'] == 'housing':
    axs[1,0].legend(['Phousing', 'Ptube' ])
if pd['COMP1'] == 'supply_tubing':
    axs[1,0].legend(['Supply Tubes', 'Housing + ET' ])
axs[1,0].set_xlabel('Time (sec)')
axs[1,0].set_ylabel('Pressure (Pa)')
axs[1,0].set_xlim(PltTMIN, PltTMAX)
axs[1,0].set_ylim(PltPrMIN, PltPrMAX)
#  plot the eversion thresholds (function of L)
axs[1,0].plot(time, pstt, linestyle='dashed', color=clrs[3])
axs[1,0].plot(time, pbat, linestyle='dashed', color=clrs[4])

# Plot 3
axs[2,0].plot(time, vol1, color=clrs[0] )
axs[2,0].plot(time, vol2, color=clrs[1] )
axs[2,0].legend(['Vhousing', 'Vtube'])
axs[2,0].set_xlabel('Time (sec)')
axs[2,0].set_ylabel('Volume (m3)')
axs[2,0].set_xlim(PltTMIN, PltTMAX)
axs[2,0].set_ylim(PltVoMIN, PltVoMAX)

axs[0,1].plot(time, f,color=clrs[0])
axs[0,1].plot(time, ft, color=clrs[1])
axs[0,1].set_xlabel('Time (sec)')
axs[0,1].set_ylabel('Flow (m3/sec)')
axs[0,1].legend(['Flow Source->C1', 'Flow C1->C2'])
axs[0,1].set_xlim(PltTMIN, PltTMAX)
axs[0,1].set_ylim(PltFlMIN, PltFlMAX)

axs[1,1].plot(time, l,color=clrs[0])
axs[1,1].plot(time, lc,color=clrs[1])
axs[1,1].set_xlabel('Time (sec)')
axs[1,1].set_ylabel('Length (m)')
axs[1,1].legend(['Tube Length', 'Crumple Length'])
axs[1,1].set_xlim(PltTMIN, PltTMAX)
axs[1,1].set_ylim(PltLeMIN, PltLeMAX)

axs[2,1].plot(time, ldot,color=clrs[0])
axs[2,1].set_xlabel('Time (sec)')
axs[2,1].set_ylabel('Tube Velocity (m/sec)')
axs[2,1].set_xlim(PltTMIN, PltTMAX)
axs[2,1].set_ylim(PltSpMIN, PltSpMAX)


#
#   Special plots for state jumps diagnosis
#

if STATEPLOT:
    # state seqence
    fig2,axs2 = plt.subplots(1)
    axs2.plot(time, state_seq)
    axs2.set_xlim(PltTMIN, PltTMAX)
    axs2.set_ylim(0, 3.5)
    #ax = plt.gca()
    stateTicks = ticker.MaxNLocator(4)
    axs2.yaxis.set_major_locator(stateTicks)
    axs2.set_ylabel('Combined State (0-3)')
    axs2.set_xlabel('Time (sec)')
    fig2.suptitle('     Eversion State Sequence\n  '+ fn.split('/')[-1] + '\n       ' + paramFileName + ', ' + compModName)

    fig2.tight_layout()

if PHASEPLOT:   # plot reel phase space (theta vs thetadot)
    dL =  (np.array(th)*pd['rReel'] - np.array(l)) # mm for plotting
    for i,d in enumerate(dL):
        dL[i] = max(0.0, d)
    fig3,axs3 = plt.subplots(1)
    axs3.plot( th,   thdot )
    axs3.legend(['th', 'thdot'])
    stateTicks = ticker.MaxNLocator(4)
    axs3.yaxis.set_major_locator(stateTicks)
    axs3.set_ylabel('thdot')
    axs3.set_xlabel('theta (rad)')
    axs3.grid()
    fig3.suptitle('   Reel Phase Plot\n  '+ fn.split('/')[-1] + '\n       ' + paramFileName + ', ' + compModName)
    fig3.tight_layout()

##########################################################33
#
#    Add Experimental Data to plot
#

if PLOT_TYPE == 'OVERLAY':
    # filename, fn, is chosen above prior to sim

    # data/simulation plot heading title super title
    fig.suptitle(fn.split('/')[-1] + '\n       ' + paramFileName + ', ' + compModName)

    # add pressure to axs[0,0] subplot
    Interval = 25
    x = ed['flow']
    y = ed['P']
    et.plot_curve_with_arrows2(x, y, ax00, Interval,color=clrs[2])
    if pd['COMP1'] == 'housing':
        ax00.legend(['Source Load Line',  'Housing', 'Everting Tube', 'Experiment'])
    if pd['COMP1'] == 'supply_tubing':
        ax00.legend(['Source Load Line',  'Supply Tubing', 'Housing + ET', 'Experiment'])


    axs[1,0].plot(np.array(ed['time']), np.array(ed['P']), '--', color=clrs[2])  # Experimental Data
    #axs[2,0].plot(ed['time'], ed['vol'], '--',color=clrs[1])  # Experimental Data
    axs[0,1].plot(ed['time'], ed['flow'], '--',color=clrs[2])  # Experimental Data
    axs[1,1].plot(ed['time'], ed['L'], '--',color=clrs[2])     # Experimental Data
    axs[2,1].plot(ed['time'], ed['Ldot'], '--',color=clrs[2])  # Experimental Data

    # Adjust layout
    fig.tight_layout()

    if False:
        print('Pressure Simulation Losses: ')
        print('Time Domain: ',  et.TD_loss(p, ed['P']))
        print('Freq Domain: ',  et.FD_loss2(p, ed['P'], dtsim, dtexp, test=True, ptitle='Pressure') )


        print('Velocity Simulation Losses: ')
        print('Time Domain: ',  et.TD_loss(ldot,  ed['Ldot']))
        print('Freq Domain: ',  et.FD_loss2(ldot, ed['Ldot'],dtsim,dtexp,test=True,ptitle='Velocity') )

    # compare dt's between experment and sim
    dtexp = ed['dtexp']
    dtsim = pd['dt']
    print('dt Sim:', dtsim, 'dt Exp:', dtexp)

    # obsolete:  changes are made outside of runtime and thus undetectable.
    #et.print_param_table2(pd,pd_orig, pu)  # print with change markers

    print('Range Report:')
    print(f"Presure: Min: {min(ed['P']):8.3E}  {max(ed['P']):8.3E}")
    print(f"Flow:    Min: {min(ed['flow']):8.3E}  {max(ed['flow']):8.3E}")

if pd['ET_RofL_mode'] != 'constant':  # if the tube shape is interesting, plot it.
    et.plot_tube_shape(pd)

print(f'\n\n           {paramFileName} was simulated by {compModName} model. \n')

pfname = input('Enter filename root for plot file (xxxxxxx.png): (cr) for none.')
if len(pfname) > 4:
    fn = 'smoothTubeWShapes/'+pfname+'.png'
    fig.savefig(fn, dpi=300)  # Save as PNG with 300 DPI resolution

plt.show()

