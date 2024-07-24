import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import simulation as sim
import et_lib as et
from et_lib import error
import sys

###     plot out data file only.

# Unit Conversions  (needed to convert to SI units from data file units)

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

paramDir = 'evtParams/'
defaultParamName = 'InitialParams.txt'
pd = et.loadParams(paramDir, defaultParamName)
pu = et.loadPUnits(paramDir, 'units_'+defaultParamName)


args = sys.argv


dataDirName = 'dataAndyMay24'

#   Get data file to compare
#
files, mdfiles = et.get_files()
for i,fn in enumerate(files):
    print(f'{i:5}', fn)

inp = input('Select file by number: ')
try:
    n = int(inp)
except:
    n = -1

if n >= 0:
    if (n>=len(files) or len(files) <=0):
        print('illegal file name index')
        quit()
    fn = files[n]   # includes datadir/

# get inertia and friction from data file name
Ji, Fric_i = et.get_inr_fric_from_name(fn)


#############################################################################
#
#   Plot the data

# Create data figure with subplots
fig, axs = plt.subplots(3, 2, figsize=(8, 10))

#vel1 = []
#pr1 = []

clrs = ['b','g','r','c','m']
#trajectory:


#  Read in the data for plotting

fig.suptitle(fn.split('/')[-1] )

print('opening: ',fn)
x=input('       ... OK?? <cr>')

ed = et.get_data_from_AL_csv(fn)
ed = et.convert_units(ed,uc)

# HACK
ed['flow'] *= 10.0


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

FPLOT = False
PltTMIN  = prd['Time'][0]
PltTMAX  = prd['Time'][1]
PltPrMIN = prd['Pressure'][0]
PltPrMAX = prd['Pressure'][1]
PltFlMIN = prd['Flow'][0]
PltFlMAX = prd['Flow'][1]
PltLeMIN = prd['Length'][0]
PltLeMAX = prd['Length'][1]
PltSpMIN = prd['Speed'][0]
PltSpMAX = prd['Speed'][1]
PltVoMIN = prd['Volume'][0]
PltVoMAX = prd['Volume'][1]

#########################################################################
#
#   Plot all the subplots
#

PltTMIN = 0.0
PltTMAX = 8.0

#   Presure vs FLow Plot (exper phase plot)
#
#

x = ed['flow']
y = ed['P']
time = ed['time']

    #HW= 0.00001
    #HL= 0.0001
Interval = 25
et.plot_curve_with_arrows2(x, y, axs[0,0], Interval,color=clrs[1])
axs[0,0].set_xlabel('Flow (m3/sec)')
axs[0,0].set_ylabel('Pressure (Pa)')
axs[0,0].set_xlim(0,0.0003)
plt.sca(axs[0,0])
plt.xticks([0.0001, 0.0002, 0.0003])

axs[0,0].set_xlim(PltFlMIN, PltFlMAX)
axs[0,0].set_ylim(PltPrMIN, PltPrMAX)

# Plot 2   # PRESSURE
axs[1,0].plot(np.array(ed['time']), np.array(ed['P']), '--', color=clrs[1])  # Experimental Data
#axs[1,0].plot(time, pc1)
#axs[1,0].plot(time, pc2)
axs[1,0].set_xlabel('Time (sec)')
axs[1,0].set_ylabel('Pressure (Pa)')


# Plot 3
axs[2,0].plot(time, ed['vol1'], time, ed['vol2'] )
axs[2,0].legend(['Vhousing', 'Vtube'])
axs[2,0].set_xlabel('Time (sec)')
axs[2,0].set_ylabel('Volume (m3)')

axs[2,0].set_xlim(PltTMIN, PltTMAX)
axs[2,0].set_ylim(PltVoMIN, PltVoMAX)


axs[0,1].plot(ed['time'], ed['flow'], '--',color=clrs[1])  # Experimental Data
#axs[0,1].plot(time, f, time, ft)
axs[0,1].set_xlabel('Time (sec)')
axs[0,1].set_ylabel('Flow (m3/sec)')
axs[0,1].legend(['Source Flow', 'Tube Flow'])

axs[0,1].set_xlim(PltTMIN, PltTMAX)
axs[0,1].set_ylim(PltFlMIN, PltFlMAX)

axs[1,1].plot(ed['time'], ed['L'], '--',color=clrs[1])  # Experimental Data
#axs[1,1].plot(time, l, time, lc)
axs[1,1].set_xlabel('Time (sec)')
axs[1,1].set_ylabel('Length (m)')

axs[1,1].set_xlim(PltTMIN, PltTMAX)
axs[1,1].set_ylim(PltLeMIN, PltLeMAX)


axs[2,1].plot(ed['time'], ed['Ldot'], '--',color=clrs[1])  # Experimental Data
#axs[2,1].plot(time, ldot)
axs[2,1].set_xlabel('Time (sec)')
axs[2,1].set_ylabel('Tube Velocity (m/sec)')

axs[2,1].set_xlim(PltTMIN, PltTMAX)
axs[2,1].set_ylim(PltSpMIN, PltSpMAX)




#if PLOT_TYPE == 'OVERLAY':

    #fig.suptitle(fn.split('/')[-1] + '\n       ' + paramFileName)

    #print('opening: ',fn)
    #x=input('       ... OK?? <cr>')

    #ed = et.get_data_from_AL_csv(fn)
    #ed = et.convert_units(ed,uc)

    ## HACK
    #ed['flow'] *= 10.0

    ## phase plot with arrows

    #HW= 0.00001
    #HL= 0.0001
    #Interval = 25

    #x = ed['flow']
    #y = ed['P']

    #Plt_ranges = [0, 0.001, pd['Patmosphere'], pd['Psource_SIu'] ]

    ##plt.figure()
    ##ax = plt.gca()
    #ax = axs[0,0]
    ##plot_curve_with_arrows(xdata, ydata, Interval, ax, arrow_scale=1.0 )
    #et.plot_curve_with_arrows2(x, y, ax, Interval,color=clrs[1])

    #axs[1,0].plot(np.array(ed['time']), np.array(ed['P']), '--', color=clrs[1])  # Experimental Data

    #dtexp = ed['dtexp']
    #dtsim = pd['dt']
    #print('dt Sim:', dtsim, 'dt Exp:', dtexp)

    ##axs[2,0].plot(ed['time'], ed['vol'], '--',color=clrs[1])  # Experimental Data
    #axs[0,1].plot(ed['time'], ed['flow'], '--',color=clrs[1])  # Experimental Data
    #axs[1,1].plot(ed['time'], ed['L'], '--',color=clrs[1])  # Experimental Data
    #axs[2,1].plot(ed['time'], ed['Ldot'], '--',color=clrs[1])  # Experimental Data

    #if False:
        #print('Pressure Simulation Losses: ')
        #print('Time Domain: ',  et.TD_loss(p, ed['P']))
        #print('Freq Domain: ',  et.FD_loss2(p, ed['P'], dtsim, dtexp, test=True, ptitle='Pressure') )


        #print('Velocity Simulation Losses: ')
        #print('Time Domain: ',  et.TD_loss(ldot,  ed['Ldot']))
        #print('Freq Domain: ',  et.FD_loss2(ldot, ed['Ldot'],dtsim,dtexp,test=True,ptitle='Velocity') )



#Adjust layout
plt.tight_layout()



if pd['ET_RofL_mode'] != 'constant':  # if the tube shape is interesting, plot it.
    et.plot_tube_shape(pd)

plt.show()

