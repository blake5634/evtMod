import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import et_lib as et
import glob, os

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


paramDir = 'evtParams/'
unitsConvfilename = 'unitConv.txt'
defaultParamName = 'InitialParams.txt'
defaultUnitsName = 'units_'+defaultParamName
paramFileName = defaultParamName
unitsConvfilename = defaultUnitsName

pd = et.loadParams(paramDir, paramFileName)
#uc = et.loadUnitConv(paramDir, unitsConvfilename)

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
m3_per_Liter =  1.0 / Liter_per_m3  # m3
Patmosphere = 101325.0    # Pascals
Psource_SIu = Patmosphere + 3.0 * Pa_per_PSI # pascals

dataDirNames = ['dataAndyMay24']

#files = ['eversion_flow-hi-inr_hi-fric_tube-1_trial-2.csv', 'eversion_flow-hi-inr_lo-fric_tube-1_trial-1.csv', 'eversion_flow-hi-inr_lo-fric_tube-1_trial-2.csv', 'eversion_flow-hi-inr_lo-fric_tube-1_trial-3.csv', 'eversion_flow-hi-inr_hi-fric_tube-2_trial-1.csv', 'eversion_flow-hi-inr_hi-fric_tube-2_trial-2.csv', 'eversion_flow-hi-inr_hi-fric_tube-3_trial-1.csv', 'eversion_flow-hi-inr_hi-fric_tube-3_trial-2.csv', 'eversion_flow-hi-inr_hi-fric_tube-3_trial-3.csv',
#]

files = []
fnRoots = []
for datadir in dataDirNames:
        hashesRemoved = set()
        tfiles = list(glob.glob(datadir + '/' + "*"))
        tfiles.sort(key=lambda x: os.path.getmtime(x),reverse=True) # newest first
        filenameroots = []
        if len(tfiles) ==0:
            cto.error('No brl_data files found')
        for f in tfiles:
            #print('found: ',f)
            if '.zip' in f:
                next
            if '.csv' in f :
                filenameroots.append(str(f).replace('.csv',''))
                files.append(f)
            if  '_meta.json' in f:
                filenameroots.append(str(f).replace('_meta.json',''))
        # dict.fromkeys better than set because preserves order (thanks google)
        filenameroots = list(dict.fromkeys(filenameroots)) # elim dupes
        #files += tfiles
        fnRoots += filenameroots

#print('file set:  ', files)
###################################################
#
#   Select Files
#
###################################################
print('Discovered Data Files: ')
for i,fn in enumerate(files):
    fn2 = fn.split('/')[1]
    print(f'{i:3}  {fn2}')

sel = str(input('Select file numbers (-1) for all: '))
nums = sel.split()
if len(nums)<1:
    nums = [-1]
fset=[] # place to collect user-selected filenames
if int(nums[0]) < 0:
    fset = range(len(files))
else:
    for n in nums:
        fset.append(int(n))


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
print('Pressure lims: ', PltPrMIN, PltPrMAX)
PltFlMIN = prd['Flow'][0]
PltFlMAX = prd['Flow'][1]
PltLeMIN = prd['Length'][0]
PltLeMAX = prd['Length'][1]
PltSpMIN = prd['Speed'][0]
PltSpMAX = prd['Speed'][1]
PltVoMIN = prd['Volume'][0]
PltVoMAX = prd['Volume'][1]


# Create a figure with subplots
fig, axs = plt.subplots(3, 2, figsize=(8, 10))

for i,fnum in enumerate(fset):
    fn = files[int(fnum)]

    #print(f'{i:5}', fn)

    #fn = dataDirName + '/' + fn

    try:
        ed = et.get_data_from_AL_csv(fn)
    except:
        print('could not open:   ', fn)
        continue
    print('       opening:   ', i, fn)

    ed = et.convert_units(ed,uc)


    # HACK  (need to get to bottom of this!!)
    ed['flow'] *= 10.0

    # compute volume via pi r**2 L
    vcomp = []
    for l in ed['L']:
        vcomp.append(np.pi*pd['ET_radius']**2*l)

    axs[0,0].plot(ed['flow'],  ed['P'])  # Experimental Data
    axs[0,0].set_xlabel('Source Flow (m3/sec)')
    axs[0,0].set_ylabel('Pressure (Pa)')

    axs[0,0].set_xlim(PltFlMIN, PltFlMAX)
    axs[0,0].set_ylim(PltPrMIN, PltPrMAX)

    plt.sca(axs[0,0])
    ax = plt.gca()
    xpressticks = ticker.MaxNLocator(3)
    ax.xaxis.set_major_locator(xpressticks)

    axs[1,0].plot(ed['time'], ed['P'])  # Experimental Data
    axs[1,0].set_xlabel('Time (sec)')
    axs[1,0].set_ylabel('Pressure (Pa)')

    axs[1,0].set_xlim(PltTMIN, PltTMAX)
    axs[1,0].set_ylim(PltPrMIN, PltPrMAX)

    #axs[2,0].plot(ed['time'], ed['vol'])  # Experimental Data
    #axs[2,0].set_xlabel('Time (sec)')
    #axs[2,0].set_ylabel('Pressure (Pa)')

    leglist = []
    for i in range(len(files)):
        leglist.append(f'file {i}')
    axs[2,0].plot(ed['time'],[ pd['Vhousing_m3'] ]*len(ed['time']))
    axs[2,0].plot(ed['time'],vcomp)
    axs[2,0].legend(leglist)
    axs[2,0].set_ylabel('Volume (m3)')
    axs[2,0].set_xlabel('Time (sec)')

    axs[2,0].set_xlim(PltTMIN, PltTMAX)
    axs[2,0].set_ylim(PltVoMIN, PltVoMAX)

    axs[2,0].set_xlim(PltTMIN, PltTMAX)
    axs[2,0].set_ylim(PltVoMIN, PltVoMAX)

    axs[0,1].plot(ed['time'], ed['flow'])  # Experimental Data
    axs[0,1].set_xlabel('Time (sec)')
    axs[0,1].set_ylabel('Flow (m3/sec)')
    axs[0,1].set_xlim(PltTMIN, PltTMAX)

    axs[0,1].set_ylim(PltFlMIN, PltFlMAX)

    axs[1,1].plot(ed['time'], ed['L'])  # Experimental Data
    axs[1,1].set_xlabel('Time (sec)')
    axs[1,1].set_ylabel('Length (m)')

    axs[1,1].set_xlim(PltTMIN, PltTMAX)
    axs[1,1].set_ylim(PltLeMIN, PltLeMAX)

    axs[2,1].plot(ed['time'], ed['Ldot'])  # Experimental Data
    axs[2,1].set_xlabel('Time (sec)')
    axs[2,1].set_ylabel('Speed (m/sec)')
    axs[2,1].set_xlim(PltTMIN, PltTMAX)

    axs[2,1].set_ylim(PltSpMIN, PltSpMAX)


fig.suptitle('Eversion Data Set, A. Lewis, May 24')

# Adjust layout
plt.tight_layout()
plt.show()
