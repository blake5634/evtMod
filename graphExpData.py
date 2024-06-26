import numpy as np
import matplotlib.pyplot as plt
import et_lib as et
import glob, os

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
        for f in files:
            #print('found: ',f)
            if '.csv' in f :
                filenameroots.append(str(f).replace('.csv',''))
            if  '_meta.json' in f:
                filenameroots.append(str(f).replace('_meta.json',''))
        # dict.fromkeys better than set because preserves order (thanks google)
        filenameroots = list(dict.fromkeys(filenameroots)) # elim dupes
        files += tfiles
        fnRoots += filenameroots

PltTMIN = 0.0
PltTMAX = 8.0

# Create a figure with subplots
fig, axs = plt.subplots(3, 2, figsize=(8, 10))

for i,fn in enumerate(files):
    #print(f'{i:5}', fn)

    #fn = dataDirName + '/' + fn

    try:
        ed = et.get_data_from_AL_csv(fn)
    except:
        print('could not open:   ', fn)
        continue
    print('       opening:   ', i, fn)

    ed = et.convert_units(ed)

    axs[0,0].plot(ed['flow'],  ed['P'])  # Experimental Data
    axs[0,0].set_xlabel('Source Flow (m3/sec)')
    axs[0,0].set_ylabel('Pressure (Pa)')
    #axs[0,0].set_xlim(0,0.01)
    axs[0,0].set_ylim(Patmosphere, Psource_SIu)


    axs[1,0].plot(ed['time'], ed['P'])  # Experimental Data
    axs[1,0].set_xlabel('Time (sec)')
    axs[1,0].set_ylabel('Pressure (Pa)')
    axs[1,0].set_xlim(PltTMIN, PltTMAX)
    axs[1,0].set_ylim(Patmosphere, Psource_SIu)

    #axs[2,0].plot(ed['time'], ed['vol'])  # Experimental Data
    #axs[2,0].set_xlabel('Time (sec)')
    #axs[2,0].set_ylabel('Pressure (Pa)')

    axs[0,1].plot(ed['time'], ed['flow'])  # Experimental Data
    axs[0,1].set_xlabel('Time (sec)')
    axs[0,1].set_ylabel('Flow (m3)')
    axs[0,1].set_xlim(PltTMIN, PltTMAX)
    axs[0,1].set_ylim(      0, 0.00020 )

    axs[1,1].plot(ed['time'], ed['L'])  # Experimental Data
    axs[1,1].set_xlabel('Time (sec)')
    axs[1,1].set_ylabel('Length (m) (reel)')
    axs[1,1].set_xlim(PltTMIN, PltTMAX)
    axs[1,1].set_ylim(      0, 0.600 )

    axs[2,1].plot(ed['time'], ed['Ldot'])  # Experimental Data
    axs[2,1].set_xlabel('Time (sec)')
    axs[2,1].set_ylabel('Speed (m/sec)')
    axs[2,1].set_xlim(PltTMIN, PltTMAX)
    axs[2,1].set_ylim(      0, 0.8 )

    leglist = []
    for i in range(len(files)):
        leglist.append(f'file {i}')
    axs[2,0].plot(ed['time'],[1.0]*len(ed['time']))
    axs[2,0].legend(leglist)
# Adjust layout
plt.tight_layout()

plt.show()
