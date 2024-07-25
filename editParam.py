import numpy as np
import matplotlib.pyplot as plt
import et_lib as et
import glob, os




#######################################################
#
#    Configure here
#
par = 'et_MPM'
newval = 0.0032
SIMULATE = False
#
#######################################################



## Unit Conversions
#uc = {
    #"sec_per_min": 60,
    #"kPa_per_Pa": 0.001,
    #"Pa_per_kPa": 1.0 / 0.001,
    #"min_per_sec": 1 / 60,
    #"Gal_per_Liter": 0.2642,
    #"Liter_per_Gal": 3.7854,
    #"Liter_per_m3": 1000.0,
    #"Liter_per_mm3": 1000.0 / 1000**3,
    #"Gal_per_mm3": (1000.0 / 1000**3) * 0.2642,
    #"mm3_per_Gal": 1.0 / ((1000.0 / 1000**3) * 0.2642),
    #"MM3perLiter": 1.0 / (1000.0 / 1000**3), # Ideal Gas Law  https://pressbooks.uiowa.edu/clonedbook/chapter/the-ideal-gas-law/
    #"m3_per_mole": 0.02241,  # m3/mol of Air
    #"moles_per_m3": 1.0 / 0.02241,
    #"Pa_per_PSI": 6894.76,
    #"atmos_Pa": 14.5 * 6894.76,
    #"et_MPM"  : 0.0032 : kg/m
    #"m3_per_Liter": 1.0 / 1000.0  # m3
#}


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

ParamDirNames = ['evtParams']

#files = ['eversion_flow-hi-inr_hi-fric_tube-1_trial-2.csv', 'eversion_flow-hi-inr_lo-fric_tube-1_trial-1.csv', 'eversion_flow-hi-inr_lo-fric_tube-1_trial-2.csv', 'eversion_flow-hi-inr_lo-fric_tube-1_trial-3.csv', 'eversion_flow-hi-inr_hi-fric_tube-2_trial-1.csv', 'eversion_flow-hi-inr_hi-fric_tube-2_trial-2.csv', 'eversion_flow-hi-inr_hi-fric_tube-3_trial-1.csv', 'eversion_flow-hi-inr_hi-fric_tube-3_trial-2.csv', 'eversion_flow-hi-inr_hi-fric_tube-3_trial-3.csv',
#]

files = []
fnRoots = []
for parDir in ParamDirNames:
        hashesRemoved = set()
        parFiles = list(glob.glob(parDir + '/' + "*"))
        parFiles.sort(key=lambda x: os.path.basename(x),reverse=False) # newest first
        filenameroots = []
        if len(parFiles) ==0:
            cto.error('No param.txt files found')
        for f in parFiles:
            if '.txt' in f and 'Set' in f:     # e.g. Set5Params.txt
                files.append(f)

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

if SIMULATE:
    print('  ...   SIMULATING  ...')
for index in fset:
    print('Param set: ', files[index])
    pd = et.loadDict('', files[index])
    print('    Loaded')
    print(f'    changing {par} from {pd[par]:12} to {newval}')
    if not SIMULATE:
        pd[par] = newval
        et.saveParams(files[index],pd)
    print('    Saved')

if SIMULATE:
    print(' ... end of SIMULATION ...')
print('done')
