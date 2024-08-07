import numpy as np
import matplotlib.pyplot as plt
import et_lib as et
import glob, os




#######################################################
#
#    Configure here
#
freeParams = ['K2drag', 'Kdrag', 'PBA_static', 'PHalt_dyn', 'Psource_SIu', 'Rsource_SIu', 'Threshold Taper']

#pars = freeParams
#newvals = [8.7556E-1,  2.0389, 1.1432E5, 1.0791E5, 1.2761E5,
           #1.0489E8, 3.1858E3 ]

pars = ['Compartments']
newvals = [2]

#paramDirList = ['evtParams/1Comp/']#, 'evtParams/2Comp/']
paramDirList = [ 'evtParams/2Comp/']

SIMULATE = False
#
#######################################################



paramDir = 'evtParams'

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

ParamDirNames = paramDirList

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
    #tmp = fn.split('/')
    fn2 = fn
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
    for i,par in enumerate(pars):
        try:
            print(f'    changing {par} from {pd[par]:12} to {newvals[i]}')
        except:
            print(f'    changing {par} from {"missing":12} to {newvals[i]}')

        if not SIMULATE:
            pd[par] = newvals[i]
    if not SIMULATE:
        et.saveParams(files[index],pd)
    print(files[index], ' ...   Saved')

if SIMULATE:
    print(' ... end of SIMULATION ...')
print('done')
