import numpy as np
import matplotlib.pyplot as plt
import et_lib as et
import glob, os
import re
import sys




#######################################################
#
#    compute the average of all "Free" params
#
#######################################################


cl_args = sys.argv

if len(cl_args) > 1:
    print('''

        compute the average of all "Free" params
        (you can select which pfiles to include in avg).

        ''')
    quit()

paramDir = 'evtParams/'

ParamDirNames = ['evtParams', 'evtParams/2Comp', 'evtParams/1Comp', 'evtParams/23-Jul-FlowData']

ParamDirNames = ['evtParams/2Comp']

unitsConvfilename = 'unitConv.txt'
defaultParamName = 'InitialParams.txt'
defaultUnitsName = 'units_'+defaultParamName
unitsConvfilename = defaultUnitsName

pRefD = et.loadParams(paramDir, defaultParamName)
pd = pRefD.copy()

print('COMP1: ', pd['COMP1'])

files = []
fnRoots = []
parFiles = []
fnums = []
filenameroots = []
for parDir in ParamDirNames:
    print (f'Checking: {parDir}')
    DparFiles = list(glob.glob(parDir + '/Set*.txt'))
    print (f'         Found:  {DparFiles}')
    for pf in DparFiles:
        fname = pf.split('/')[-1]
        try:
            SetNo = int(re.search(r'\d+', fname).group())
        except:
            et.error('No match to int in file name: '+pf)
        fnums.append(SetNo)
    tmp = list(zip(fnums, DparFiles))
    tmp2 = sorted(tmp, key=lambda x: x[0] ) # numerical order
    for n, fn in tmp2:
        parFiles.append(fn)
if len(parFiles) ==0:
    et.error('No SetXXparam.txt files found')
files = parFiles

#print('file set:  ', files)
###################################################
#
#   Select Files
#
###################################################
print('\n Discovered Data Files: ')
for i,fn in enumerate(files):
   # parts = fn.split('/')
   # fn2 = parts[-1]
   # fnRoot = '/'.join(parts[:-1])
   # print(f'{i:3}  {fnRoot} {fn2}')
   print(f'{i:3}  {fn}')

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

print('File number set: ', fset)

freeParams = [ 'K2drag', 'Kdrag', 'PBA_static', 'PHalt_dyn', 'Psource_SIu', 'Rsource_SIu',
    'Threshold Taper', 'Tau_coulomb', 'Tau_coulomb', ]

n = len(files)
fpvals = {}
for fp in freeParams:
    fpvals[fp] = np.zeros(len(fset))   # initialize (store all param vals for computing mu,sig)

pu = et.loadPUnits('evtParams','units_InitialParams.txt')

for index in fset:  # load the pfiles we are averaging and read each par.
    #print('Param set: ', files[index])
    pd = et.loadParams('', files[index])
    print(index, ' COMP1: ', pd['COMP1'], files[index])
    for fp in freeParams:
        fpvals[fp][index] = pd[fp]

# compute avgs
fpAvgs = {}
fpStds = {}
for fp in freeParams:
    if fp == 'Compartments':
        print('Compartments data:')
        print(fpvals[fp])
    fpAvgs[fp] = np.mean(fpvals[fp])
    fpStds[fp] = np.std(fpvals[fp])

print('\n---------------------------------------------------')
print('  Parameter set including Averaged Free Parameters (n=',len(fset),')' )

## fill in any missing values
#for k in pRefD.keys():
    #try:   # defaults for missing values
        #x = pd[k]
    #except:
        #print('missing key: ', k)
        #pd[k] = pRefD[k]
#print('      ...   fixed: ', pd['COMP1'])

for fp in freeParams:
    try:
        x=pd[fp]
    except:
        et.error('updating params with avgs: something wrong: ', fp)
    pd[fp] = fpAvgs[fp]

et.print_param_table(pd, pu)

print('\n------Standard Deviations (free Params)---------------')
for k in pRefD.keys():
    if type(pd[k]) == type('hello'):
        print(f'{k:16}: {pd[k]:10}')
    else:
        if k in freeParams:
            if abs(pd[k]) > 0.0000001: # check for divide by zero
                pct = 100.0 * fpStds[k]/pd[k]
                print(f'{k:16}: {pd[k]:10.4E} +/- {fpStds[k]:10.4E} ({pct:3.1f}% stdev)')
            else:
                pct = 100.0 * fpStds[k] #scale by 1.0 arbitrarily
                print(f'{k:16}: {pd[k]:10.4E} +/- {fpStds[k]:10.4E} ')

print('---------------------------\n\n')
#for paramName in pd.keys(l|c|l):
    #try:
        #print(f'{paramName:16}  : {fpAvgs[paramName]:10.5E}')
    #except:
        #print(f'{paramName:16}  : {pd[paramName]:10.5E}')

print('-------------- Latex Table Output ----------------')
et.print_param_table_latex(pd, pu)
print('---------------------------\n\n')

print('-------------- Complete editParam.py config code output ----------------')
pdEdits = pRefD.copy()
del pdEdits['Compartments']  # these are fixed pars
del pdEdits['DataFile']
del pdEdits['COMP1']
del pdEdits['ET_RofL_mode']
del pdEdits['Patmosphere']
del pdEdits['RT']
del pdEdits['T']
del pdEdits['ET_radius']
del pdEdits['dt']
del pdEdits['et_MPM']
del pdEdits['rReel']

newvals = '['
for k in pdEdits.keys():
    newvals += f"{pd[k]}, "
newvals += ']\n'
parlist = '['
for k in pdEdits.keys():
    parlist += f"'{k}', "
parlist += ']\n'

print('pars = ',parlist)
print('newvals = ', newvals)
print('---------------------------\n\n')



print('done')
