import numpy as np
import matplotlib.pyplot as plt
import et_lib as et
import glob, os
import re




#######################################################
#
#    Configure here
#

ParamDirNames = ['evtParams/2Comp/']

# example of multiple changes:
#freeParams = ['K2drag', 'Kdrag', 'PBA_static', 'PHalt_dyn', 'Psource_SIu', 'Rsource_SIu', 'Threshold Taper']
#pars = freeParams
#newvals = [8.7556E-1,  2.0389, 1.1432E5, 1.0791E5, 1.2761E5,
           #1.0489E8, 3.1858E3 ]

pars = ['K2drag']
newvals = [0.0]

pars =  ['J', 'K2drag', 'Kdrag', 'Lmax', 'PBA_static', 'PHalt_dyn', 'Psource_SIu', 'Rsource_SIu', 'Tau_coulomb', 'Threshold Taper'  ]

newvals =  [0.000262, 5.555555555555555, 2.111111111111111, 0.7, 112943.33333333333, 108844.44444444444, 125822.22222222222, 135000000.0, 0.0032222222222222222, 0.0  ]




######
SIMULATE = False       # don't do the edit for real: only test
#
#
##################################################################


paramDir = ParamDirNames[0]
defaultParamDir = 'evtParams/'
unitsConvfilename = 'unitConv.txt'
defaultParamName = 'InitialParams.txt'
defaultUnitsName = 'units_'+defaultParamName
unitsConvfilename = defaultUnitsName

fnums = []
parFiles = []
DparFiles = list(glob.glob(paramDir + '/Set*.txt'))
if len(DparFiles) ==0:
    et.error('No Setxxparam.txt files found')
for pf in DparFiles:
    fname = pf.split('/')[-1]
    try:
        SetNo = int(re.search(r'\d+', fname).group())
    except:
        et.error('No match to int in file name: '+str(fname) + ' fullpath:'+pf)
    fnums.append(SetNo)
tmp = list(zip(fnums, DparFiles))
tmp2 = sorted(tmp, key=lambda x: x[0] ) # numerical order
for n, fn in tmp2:
    parFiles.append(fn)

files = parFiles
#print('file set:  ', files)
###################################################
#
#   Select Files
#
###################################################
print('Discovered Data Files: ')
for i,fn in enumerate(parFiles):
    fn2 = '/'.join(fn.split('/')[1:])
    print(f'{i:3}  {fn2}')

if SIMULATE:
    print('\n            Simulating ...\n')
else:
    x = input('NOT Simulating... are you sure? <cr>')

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
    print('  ...   SIMULATING (no file changes) ...')
for index in fset:
    print('Param set: ', files[index])
    pd = et.loadDict('', files[index])
    for i,par in enumerate(pars):
        print(f'    changing {par:17} from:{pd[par]:8.2}   to:{newvals[i]:8.2}')
        if not SIMULATE:
            pd[par] = newvals[i]
    if not SIMULATE:
        et.saveParams(files[index],pd)
        print(files[index], ' ...   Saved')
    else:
        print('Simulating: ... save modified params to ', files[index])

if SIMULATE:
    print(' ... end of SIMULATION ...')
print('done')
