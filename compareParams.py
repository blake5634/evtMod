import numpy as np
import matplotlib.pyplot as plt
import et_lib as et
import glob, os
import sys
import re




#######################################################
#
#    Configure here
#
freeParams = ['K2drag', 'Kdrag', 'PBA_static', 'PHalt_dyn', 'Psource_SIu', 'Rsource_SIu', 'Threshold Taper']

#pars = freeParams
#newvals = [8.7556E-1,  2.0389, 1.1432E5, 1.0791E5, 1.2761E5,
           #1.0489E8, 3.1858E3 ]

pars = ['K2drag']
newvals = [0.0]

SIMULATE = False
#
#######################################################

paramDir = 'evtParams/'
unitsConvfilename = 'unitConv.txt'
defaultParamName = 'InitialParams.txt'
defaultUnitsName = 'units_'+defaultParamName
paramFileName = defaultParamName
unitsConvfilename = defaultUnitsName

pd = et.loadParams(paramDir, paramFileName)
#uc = et.loadUnitConv(paramDir, unitsConvfilename)

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

groupsfile = 'evtParams/ParamGroups.txt'


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
        et.error('No match to int in file name: '+pf)
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

sel = str(input('Select file numbers to compare (-1) for all: '))
nums = sel.split()
#if len(nums)<1:
    #nums = [-1]
if len(nums) != 2:
    et.error(' I can only compare 2 files at this time...')
fset=[] # place to collect user-selected filenames
if int(nums[0]) < 0:
    fset = range(len(files))
else:
    for n in nums:
        fset.append(int(n))

fnum1  =fset[0]
fnum2 = fset[1]
pd1 = et.loadDict('', files[fnum1])
pd2 = et.loadDict('', files[fnum2])

print(f'Comparing {files[fnum1]} to {files[fnum2]}')
for k in pd1.keys():
    star = ''
    try:
        pd2[k]
    except:
        pd2[k] = 'None'
    if type(pd1[k])==type(5.78):
        if abs(pd1[k] - pd2[k]) > pd1[k]*0.02:
            star = '*'
        print(f'{k:16}: {pd1[k]:10.4E} ... {pd2[k]:10.4E} {star}')
    else:
        if pd1[k] != pd2[k]:
            star = '*'
        print(f'{k:16}: {pd1[k]:10} ... {pd2[k]:10} {star}')



print('done')
