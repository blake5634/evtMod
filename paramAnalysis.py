import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import simulation as sim
import et_lib as et
#from et_lib import error
import glob
import sys

def print_params_by_group(pnames, pd, pg, tstr='\n      Parameter Table'):
    print(tstr)
    # sort param names by their group
    pngs = []
    for pn in pnames:
        pngs.append(pg[pn])
    parByGroup = sorted(zip(pnames, pngs), key= lambda x: x[1])
    for k2 in parByGroup:
        k = k2[0]
        val = pd[k]
        if type(val) != type('x'):
            print(f'{k:18}  {pd[k]:8.4E}  {pu[k]:15}   {pg[k]}')
        else:
            print(f'{k:18}  {pd[k]:10}  {pu[k]:15}  (string????)')

    print('')

#
#  load parameter classes
#

paramDir = 'evtParams/'
paramGFile = 'ParamGroups.txt'   # conceptual groupings of parms

pg = et.loadDict(paramDir, paramGFile)
pu = et.loadPUnits(paramDir, 'units_InitialParams.txt')

#files, mdfiles = et.get_files()
#print('Parameter File Analysis')
pfiles = list(glob.glob(paramDir + "Set*.txt"))
pfiles.sort(key=lambda x: x[3]) # Set0, Set1, ...
if len(pfiles) == 0:
    et.error('No eversion parameter files found')
pdicts = []
for i,fn in enumerate(pfiles):
    print(f'{i:5}  Loading: {fn}')
    pdicts.append(et.loadParams('',fn))

# remove all params in group "Admin"
for pd in pdicts:
    parlist = list(pd.keys())
    for par in parlist:
        try:
            pgroup = pd[par]
            #print('checking par: ', par, pgroup)
            if pg[par] == 'Admin':    # check against parameter groups
                #print('             deleting key: ',par)
                del pd[par]
            if pg[par] == 'Physics':    # check against parameter groups
                #print('             deleting key: ',par)
                del pd[par]
        except:
            pass

pnames = pdicts[0].keys()  # assume they're the same

parAvg = {}
for k in pnames:
    parAvg[k] = 0.0

n = 0
for f in pdicts:
    n += 1
    for k in pnames:
        parAvg[k] += f[k]

print('Average parameters')
for k in pnames:
    parAvg[k] /= n

print_params_by_group(pnames, parAvg, pg, tstr='\n     Average Parameters by Group')


