import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import simulation as sim
import et_lib as et
#from et_lib import error
import glob
import sys


def compare_params(pavg, pcomps):
    nscale = 20
    changes = {}
    #pavg is the dict of avg params
    #pcomp is a list of pdicts
    for pk in pavg.keys():
        diffs = []
        for pd in pcomps:
            diffs.append(pdiff(pavg[pk],pd[pk],nscale))
        changes[pk] = change_str(diffs,nscale)
    return changes


def change_str(diffs,nscale):
    left = '*0.1'
    right = '*10'
    scale = ['.']*((nscale-1)*2 + 1)
    NoDiff = True
    for d in diffs:
        #print(scale, d)
        if abs(d) > 0.1:
            NoDiff = False
            dIndex = int(len(scale)/2 + 0.5 + d)
            if dIndex > -nscale:
                scale[dIndex] = 'X'
    scale[nscale] = '_'
    if NoDiff:
        return ''
    else:
        return left+str(''.join(scale))+right

def pdiff(a,b,nscale):
    n = int(nscale*(np.log10(b/a)))
    n = max(n, -nscale)
    n = min(n, nscale)
    return n


def print_params_by_group(pnames, pd, pg, pch=None, tstr='\n      Parameter Table'):
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
            if pch is None:
                print(f'{k:18}  {pd[k]:8.4E}  {pu[k]:15}   {pg[k]}')
            else:
                print(f'{k:18} {pd[k]:8.4E} {pu[k]:15} {pg[k]:12} {pch[k]}')
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


print('Parameter Changes vs Set5')

pchanges = compare_params(parAvg, [ pdicts[1],  pdicts[6],  pdicts[7]  ])

print_params_by_group(pnames, parAvg, pg, pch=pchanges, tstr='\n    Compare Avg w lo fric')

pchanges = compare_params(parAvg, [ pdicts[0],  pdicts[2],  pdicts[3], pdicts[4],  pdicts[5],  pdicts[8]  ])

print_params_by_group(pnames, parAvg, pg, pch=pchanges, tstr='\n    Compare Avg w HI fric')


