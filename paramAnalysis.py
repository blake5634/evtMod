import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import simulation as sim
import et_lib as et
#from et_lib import error
import glob
import sys

def compare_param_lists(pdl1,pdl2,nscale):
    if len(pdl1) < len(pdl2):
        minList = pdl1
    else:
        minList = pdl2
    for d in pdl1:


def compare_param_list_avgs(pdl1, pdl2,nscale):
    avd1 = {}
    avd2 = {}
    n=0
    for k in pdl1[0].keys():
        avd1[k] = 0.0
        avd2[k] = 0.0
    for d in pdl1:
        for k in d.keys():
            avd1[k] += d[k]
        n += 1
    n=0
    for d in pdl2:
        for k in d.keys():
            avd2[k] += d[k]
        n += 1

    for k in pdl1[0].keys():
        avd1[k] /= n
        avd2[k] /= n
    changes = compare_param_dicts(avd1, avd2, nscale)
    return changes

# compare two param dicts
def compare_param_dicts(pd1, pd2, nscale):
    changes = {}
    #pd1 and pd2 are param dicts
    for pk in pd1.keys():
        diff= [ pdiff(pd1[pk],pd2[pk],nscale) ]
        #print('comparing: ', pd1[pk], pd2[pk], pdiff(pd1[pk],pd2[pk],nscale))
        changes[pk] = change_str(diff,nscale)
    return changes


def change_str(diffs,nscale):
    left = '*0.1'
    right = '*10'
    scale = ['.']*((nscale-1)*2 + 1)
    NoDiff = True
    for d in diffs:
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


def print_params_by_group(pg1, pGroupsD, pu, pg2=None, pch=None, tstr='\n      Parameter Table'):
    print(tstr)
    pnames = pg1.keys()
    # sort param names by their group
    pngs = []
    for pn in pnames: # parameter names
        pngs.append(pGroupsD[pn])  # list of groups
    parByGroup = sorted(zip(pnames, pngs), key= lambda x: x[1])
    for k2 in parByGroup: # sorted 2xtuples
        k =    k2[0]
        val1 = pg1[k]
        if pg2:
            val2 = pg2[k]
        if type(val1) != type('x'):
            if pg2 is None:
                print(f'{k:18} {pg1[k]:8.4E} {pu[k]:15} ')
            else:
                print(f'{k:18} {pg1[k]:8.4E}-->{pg2[k]:8.4E} {pu[k]:15}   {pch[k]}')
        else:
            print(    f'{k:18} {pg1[k]:10} {pu[k]:15}  (string????)')

    print('')

#
#  load parameter classes
#
def analyze_params():
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

    print_params_by_group(parAvg, pg, pu, tstr='\n     Average Parameters by Group')


    print('Parameter Changes vs Set5')

    loFricParms = [ pdicts[1],  pdicts[6],  pdicts[7]  ]
    hiFricParms = [ pdicts[0],  pdicts[2],  pdicts[3], pdicts[4],  pdicts[5],  pdicts[8]  ]

    nscale = 20
    pchanges = compare_param_lists(loFricParms, hiFricParms,nscale)

    print_params_by_group(pdicts[6], pg, pu, pg2=pdicts[2], pch=pchanges, tstr='\n    Compare Lo to Hi fric')

    #pchanges = compare_params(parAvg, hiPressureParms)
    #print_params_by_group(parAvg, pg, pch=pchanges, tstr='\n    Compare Avg w HI fric')


if __name__ == '__main__':

    test = True

    if not test:
        analyze_params()

    else:

        print('TESTING:  change_str()')
        # test change string "bar graph"
        diffs = [-12, -6, 0, 6, 12]
        npts = 20
        ch_str = change_str(diffs, npts)
        print('--------------')
        print(ch_str)
        print('--------------')
        assert ch_str == '*0.1........X.....X....._.....X.....X......*10', 'failed change string test'

        print('                             change_str()    PASS\n')


        f = np.sqrt(10)
        d1 = {'A':1.2,   'B':2.2}
        d2 = {'A':1.2/f, 'B':2.2*f}   # changes from d1
        d3 = {'A':2.2,   'B':3.2}
        d4 = {'A':2.2/(2*f), 'B':3.2*(2*f)}   # changes from d1
        d5 = {'A':3.2,   'B':4.2}
        d6 = {'A':3.2/(3*f), 'B':4.2*(3*f)}   # changes from d1
        dg = {'A':'Reel', 'B':'Tube'}
        pu = {'A':'m' , 'B':'kg'}

        print('TESTING:  compare_param_dicts')
        ch_d = compare_param_dicts(d1, d2, npts)
        print('--------------')
        print(ch_d)
        print('--------------')
        assert ch_d['A'] == '*0.1..........X........._..................*10'
        assert ch_d['B'] == '*0.1...................._.........X........*10'

        print('                             compare_param_dicts()    PASS\n')


        print('\nTESTING: compare_param_lists')
        ld1 = [d1,d3,d5]
        ld2 = [d2,d4,d6]
        ch_ds = compare_param_lists(ld1,ld2,npts)

        print('--------------')
        print(ch_ds)
        print('--------------')



        #pchanges = compare_param_lists([d1],[d2],npts)
        #print_params_by_group(d1,d2,pchanges,tstr='\n\n       TEST \n\n')


