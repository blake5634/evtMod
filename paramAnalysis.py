import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import simulation as sim
import et_lib as et
#from et_lib import error
import glob
import sys

def avg_dict_list(dl):
    avgD = {}
    for k in dl[0].keys():
        avgD[k] = 0.0
    # caution: assumes all d's have same keys
    for d in dl:
        for k in d.keys():
            if type(d[k]) == type(5.7):
                avgD[k] += d[k]
            else:
                avgD[k] = None
    for k in d.keys():
        if avgD[k] is not None:
            avgD[k] /= len(dl)
    return avgD

def compare_param_list_avgs(pdl1, pdl2,nscale):
    avd1 = {}
    avd2 = {}
    n1=0
    for k in pdl1[0].keys():
        avd1[k] = 0.0
        avd2[k] = 0.0
    for d in pdl1:
        for k in d.keys():
            try:
                avd1[k] += d[k]
            except:
                pass
        n1 += 1
    n2=0
    for d in pdl2:
        for k in d.keys():
            try:
                avd2[k] += d[k]
            except:
                pass
        n2 += 1

    for k in pdl1[0].keys():
        avd1[k] /= n1
        avd2[k] /= n2
    #print('lists: ', pdl1, pdl2)
    #print('averaged dicts: \n', avd1, '\n', avd2)
    changes = compare_param_dicts(avd1, avd2, nscale)
    return changes

# compare two param dicts
def compare_param_dicts(pd1, pd2, nscale):
    changes = {}
    #pd1 and pd2 are param dicts
    for pk in pd1.keys():
        try:
            diff= [ pdiff(pd1[pk],pd2[pk],nscale) ]
            #print('comparing param dicts: ', pd1[pk], pd2[pk], pdiff(pd1[pk],pd2[pk],nscale))
            changes[pk] = change_str(diff,nscale)
        except:
            changes[pk] = 'None'
    return changes


def change_str(diffs,nscale):
    left = '*0.1 '
    right = ' *10'
    scale = ['.']*((nscale)*2 + 1)
    scale[nscale] = '_'
    NoDiff = True
    for d in diffs:
        if abs(d) > 0.000025:
            NoDiff = False
            dIndex = int(len(scale)/2 + 0.5 + d)
            dIndex = max(0, min( nscale*2 , dIndex))
            scale[dIndex] = 'X'
    if NoDiff:
        return ''
    else:
        return left+str(''.join(scale))+right

def pdiff(a,b,nscale):
    n = (nscale*(np.log10(b/a)))
    #n = max(n, -nscale)
    #n = min(n, nscale)
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
    pfiles.sort() # Set0, Set1, ...
    print('Testing: \n files in order: ')
    for i in range(len(pfiles)):
        print(i, pfiles[i])

    if len(pfiles) == 0:
        et.error('No eversion parameter files found')
    pdicts = []
    for i,fn in enumerate(pfiles):
        print(f'{i:5}  Loading: {fn}')
        pdicts.append(et.loadParams('',fn))

    # remove all params in group "Admin" and "text"
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
                if pg[par] == 'text':    # check against parameter groups
                    #print('             deleting key: ',par)
                    del pd[par]
            except:
                pass

    pnames = pdicts[0].keys()  # assume they're the same

    parAvg = {}
    parMaxMin = {}

    for k in pnames:
        parAvg[k] = 0.0
        parMaxMin[k] = (-999999.999E+50, 999999.9999E+50)
    n = 0
    for f in pdicts:
        n += 1
        for k in pnames:
            try:
                val = float(f[k])
                parAvg[k] += val
                mx, mi = parMaxMin[k]
                if val > mx:
                    mx = val
                if val < mi:
                    mi = val
                parMaxMin[k] = (mx,mi)
            except:
                parAvg[k] = 0.0

    print('Average parameters')
    for k in pnames:
        parAvg[k] /= n

    print_params_by_group(parAvg, pg, pu, tstr='\n     Average Parameters by Group')

    print('Range for each Numerical Param')
    for k in pnames:
        mx, mi = parMaxMin[k]
        print(f'{k:20}  {mi} --- {mx}')

    print('')

    print('Parameter Changes vs Set5')

    loFricParms = [ pdicts[1],  pdicts[6],  pdicts[7]  ]
    hiFricParms = [ pdicts[0],  pdicts[2],  pdicts[3], pdicts[4],  pdicts[5],  pdicts[8]  ]

    print('testing: ')
    print('   TauFric d6:', pdicts[6]['Tau_coulomb'])
    print('   TauFric d2:', pdicts[2]['Tau_coulomb'])
    nscale = 20

    pchanges = compare_param_list_avgs(loFricParms, hiFricParms,nscale)

    print_params_by_group(pdicts[6], pg, pu, pg2=pdicts[2], pch=pchanges, tstr='\n    Compare Lo to Hi fric')


    # compare trail #1 w trial #3 on each tube
    T1parms = [pdicts[1], pdicts[4], pdicts[6]]
    T3parms = [pdicts[3], pdicts[7]]

    pchanges = compare_param_list_avgs(loFricParms, hiFricParms,nscale)

    print_params_by_group(pdicts[1], pg, pu, pg2=pdicts[3], pch=pchanges, tstr='\n    Compare trial 1 to trial 3')

    avgTrial1 = avg_dict_list(T1parms)
    avgTrial2 = avg_dict_list(T3parms)

    pchanges = compare_param_dicts(avgTrial1, avgTrial2,20)
    print_params_by_group(avgTrial1, pg, p2=avgTrial2, pch=pchanges, tstr='\n    Compare Avg Trial 1 to Avg trial 2')


if __name__ == '__main__':

    test = False
    #test = True

    if not test:

        args = sys.argv
        if len(args) == 3: # compare two param sets
            set1 = int(args[1])
            set2 = int(args[2])
        else:
            analyze_params()
            quit()

        paramDir = 'evtParams/'
        paramDirL2 = 'evtParams/1Comp/'
        paramGFile = 'ParamGroups.txt'   # conceptual groupings of parms

        pg = et.loadDict(paramDir, paramGFile)
        pu = et.loadPUnits(paramDir, 'units_InitialParams.txt')

        #files, mdfiles = et.get_files()
        #print('Parameter File Analysis')
        pfiles = list(glob.glob(paramDirL2 + "Set*.txt"))
        pfiles.sort() # Set0, Set1, ...
        print('Pfiles:',pfiles)
        for pf in pfiles:
            print('testing: ',pf)
            if f'Set{str(set1)}Params.txt' in pf:
                file1 = pf
            if f'Set{str(set2)}Params.txt' in pf:
                file2 = pf
        pd1 = et.loadParams('', file1)
        pd2 = et.loadParams('', file2)

        pchanges = compare_param_dicts(pd1,pd2,18)

        print_params_by_group(pd1, pg, pu, pg2=pd2, pch=pchanges, tstr=f'\n {file1} to {file2} ')


    else:

        print('TESTING: avg_dict_list()')
        d1 = {'A':1.0,   'B':2.0}
        d2 = {'A':1.1,   'B':2.2}   # changes from d1
        d3 = {'A':0.9,   'B':1.8}
        avgd = avg_dict_list([d1,d2,d3])
        assert avgd['A'] == 1.0
        assert avgd['B'] == 2.0

        print('\n                          avg_dict_list()             PASS\n')


        print('TESTING:  change_str()')
        # test change string "bar graph"
        diffs = [-12, -6, 0, 6, 12]
        npts = 20
        ch_str = change_str(diffs, npts)
        #print('--------------')
        #print(diffs)
        #print(ch_str)
        #print('--------------')
        assert ch_str == '*0.1 .........X.....X...._......X.....X....... *10', 'failed change string test'

        print('\n                             change_str()             PASS\n')


        f = np.sqrt(10)
        #f = 20
        d1 = {'A':1.2,   'B':2.2}
        d2 = {'A':1.2/f, 'B':2.2*f}   # changes from d1
        d3 = {'A':2.2,   'B':3.2}
        d4 = {'A':2.2/(2*f), 'B':3.2*(2*f)}   # changes from d1
        d5 = {'A':3.2,   'B':4.2}
        d6 = {'A':3.2/(3*f), 'B':4.2*(3*f)}   # changes from d1
        dg = {'A':'Reel', 'B':'Tube'}
        pu = {'A':'m' , 'B':'kg'}

        print('TESTING:  compare_param_dicts')
        ch_ds = compare_param_dicts(d1, d2, npts)
        #print('--------------')
        #print('d1,d2: ', d1, d2)
        #print(ch_ds)
        #print('--------------')
        assert ch_ds['A'] == '*0.1 ...........X........_.................... *10'
        assert ch_ds['B'] == '*0.1 ...................._..........X......... *10'

        print('\n                             compare_param_dicts()    PASS\n')


        print('\nTESTING: compare_param_lists')
        ld1 = [d1,d3,d5]
        ld2 = [d2,d4,d6]
        ch_ds = compare_param_list_avgs(ld1,ld2,npts)

        #print('--------------')
        #print(ch_ds)
        #print('--------------')
        assert ch_ds['A'] == '*0.1 .....X.............._.................... *10'
        assert ch_ds['B'] == '*0.1 ...................._................X... *10'

        print('\n                             compare_param_lists()    PASS\n')


        #pchanges = compare_param_lists([d1],[d2],npts)
        #print_params_by_group(d1,d2,pchanges,tstr='\n\n       TEST \n\n')


