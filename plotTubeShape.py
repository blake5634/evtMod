import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import simulation as sim
import glob
import et_lib as et
from et_lib import error
import sys



paramDir = 'evtParams/'
defaultParamName = 'InitialParams.txt'

pu = et.loadPUnits(paramDir, 'units_'+defaultParamName)

paramDir += '2Comp/'

args = sys.argv

#get parameter files of interest
pfiles = list(glob.glob(paramDir + '/' + "Set*.txt"))

#pfiles.sort(key=lambda x: os.path.getmtime(x),reverse=True) # newest first
if len(pfiles) == 0:
    cto.error('No brl_data files found')

for i,f in enumerate(pfiles):
    print(f'{i:3d}  {f}')


inp = input('Select file by number: ')
try:
    inp = int(inp)
except:
    et.error('bad input: ', inp)

if (inp>=len(pfiles) or inp <=0):
    et.error('illegal file name index')

pd = et.loadParams('',pfiles[inp])


#############################################################################
#
#   Plot the Tube Shape
#

et.plot_tube_shape(pd)

plt.show()
