import numpy as np
import glob,os

import matplotlib.pyplot as plt


#
#   2Compartment model functions
#

# flow resistance of everting tube as a function of Length
def ResET(L,pd):
    # Rmin = 1.0E4 # 10% of Rsource
    # r =  pd['ET_Res_per_m'] * L / pd['Lmax'] + Rmin
    #print(f'                    ET Resistance: L:{L:7.4f} r: {r:6.3E} Rsource:{pd["Rsource_SIu"]:6.3E}')
    rhack = pd['ET_Res_ratio']*pd['Rsource_SIu']
    #print(f'      FIXED         ET Resistance: rhack: {rhack:6.3E}')
    return rhack

 # volume of ET
tubeLinit = 0.0025
def Vet(L,pd):
     # LP =  previous L value
    if Vet.LP < 0:  # init condition flag (in loop 'cause need pd)
        Vet.LP = tubeLinit  # initial non-zero length
        Vet.et_vol = Vet.LP * Ret(L,pd)**2 * np.pi
        return
    else:
        Vet.dL =  max(0.0, L-Vet.LP)
        Vet.et_dVol_dL = np.pi* Ret(L,pd)**2    #dV/dL
        Vet.et_vol +=   Vet.et_dVol_dL * Vet.dL
        Vet.LP = L
        return Vet.et_vol
# initialize 'static' variable (func attribs)
# flag for initial condition setting
Vet.LP= -1    # flag for init. cond. setting

#
#    functions to support variable diameter with length ( V(L) )
#
# radius of everting tube as function of Length
def Ret(L,pd):
    mode = pd['ET_RofL_mode']
    if mode == 'constant':
        return pd['ET_radius']
    elif mode == 'box':
        return etr_box(L,pd)
    elif mode == 'constrict':
        return etr_constrict(L,pd)
    elif mode == 'gap':
        return etr_gap(L,pd)
    elif mode == 'ramp':
        return etr_ramp(L,pd)
    else:
        error('R of (L): unknown radius mode: '+pd['ET_RofL_mode'] )

def etr_box(L,pd):
    if L < 0.15:
        return pd['ET_radius']
    elif L >= 0.15 and L < 0.3:
        return 1.25*pd['ET_radius']
    else:
        return pd['ET_radius']

def etr_gap(L,pd):
    if L < 0.15:
        return pd['ET_radius']
    elif L >= 0.15 and L < 0.25:
        return 0.25*pd['ET_radius']
    else:
        return pd['ET_radius']

def etr_constrict(L,pd):
    if L > 0.20:
        return 0.5*pd['ET_radius']
    else:
        return pd['ET_radius']

def etr_ramp(L,pd):
    drdL = (0.5 * pd['ET_radius']) / (pd['Lmax'] - 0.2)
    if L < 0.2:
        return pd['ET_radius']
    else:
        return pd['ET_radius'] + (L-0.2)*drdL

def plot_tube_shape(pd):
    l = np.linspace(0.0, pd['Lmax'], 150)
    fig, axs = plt.subplots(1,1, figsize=(10,2))
    y1 = []
    y2 = []
    for L in l:
        y1.append(1000.0*Ret(L,pd))
        y2.append(-1000.0*Ret(L,pd))
    axs.plot(l,y1,l,y2)
    axs.set_title('Tube Geometry: '+pd['ET_RofL_mode'])
    axs.set_xlabel('Length (m)')
    axs.set_ylim([-20,20])
    axs.set_ylabel('tube width (mm)')
    # Adjust layout
    plt.tight_layout()
    return    # must then call plot.show()


# Unit Conversions (only here to generate initial param units file)

uc = {
    "sec_per_min": 60,
    "kPa_per_Pa": 0.001,
    "Pa_per_kPa": 1.0 / 0.001,
    "min_per_sec": 1 / 60,
    "Gal_per_Liter": 0.2642,
    "Liter_per_Gal": 3.7854,
    "Liter_per_m3": 1000.0,
    "Liter_per_mm3": 1000.0 / 1000**3,
    "Gal_per_mm3": (1000.0 / 1000**3) * 0.2642,
    "Gal_per_m3" : 264.17205,
    "mm3_per_Gal": 1.0 / ((1000.0 / 1000**3) * 0.2642),
    "m3_per_Gal" : 1.0 / 264.17205,
    "m3_per_Liter": 1.0 / 1000.0,  # m3
    "MM3perLiter": 1.0 / (1000.0 / 1000**3), # Ideal Gas Law  https://pressbooks.uiowa.edu/clonedbook/chapter/the-ideal-gas-law/
    "m3_per_mole": 0.02241,  # m3/mol of Air
    "moles_per_m3": 1.0 / 0.02241,
    "Pa_per_PSI": 6894.76,
    "atmos_Pa": 14.5 * 6894.76,
}

def get_param_files():
    paramDirNames = ['evtParams']
    #  Collect existing data files for overplotting
    files = []
    dirfiles = []
    for pdir in paramDirNames:
            # these files *include* the path relative to working dir.
            pfiles = list(glob.glob(pdir + '/' + "*arams.txt"))
            if len(pfiles) == 0:
                cto.error('No brl_data files found')
            dirfiles = []
            for f in pfiles:
                dirfiles.append(f)

    print ('get_param_files: ', dirfiles)
    return dirfiles

def get_files():   # load the data files
    dataDirNames = ['dataAndyMay24']
    #  Collect existing data files for overplotting
    files = []
    mdfiles = []  # markdown files for metadata
    for datadir in dataDirNames:
            # these files *include* the path relative to working dir.
            #  e.g. 'datadir/filename.csv'
            pfiles = list(glob.glob(datadir + '/' + "*"))
            pfiles.sort(key=lambda x: os.path.getmtime(x),reverse=True) # newest first
            if len(pfiles) == 0:
                cto.error('No brl_data files found')
            dirfiles = []
            dirmdfiles = []
            for f in pfiles:
                #print('found: ',f)
                if f[-4:] == '.csv':          # data
                    dirfiles.append(f)
                elif f[-10:]==  '_meta.json': # metadata
                    dirmdfiles.append(f)

            files += dirfiles
            mdfiles += dirmdfiles
    print ('get_files: ', files)
    return files, mdfiles

#
#   Inertia value and Friction Value is encoded in some filenames
#
def get_inr_fric_from_name(fn):
    I = 1
    F = 1
    tags = ['hi-inr','lo-iner','hi-fric','lo-fric']
    if 'hi-inr' in fn:
        I = 2
        print(' Inertia: HIGH')
    if 'lo-inr' in fn:
        I = 0
        print(' Inertia: LO')
    if 'hi-fric' in fn:
        F = 2
        print('  Fric: HIGH')
    if 'lo-fric' in fn:
        F = 0
        print('  Fric: LO')
    return I, F

#  read data from Andy Lewis' data files
def get_data_from_AL_csv(fn):
    #read data
    header = [
    'time_sec', 'tof_pos_mm_1', 'tof_pos_mm_2', 'tof_pos_mm_3',
    'reel_radians_1', 'reel_radians_2', 'reel_radians_3',
    'reel_lin_mm_1', 'reel_lin_mm_2', 'reel_lin_mm_3',
    'press_psi_1', 'press_psi_2', 'press_psi_3',
    'flow_lpm_1', 'flow_lpm_2', 'flow_lpm_3',
    'fric_level', 'inertia_level', 'hash'
    ]
    hd = {}
    for i in range(len(header)):
        hd[header[i]] = i

    # read entire CSV into an array - neat (ignores header too).
    data = np.genfromtxt(fn, delimiter=',')

    #fill ed[] dictionary
    # experiment data
    ed = {}
    idfric = hd['fric_level']
    idIntertia = hd['inertia_level']
    idxflow = hd['flow_lpm_1']
    idxpress = hd['press_psi_1']
    idxL = hd['reel_lin_mm_1']
    idxLdot =  hd['reel_lin_mm_2']

    # time index = 0
    ed['time'] = data[1:,0]       # col 0, skip header
    #ed['dtexp'] = data[2,0]-data[1,0] # delta t
    d2 = data[1:,0][:-1]
    ed['dtexp'] = np.mean(data[2:,0]-d2) # average dt
    ed['flow'] = data[1:,idxflow]
    ed['P'] = data[1:,idxpress]
    ed['L'] = data[1:,idxL]
    ed['Ldot'] = data[1:,idxLdot]


    #  ed['f'] = flow data
    #  ed['time'] = time axis
    #   etc.
    return ed

#
#  do some unit conversions in the exper data to match SI units
def convert_units(ed,uc):
    ed['flow'] *= uc['m3_per_Liter'] / uc['sec_per_min']  # L/min --> m3/sec
    ed['P'] = uc['atmos_Pa'] + ed['P']*uc['Pa_per_PSI']
    ed['L'] *= 0.001  # m/mm
    ed['Ldot'] *= 0.001  # mm/sec -> m/sec
    return ed

#
#   Time domain match betweeen sim and exper
def TD_loss(sim, exper):
    sse = 0.0
    for i,s in enumerate(sim):
        if i<len(exper):
            sse += (s-exper[i])**2
        else:
            break
    sse /= len(sim)
    rmse = np.sqrt(sse)
    return rmse

#
# Frequency domain match between sim and exper
def FD_loss2(sim,exper,dtsim,dtexp,test=False,ptitle='Data vs Simulation'):    # loss = | peak shift |
    decimate =  15
    dtsim *= decimate
    dtexp *= decimate
    sim = sim[::decimate]
    exper = exper[::decimate]
    npts = np.max([len(sim), len(exper)])
    if test:
        print('  FD_loss2():')
        print('    Npts: ',npts)
    fsim = np.absolute(np.fft.rfft(sim, npts))
    fexp = np.absolute(np.fft.rfft(exper,npts))

    if test:
        print('    shapes: ',fsim.shape, fexp.shape)
        # Create simulation output figure with subplots
        fig, axs = plt.subplots(2,1, figsize=(8, 10))
        frangeS = np.arange( len(fsim))/(npts*dtsim)
        frangeE = np.arange( len(fexp))/(npts*dtexp) # convert to Freq.(Hz)
        FminCutoff = 0.250  #Hz
        DCfiltSim = int(FminCutoff*(npts*dtsim))  # clear below 1Hz
        DCfiltExp = int(FminCutoff*(npts*dtexp))  # clear below 1Hz
        fsim[0:DCfiltSim] = 0
        fexp[0:DCfiltSim] = 0
        # freq domain
        axs[0].plot(frangeS, fsim, frangeE,fexp)
        axs[0].set_xlabel('Freq. (Hz)')
        axs[0].set_title('Frequency Domain: '+ptitle)

        #time domain
        simtime = dtsim*np.array(range(  len(sim)))
        exptime = dtexp*np.array(range(len(exper)))
        axs[1].plot(simtime, sim, exptime, exper)
        axs[1].set_xlabel('Time (sec)')
        axs[1].set_ylabel('Amplitude')
        axs[1].set_title('T.D. Comparison: '+ptitle)
    pksim = np.argmax(fsim[3:]) # eliminate big DC peak
    pkexp = np.argmax(fexp[3:])

    if test:
        print('    peak indeces: ', pksim, pkexp)
    rmse = np.sqrt((pksim-pkexp)**2)
    return rmse


def error(str):
    print('Error: ',str)
    print('Quitting')
    quit()

def plot_curve_with_arrows2(x, y, axis, interval, color='b',HW =0.1, arrow_scale=1.0, arrow_color='red'):
    """
    Plot a curve with arrows indicating the direction of the curve.

    Parameters:
        y (array-like): y-coordinates of the data points.
        axis (matplotlib axis): The axis on which to plot the curve and arrows.
        interval (int): Spacing between arrows along the curve.
        arrow_scale (float): Scale factor for the arrow lengths.
        arrow_color (str): Color of the arrows.
    """
    # scale factors for x and y
    #if max(ranges) == 0:
    xsf = 0.03* arrow_scale / np.abs(np.max(x)-np.min(x))
    ysf = 0.03* arrow_scale / np.abs(np.max(y)-np.min(y))
    #else:
        #xsf = 1/np.abs(ranges[1]-ranges[0])
        #ysf = 1/np.abs(ranges[3]-ranges[2])

    # Compute the direction of the curve (tangent vectors)
    dx = np.gradient(x)

    dy = np.gradient(y)
    th = np.arctan2(dy, dx)

    #  get unit vectors
    dx_unit = np.cos(th)
    dy_unit = np.sin(th)

    # Plot the curve
    axis.plot(x, y,color=color)

    xsf = .01
    ysf = .01

    xq = []
    yq = []
    uq = []
    vq = []
    ## Plot arrows along the curve at specific intervals
    for i in range(0, len(x), interval):

        #print('plotting arrow: ', x[i], y[i], dx_unit[i]*xsf, dy_unit[i]*xsf)
        #axis.arrow(x[i], y[i], dx_unit[i] * xsf, dy_unit[i] * ysf)
                #head_width=0.1, head_length=0.1, fc=arrow_color, ec=arrow_color)
        dx = dx_unit[i]
        dy = dy_unit[i]
        #start = (x[i],y[i])
        #end   = (x[i]+dx*xsf, y[i]+dy*ysf)
        #end   = (x[i]+dx_unit[i], y[i]+dy_unit[i])

        w = 1

        xsf=0.002
        ysf=0.002

        # fill the arrow quiver starts
        xq.append(x[i])
        yq.append(y[i])
        uq.append(dx * xsf)
        vq.append(dy * ysf)
    #axis.quiver(xq,yq,uq,vq, scale_units='xy', scale=1, width=0.01)
    axis.quiver(xq,yq,uq,vq, angles='xy',color=color, scale=0.1)

    return




###############################################################################
#
#    Manage System Parameters
#

# set to nominal values
def setup_params():
    pd = {}
    pu = {}


    # Physical Constants
    R = 8.314 # wikipedia Gas Constant
    T= 295.4  # 72F in deg K

    pd['T'] = T
    pu['T'] = 'K'
    pd['RT'] = R*T
    pu['RT']= '??'

    # atmospheric pressure
    Patmosphere = 101325.0,  # Pascals
    pd['Patmosphere'] = 101325.0
    pu['Patmosphere'] = 'Pascals'

    # Reel Params
    J = 1.0E-04 * [4.67, 5.10, 5.64][1]  # kg/m2: Reel rot. Inertia  Table I
    rReel_SIu = 10.9 / 1000  # meters
    pd['J'] = J
    pu['J'] = 'kg/m2'
    pd['rReel'] = rReel_SIu
    pu['rReel'] = 'm'

    #  Reel friction (Coulomb type)
    Tau_coulomb= [ 0.0029, 0.0174, 0.0694][1] # N m :   Table II
    pd['Tau_coulomb'] = Tau_coulomb
    pu['Tau_coulomb'] = 'Nm'

    th_dot_minCoulomb = 0.05 # rad/sec  (no coulomb fric below this vel)
    pd['th_dot_minCoulomb'] = th_dot_minCoulomb
    pu['th_dot_minCoulomb'] = 'rad/sec'

    # housing
    r_hous = 0.075 / 2 # 75mm diam.
    px_p_m = (441-322)/(5*2.54/100.0)
    l_hous = (328-18)/px_p_m       # m - measured off photo
    print(f'Housing r: {r_hous:5.3f} Lh: {l_hous:5.3f}')
    Vhousing_m3 = np.pi * r_hous**2 * l_hous
    pd['Vhousing_m3'] = Vhousing_m3
    pu['Vhousing_m3'] = 'm3'

    #Tube diameter
    pd['ET_RofL_mode'] = 'constant'  # ['constant', 'box', 'ramp', 'constrict' ]
    pu['ET_RofL_mode'] = 'text'
    Diam_MM = 25
    #pd['ET_diam'] = Diam_MM/1000
    #pu['ET_diam'] = 'm'
    pd['ET_radius'] = 0.001* Diam_MM/2
    pu['ET_radius'] = 'm'

    #area_m2  = np.pi * pd['ET_radius']**2
    #pd['ET_area'] = area_m2
    #pu['ET_area'] = 'm2'

    # Tube material
    tubing_mass_per_meter = 0.100 # kg/meter
    pd['et_MPM'] = tubing_mass_per_meter
    pu['et_MPM'] = 'kg/m'


    Kdrag = 8 # N / m / sec2  (pure speculation!!)
    pd['Kdrag'] = Kdrag           # velocity coeff
    pu['Kdrag'] = 'N / m2 / sec'
    K2drag = 8.0/0.8    # N m    (drag increases w/ len)
    pd['K2drag'] = K2drag         # inverse Length coeff
    pu['K2drag'] = 'Nm'

    #
    #  Tube airflow resistance
    #
    pd['ET_Res_per_m'] = 1.0E7      # airflow resistance per meter
    pu['ET_Res_per_m'] = 'Pa/m2/sec'

    #
    #  FLOW load line for source
    #

    #  Source Pressure (P intercept)
    print(pd)
    Psource_SIu = pd['Patmosphere'] + 3.0 * uc['Pa_per_PSI']    # pascals
    pd['Psource_SIu'] = Psource_SIu
    pu['Psource_SIu'] = 'Pascals'

    # Flow intercept
    LLine_meas = 0.6   # L/ min / kPa  (Lewis paper, page 5 col 2)
    LLine_SIu = LLine_meas * uc['m3_per_Liter'] /(uc['sec_per_min'] * uc['Pa_per_kPa'])  # m3 / sec / Pa
    Fopen = LLine_SIu * pd['Psource_SIu']  # open cavity flow (no ev tube)
    pd['Rsource_SIu'] = Psource_SIu / Fopen
    pu['Rsource_SIu'] = 'Pa/m3/sec'
    pd['LLine_SIu']  = LLine_SIu
    pu['LLine_SIu']  = 'm3/sec / Pascal'


    ## these seem not used???
    #FlowMax_SIu = (pd['Psource_SIu'] - pd['Patmosphere'])/pd['Rsource_SIu']
    #pd['Fintercept'] =  FlowMax_SIu
    #pu['Fintercept'] = 'm3/sec'

    ## VELOCITY load line
    #Vintercept = pd['Fintercept'] / area_m2
    #pd['Vintercept'] =  Vintercept
    #pu['Vintercept'] = 'm/sec'


    #
    #    Eversion Thresholds
    #

    PBA_static = pd['Patmosphere'] + 2.0 * uc['Pa_per_PSI']  # Pascals  Static Breakaway Pressure
    pd['PBA_static'] = PBA_static
    pu['PBA_static'] = 'Pascals'

    PHalt_dyn = pd['Patmosphere'] + 1.25 * uc['Pa_per_PSI']  # Pascals   dynamic stopping pressure
    pd['PHalt_dyn'] = PHalt_dyn
    pu['PHalt_dyn'] = 'Pascals'

    #PBA_static = Patmosphere + 2.0  * Pa_per_PSI   #  Pascals  Static Breakaway Pressure
    #PHalt_dyn  = Patmosphere + 1.25 * Pa_per_PSI  #  Pascals   dynamic stopping pressure
    #P1 = PHalt_dyn
    #P2 = PBA_static


    #  TESTING: taper the thresholds together depending time
    dP1dL = 0.5*(pd['PBA_static'] - pd['PHalt_dyn'])/0.8  # pa/m  (0.8m is fully everted)
    dP2dL = -dP1dL
    pd['Threshold Taper'] = dP1dL
    pu['Threshold Taper'] = 'Pa /m'

    # force threshold taper
    pd['dF1dL'] =  0.5*(26-18)/1.1
    pu['dF1dL'] = 'N/m??' #??


    # simulation details
    pd['dt'] = tubeLinit5
    pu['dt'] = 'sec'


    return (pd,pu)

def saveParams(fname, pd):
    saveDict(fname, pd)

def savePUnits(fname, pu):
    saveDict(fname, pu)

def saveUnitConv(fname, uc):
    saveDict(fname, uc)

def saveDict(fname, d):
    f = open(fname, 'w')
    skeys = sorted(d.keys())
    for k in skeys:
        value = d[k]
        if type(value) == type(3.1416):
            print(f'{k:20}: {value:10.7E}',file=f)
        else:
            print(f'{k:20}: {value}',file=f)
    f.close()
    return

def loadParams(dir, fname):
    if len(dir) > 0:
        return loadDict(dir, fname)
    else:
        return loadDict('', fname)

def loadPUnits(dir, fname):
    return loadDict(dir,fname)

def loadUnitConv(dir,fname):
    return loadDict(dir,fname)

def loadPlotRanges(dir, fname):
    return loadDict(dir,)

def loadDict(folder, fname):
    if len(folder) > 0 and folder != '/':
        if '/' == folder[-1]:
            f = open(folder+fname, 'r')
        else:
            f = open(folder+'/'+fname, 'r')
    else:
        f = open(fname, 'r')

    d = {}
    for line in f:
        k,v = line.split(':')
        try:
            d[k.strip()] = float(v)
        except:
            d[k.strip()] = v.strip()
        if '.' not in v: # detect hash that looks like float:  123e4567
            d[k.strip()] = v.strip()  # keep it as string
    return d

def print_param_table(pd,pu):
    print('\n Parameter Table ')
    for k in sorted(pd.keys()):
        val = pd[k]
        #print('        --- ',k,pd[k])
        if type(val) == type(['x','y']):  # list param
            print(f'{k:18} \n    {pd[k]}      {pu[k]:15}')
        elif type(val) == type('x'):
            print(f'{k:18}  {pd[k]:10}  {pu[k]:15})')  # string params
        else: #  ints and floats
            print(f'{k:18}  {pd[k]:8.4E}  {pu[k]:15}')

    print('')

def print_param_table2(pd,pdo,pu):  # highlight changes in params due to hacks
    print('\n Parameter Table ')
    for k in sorted(pd.keys()):
        flg = ' '
        try:
            if pd[k] != pdo[k]:
                flg = '*'  # flag the changes
        except:
            pass

        #print('        --- ',k,pd[k])
        val = pd[k]
        if type(val) == type(['x','y']):  # list param
            print(f'{k:18} \n    {flg} {pd[k]}      {pu[k]}')
        elif type(val) == type('x'):
            print(f'{k:18}  {flg} {pd[k]:10}  {pu[k]}')  # string params
        else: #  ints and floats
            print(f'{k:18}  {flg} {pd[k]:8.4E}  {pu[k]}')

    print('')



if __name__=='__main__':

    if True:    # doing testing

        ###  Test loss functions in et_lib

        dt = 0.005  # sec
        test_data_len = 200
        data1 = np.zeros(test_data_len)
        period1 = int(1.0/20.0/dt)  # 20Hz
        period2 = int(1.0/12.0/dt)  # 12Hz
        data1 = np.sin(np.array(range(test_data_len))*2*np.pi/period1)
        data2 = data1 * 0.95
        data3 = np.sin(np.array(range(test_data_len))*2*np.pi/period2)


        print('Time Domain Loss: ',  TD_loss(data1,data2))
        print('Freq Domain Loss2: ', FD_loss2(data1,data3,dt,dt,test=True), 'Hz')
        plt.show()
