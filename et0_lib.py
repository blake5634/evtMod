import numpy as np
import glob,os

import matplotlib.pyplot as plt


# unit constants
sec_per_min = 60
kPa_per_Pa = 0.001
Pa_per_kPa = 1.0/kPa_per_Pa
min_per_sec = 1/sec_per_min
Gal_per_Liter = 0.2642
Liter_per_Gal = 3.7854
Liter_per_m3  = 1000.0
m3_per_Liter =  1.0 / 1000.0  # m3
Liter_per_mm3 = Liter_per_m3 / 1000**3
Gal_per_mm3 = Liter_per_mm3 *Gal_per_Liter
mm3_per_Gal = 1.0/Gal_per_mm3
MM3perLiter = 1.0 / Liter_per_mm3
# Ideal Gas Law  https://pressbooks.uiowa.edu/clonedbook/chapter/the-ideal-gas-law/
m3_per_mole = 0.02241 # m3/mol of Air
moles_per_m3 = 1.0/m3_per_mole
Pa_per_PSI  = 6894.76
atmos_Pa = 14.5 * Pa_per_PSI
Patmosphere = 101325.0    # Pascals

def get_files():
    dataDirNames = ['dataAndyMay24']
    #  Collect existing data files for overplotting
    files = []
    mdfiles = []  # markdown files for metadata
    for datadir in dataDirNames:
            tfiles = list(glob.glob(datadir + '/' + "*"))
            tfiles.sort(key=lambda x: os.path.getmtime(x),reverse=True) # newest first
            if len(tfiles) ==0:
                cto.error('No brl_data files found')
            dirfiles = []
            dirmdfiles = []
            for f in tfiles:
                #print('found: ',f)
                if '.csv' in f :
                    dirfiles.append(f)
                elif  '_meta.json' in f:
                    dirmdfiles.append(f)

            # dict.fromkeys better than set because preserves order (thanks google)
            #filenameroots = list(dict.fromkeys(filenameroots)) # elim dupes
            files += dirfiles
            mdfiles += dirmdfiles
    return files, mdfiles

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


    ed = {}
    idfric = hd['fric_level']
    idIntertia = hd['inertia_level']
    idxflow = hd['flow_lpm_1']
    idxpress = hd['press_psi_1']
    idxL = hd['reel_lin_mm_1']
    idxLdot =  hd['reel_lin_mm_2']

    #ed['f'] =
    # time index = 0
    ed['time'] = data[1:,0]       # col 0, skip header
    ed['dtexp'] = data[2,0]-data[1,0] # delta t
    d2 = data[1:,0][:-1]
    ed['dtexp'] = np.mean(data[2:,0]-d2)
    ed['flow'] = data[1:,idxflow]
    ed['P'] = data[1:,idxpress]
    ed['L'] = data[1:,idxL]
    ed['Ldot'] = data[1:,idxLdot]


    #  ed['f'] = flow data
    #  ed['time'] = time axis
    #   etc.
    return ed


def convert_units(ed):
    ed['flow'] *= m3_per_Liter / sec_per_min  # L/min --> m3/sec
    ed['P'] = atmos_Pa + ed['P']*Pa_per_PSI
    ed['L'] *= 0.001  # m/mm
    ed['Ldot'] *= 0.001  # mm/sec -> m/sec
    return ed

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

#def FD_loss1(sim,exper):    # primitive!  should look for peak shifts instead
    #npts = np.max([len(sim), len(exper)])
    #fsim = np.absolute(np.fft.rfft(sim, npts))
    #fexp = np.absolute(np.fft.rfft(exper,npts))
    #sse=0.0
    #for i,msim in enumerate(fsim):
        #sse += (msim-fexp[i])**2
    #sse /= len(sim)
    #rmse = np.sqrt(sse)
    #return rmse

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
        print('Freq Domain Loss2: ', FD_loss2(data1,data3,dt,test=True), 'Hz')
        plt.show()
