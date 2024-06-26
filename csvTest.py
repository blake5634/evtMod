
import numpy as np
import matplotlib.pyplot as plt

# unit constants
sec_per_min = 60
kPa_per_Pa = 0.001
Pa_per_kPa = 1.0/kPa_per_Pa
min_per_sec = 1/sec_per_min
Gal_per_Liter = 0.2642
Liter_per_Gal = 3.7854
Liter_per_m3  = 1000.0
Liter_per_mm3 = Liter_per_m3 / 1000**3
Gal_per_mm3 = Liter_per_mm3 *Gal_per_Liter
mm3_per_Gal = 1.0/Gal_per_mm3
MM3perLiter = 1.0 / Liter_per_mm3
# Ideal Gas Law  https://pressbooks.uiowa.edu/clonedbook/chapter/the-ideal-gas-law/
m3_per_mole = 0.02241 # m3/mol of Air
moles_per_m3 = 1.0/m3_per_mole
Pa_per_PSI  = 6894.76
m3_per_Liter =  1.0 / Liter_per_m3  # m3
atmos_Pa = 14.5 * Pa_per_PSI

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

    data = np.genfromtxt(fn, delimiter=',')

    #fill ed[] dictionary


    ed = {}
    idfric = hd['fric_level']
    idIntertia = hd['inertia_level']
    idxflow = hd['flow_lpm_1']
    idxpress = hd['press_psi_1']
    idxL = hd['reel_lin_mm_1']

    #ed['f'] =
    ed['time'] = data[:,0]
    ed['flow'] = data[:,idxflow]
    ed['P'] = data[:,idxpress]
    ed['L'] = data[:,idxL]
    #  ed['f'] = flow data
    #  ed['time'] = time axis
    #   etc.
    return ed


def convert_units(ed):
    ed['flow'] *= m3_per_Liter
    ed['P'] = atmos_Pa + ed['P']*Pa_per_PSI
    ed['L'] *= 0.001  # m/mm
    return ed

if __name__=='__main__':

    dirname = 'dataAndyMay24'

    #open file
    files = ['eversion_flow-hi-inr_hi-fric_tube-1_trial-2.csv',
'eversion_flow-hi-inr_lo-fric_tube-1_trial-1.csv', 'eversion_flow-hi-inr_lo-fric_tube-1_trial-2.csv',
'eversion_flow-hi-inr_lo-fric_tube-1_trial-3.csv' ]
    for i,fn in enumerate(files):
        print(f'{i:5}', fn)
    n = int(input('Select file by number: '))
    if (n>=len(files) or len(files) <0):
        print('illegal file name index')
        quit()
    fn = dirname + '/' + files[n]
    ed = get_data_from_AL_csv(fn)
    ed = convert_units(ed)

    plt.figure()
    plt.plot(ed['flow'],ed['P'])


    plt.figure() #FLOW   m3/sec
    plt.plot(ed['time'],ed['flow'])

    plt.figure() # Pressure  Pa
    plt.plot(ed['time'],ed['P'])

    plt.figure() # Length m
    plt.plot(ed['time'], ed['L'])
    plt.show()
