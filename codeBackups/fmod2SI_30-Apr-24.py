import numpy as np
import matplotlib.pyplot as plt

pd = {}  # let's keep track of our params and their units
pu = {}
def print_param_table(pd,pu):
    print('\n Parameter Table ')
    for k in pd.keys():
        print(f'{k:12}  {pd[k]:8.4E}  {pu[k]:15}')
    print('')

def print_param_table2(pd,pdo,pu):  # highlight changes in params due to hacks
    print('\n Parameter Table ')
    for k in pd.keys():
        f = ' '
        if pd[k] != pdo[k]:
            f = '*'
        print(f'{k:12} {f} {pd[k]:8.4E}  {pu[k]:15}')
    print('')


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
Pa_per_PSI  = 6894.7

m3_per_Liter =  1.0 / Liter_per_m3  # m3
pd['moles_per_m3'] = moles_per_m3
pu['moles_per_m3'] = 'moles / m3 of Air'

#  Data
Patmosphere = 14.5 * Pa_per_PSI
Psource_SIu = Patmosphere + 3.0 * Pa_per_PSI # pascals

LLine_meas = 0.6   # L/ min / kPa  (page 5 col 2)
LLine_SIu = LLine_meas * m3_per_Liter /(sec_per_min * Pa_per_kPa)  # m3 / sec / Pa

Fopen = LLine_SIu * Psource_SIu

Rsource_SIu = Psource_SIu / Fopen

pd['Psource_SIu'] = Psource_SIu
pu['Psource_SIu'] = 'Pascals'
pd['LLine_SIu']  = LLine_SIu
pu['LLine_SIu']  = 'm3/sec / Pascal'
pd['Rsource_SIu'] = Rsource_SIu
pu['Rsource_SIu'] = 'Pa/m3/sec'


# Reel Params
J = 1.0E-04 * [4.67, 5.10, 5.64][1]  # kg/m2: Reel rot. Inertia  Table I
rReel_SIu = 10.9 / 1000  # meters
#  Reel friction (Coulomb type)
Tau_coulomb = [0.0029, 0.0174, 0.0694][1] # N m :   Table II

#Vhousing_m3 = 0.5 * m3_per_Liter  # 0.5L in m3

r_hous = 0.075 / 2 # 75mm diam.
px_p_m = (441-322)/(5*2.54/100.0)
l_hous = (328-18)/px_p_m       # m - measured off photo
print(f'Housing r: {r_hous:5.3f} Lh: {l_hous:5.3f}')
Vhousing_m3 = np.pi * r_hous**2 * l_hous

pd['Vhousing_m3'] = Vhousing_m3
pu['Vhousing_m3'] = 'm3'


#Tube
Diam_MM = 25
Radius_MM = Diam_MM/2
Radius_m = Radius_MM/1000.0
#area_mm2 = np.pi * Radius_MM**2
area_m2  = np.pi * Radius_m**2

Kdrag = .3 # N / m / sec2  (pure speculation!!)

pd['Kdrag'] = Kdrag
pu['Kdrag'] = 'N / m2 / sec'
pd['area_m2'] = area_m2
pu['area_m2'] = 'm2'

#
#  FLOW load line for source
#

# Presure intercept:
Pintercept = Psource_SIu
pd['Pintercept'] = Pintercept
pu['Pintercept'] = 'Pascals'

# Velocity intercept:
FlowMax_SIu = (Psource_SIu - Patmosphere)/Rsource_SIu
Fintercept = FlowMax_SIu
pd['Fintercept'] =  FlowMax_SIu
pu['Fintercept'] = 'm3/sec'



# source VELOCITY load line
Vintercept = FlowMax_SIu / area_m2
pd['Vintercept'] =  Vintercept
pu['Vintercept'] = 'm/sec'


x1 = 0
y1 = Pintercept
x2 = Fintercept
y2 = Patmosphere
#
#def P2V_SourceLL_SIu(P):  # get evert Vel from Pressure via Source LL
    #dydxLL = (y2-y1) / (x2-x1)
    #return dydxLL

print('Pressure and velocity intercept points (SI units)')
print(Pintercept, 'Pascals')
print(Vintercept, 'meters/sec')


#
#    Eversion Thresholds
#
PBA_static = Patmosphere + 2.5 * Pa_per_PSI   #  Pascals  Static Breakaway Pressure
PHalt_dyn  = Patmosphere + 1.5 * Pa_per_PSI  #  Pascals   dynamic stopping pressure
P1 = PHalt_dyn
P2 = PBA_static

dP1dt = 2500/5
dP2dt = -dP1dt


pd['PBA_static'] = PBA_static
pu['PBA_static'] = 'Pascals'
pd['PHalt_dyn'] = PHalt_dyn
pu['PHalt_dyn'] = 'Pascals'

# States
PRESSURE_TEST = 2
GROWING = 1
STUCK = 0

# Constant Eversion forces
f_Brake_SIu = 1.0  # Newtons (positive opposes motion)
fEverForce_SIu = 1.0 # Newtons (could be prop to L!)


# PV = NRT !  # need SI units here!
# https://en.wikipedia.org/wiki/Ideal_gas_law
#  P  Pascals
#  V  m3
#  N  moles
#  T  degK
#  R = 8.314 (ideal gas constant)

# Physical Constants
R = 8.314 # wikipedia Gas Constant
T= 295.4  # 72F in deg K

print_param_table(pd,pu)  # before hacks
pd_orig = pd.copy()

#################################################
#
##  Initial Conditions
#

P = Patmosphere  # 1 atmosphere in Pa
#

# Starting STATE
state = STUCK
#state = PRESSURE_TEST

# state var initial values

N = P*Vhousing_m3/(R*T)    #   N = PV/(RT), ideal gas law
L = 0       # ET length (m)
Ldot  = 0
Lddot = 0

# declare data storage
p = []      # Pressure (Pa)
vol = []    # Volume (m3)
f = []      # SourceFlow (m3/sec)
fet = []    # Flow out due to everting
l = []      # L (m)
ldot = []   # velocity (m/sec)

#######################################################3
#
# model the burst and stick behavior
#
##   Simulation Parameters
dt = 0.005  # sec

TMIN = 0.0
TMAX = 5.0

# Time range for plots
PltTMIN = TMIN
PltTMAX = TMAX

#PltTMIN = 0
#PltTMAX = 0.5


##########################################################
#
#  Parameter hacks needed to match results (Fig 3)
#
#   TESTING  HACK

#Vhousing_m3 *= 0.8
#pd['Vhousing_m3'] = Vhousing_m3

time = np.arange(TMIN,TMAX,dt)

for t in time:

    # Eq 1
    Voltot = Vhousing_m3 + L*area_m2

    # Eq 2
    P = N*R*T / Voltot

    # Eq 5.1
    Fcoulomb = Tau_coulomb / rReel_SIu

    # Eq 5.2
    Fdrag = L * Kdrag * 2.0 * Ldot  # material speed is 2x Ldot

    # Eq 6

    if state == GROWING:
        #  0.5 due to eversion kinematics
        Lddot = max(0, (rReel_SIu**2/J) * ( 0.5 * P*area_m2 -Fdrag -Fcoulomb ) )
    else: # smooth slow down
        alpha = 2000 * dt  # empirical fit
        Lddot = -1 * max(0, alpha * Ldot)
        #Ldot  = 0


    # Eq 2.5
    sourceFlowIn = (Psource_SIu - P)/Rsource_SIu

    # Eq 3
    #ETGrowthFlow = Ldot*area_m2    # vol is growing so N is growing
    Ndot =+ (sourceFlowIn) * moles_per_m3

    # Eq 3.5
    #  thresholds seem to grow closer together with eversion
    P1 += dP1dt * dt
    P2 += dP2dt * dt

    # Eq 4
    if P > P2:
        state = GROWING
    if P < P1:
        state = STUCK

    # Eq 7
    Ldot += Lddot * dt
    L    += Ldot  * dt
    N    += Ndot  * dt



    # Record data
    # tube length
    l.append(L)
    # tube eversion velocity
    ldot.append(Ldot)
    # flow from source
    f.append(sourceFlowIn)
    fet.append(Ldot * area_m2)  #flow out into the tube
    # Pressure
    p.append(P) # Pa
    # Volume
    vol.append(Voltot)




# Create a figure with subplots
fig, axs = plt.subplots(3, 2, figsize=(8, 10))

# Plot 1
#axs[0].plot(vel, p)
# "ideal" loadline
axs[0,0].plot([x1,x2], [y1,y2])  # x1,..y2 defined above

vel1 = []
pr1 = []


if state==PRESSURE_TEST:
    axs[0].text(0.5, 10000, 'PRESSURE_TEST (no eversion)')
#trajectory:
axs[0,0].plot(fet,p) #pressure is Pa
axs[0,0].legend(['Source Load Line',  'Trajectory'])
axs[0,0].set_xlabel('Source Flow (m3/sec)')
axs[0,0].set_ylabel('Pressure (Pa)')
axs[0,0].set_xlim(0,0.001)
plt.sca(axs[0,0])
#plt.xticks([0.0001, 0.0002, 0.0003])
axs[0,0].set_ylim(Patmosphere, Psource_SIu)

# Plot 2   # PRESSURE
axs[1,0].plot(time, p)
axs[1,0].legend(['Pressure (Pa)' ])
axs[1,0].set_xlabel('Time (sec)')
axs[1,0].set_ylabel('Pressure (Pa)')
axs[1,0].set_xlim(PltTMIN, PltTMAX)
axs[1,0].set_ylim(Patmosphere, Psource_SIu)
axs[1,0].plot([0, PltTMAX] , [PHalt_dyn, P1] ,'g', linestyle='dashed' )
axs[1,0].plot([0, PltTMAX] , [PBA_static, P2] ,'r', linestyle='dashed' )


# Plot 3
axs[2,0].plot(time, vol )
axs[2,0].legend(['Vol (m3)'])
axs[2,0].set_xlabel('Time (sec)')
axs[2,0].set_ylabel('Volume (m3)')
axs[2,0].set_xlim(PltTMIN, PltTMAX)
axs[2,0].set_ylim( 0.0010, 0.003 )

axs[0,1].plot(time, f)
axs[0,1].set_xlabel('Time (sec)')
axs[0,1].set_ylabel('Flow (m3/sec)')
axs[0,1].set_xlim(PltTMIN, PltTMAX)
axs[1,1].set_ylim(      0, 0.010 )

axs[1,1].plot(time, l)
axs[1,1].set_xlabel('Time (sec)')
axs[1,1].set_ylabel('Tube Len (m)')
axs[1,1].set_xlim(PltTMIN, PltTMAX)
axs[1,1].set_ylim(      0, 0.800 )

axs[2,1].plot(time, ldot)
axs[2,1].set_xlabel('Time (sec)')
axs[2,1].set_ylabel('Tube Velocity (m/sec)')
axs[2,1].set_xlim(PltTMIN, PltTMAX)
axs[2,1].set_ylim(      0, 0.8 )


# Adjust layout
plt.tight_layout()

print_param_table2(pd,pd_orig, pu)

plt.show()

