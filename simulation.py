#
#
#   Simulation code for Everting tubes
#
import numpy as np
import et_lib as et

# States
statenames = ['STUCK','GROWING', 'PRESSURE_TEST']
PRESSURE_TEST = 2
GROWING = 1
STUCK = 0



# some force / torque physics
def F_drag(L,Ldot, pd):
    #return Kdrag * v + K2drag / max(L, 0.100)
    # 2 comes from eversion kinematics
    return  2 * (pd['Kdrag'] +  pd['K2drag'] * L) * Ldot


def Tau_brake(Tc, th_dot, pd):
    # eliminate braking torque near zero velocity
    if abs(th_dot) < pd['th_dot_minCoulomb']:
        return Tc  # HACK!  Disable this feature
    else:
        return Tc


####################################################################
#
#   Perform Simulation

def simulate(pd,uc,tmin=0,tmax=8.0):
    error_count = 0
    #########
    #
    ##  Initial Conditions
    #

    P = pd['Patmosphere']  # 1 atmosphere in Pa
    #

    # Starting STATE
    state = STUCK
    #state = PRESSURE_TEST

    # state var initial values

    N = P*pd['Vhousing_m3']/pd['RT']    #   N = PV/(RT), ideal gas law
    L = 0       # ET length (m)
    Lc = 0  #  length of tube crumpled in the housing
    Ldot  = 0
    Lddot = 0
    theta = 0
    th_dot = 0
    th_ddot = 0

    # declare data storage
    p = []      # Pressure (Pa)
    vol = []    # Volume (m3)
    f = []      # SourceFlow (m3/sec)  (lower case f = flow)
    fet = []    # Flow out due to everting
    l = []      # L (m)
    lc = []     # Len of crumple material (m)
    ldot = []   # velocity (m/sec)
    pstt = []     # stuck threshold
    pbat = []    # break away thresh

    F_e = []   # eversion force   (upper case F = force)
    F_c = []   # coulomb force
    F_d = []   # drag force
    F_j = []   # reel inertia force

    #######################################################3
    #
    # model the burst and stick behavior
    #
    ##   Simulation Parameters
    dt = pd['dt']  # sec

    TMIN = tmin
    TMAX = tmax

    # Time range for plots
    PltTMIN = TMIN
    PltTMAX = TMAX

    # Detect thresholds too close
    CONVERGED = False

    time = np.arange(TMIN,TMAX,dt)
    tdata = []
    for t in time:
        tdata.append(t)  # in case of error exit
        #breakpoint()

        # Eq 1
        Voltot = pd['Vhousing_m3'] + L*pd['ET_area']

        # Eq 2
        P = N*pd['RT'] / Voltot

        # Eq 5.1
        Fcoulomb = Tau_brake(pd['Tau_coulomb'], th_dot, pd)
        F_c.append(Fcoulomb)


        #
        # F_dr
        F_dr = F_drag(L,Ldot,pd)     # material speed is 2x Ldot
        F_d.append(-F_drag(L,Ldot,pd))
        F_j.append(-Lddot*pd['J']/pd['rReel']**2)

        # F_ever

        F_ever = max(0, 0.5 * P*pd['ET_area'])
        F_e.append(F_ever)   # eversion force

        #  Crumple length
        Lc = max(0, pd['rReel']*theta - L)
        lc.append(Lc)


        #tubing mass
        Mt =  (L+0.3) *  pd['et_MPM']
        #Bt = 10000.0     # N/m/sec  crumple damping
        #Mt_epsilon = 100 * pd['rReelpd']*pd['et_MPM']

        #Lc = 0
        if   state == GROWING and Lc > 0: # CRUMPLE zone active
            Lddot = (F_ever - F_drag(L,Ldot,pd) - Fcoulomb) / Mt
            th_ddot = -1 * Tau_brake(pd['Tau_coulomb'], th_dot,pd)/pd['J']     # damping out overspin

        elif state == GROWING and Lc <=  0:  # no crumple: TAUGHT
            Lddot = (F_ever - F_drag(L,Ldot,pd) - Fcoulomb ) / (Mt + (pd['J']/pd['rReel']**2) )
            th_ddot = Lddot/pd['rReel']   # kinematic relation

        elif state == STUCK:
            # smooth slow down
            alpha = 2000 * pd['dt']  # empirical fit
            Lddot =   -1 * max(0, alpha * Ldot)

            th_ddot = -1 * Tau_brake(pd['Tau_coulomb'],th_dot,pd)/pd['J']

        else:
            error('Invalid State: '+state)

        # Eq 2.5
        sourceFlowIn = (pd['Psource_SIu'] - P)/pd['Rsource_SIu']

        # Eq 3
        #ETGrowthFlow = Ldot*area_m2    # vol is growing so N is growing
        Ndot = sourceFlowIn * uc['moles_per_m3']  # w

        # Eq 3.5
        #  thresholds seem to grow closer together with eversion
        # lower threshold (Halting)
        P1 = pd['PHalt_dyn']  + pd['Threshold Taper'] * L
        # upper threshold (break-away)
        P2 = pd['PBA_static'] - pd['Threshold Taper'] * L

        PdTol = 0.0075
        midpoint = (pd['PBA_static'] + pd['PHalt_dyn'])/2.0
        PdMax = PdTol * midpoint
        if abs(P2-P1) < PdMax:  # prevent crossover (optional)
            CONVERGED=True  # lock this state
        if CONVERGED:
            P1 = (1.0-0.5*PdTol) * midpoint
            P2 = (1.0+0.5*PdTol) * midpoint
        pstt.append(P1)
        pbat.append(P2)

        # if used in Force break mode
        F1 = 18  + pd['dF1dL'] * L * 0.5 # Newtons
        F2 = 26  - pd['dF1dL'] * L * 0.5

        PRESSURE_BREAK = True  # False = Force breakaway
        # Eq 4

        if t > 4.570 and t < 4.600:
            print('')
            print(f't: {t:8.4f}   State: {statenames[state]}     Pressure: {P:8.1f}')
            print(f'SourceFlow: {sourceFlowIn:8.3f}   Rsource_SIu: {pd["Rsource_SIu"]:8.3E}')
            #Lddot = (F_ever - F_drag(L,Ldot,pd) - Fcoulomb ) / (Mt + (pd['J']/pd['rReel']**2) )
            print(f'F_ever {F_ever:8.3f}  F_drag: {F_drag(L,Ldot,pd):8.3f}  Fcoul: {Fcoulomb:8.3f}  Mt: {Mt:8.3f}')
            print(f'Ndot: {Ndot:8.3f}   N: {N:8.3f},  Voltot: {Voltot:8.3f}')
            print(f'Lddot: {Lddot:8.3f}  Ldot: {Ldot:8.3f}, Mt: {Mt:8.3f} L: {L:8.3f} ')
        #if t > 4.5820-pd['dt']:
            #breakpoint()

            #print('')
            #print('Thresholds Crossed!')
            #print(f'midpoint: {midpoint:8.3f}')
            #print(f't: {t:8.4f}   State: {statenames[state]}     Pressure: {P:8.1f}')
            #print(f'Thresholds: Break: {P2:8.1f}    Halt: {P1:8.1f}')
            #print(f'SourceFlow: {sourceFlowIn:8.3f}   Rsource_SIu: {pd["Rsource_SIu"]:8.3E}')
            #print(f'Ndot: {Ndot:8.3f}   N: {N:8.3f},  Voltot: {Voltot:8.3f}')
            #print(f'Lddot: {Lddot:8.3f}  Ldot: {Ldot:8.3f},  L: {L:8.3f} ')

        if PRESSURE_BREAK:
            if P > P2:
                state = GROWING
            if P < P1:
                state = STUCK
        else:                 # switch states based on Force instead of pressure
            if F_ever > F2:
                state = GROWING
            if F_ever < F1:
                state = STUCK

        #if t> 4.50:
                #print('')
                #print(f't: {t:8.4f}   State: {statenames[state]}     Pressure: {P:8.1f}')
                #print(f'Thresholds: Break: {P2:8.1f}    Halt: {P1:8.1f}')
                #print(f'SourceFlow: {sourceFlowIn:8.3f}   Rsource_SIu: {pd["Rsource_SIu"]:8.3E}')
                #print(f'Ndot: {Ndot:8.3f}   N: {N:8.3f},  Voltot: {Voltot:8.3f}')
                #print(f'Lddot: {Lddot:8.3f}  Ldot: {Ldot:8.3f},  L: {L:8.3f} ')
                #breakpoint()

        # Integrate the state variables.
        Ldot   += Lddot   * dt
        L      += Ldot    * dt
        th_dot += th_ddot * dt
        theta  += th_dot  * dt
        N      += Ndot    * dt

        # Record data
        # tube length
        l.append(L)
        # tube eversion velocity
        ldot.append(Ldot)
        # flow from source
        f.append(sourceFlowIn)
        fet.append(Ldot * pd['ET_area'])  #flow out into the tube
        # Pressure
        p.append(P) # Pa
        # Volume
        vol.append(Voltot)

        # results to be returned:
            ## declare data storage
            #p = []      # Pressure (Pa)
            #vol = []    # Volume (m3)
            #f = []      # SourceFlow (m3/sec)  (lower case f = flow)
            #fet = []    # Flow out due to everting
            #l = []      # L (m)
            #lc = []     # Len of crumple material (m)
            #ldot = []   # velocity (m/sec)

            #F_e = []   # eversion force   (upper case F = force)
            #F_c = []   # coulomb force
            #F_d = []   # drag force
            #F_j = []   # reel inertia force
        if error_count > 10:
            break
    return (tdata,l,lc,ldot,f, fet, p, pbat, pstt, vol, F_e,F_c,F_d,F_j)  # return the simulation results
