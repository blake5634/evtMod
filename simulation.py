#
#
#   Simulation code for Everting tubes
#
import numpy as np
import et_lib as et

# States
PRESSURE_TEST = 2
GROWING = 1
STUCK = 0


####################################################################
#
#   Perform Simulation

def simulate(pd,uc,tmin=0,tmax=8.0):
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

    time = np.arange(TMIN,TMAX,dt)

    for t in time:
        #breakpoint()

        # Eq 1
        Voltot = pd['Vhousing_m3'] + L*pd['ET_area']

        # Eq 2
        P = N*pd['RT'] / Voltot

        # Eq 5.1
        Fcoulomb = et.Tau_brake(pd['Tau_coulomb'], th_dot, pd)
        F_c.append(Fcoulomb)


        #
        # F_dr
        F_dr = et.F_drag(L,Ldot,pd)     # material speed is 2x Ldot
        F_d.append(-et.F_drag(L,Ldot,pd))
        F_j.append(-Lddot*pd['J']/pd['rReel']**2)

        # F_ever

        F_ever = max(0, 0.5 * P*pd['ET_area'])
        F_e.append(F_ever)   # eversion force

        #  Crumple length
        Lc = max(0, pd['rReel']*theta - L)
        lc.append(Lc)


        Mt =  (L+0.3) * pd['rReel']*pd['et_MPM']
        #Bt = 10000.0     # N/m/sec  crumple damping
        #Mt_epsilon = 100 * pd['rReelpd']*pd['et_MPM']

        #Lc = 0
        if   state == GROWING and Lc > 0: # CRUMPLE zone active
            Lddot = (F_ever - et.F_drag(L,Ldot,pd) - Fcoulomb) / Mt
            th_ddot = -1 * et.Tau_brake(pd['Tau_coulomb'], th_dot,pd)/pd['J']     # damping out overspin

        elif state == GROWING and Lc <=  0:  # no crumple: TAUGHT
            Lddot = (F_ever - et.F_drag(L,Ldot,pd) - Fcoulomb ) / (Mt + (pd['J']/pd['rReel']**2) )
            th_ddot = Lddot/pd['rReel']   # kinematic relation

        elif state == STUCK:
            # smooth slow down
            alpha = 2000 * dt  # empirical fit
            Lddot = -1 * max(0, alpha * Ldot)

            th_ddot = -1* et.Tau_brake(pd['Tau_coulomb'],th_dot,pd)/pd['J']

        else:
            error('Invalid State: '+state)


        # Eq 2.5
        sourceFlowIn = (pd['Psource_SIu'] - P)/pd['Rsource_SIu']

        # Eq 3
        #ETGrowthFlow = Ldot*area_m2    # vol is growing so N is growing
        Ndot = sourceFlowIn * uc['moles_per_m3']  # w



        # Eq 3.5
        #  thresholds seem to grow closer together with eversion
        P1 = pd['PHalt_dyn']  + pd['Threshold Taper'] * L
        P2 = pd['PBA_static'] - pd['Threshold Taper'] * L

        # if used in Force break mode
        F1 = 18  + pd['dF1dL'] * L * 0.5 # Newtons
        F2 = 26  - pd['dF1dL'] * L * 0.5

        PRESSURE_BREAK = True  # False = Force breakaway

        # Eq 4
        if PRESSURE_BREAK:
            if P > P2:
                state = GROWING
            if P < P1:
                state = STUCK
        else:
            if F_ever > F2:
                state = GROWING
            if F_ever < F1:
                state = STUCK

        print(f'SourceFlow: {sourceFlowIn:8.3f}   Rsource_SIu: {pd["Rsource_SIu"]:8.3E}')
        print(f'Ndot: {Ndot:8.3f}   N: {N:8.3f},  Voltot: {Voltot:8.3f}')
        print(f'Ldot: {Ldot:8.3f},  L: {L:8.3f}  P: {P:8.3f}')
        breakpoint()

        # Integrate the state variables.
        Ldot   += Lddot * dt
        L      += Ldot  * dt
        th_dot += th_ddot * dt
        theta  += th_dot * dt
        N      += Ndot  * dt

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

    return (time,l,lc,ldot,f, fet, p, vol, F_e,F_c,F_d,F_j)  # return the simulation results
