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
    everDrag =   2 * (pd['Kdrag'] +  pd['K2drag'] * L) * Ldot
    #return everDrag

    #  increase drag at end of tubing (max eversion length)
    Fsteadystate = pd['Psource_SIu']*np.pi*et.Ret(L,pd)**2
    F_limit = (L/pd['Lmax'])**7 * Fsteadystate
    #return max(everDrag, min(Fsteadystate, F_limit))
    return max(everDrag,  F_limit)

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


    # declare data storage
    #2Compartment
    Phous = []      # Housing Pressure (Pa)
    Ptube = []      # tubing tip pressure
    vol1 = []    # Housing Volume (m3)
    vol2 = []   # tubing vol
    f = []      # SourceFlow (m3/sec)  (lower case f = flow)
    fint = []   # internal flow to housing
    ft = []     # Flow out due to everting
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

    #########################################################
    #
    #   Initial Conditions
    #

    # Starting STATE
    state = STUCK
    #state = PRESSURE_TEST

    # state var initial values (initial conditions)
    PC1 = pd['Patmosphere']  # 1 atmosphere in Pa
    PC2 = PC1
    #
    et.Vet(0.0,pd)  # initialize everting tube computation (L param not used)
    N1 = PC1*pd['Vhousing_m3']/pd['RT']    #   N = PV/(RT), ideal gas law
    N2 = PC2*et.Vet.et_vol / pd['RT']
    #print(f'Initial Cond:  PC1: {PC1}, PC2: {PC2}')
    #print(f'Initial Cond:  N1: {N1},   N2: {N2}')
    #x = input('pause: (cr)')
    fint = 0
    fT = 0
    L = et.tubeLinit       # ET length (m)
    Lc = 0  #  length of tube crumpled in the housing
    Ldot  = 0
    Lddot = 0
    N1 = PC1*pd['Vhousing_m3']/pd['RT']
    N2 = PC2* et.Vet.et_vol/pd['RT']
    N1dot = 0
    N2dot = 0
    theta = 0
    th_dot = 0
    th_ddot = 0
    time = np.arange(TMIN,TMAX,dt)
    tdata = []
    for t in time:
        tdata.append(t)  # in case of error exit
        #breakpoint()


        # Solve Pressures
        #2Compartment
        PC1 = N1*pd['RT'] / pd['Vhousing_m3']
        PC2 = N2*pd['RT'] / et.Vet.et_vol

        # Eq 2.5
        fsource = (pd['Psource_SIu'] - PC1)/pd['Rsource_SIu']
        # if True:
        #     if PC2>=PC1:
        #     #if False:
        #         fT = et.Vet.et_vol * Ldot # Ldot = constant 0.05e.g.
        #     else:
        #         fT = (PC1-PC2)/et.ResET(L,pd)
        #         #fT = (PC1-PC2)/pd['Rsource_SIu']
        # else:
        #     fT = 0.0000002  # HACK
        #
        fT = (PC1-PC2)/et.ResET(L,pd)
        fint = fsource - fT

        # 2Compartment:
        N1dot = fint * uc['moles_per_m3']
        N2dot = fT   * uc['moles_per_m3']

        #print(f'{t:6.4f}  Ptube: {PC2}  Phouse: {PC1} etVol: {et.Vet.et_vol:5.2E}')
        #print(f"        tube flow: {(PC1-PC2)/et.ResET(L,pd):10.5f}, source flow: { (pd['Psource_SIu'] - PC1)/pd['Rsource_SIu']:10.5f}")
        #print(f'        N2dot: {N2dot:10.5f}      N1dot: {N1dot:10.5f}')
        #print(f'           N2: {N2:10.5f}         N1: {N1:10.5f}')

        #breakpoint()

        # Eq 5.1
        Fcoulomb = Tau_brake(pd['Tau_coulomb'], th_dot, pd)
        F_c.append(Fcoulomb)

        #
        # F_dr
        F_dr = F_drag(L,Ldot,pd)     # material speed is 2x Ldot
        F_d.append(-F_drag(L,Ldot,pd))
        F_j.append(-Lddot*pd['J']/pd['rReel']**2)

        # F_ever
        # 2Compartment
        F_ever = max(0, 0.5 * (PC2-pd['Patmosphere']) * np.pi*et.Ret(L,pd)**2) # C2 pres. applied to tip
        F_e.append(F_ever)   # eversion force

        #  Crumple length
        Lc = max(0, pd['rReel']*theta - L)
        lc.append(Lc)

        #tubing mass depends on length
        Mt =  (L+0.3) *  pd['et_MPM']
        #Bt = 10000.0     # N/m/sec  crumple damping
        #Mt_epsilon = 100 * pd['rReelpd']*pd['et_MPM']

        #Lc = 0
        if   state == GROWING and Lc > 0: # CRUMPLE zone active
            # note pulled tubing accelerates at Lddot/2
            Lddot = (F_ever - F_drag(L,Ldot,pd) - Fcoulomb) / (Mt/2)
            th_ddot = -1 * Tau_brake(pd['Tau_coulomb'], th_dot,pd)/pd['J']     # damping out overspin

        elif state == GROWING and Lc <=  0.01:  # no crumple: TAUGHT
            Lddot = (F_ever - F_drag(L,Ldot,pd) - Fcoulomb ) / (Mt/2 + (pd['J']/pd['rReel']**2) )
            th_ddot = Lddot/pd['rReel']   # kinematic relation

        elif state == STUCK:
            # smooth slow down
            alpha = 20 / pd['dt']  # empirical fit
            Lddot =   -1 * max(0, alpha * Ldot)

            th_ddot = -1 * Tau_brake(pd['Tau_coulomb'],th_dot,pd)/pd['J']

        else:
            error('Invalid State: '+state)

        # Eq 3.5
        #  thresholds seem to grow closer together with eversion
        # lower threshold (Halting)
        Pth1 = pd['PHalt_dyn'] # + pd['Threshold Taper'] * L
        # upper threshold (break-away)
        Pth2 = pd['PBA_static'] #- pd['Threshold Taper'] * L

        PdTol = 0.0075
        midpoint = (pd['PBA_static'] + pd['PHalt_dyn'])/2.0
        PdMax = PdTol * midpoint
        if abs(Pth2-Pth1) < PdMax:  # prevent crossover (optional)
            CONVERGED=True  # lock this state
        if CONVERGED:
            Pth1 = (1.0-0.5*PdTol) * midpoint
            Pth2 = (1.0+0.5*PdTol) * midpoint
        pstt.append(Pth1)
        pbat.append(Pth2)

        ## if used in Force break mode
        F1 = 18  + pd['dF1dL'] * L * 0.5 # Newtons
        F2 = 26  - pd['dF1dL'] * L * 0.5

        PRESSURE_BREAK = True  # False = Force breakaway
        # Eq 4

        if PRESSURE_BREAK:
            if PC2 > Pth2:       # use tube (PC2) OR housing (PC1) pressure(!)
                state = GROWING
            if PC2 < Pth1:
                state = STUCK
        else:                 # switch states based on ET Force instead of pressure
            if F_ever > F2:
                state = GROWING
            if F_ever < F1:
                state = STUCK

        # Record data
        l.append(L)                      # tube length
        ldot.append(Ldot)                # tube eversion velocity
        f.append(fsource)                # flow from source

        #ft.append(Ldot * np.pi*et.Ret(L,pd)**2)  #flow out into the tube

        # 2Compartment
        ft.append(fT)                    # tubing airflow
        # Pressure
        Phous.append(PC1)                # Housing pressure Pa
        Ptube.append(PC2)                # Tubing Pressure
        vol1.append(pd['Vhousing_m3'])   # right now housing vol is constant(!)
        vol2.append(et.Vet.et_vol)       # query the Etube volume

        # Integrate the state variables.
        Ldot   += Lddot   * dt
        Ldot   =  max(0,Ldot)    # Ldot can never go negative
        L      += Ldot    * dt
        th_dot += th_ddot * dt
        theta  += th_dot  * dt
        # 2Compartment state integration
        N1     += N1dot   * dt
        N2     += N2dot   * dt
        et.Vet(L,pd)             #integrate et volume

    return (tdata,l,lc,ldot,f, ft, Phous, Ptube, pbat, pstt, vol1, vol2, F_e, F_c, F_d, F_j)  # return the simulation results
