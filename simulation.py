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

TAUT = 3
SLACK = 4



# some force / torque physics
def F_drag(L,Ldot, pd):
    # 2Ldot comes from eversion kinematics
    everDrag =    (pd['Kdrag'] +  pd['K2drag'] * L) * 2 * Ldot
    #
    #   28-Jul:   K2drag seems to have little effect: set to 0 in params?

    #  increase drag at end of tubing (max eversion length)
    #     i.e. when the tubing runs out and is stuck to reel.
    #     0.5*F comes from eversion kinematics
    Fsteadystate = 0.5 *  pd['Psource_SIu']*np.pi*et.Ret(L,pd)**2
    F_limit = (L/pd['Lmax'])**7 * Fsteadystate  # a very steep increase
    return max(everDrag,  F_limit)

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
    state_seq = []  # sequence of states
    th = []     # reel pos
    thdot = []  # reel speed

    def stateSymbol(state, tsk):
        if state==GROWING and tsk==TAUT:
            return 0
        if state==GROWING and tsk==SLACK:
            return 1
        if state==STUCK and tsk==TAUT:
            return 2
        if state==STUCK and tsk==SLACK:
            return 3

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
    state = STUCK   # growing vs stuck
    state2 = TAUT   # spool crumple (SLACK) or TAUT

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
    L = et.tubeLinit       # ET nominal non-zero init length (m)
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

        # Eq 2.5   HACK to fit saturation
        #fsource = min(0.000332, (pd['Psource_SIu'] - PC1)/pd['Rsource_SIu'])
        fsource = (pd['Psource_SIu'] - PC1)/pd['Rsource_SIu']
        # if True:
        #     if PC2>=PC1:
        #     #if False:
        #         fT = et.Vet.et_vol * Ldot # Ldot = constant 0.05e.g.
        #     else:
        #         fT = (PC1-PC2)/et.ResET(L,pd)
        #         #fT = (PC1-PC2)/pd['Rsource_SIu']
        # else:
        #     fT = 0.0000002
        #
        #fT = (PC1-PC2)/et.ResET(L,pd)
        rtube = pd['Rsource_SIu'] * pd['ET_Res_ratio']
        fT = (PC1-PC2)/rtube
        fint = fsource - fT

        # 2Compartment:

        N1dot = fint * uc['moles_per_m3']
        #N1dot = min(fint, 0.000332)* uc['moles_per_m3']
        N2dot = fT   * uc['moles_per_m3']

        #print(f'{t:6.4f}  Ptube: {PC2}  Phouse: {PC1} etVol: {et.Vet.et_vol:5.2E}')
        #print(f"        tube flow: {(PC1-PC2)/et.ResET(L,pd):10.5f}, source flow: { (pd['Psource_SIu'] - PC1)/pd['Rsource_SIu']:10.5f}")
        #print(f'        N2dot: {N2dot:10.5f}      N1dot: {N1dot:10.5f}')
        #print(f'           N2: {N2:10.5f}         N1: {N1:10.5f}')

        #breakpoint()

        # Eq 5.1
        # this is coulomb fric due to reel brake (TAUT state only)
        Fcoulomb = pd['Tau_coulomb'] / pd['rReel']
        F_c.append(Fcoulomb)

        #
        # F_dr
        F_dr = F_drag(L,Ldot,pd)     # material speed is 2x Ldot
        F_d.append(-F_drag(L,Ldot,pd))
        F_j.append(-Lddot*pd['J']/pd['rReel']**2)

        # F_ever
        F_ever = max(0, 0.5 * (PC2-pd['Patmosphere']) * np.pi*et.Ret(L,pd)**2) # 1/2 pres. applied to tip
        F_e.append(F_ever)   # eversion force

        #  Crumple length
        Lc = pd['rReel']*theta - 2*L
        Lc = max(0.0, Lc)  # can't go negative
        #print('.... testing: r*theta, Lc, L: ',pd['rReel']*theta, Lc, 2*L)
        lc.append(Lc)

        #tubing mass depends on length
        Mt =  (L+0.1) *  pd['et_MPM']

        LCmin = 0.001

             # COMB STATE 1
        #if  state == GROWING and Lc >= LCmin: # CRUMPLE zone active
        if  state == GROWING and state2==SLACK: # CRUMPLE zone active
            #state2= SLACK  # stateseq == 1
            # note pulled tubing accelerates at 2*Lddot
            #   no Fcoulomb in TAUT state (29-Jul)
            Lddot = (F_ever - F_drag(L,Ldot,pd)) / (Mt*2)
            th_ddot = -1 * pd['Tau_coulomb']/pd['J']   # damping out overspin

             # COMB STATE 0
        #elif state == GROWING and Lc < LCmin:  # no crumple: TAUT
        elif state == GROWING and state2==TAUT:  # no crumple: TAUT
            #state2 = TAUT  # stateseq == 0
            #if t>1.0:
                #breakpoint()
            Ldot = max(0.0,Ldot)  # enforce that L only grows

            Lddot = (F_ever - F_drag(L,Ldot,pd) - Fcoulomb ) / (2*(Mt + (pd['J']/(pd['rReel']**2))))

            # corrected 28-Jul: (2x)
            th_ddot = 2*Lddot/pd['rReel']   # kinematic relation

        elif state == STUCK:  # can only change state2 in STUCK
            if Lc > LCmin:
                state2 = SLACK  # stateseq == 3
            else:
                state2 = TAUT   # stateseq == 2
            # smooth slow down
            alpha = 100000 * pd['dt']  # empirical fit
            Lddot =   -1 * max(0, alpha * Ldot)
            if abs (th_dot) > 0.0005:
                th_ddot = -1 * pd['Tau_coulomb'] / pd['J']
            else:
                th_ddot = 0.0

        else:
            error('Invalid State: '+state)

        prevState2 = state2

        # Eq 3.5
        #  thresholds seem to grow closer together with eversion
        # lower threshold (Halting)
        Pth1 = pd['PHalt_dyn']  + pd['Threshold Taper'] * L
        # upper threshold (break-away)
        Pth2 = pd['PBA_static'] - pd['Threshold Taper'] * L

        PdTol = 0.0075
        midpoint = (pd['PBA_static'] + pd['PHalt_dyn'])/2.0
        PdMax = PdTol * midpoint
        if abs(Pth2-Pth1) < PdMax:  # prevent crossover (optional)
            CONVERGED=True  # lock this state
        if CONVERGED:
            Pth1 = (1.0-0.5*PdTol) * midpoint
            Pth2 = (1.0+0.5*PdTol) * midpoint

        ###  HACK disable Threshold Taper Completely
        Pth1 = pd['PHalt_dyn']
        Pth2 = pd['PBA_static']

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
        ft.append(fT)                    # tubing airflow
        state_seq.append(stateSymbol(state,state2))  # record state as 0-3 int
        th.append(theta)
        thdot.append(th_dot)
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
        th_dot =  max(0,th_dot)  # th_dot can never go neg.
        theta  += th_dot  * dt

        # 2Compartment state integration
        N1     += N1dot   * dt
        N2     += N2dot   * dt
        et.Vet(L,pd)             #integrate et volume

    return (tdata,th, thdot, state_seq, l,lc,ldot,f, ft, Phous, Ptube, pbat, pstt, vol1, vol2, F_e, F_c, F_d, F_j)  # return the simulation results
