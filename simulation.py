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
    #   Kdrag is viscous drag,   K2drag is lenght-proportional viscous drag

    #  increase drag at end of tubing (max eversion length)
    #     i.e. when the tubing runs out and is stuck to reel.
    #     0.5*F comes from eversion kinematics
    Fsteadystate = 0.5 *  pd['Psource_SIu']*np.pi*et.Ret(L,pd)**2
    F_limit = (L/pd['Lmax'])**17 * Fsteadystate  # a very steep increase
    return max(everDrag,  F_limit)

####################################################################
#
#   Perform Simulation

def simulate(pd,uc,tmin=0,tmax=8.0):
    error_count = 0
    # normally we are doing two compartment
    ONECOMPARTMENT = False
    try:
        if pd['Compartments'] < 2:
            ONECOMPARTMENT = True
    except:
        pass
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
    PC2 = PC1    # we will use PC2 only in two-compartment
    #
    et.Vet(0.0,pd)  # initialize everting tube computation (L param not used)
    fint = 0
    fT = 0
    L = et.tubeLinit       # ET nominal non-zero init length (m)
    Lc = 0  #  length of tube crumpled in the housing
    Ldot  = 0
    Lddot = 0

    try:
        COMP1 = pd['COMP1']
    except:
        # define what is compartment 1 (also determines comp2)
        #COMP1 = 'housing'       # comp2 = et vol
        COMP1 = 'supply_tubing'  # comp2 = housing+et vol

    print('\n Compartment 1 represents: ', COMP1, '\n')
    pd['COMP1'] = COMP1  # store it.

    VsupTube = 0.1*pd['Vhousing_m3']
    if ONECOMPARTMENT:
        N1 = PC1*(pd['Vhousing_m3']+et.Vet.et_vol)/pd['RT']
        N2 = 0  # not used
    else:
        if COMP1 == 'housing':
            N1 = PC1*pd['Vhousing_m3']/pd['RT']
            N2 = PC2* et.Vet.et_vol/pd['RT']
        elif COMP1 == 'supply_tubing':
            N1 = PC1*VsupTube/pd['RT']
            N2 = PC2*(pd['Vhousing_m3']+et.Vet.et_vol)/pd['RT']
        else:
            et.error('Illegal Compartment One def: ', COMP1)

    N1dot = 0
    N2dot = 0
    theta = 0
    th_dot = 0
    th_ddot = 0
    time = np.arange(TMIN,TMAX,dt)
    tdata = []

    #
    #  constrain flow resistances (experimental 5-Aug-24)
    #
    r_source_C1, r_C1_C2 = et.constrainR(pd)
    it = 0
    for t in time:
        it +=1
        tdata.append(t)  # in case of error exit
        #breakpoint()


        # Solve Pressures
        #2Compartment

        # Experimental: Comp1 = air tubing, Comp2 = housing+ET
        #VC1 = 0.2*pd['Vhousing_m3'] # estimate tubing vol
        #VC2 = (pd['Vhousing_m3'] + et.Vet.et_vol)
        # original model

        if COMP1 == 'housing':
            VC1 = pd['Vhousing_m3']  # fixed
            VC2 = et.Vet.et_vol      # changes w/ length
        elif COMP1 == 'supply_tubing':
            VC1 = VsupTube           # defined in ICs above
            VC2 = pd['Vhousing_m3'] + et.Vet.et_vol
        else:
            et.error('Illegal Compartment definition: ', COMP1)

        if ONECOMPARTMENT:
            PC1 = N1 * pd['RT'] / (VC1 + VC2)
            PC2 = 0
        else:  # Two compartment
            PC1 = N1*pd['RT'] / VC1
            PC2 = N2*pd['RT'] / VC2

        # source flow
        fsource = (pd['Psource_SIu'] - PC1) / r_source_C1
        # other flows for 2comp
        if not ONECOMPARTMENT:
            # Two Comp
            fT = (PC1-PC2)/r_C1_C2
            fint = fsource - fT
            N1dot = fint * uc['moles_per_m3']
            N2dot = fT   * uc['moles_per_m3']
        else:  # ONE Comp.
            N1dot = fsource * uc['moles_per_m3']
            N2dot = 0.0  # not used

        #
        #  Experimental: sensed pressure middle of Rsource
        ##
        #beta = 0.5
        #Rs1 = beta*pd['Rsource_SIu']
        #Rs2 = (1.0-beta)*pd['Rsource_SIu']
        #Psens.append(pd['Psource_SIu'] - fsource*Rs1)

        # Eq 5.1
        # this is coulomb fric due to reel brake (TAUT state only)
        Fcoulomb = pd['Tau_coulomb'] / pd['rReel']
        F_c.append(Fcoulomb) # not scaled because constant w.r.t. |vel|

        #
        # F_drag
        #
        F_dr = F_drag(L,Ldot,pd)     # material speed is 2x Ldot
        F_d.append(F_dr)

        # F_ever
        if ONECOMPARTMENT:
            drivepress = PC1
        else:
            drivepress = PC2
        F_ever = max(0, 0.5 * (drivepress-pd['Patmosphere']) * np.pi*et.Ret(L,pd)**2) # 1/2 pres. applied to tip
        F_e.append(F_ever)   # eversion force

        # Net Force
        #F_j.append(-2*Lddot*pd['J']/pd['rReel']**2)
        F_j.append(F_ever-F_dr-Fcoulomb)  # now 'Net Force'


        #tubing mass depends on length
        Mt =  (L+0.4) *  pd['et_MPM']

        #
        #    Determine Stuck or Growing
        #

        PRESSURE_BREAK = True  # False = Force breakaway
        # Eq 4

        # Eq 3.5
        #  thresholds seem to grow closer together with eversion
        # lower threshold (Halting)
        Pth1 = pd['PHalt_dyn']  + pd['Threshold Taper'] * L
        # upper threshold (break-away)
        Pth2 = pd['PBA_static'] - pd['Threshold Taper'] * L

        PdTol = 0.0075  # minimum gap between Pth1 and Pth2
                        # normalized to midpoint
        midpoint = (pd['PBA_static'] + pd['PHalt_dyn'])/2.0
        PdMax = PdTol * midpoint
        if abs(Pth2-Pth1) < PdMax:  # prevent crossover (optional)
            CONVERGED=True  # lock this state
        if CONVERGED:
            Pth1 = (1.0-0.5*PdTol) * midpoint
            Pth2 = (1.0+0.5*PdTol) * midpoint

        ####  Option: disable Threshold Taper Completely
        #Pth1 = pd['PHalt_dyn']
        #Pth2 = pd['PBA_static']

        pstt.append(Pth1)   # store thresholds
        pbat.append(Pth2)

        ## if used in Force break mode
        F1 = 18  + pd['dF1dL'] * L * 0.5 # Newtons
        F2 = 26  - pd['dF1dL'] * L * 0.5

        if PRESSURE_BREAK:
            # drivepress is tube (PC2) OR housing (PC1) pressure(!)
            if drivepress > Pth2:
                state = GROWING
                # don't change state2 until slack is taken up
            if drivepress < Pth1:
                state = STUCK
                #state2 = SLACK
        else:                 # switch states based on ET Force instead of pressure
            if F_ever > F2:
                state = GROWING
            if F_ever < F1:
                state = STUCK

        #  Crumple length
        Lc = pd['rReel']*theta - 2*L  # reel must supply 2x length
        if Lc < -0.050005 and it%500==0:  # for some reason frequently close to -5mm
            print(f't:{t:5.3f} somethings wrong with Lc: {Lc:5.3f}')
        Lc = max(0.0, Lc)  # can't go negative
        lc.append(Lc)

        LCmin = 0.001
        if Lc > LCmin:
            state2 = SLACK  # stateseq == 3
        else:
            state2 = TAUT   # stateseq == 2

        #
        #     Compute accelerations via dynamic equations
        #
             # COMB STATE 0
        if state == GROWING and state2==TAUT:  # no crumple: TAUT
            #state2 = TAUT  # stateseq == 0
            #if t>1.0:
                #breakpoint()
            Ldot = max(0.0,Ldot)  # enforce that L only grows
            Lddot = (F_ever - F_dr - Fcoulomb ) / (2*(Mt + (pd['J']/(pd['rReel']**2) ) ) )
            # corrected 28-Jul: (2x)
            th_ddot = 2*Lddot/pd['rReel']   # kinematic relation
            #
            # further constrain theta
            th_dot = 2*Ldot/pd['rReel']
            theta = 2*L/pd['rReel']

             # COMB STATE 1
        elif  state == GROWING and state2==SLACK: # CRUMPLE zone active
            # note pulled tubing accelerates at 2*Lddot
            #   no Fcoulomb in TAUT state (29-Jul)
            Lddot = (F_ever - F_dr) / (Mt*2)
            th_ddot = -1 * pd['Tau_coulomb']/pd['J']   # damping out overspin

        elif state == STUCK:
            alpha = 100000 * pd['dt']  # empirical fit
            Lddot =   -1 * max(0, alpha * Ldot)
            #disable smooth slow down
            #Lddot=0
            #Ldot = 0
            if abs (th_dot) > 0.0005:
                th_ddot = -1 * pd['Tau_coulomb'] / pd['J']
            else:
                th_ddot = 0.0

        else:
            error('Invalid State: '+state)

        prevState2 = state2



        # Record data
        l.append(L)                      # tube length
        ldot.append(Ldot)                # tube eversion velocity
        f.append(fsource)                # flow from source
        ft.append(fT)                    # tubing airflow
        state_seq.append(stateSymbol(state,state2))  # record state as 0-3 int
        th.append(theta)
        thdot.append(th_dot)
        # Pressure
        Phous.append(PC1)                # Housing pressure Pa (or 1comp pressure)
        Ptube.append(PC2)                # Tubing Pressure
        vol1.append(pd['Vhousing_m3'])   # right now housing vol is constant(!)
        vol2.append(et.Vet.et_vol)       # query the Etube volume

        # Integrate the state variables.
        Ldot   += Lddot   * dt
        Ldot   =  max(0,Ldot)    # Ldot can never go negative
        LdotMAX = 0.5  # m/sec
        #Ldot   =  min(Ldot, LdotMAX)  # or exceed a max
        L      += Ldot    * dt
        if state2 != TAUT:  # don't integrate when constrained
            th_dot += th_ddot * dt
            th_dot =  max(0,th_dot)  # th_dot can never go neg.
            theta  += th_dot  * dt

        # 2Compartment state integration
        N1     += N1dot   * dt
        N2     += N2dot   * dt

        et.Vet(L,pd)             #integrate et volume

    return (tdata,th, thdot, state_seq, l,lc,ldot,f, ft, Phous, Ptube, pbat, pstt, vol1, vol2, F_e, F_c, F_d, F_j)  # return the simulation results
