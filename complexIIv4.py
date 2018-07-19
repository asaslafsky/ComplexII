# import the pysb module and all its methods and functions
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
import seaborn as sns
from scipy import stats
from pysb.simulator import StochKitSimulator
from pysb.simulator import ScipyOdeSimulator
from pysb import *
#import random

# Definitions
NUM_SSA_RUNS = 10000 #How many times SSA will be ran

# instantiate a model
Model()

#Monomers: Name, binding site (named according to what it binds to or if the site gets modified=mod), optional (the states the binding site can be in)
Monomer('TNF', ['tnfr'])
Monomer('TNFR', ['tnf', 'traddrip1'])
Monomer('TRADD', ['tnfr', 'rip1'])
#ub=ubiquitinated, db=deubiquitinated, p=phosphorylated
Monomer('RIP1', ['tnfr', 'tradd', 'a20', 'fadd', 'rip3c8', 'mod'], {'mod': ['ub', 'db', 'p']})
Monomer('A20', ['rip1'])
Monomer('FADD', ['rip1', 'c8', 'flip'])
#unmod=unmodified p=phosphorylated
Monomer('RIP3', ['rip1', 'mlkl', 'mod'], {'mod': ['unmod', 'p']})
#i=inactive, p=phosphorylated
Monomer('MLKL', ['rip3', 'mod'], {'mod': ['i', 'p']})
#i=inactive, a=active
Monomer('C8', ['fadd', 'flip', 'bid', 'mod'], {'mod': ['i', 'a']})
Monomer('FLIP', ['c8', 'fadd'])
Monomer('BID', ['c8', 'mod'], {'mod': ['unmod', 'trunc']})

#Parameters: Rates of reactions
Parameter('kc_degrade_TNF', 1.00e-03) #1.00e-03
Parameter('kf_TNF_TNFR_bind', 1e-6) #2.26e-06
Parameter('kr_TNF_TNFR_bind', 1e-3) #6.00e-03 **
Parameter('kf_TNF_TNFR_bind_TRADD', 1.00e-06) #1.00e-06
Parameter('kr_TNF_TNFR_bind_TRADD', 6.02e-02) #1.00e-03
Parameter('kf_TNF_TNFR_bind_RIP1', 1.00e-06) #1.00e-06
Parameter('kr_TNF_TNFR_bind_RIP1', 1e-3) #1.44e-01
Parameter('kf_TNF_TNFR_RIP1_bind_TRADD', 1.00e-06) #1.00e-06
Parameter('kr_TNF_TNFR_RIP1_bind_TRADD', 1e-3) #2.23e-03
Parameter('kf_TNF_TNFR_RIP1_TRADD_bind_A20', 1.00e-06) #1.00e-06
Parameter('kr_TNF_TNFR_RIP1_TRADD_bind_A20', 1.00e-03) #1.00e-03
Parameter('kc_TNF_TNFR_RIP1_TRADD_A20_disassociate', .1) #2.23e-03
Parameter('kf_RIP1_TRADD_bind_FADD', .1) #
Parameter('kr_RIP1_TRADD_bind_FADD', 3.11e-7) #
Parameter('kf_RIP1_TRADD_FADD_bind_RIP3', 3.27e-06) #8.72e-05 #same binding rate as C8
Parameter('kr_RIP1_TRADD_FADD_bind_RIP3', 0.018) #1.08 #same binding rates as C8
Parameter('kc_RIP1_TRADD_FADD_RIP3_disassociate', .1) #6.00
Parameter('kc_RIP3_phosphorylated', 1e-2)
Parameter('kc_RIP1_phosphorylated', 1e-3)
Parameter('kf_RIP1_RIP3_bind_MLKL', 1e-3)
Parameter('kr_RIP1_RIP3_bind_MLKL', 1e-6)
Parameter('kc_RIP1_RIP3_activate_MLKL', 1)
Parameter('kf_RIP1_TRADD_FADD_bind_C8', 3.27e-06) #8.72e-05
Parameter('kr_RIP1_TRADD_FADD_bind_C8', 0.018) #1.08
Parameter('kf_RIP1_TRADD_FADD_C8_bind_FLIP', 3.27e-06) #8.72e-05
Parameter('kr_RIP1_TRADD_FADD_C8_bind_FLIP', 0.018) #1.08
Parameter('kc_FLIP_activates_C8', 3.27e-02) #6
Parameter('kc_RIP1_TRADD_FADD_C8_FLIP_disassociate', .1) #6.00
Parameter('kf_C8a_binds_BID', 1e-6)
Parameter('kr_C8a_binds_BID', 1e-3)
Parameter('kc_C8a_truncates_BID', 1)

#TNF INTERACTIONS
#Unbound TNF is degraded
Rule('degrade_TNF', TNF(tnfr=None) >> None, kc_degrade_TNF)

#Unbound TNF and unbound TNFR reversibly bind to create the TNF/TNFR complex
Rule('TNF_TNFR_bind', TNF(tnfr=None) + TNFR(tnf=None, traddrip1=None) | TNF(tnfr=1) % TNFR(tnf=1, traddrip1=None), kf_TNF_TNFR_bind, kr_TNF_TNFR_bind)


#TRADD BINDING FIRST: DOES NOTHING
#TNF/TNFR complex and unbound TRADD reversibly bind to create the TNF/TNFR/TRADD complex
Rule('TNF_TNFR_bind_TRADD', TNF(tnfr=1) % TNFR(tnf=1, traddrip1=None) + TRADD(tnfr=None, rip1=None) | TNF(tnfr=1) % TNFR(tnf=1, traddrip1=2) % TRADD(tnfr=2, rip1=None),
     kf_TNF_TNFR_bind_TRADD, kr_TNF_TNFR_bind_TRADD)


#RIP1 BINDING FIRST: GOES TO RIPP/TRADD/FADD COMPLEX
#TNF/TNFR complex and unbound RIP1 reversibly bind to create the TNF/TNFR/RIP1 complex
Rule('TNF_TNFR_bind_RIP1', TNF(tnfr=1) % TNFR(tnf=1, traddrip1=None) + RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=None, mod='ub')
     | TNF(tnfr=1) % TNFR(tnf=1, traddrip1=2) % RIP1(tnfr=2, tradd=None, a20=None, fadd=None, rip3c8=None, mod='ub'),
     kf_TNF_TNFR_bind_RIP1, kr_TNF_TNFR_bind_RIP1)

#TNF/TNFR/RIP1 complex and TRADD reversibly bind to create the TNF/TNFR/RIP1/TRADD complex. TRADD is only bound to RIP1.
Rule('TNF_TNFR_RIP1_bind_TRADD',
     TNF(tnfr=1) % TNFR(tnf=1, traddrip1=2) % RIP1(tnfr=2, tradd=None, a20=None, fadd=None, rip3c8=None, mod='ub') + TRADD(tnfr=None, rip1=None)
     | TNF(tnfr=1) % TNFR(tnf=1, traddrip1=2) % RIP1(tnfr=2, tradd=3, a20=None, fadd=None, rip3c8=None, mod='ub') % TRADD(tnfr=None, rip1=3),
     kf_TNF_TNFR_RIP1_bind_TRADD, kr_TNF_TNFR_RIP1_bind_TRADD)

#TNF/TNFR/RIP1/TRADD complex and A20 reversibly bind to create the TNF/TNFR/RIP1/TRADD/A20 complex
Rule('TNF_TNFR_RIP1_TRADD_bind_A20',
     TNF(tnfr=1) % TNFR(tnf=1, traddrip1=2) % RIP1(tnfr=2, tradd=3, a20=None, fadd=None, rip3c8=None, mod='ub') % TRADD(tnfr=None, rip1=3) + A20(rip1=None)
    | TNF(tnfr=1) % TNFR(tnf=1, traddrip1=2) % RIP1(tnfr=2, tradd=3, a20=4, fadd=None, rip3c8=None, mod='ub') % TRADD(tnfr=None, rip1=3) % A20(rip1=4),
     kf_TNF_TNFR_RIP1_TRADD_bind_A20, kr_TNF_TNFR_RIP1_TRADD_bind_A20)

#TNF/TNFR/RIP1/TRADD/A20 complex nonreversibly disassociates leaving only the TRADD/RIP1 complex. RIP1 is deubiquitinated by A20 in the process.
Rule('TNF_TNFR_RIP1_TRADD_A20_disassociate',
     TNF(tnfr=1) % TNFR(tnf=1, traddrip1=2) % RIP1(tnfr=2, tradd=3, a20=4, fadd=None, rip3c8=None, mod='ub') % TRADD(tnfr=None, rip1=3) % A20(rip1=4)
     >> TNF(tnfr=None) + TNFR(tnf=None, traddrip1=None) + TRADD(tnfr=None, rip1=3) % RIP1(tnfr=None, tradd=3, a20=None, fadd=None, rip3c8=None, mod='db') + A20(rip1=None),
     kc_TNF_TNFR_RIP1_TRADD_A20_disassociate)

#RIP1/TRADD and FADD reversibly bind to create the RIP1/TRADD/FADD complex
Rule('RIP1_TRADD_bind_FADD', TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=None, rip3c8=None, mod='db') + FADD(rip1=None, c8=None, flip=None)
     | TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=None, flip=None),
     kf_RIP1_TRADD_bind_FADD, kr_RIP1_TRADD_bind_FADD)


#RIP/TRADD/FADD COMPLEX: GOES TO NECROPTOSIS
#RIP1/TRADD/FADD complex and RIP3 reversibly bind to create the RIP1/TRADD/FADD/RIP3 complex
Rule('RIP1_TRADD_FADD_bind_RIP3',
     TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=None, flip=None) + RIP3(rip1=None, mlkl=None, mod='unmod')
     | TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=6, mod='db') % FADD(rip1=5, c8=None, flip=None) % RIP3(rip1=6, mlkl=None, mod='unmod'),
     kf_RIP1_TRADD_FADD_bind_RIP3, kr_RIP1_TRADD_FADD_bind_RIP3)

#RIP1/TRADD/FADD/RIP3 nonreversibly disassociate leaving only the RIP1/RIP3 complex
Rule('RIP1_TRADD_FADD_RIP3_disassociate',
     TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=6, mod='db') % FADD(rip1=5, c8=None, flip=None) % RIP3(rip1=6, mlkl=None, mod='unmod')
     >> TRADD(tnfr=None, rip1=None) + RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='db') % RIP3(rip1=6, mlkl=None, mod='unmod') + FADD(rip1=None, c8=None, flip=None),
     kc_RIP1_TRADD_FADD_RIP3_disassociate)

#RIP1/RIP3 nonreversibly phosphorylates RIP3
Rule('RIP3_phosphorylated', RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='db') % RIP3(rip1=6, mlkl=None, mod='unmod')
     >> RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='db') % RIP3(rip1=6, mlkl=None, mod='p'),
     kc_RIP3_phosphorylated)

#RIP1/RIP3p nonreversibly phosphorylates RIP1
Rule('RIP1_phosphorylated', RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='db') % RIP3(rip1=6, mlkl=None, mod='p')
     >> RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='p') % RIP3(rip1=6, mlkl=None, mod='p'),
     kc_RIP1_phosphorylated)

#RIP1p/RIP3p reversibly binds MLKL
Rule('RIP1_RIP3_bind_MLKL', RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='p') % RIP3(rip1=6, mlkl=None, mod='p') + MLKL(rip3=None, mod='i')
     | RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='p') % RIP3(rip1=6, mlkl=8, mod='p') % MLKL(rip3=8, mod='i'),
     kf_RIP1_RIP3_bind_MLKL, kr_RIP1_RIP3_bind_MLKL)

#RIP1p/RIP3p nonreversibly unbinds MLKL and activates/phosphorylates MLKL
Rule('RIP1_RIP3_activate_MLKL', RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='p') % RIP3(rip1=6, mlkl=8, mod='p') % MLKL(rip3=8, mod='i')
     >> RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='p') % RIP3(rip1=6, mlkl=None, mod='p') + MLKL(rip3=None, mod='p'),
     kc_RIP1_RIP3_activate_MLKL)


#RIP/TRADD/FADD COMPLEX: GOES TO APOPTOSIS
#RIP1/TRADD/FADD complex and C8 reversibly bind to create the RIP1/TRADD/FADD/C8 complex
Rule('RIP1_TRADD_FADD_bind_C8',
     TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=None, flip=None) + C8(fadd=None, flip=None, bid=None, mod='i')
     | TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=6, flip=None) % C8(fadd=6, flip=None, bid=None, mod='i'),
     kf_RIP1_TRADD_FADD_bind_C8, kr_RIP1_TRADD_FADD_bind_C8)

#RIP1/TRADD/FADD/C8 complex and FLIP reversibly bind to create the RIP1/TRADD/FADD/C8/FLIP complex
Rule('RIP1_TRADD_FADD_C8_bind_FLIP',
     TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=6, flip=None) % C8(fadd=6, flip=None, bid=None, mod='i') + FLIP(c8=None, fadd=None)
     | TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=6, flip=9) % C8(fadd=6, flip=7, bid=None, mod='i') % FLIP(c8=7, fadd=9),
     kf_RIP1_TRADD_FADD_bind_C8, kr_RIP1_TRADD_FADD_bind_C8)

#RIP1/TRADD/FADD/C8/FLIP complex, FLIP nonreversibly activates C8
Rule('FLIP_activates_C8',
     TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=6, flip=9) % C8(fadd=6, flip=7, bid=None, mod='i') % FLIP(c8=7, fadd=9)
     >> TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=6, flip=9) % C8(fadd=6, flip=7, bid=None, mod='a') % FLIP(c8=7, fadd=9),
     kc_FLIP_activates_C8)

#RIP1/TRADD/FADD/C8_activated/FLIP complex nonreversibly disassociates
Rule('RIP1_TRADD_FADD_C8_FLIP_disassociate',
     TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=6, flip=9) % C8(fadd=6, flip=7, bid=None, mod='a') % FLIP(c8=7, fadd=9)
     >> TRADD(tnfr=None, rip1=None) + RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=None, mod='db') + FADD(rip1=None, c8=None, flip=None) + C8(fadd=None, flip=None, bid=None, mod='a') + FLIP(c8=None, fadd=None),
     kc_RIP1_TRADD_FADD_C8_FLIP_disassociate)

#Activated C8 reversibly binds BID
Rule('C8a_binds_BID', C8(fadd=None, flip=None, bid=None, mod='a') + BID(c8=None, mod='unmod') | C8(fadd=None, flip=None, bid=10, mod='a') % BID(c8=10, mod='unmod'),
     kf_C8a_binds_BID, kr_C8a_binds_BID)

#C8/Bid nonreversibly disassociates, Bid is truncated into tBID
Rule('C8a_truncates_BID', C8(fadd=None, flip=None, bid=10, mod='a') % BID(c8=10, mod='unmod') >> C8(fadd=None, flip=None, bid=None, mod='a') + BID(c8=None, mod='trunc'),
     kc_C8a_truncates_BID)


#Initial Conditions
Parameter('TNF_0', 2326)#2326     #698 is 30ng/ml of TNF
Parameter('TNFR_0', 4800)#4800    #0.00246
Parameter('TRADD_0', 9000)#9000
Parameter('RIP1_0', 40000)#random.randint(400, 400000)  #47000 0.04
Parameter('A20_0', 9000)#random.randint(90, 90000)     #2256
Parameter('FADD_0', 8030)#8030    #0.0033
Parameter('RIP3_0', 40000)#40000  #20000
Parameter('MLKL_0', 10000)
Parameter('C8_0', 9000)#random.randint(90, 90000)      #0.004 #0.09
Parameter('FLIP_0', 3900)#3900    #0.004 #0.09
Parameter('BID_0', 10000)
Initial(TNF(tnfr=None), TNF_0)
Initial(TNFR(tnf=None, traddrip1=None), TNFR_0)
Initial(TRADD(tnfr=None, rip1=None), TRADD_0)
Initial(RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=None, mod='ub'), RIP1_0)
Initial(A20(rip1=None), A20_0)
Initial(FADD(rip1=None, c8=None, flip=None), FADD_0)
Initial(RIP3(rip1=None, mlkl=None, mod='unmod'), RIP3_0)
Initial(MLKL(rip3=None, mod='i'), MLKL_0)
Initial(C8(fadd=None, flip=None, bid=None, mod='i'), C8_0)
Initial(FLIP(c8=None, fadd=None), FLIP_0)
Initial(BID(c8=None, mod='unmod'), BID_0)

#Observables
Observable('obsComplexI', TNF(tnfr=1) % TNFR(tnf=1, traddrip1=2) % RIP1(tnfr=2, tradd=3, a20=4, fadd=None, rip3c8=None, mod='ub') % TRADD(tnfr=None, rip1=3) % A20(rip1=4))
Observable('obsComplexIIa', TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=6, flip=None) % C8(fadd=6, flip=None, bid=None, mod='i'))
Observable('obsComplexIIb', TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=6, mod='db') % FADD(rip1=5, c8=None, flip=None) % RIP3(rip1=6, mlkl=None, mod='unmod'))
Observable('obsMLKLp', MLKL(rip3=None, mod='p'))
Observable('obstBID', BID(c8=None, mod='trunc'))


tspan = np.linspace(0, 2160, 2161) #Length of simulation
plt.ioff() #turn off graph showing up
path = '/home/asasla/main/ComplexII/' #Set Path to save figures
TNF_LOOP = [('100 ng/ml TNF', 2326), ('10 ng/ml TNF', 232), ('1 ng/ml TNF', 23), ('.1 ng/ml TNF', 2)]
all_times = [420, 600, 1440, 2160] #Array of time points of interest to create a probability density function for, in minutes

# path = '/Users/ariella/PycharmProjects/ComplexII/'
# obs_plot_num = [('obsComplexI', 511), ('obsComplexIIa', 512), ('obsComplexIIb', 513), ('obsMLKLp', 514), ('obstBID', 515)] #Observable name, number for figure subplot

#RUN THROUGH EACH AMOUNT OF TNF: HOW DOES THAT AFFECT SSA VS ODE
for tnf_title, dose in TNF_LOOP:

    #RUN STOCHASTIC SIMULATION ALGORITHM (SSA)
    ssa_sim = StochKitSimulator(model, tspan=tspan, verbose=True)
    ssa_sim_res = ssa_sim.run(initials={TNF(tnfr=None): dose}, n_runs=NUM_SSA_RUNS)
    df = ssa_sim_res.dataframe

    #FOR EACH OBSERVABLE AVERAGE THE SSA RUNS AT EACH TIME POINT
    avg = df.groupby(level='time').mean()

    #RUN ODE SIMULATION
    ode_sim = ScipyOdeSimulator(model, tspan=tspan)
    ode_sim_res = ode_sim.run(initials={TNF(tnfr=None): dose})

    #PLOT STOCHASTIC SIMULATION ALGORITHM (SSA) WITH ODE
    for obs in model.observables:
        plt.figure()
        for _, run in df.groupby('simulation'):
                plt.plot(tspan / 60, run.loc[:, obs.name])
        plt.plot(tspan / 60, ode_sim_res.observables[obs.name], 'black', linewidth=3)
        plt.xlabel("Time (in hr)", fontsize=15)
        plt.ylabel("Molecules/Cell", fontsize=15)
        plt.title('%s Trajectories' % obs.name, fontsize=18)
        ssa_name = path + '%d_SSA_%s.png' % (dose, obs.name)
        plt.savefig(ssa_name, bbox_inches='tight')

    # FOR EACH OBSERVABLE: AVERAGE THE SSA RUNS AT EACH TIME POINT AND PLOT VS ODE
    avg = df.groupby(level='time').mean()

    for obs in model.observables:
        plt.figure()
        plt.plot(tspan / 60, avg.loc[:, obs.name], 'blue')
        plt.plot(tspan / 60, ode_sim_res.observables[obs.name], 'black')
        plt.xlabel("Time (in hr)", fontsize=15)
        plt.ylabel("Molecules/Cell", fontsize=15)
        plt.title('%s Trajectories' % obs.name)
        avg_ssa = path + '%d_Avg_%s' % (dose, obs.name)
        plt.savefig(avg_ssa, bbox_inches='tight')

    #AT HIGH VARIABILITY AND END TIMEPOINTS FOR EACH OBSERVABLE: PLOT ALL SSA RUNS OF THAT TIME POINT WITH A DENSITY PLOT
    # Create dataframe
    idx = pd.MultiIndex.from_product([all_times, (range(0, NUM_SSA_RUNS))], names=['timepoint', 'simulation'])
    col = ['obsComplexI', 'obsComplexIIa', 'obsComplexIIb', 'obsMLKLp', 'obstBID']
    df_dens_plot = pd.DataFrame(np.nan, idx, col)

    # Fill in dataframe to plot
    for sim_num, run in df.groupby('simulation'):
        for t_point in all_times:
            for obs in model.observables:
                run_slice = run.loc[[(sim_num, t_point)], [obs.name]]
                conc = run_slice.values[0]
                df_dens_plot.loc[[(t_point, sim_num)], [obs.name]] = conc

    # Plot dataframe
    for t_point in all_times:
        for obs in model.observables:
            plt.figure()
            sns.distplot(df_dens_plot.loc[[t_point], [obs.name]], kde=True)
            plt.xlabel("Molecules/Cell", fontsize=15)
            plt.ylabel("Density", fontsize=15)
            plt.title('%s Kernel Density Estimation at %d Hours' % (obs.name, t_point / 60), fontsize=18)
            pdf_name = path + '%d_KDE_%dhrs_%s.png' % (dose, t_point / 60, obs.name)
            plt.savefig(pdf_name, bbox_inches='tight')

# #DETERMINE ODE VALUE OF VARIABLE TIME POINT USED ABOVE
# for t_point in all_times:
#      for obs in model.observables:
#           print("At %d hour(s) %s amount is %d molecules/cell" % (t_point/60, obs.name, ode_sim_res.observables[t_point][obs.name]))

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
# plt.figure()
# fig, axs = plt.subplots(nrows=5, ncols=1, sharex=True)
# for plot_num, obs in enumerate(model.observables):
#     # for obs, plot_num in obs_plot_num:
#     axs[plot_num].plot(tspan / 60, df.loc[:, obs.name])
#     plt.xlabel("Time (in hr)", fontsize=15)
#     plt.ylabel("%s [Molecules/Cell]" % obs.name, fontsize=15)
#     plt.title(obs.name)
# plt.suptitle('SSA Trajectories: %s' % tnf_title, size=20)
# plt.tight_layout()
# plt.show()
   # #PLOT STOCHASTIC SIMULATION ALGORITHM (SSA)
   #  plt.figure()
   #  fig, axs = plt.subplots(nrows=5, ncols=1, sharex=True)
   #  for plot_num, obs in enumerate(model.observables):
   #  # for obs, plot_num in obs_plot_num:
   #      axs[plot_num].plot( for _, run in df.groupby('simulation'): (run.loc[:, obs.name])
   #      plt.xlabel("Time (in hr)", fontsize=15)
   #      plt.ylabel("%s [Molecules/Cell]" % obs.name, fontsize=15)
   #      plt.title(obs.name)
   #  plt.suptitle('SSA Trajectories: %s' % tnf_title, size=20)
   #  plt.tight_layout()
   #  plt.show()
   #  # ssa_name = path + '%d_SSA_%s.png' % (dose, obs)
   #  # plt.savefig(ssa_name, bbox_inches='tight')