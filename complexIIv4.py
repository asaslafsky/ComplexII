#hello

# pysb.pathfinder.set_path(ariella/anaconda3/envs/pysb/lib/stochkit)

# import the pysb module and all its methods and functions
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
from pysb.simulator import StochKitSimulator
from pysb.simulator import ScipyOdeSimulator
from pysb import *
#import random

# Definitions
RUN_TYPE = 'Test 10,000 runs' #Gives Titles to Saved Graphs: name according to what is changed i.e. protein fluctuations
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
Monomer('MLKL', ['rip3', 'mod'], {'mod': ['i', 'a']})
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
Parameter('kf_RIP1_TRADD_bind_FADD', 1.25e-5) #
Parameter('kr_RIP1_TRADD_bind_FADD', 10) #
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

#RIP1p/RIP3p nonreversibly unbinds MLKL and activates MLKL
Rule('RIP1_RIP3_activate_MLKL', RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='p') % RIP3(rip1=6, mlkl=8, mod='p') % MLKL(rip3=8, mod='i')
     >> RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='p') % RIP3(rip1=6, mlkl=None, mod='p') + MLKL(rip3=None, mod='a'),
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
Observable('obsMLKLa', MLKL(rip3=None, mod='a'))
Observable('obstBID', BID(c8=None, mod='trunc'))

#LENGTH OF SIMULATION
tspan = np.linspace(0, 2160, 2161) #Length of sim (start, stop, number of samples to generate). 4320 min, 72 hours.
plt.ioff()

#RUN STOCHASTIC SIMULATION ALGORITHM (SSA)
ssa_sim = StochKitSimulator(model, tspan=tspan, verbose=True)
ssa_sim_res = ssa_sim.run(n_runs=NUM_SSA_RUNS)
df = ssa_sim_res.dataframe
print(df.head())
print(df.index.names)

for obs in model.observables:
    plt.figure()
    for _, run in df.groupby('simulation'):  #groupby returns tuple: index and run
            plt.plot(tspan / 60, run.loc[:, obs.name])
    plt.xlabel("Time (in hr)", fontsize=15)
    plt.ylabel("%s [Molecules/Cell]" % obs.name, fontsize=15)
    plt.title('%s Trajectories' % obs.name)
    plt.savefig('%s: SSA %s' % (RUN_TYPE, obs.name), bbox_inches='tight')

#RUN ODE SIMULATION
ode_sim = ScipyOdeSimulator(model, tspan=tspan)
ode_sim_res = ode_sim.run()

for obs in model.observables:
    plt.figure()
    plt.plot(tspan / 60, ode_sim_res.observables[obs.name], label='100 ng/ml TNF')
    plt.xlabel("Time (in hr)", fontsize=15)
    plt.ylabel("%s [Molecules/Cell]" % obs.name, fontsize=15)
    plt.title('%s Trajectories' % obs.name)
    plt.legend(loc='best')
    plt.savefig('%s: ODE %s' % (RUN_TYPE, obs.name), bbox_inches='tight')

# FOR EACH OBSERVABLE: AVERAGE THE SSA RUNS AT EACH TIME POINT AND PLOT
avg = df.groupby(level='time').mean()

for obs in model.observables:
    plt.figure()
    plt.plot(tspan / 60, avg.loc[:, obs.name])
    plt.xlabel("Time (in hr)", fontsize=15)
    plt.ylabel("%s [Molecules/Cell]" % obs.name, fontsize=15)
    plt.title('Average SSA %s Trajectories' % obs.name)
    plt.legend(loc='best')
    plt.savefig('%s: Average SSA %s' % (RUN_TYPE, obs.name), bbox_inches='tight')

# #AT HIGH VARIABILITY AND END TIMEPOINTS FOR EACH OBSERVABLE: PLOT ALL SSA RUNS OF THAT TIME POINT WITH A DENSITY PLOT
# all_times = [720, 1440, 2880] #Array of time points of interest: 24, 48, 72
#
# #Create dataframe
# idx = MultiIndex.from_product([all_times, range(1, NUM_SSA_RUNS + 1)], names=['timepoint', 'simulation'])
# col = ['obsComplexI', 'obsComplexIIa', 'obsComplexIIb']
# df_dens_plot = pd.DataFrame('', idx, col)
#
# #Fill in dataframe to plot
# for sim_num, run in df.groupby('simulation'):
#     run['Index'] = range(0, len(run))
#     for t_point in all_times:
#         obs_row = run.loc[run['Index'] == t_point]
#         for obs in model.observables:
#             df_dens_plot.add(obs_row.loc[:, obs.name], index=[t_point, sim_num], columns=obs.name)
#
# print(df_dens_plot)

#Plot dataframe


# for t_point in all_times:
#      for obs in model.observables:
#           plt.figure()
#           for _, run in df.groupby('simulation'):
#                run['Index'] = range(0, len(run))
#                time = run.loc[run['Index'] == t_point]
#                plt.hist(time.loc[:, obs.name], density=True, histtype='stepfilled')
#           plt.xlabel("%s [Molecules/Cell]" % obs.name, fontsize=15)
#           plt.ylabel("Frequency", fontsize=15)
#           plt.title('Density Plot for %s at %d Hour(s)' % (obs.name, t_point/60))
#           plt.legend(loc='best')
#           plt.savefig('%s: SSA variable timepoint %d for %s' % (RUN_TYPE, t_point, obs.name), bbox_inches='tight')



# #AT HIGH VARIABILITY AND END TIMEPOINTS FOR EACH OBSERVABLE: PLOT ALL SSA RUNS OF THAT TIME POINT WITH A DENSITY PLOT
# all_times = [1440, 2880, 4320] #Array of time points of interest: 24, 48, 72
#
# for t_point in all_times:
#      for obs in model.observables:
#           plt.figure()
#           for _, run in df.groupby('simulation'):
#                run['Index'] = range(0, len(run))
#                time = run.loc[run['Index'] == t_point]
#                plt.hist(time.loc[:, obs.name], density=True, histtype='stepfilled')
#           plt.xlabel("%s [Molecules/Cell]" % obs.name, fontsize=15)
#           plt.ylabel("Frequency", fontsize=15)
#           plt.title('Density Plot for %s at %d Hour(s)' % (obs.name, t_point/60))
#           plt.legend(loc='best')
#           plt.savefig('%s: SSA variable timepoint %d for %s' % (RUN_TYPE, t_point, obs.name), bbox_inches='tight')

# #DETERMINE ODE VALUE OF VARIABLE TIME POINT USED ABOVE
# for t_point in all_times:
#      for obs in model.observables:
#           print("At %d hour(s) %s amount is %d molecules/cell" % (t_point/60, obs.name, ode_sim_res.observables[obs.name][t_point]))

#plt.show()

