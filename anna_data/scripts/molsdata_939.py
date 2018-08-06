# for reading and plotting mols data in smoldyn simulation


import matplotlib.pyplot as plt
import sys
import numpy as np
from pysb.simulator import ScipyOdeSimulator
from pysb import *
from pysb.simulator import StochKitSimulator

plt.switch_backend('agg')

#from complexIIv4_normalization import model
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
Parameter('kr_TNF_TNFR_bind_TRADD', 0.001) #6.02e-02) #1.00e-03
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
Parameter('TNF_0', 9390)#2326     #698 is 30ng/ml of TNF
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

##Observables
Observable('O_1', TNF(tnfr=None))
Observable('O_2', TNFR(tnf=None,traddrip1=None))
Observable('O_3', TRADD(tnfr=None, rip1=None))

Observable('O_4',RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=None, mod='ub'))
Observable('O_5',A20(rip1=None))
Observable('O_6',FADD(rip1=None, c8=None, flip=None))
Observable('O_7',FLIP(c8=None, fadd=None))
Observable('O_8',C8(fadd=None, flip=None, bid=None, mod='i'))
Observable('O_9',RIP3(rip1=None, mlkl=None, mod='unmod'))
Observable('O_10',TNF(tnfr=1) % TNFR(tnf=1, traddrip1=None))
Observable('O_11',TNF(tnfr=1) % TNFR(tnf=1, traddrip1=2) % TRADD(tnfr=2, rip1=None))
Observable('O_12',TNF(tnfr=1) % TNFR(tnf=1, traddrip1=2) % RIP1(tnfr=2, tradd=None, a20=None, fadd=None, rip3c8=None, mod='ub'))
Observable('O_13',TNF(tnfr=1) % TNFR(tnf=1, traddrip1=2) % RIP1(tnfr=2, tradd=3, a20=None, fadd=None, rip3c8=None, mod='ub') % TRADD(tnfr=None, rip1=3))
Observable('O_14',TNF(tnfr=1) % TNFR(tnf=1, traddrip1=2) % RIP1(tnfr=2, tradd=3, a20=4, fadd=None, rip3c8=None, mod='ub') % TRADD(tnfr=None, rip1=3) % A20(rip1=4))
Observable('O_15',TRADD(tnfr=None, rip1=3) % RIP1(tnfr=None, tradd=3, a20=None, fadd=None, rip3c8=None, mod='db'))
Observable('O_16',TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=None, flip=None))
Observable('O_17',TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=6, flip=None) % C8(fadd=6, flip=None, bid=None, mod='i'))

Observable('O_18',TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=6, mod='db') % FADD(rip1=5, c8=None, flip=None) % RIP3(rip1=6, mlkl=None, mod='unmod'))
Observable('O_19',TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=6, flip=9) % C8(fadd=6, flip=7, bid=None, mod='i') % FLIP(c8=7, fadd=9))

Observable('O_20',TRADD(tnfr=None, rip1=4) % RIP1(tnfr=None, tradd=4, a20=None, fadd=5, rip3c8=None, mod='db') % FADD(rip1=5, c8=6, flip=9) % C8(fadd=6, flip=7, bid=None, mod='a') % FLIP(c8=7, fadd=9))

Observable('O_21',C8(fadd=None, flip=None, bid=None, mod='a'))

Observable('O_22',RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='db') % RIP3(rip1=6, mlkl=None, mod='unmod'))
Observable('O_23',RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='db') % RIP3(rip1=6, mlkl=None, mod='p'))
Observable('O_24',RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='p') % RIP3(rip1=6, mlkl=None, mod='p'))
Observable('O_25',RIP1(tnfr=None, tradd=None, a20=None, fadd=None, rip3c8=6, mod='p') % RIP3(rip1=6, mlkl=8, mod='p') % MLKL(rip3=8, mod='i'))
Observable('O_26',MLKL(rip3=None, mod='i'))

Observable('O_27',MLKL(rip3=None, mod='p'))
#Observable('O_28',BID(c8=None, mod='unmod'))
Observable('O_29',BID(c8=None, mod='trunc'))
#Observable('O_30',C8(fadd=None, flip=None, bid=10, mod='a') % BID(c8=10, mod='unmod'))



#from complexIIv4_normalization import model
NUM_SSA_RUNS = 10000

datafiles = []
if sys.argv[1] == 'file':
    datafiles = sys.argv[2:]
else:
    datalist_f = sys.argv[1]
    dose = int(sys.argv[2])
    with open(datalist_f,'r') as f:
        for line in f.readlines():
            datafiles.append(line.replace('\n',''))
        f.close()

# list of species to make average plots for with ODE
avgs = ['C8.1.FADD.1.RIP1.1.TRADD.1.0','FADD.1.RIP1.1.RIP3.1.TRADD.1.0','A20.1.RIP1.1.TNF.1.TNFR.1.TRADD.1.0','MLKL.1.1','Bid.1.1','C8.1.1']

# plot titles
plotnames = {'C8.1.FADD.1.RIP1.1.TRADD.1.0':'Complex IIa Formation',
            'FADD.1.RIP1.1.RIP3.1.TRADD.1.0':'Complex IIb Formation',
            'A20.1.RIP1.1.TNF.1.TNFR.1.TRADD.1.0':'Complex I Formation',
            'MLKL.1.1':'MLKL Phosphorylation',
            'Bid.1.1':'Bid Activation',
            'C8.1.1':'C8 Activation'}

# totals to normalize against (calculated in norm.xlsx)
norms = {'C8.1.FADD.1.RIP1.1.TRADD.1.0':335,
        'FADD.1.RIP1.1.RIP3.1.TRADD.1.0':335,
        'A20.1.RIP1.1.TNF.1.TNFR.1.TRADD.1.0':200,
        'MLKL.1.1':417,
        'Bid.1.1':417,
        'C8.1.1':375}

ariella_norms = {'C8.1.FADD.1.RIP1.1.TRADD.1.0':8030,
                'FADD.1.RIP1.1.RIP3.1.TRADD.1.0':8030,
                'A20.1.RIP1.1.TNF.1.TNFR.1.TRADD.1.0':4800,
                'MLKL.1.1':10000,
                'Bid.1.1':10000,
                'C8.1.1':9000}

# mod is the step (which data points to use)
mod = 1
datalist = []
keys = []

plt.rc('legend', fontsize=12)

def readlines(datafile):
    cont_lns_list= []
    f = open(datafile, 'r')
    if f.mode == 'r':
        contents = f.read()
    f.close()
    cont_lns = contents.split('\n')
    cont_lns.pop() #remove last entry which is empty line
    for ln in cont_lns:
        line = ln.split()
        cont_lns_list.append(line)
    return cont_lns_list

def popdict(cont_lns_list,pop):
    data = {}
    for entry in cont_lns_list[0]:
        data[entry] = []
        if not pop:
            keys.append(entry)
    data['Average'] = []
    if not pop:
        keys.append('Average')
    datalist.append(data)
    return True

def enterdata(cont_lns_list,index):
    data = datalist[index]
    for dataset in cont_lns_list[1:]:
        i = 0
        tot = 0
        for entry in dataset:
            if i % mod == 0:
                if i == 0:
                    data[keys[0]].append(float(entry) * 0.000277778)           # convert seconds to hours
                else:
                    data[keys[i]].append(int(entry))
                    tot += int(entry)
            i += 1
        data['Average'].append(tot / (len(keys)-1))

def normalize(num,key):
    return num/norms[key]

i = 0
pop = False
for f in datafiles:
    print(f)
    cont = readlines(f)
    pop = popdict(cont,pop)
    enterdata(cont,i)
    i += 1

 #Length of sim
tspan = np.linspace(0, 1440, 1441)

#Run ODE
print('simulating')
ode_sim = ScipyOdeSimulator(model, tspan=tspan)
ode_sim_res = ode_sim.run()
yout = ode_sim_res.all

#SSA Simulation
ssa_sim = StochKitSimulator(model, tspan=tspan, verbose=True)
ssa_sim_res = ssa_sim.run(n_runs=NUM_SSA_RUNS)
#ssa_sim_res = ssa_sim.run(initials={TNF(tnfr=None): dose}, n_runs=NUM_SSA_RUNS)
df = ssa_sim_res.dataframe

#FOR EACH OBSERVABLE AVERAGE THE SSA RUNS AT EACH TIME POINT
ssa_avg = df.groupby(level='time').mean()

for key in keys:
    print(key)
    if key in avgs:
    #if key != 'time':
        color = iter(plt.cm.rainbow(np.linspace(0,1,len(datafiles)+1)))
        f = plt.figure()
        ax = f.add_subplot(1,1,1)
        '''
        for i in range(len(datafiles)):
            data = datalist[i]
            c = next(color)
            ax.plot(data['time'],data[key],color=c,linestyle='-',marker='')
        '''
        if key in avgs:
            #grpnames = ['~100 ng/mL TNF','~10 ng/mL TNF','~0.1 ng/mL TNF']
            grpnames = ['Smol, Membrane','Smol, No membrane']
            cols = ['tomato','cornflowerblue']
            #cols = ['navy', 'indigo', 'darkorange']
            grps = [60,10]
            index = 0
            for n in range(len(grps)):
                avg_list = []
                ster = []
                for i in range(len(datalist[0]['time'])):
                    tot = 0
                    count = 0
                    for j in range(grps[n]):
                        tot += datalist[j+index][key][i]
                        count += 1
                    avg = float(tot) / count
                    avg_list.append(avg)
                loops = 0
                for i in range(len(avg_list)):
                    tot = 0
                    if True:
                        for j in range(grps[n]):
                            tot += (avg_list[i] - datalist[j+index][key][i])**2
                        #ster.append(((tot/grps[n])**0.5)/grps[n]**0.5)
                        ster.append(normalize(((tot / grps[n])**0.5) / grps[n]**0.5,key))          #for standard error
                        # NORMALIZE STANDARD ERROR HERE (above)
                        #ster.append((tot/grps[n])**0.5)                            #for standard deviation
                    loops += 1

                # NORMALIZE AVG_LIST HERE
                print(avg_list)
                for i in range(len(avg_list)):
                    avg_list[i] = normalize(avg_list[i],key)
                print(avg_list)

                # erlow=[]
                # for i in range(len(ster)):
                #     erlow.append(avg_list[i]-ster[i])
                # erhigh=[]
                # for i in range(len(ster)):
                #     erhigh.append(avg_list[i]+ster[i])
                        #ax.fill_between(datalist[0]['time'][i], avg_list[i]-ster, avg_list[i]+ster,alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF') #, antialiased=True)
                        #ax.errorbar(datalist[0]['time'][i], avg_list[i], yerr=((tot / grps[n])**0.5) / grps[n]**0.5,color='red',capsize=5)

                c = next(color)
                ax.plot(datalist[0]['time'],avg_list,color=cols[n],linestyle='-',marker='',label=grpnames[n],linewidth=1.5)
                #ax.plot(datalist[0]['time'],avg_list,color=c,linestyle='-',marker='',label=grpnames[n],linewidth=1)
                index += grps[n]
                #ax.fill_between(datalist[0]['time'],erlow,erhigh,facecolor=cols[n],alpha=0.3)

        for obs in model.observables:
            # Crappy quick fix, but observables must be in same order in model as in data files
            if key == keys[int(obs.name.replace('O_',''))]:
                ax.plot(tspan/60, ssa_avg.loc[:, obs.name]/ariella_norms[key], color='gold', label='SSA Avg',linewidth=1.5)
                ax.plot(tspan/60, ode_sim_res.observables[obs.name]/ariella_norms[key], color = 'black',linestyle=':',label='ODE',linewidth=3)
                break

        ax.set_xlabel('Time (hrs)')
        ax.set_ylabel('Proportion of Total Possible Species')
        ax.xaxis.label.set_fontsize(16)
        ax.yaxis.label.set_fontsize(16)
        if key in plotnames.keys():
            ax.set_title(plotnames[key])
            ax.title.set_fontsize(18)
        else:
            ax.set_title(key)
            ax.title.set_fontsize(18)
        ax.legend(loc='best', frameon=False)
        plotname = key + '.png'
        f.savefig(plotname, bbox_inches='tight')
        plt.close(f)

'''
def compoundplot(species):
    print(species)
    cf = plt.figure()
    ax = cf.add_subplot(1,1,1)
    for entry in species:
        ax.plot(data['time'],data[entry], label=entry)
    ax.set_xlabel('Time (hrs)')
    ax.set_ylabel('Number of molecules')
    ax.legend(loc='upper left', frameon=False)
    return cf

# add which compound plot you want here
#plot = compoundplot(['C8.1.0','C8.1.1'])
#plt.show(plot)
#plot.savefig('cplot', bbox_inches='tight')
'''
