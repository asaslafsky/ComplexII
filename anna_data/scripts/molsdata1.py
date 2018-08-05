# for reading and plotting mols data in smoldyn simulation


import matplotlib.pyplot as plt
import sys
import numpy as np
from pysb.simulator import ScipyOdeSimulator
from pysb import *
import pandas as pd
import seaborn as sns
from scipy import stats
from pysb.simulator import StochKitSimulator

plt.switch_backend('agg')

from complexIIv4_normalization import model

#from complexIIv4_normalization import model
NUM_SSA_RUNS = 10

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
# ssa_sim_res = ssa_sim.run(initials={TNF(tnfr=None): dose}, n_runs=NUM_SSA_RUNS)
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
            grpnames = ['Membrane','No membrane']
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
                    avg = tot / count
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
                for i in range(len(avg_list)):
                    avg_list[i] = normalize(avg_list[i],key)

                erlow=[]
                for i in range(len(ster)):
                    erlow.append(avg_list[i]-ster[i])
                erhigh=[]
                for i in range(len(ster)):
                    erhigh.append(avg_list[i]+ster[i])
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
                ax.plot(tspan / 60, ssa_avg.loc[:, obs.name]/ariella_norms[key], color='mediumpurple', label='SSA Avg',linewidth=1.5)
                ax.plot(tspan/60, ode_sim_res.observables[obs.name]/ariella_norms[key], color = 'black',linestyle=':',label='ODE',linewidth=3)
                break

        ax.set_xlabel('Time (hrs)')
        ax.set_ylabel('Proportion of total possible')
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
