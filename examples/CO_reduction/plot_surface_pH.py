import math
import numpy as np
import matplotlib.pyplot as plt
from catint.plot import Plot
from catint.transport import Transport
from glob import glob
from catint.io import read_all
from itertools import cycle
import sys
from units import *
import re

colors=cycle(['orange','b']) #,'k','r','y'])
#colors=cycle(['orange','b','k','r','y'])
styles=cycle(['-','--']) #,'-.',':','o','d','x'])
folders=sys.argv[1:] #glob('test') #calc_std_settings*') #obj') #calc_std*')
tp=Transport()
oldlabel=''
ax1=plt.subplot('211')
ax2=plt.subplot('212')

all_coverages=[]

for f in folders:
    read_all(tp,f,only=['all_data','species','descriptors','comsol_outputs','comsol_outputs_data'])
    for isp,sp in enumerate(tp.species):
        if sp=='OH-':
            index=isp

    pHs=[]
    desc=[]
    for key in tp.all_data:
        for key2 in tp.all_data[key]:
            #print tp.all_data[key][key2]['system']['cout']
#            print tp.all_data[key][key2]['species']['OH-']['concentration']
            #plt.plot(tp.all_data[key][key2]['species']['OH-']['concentration'])

            #pHs.append(14+np.log(tp.all_data[key][key2]['species']['OH-']['concentration'][0]/1000.))
            pHs.append([float(key),14+np.log10(tp.all_data[key][key2]['system']['cout'][-1][index]/1000.)])
            desc.append(float(key))
    pHs=np.array(pHs)
    pHs.view(pHs.dtype.str+','+pHs.dtype.str).sort(order=['f0'], axis=0)
    #print pHs
    ax1.plot(pHs[:,0],pHs[:,1],'-',label=f.replace('calc_std_settings','std'))

    #2nd the coverages:
    cov=tp.comsol_outputs_data['coverage']
    coverages=[]
    for key1,key2 in cov:
        if not math.isnan(cov[(key1,key2)][0]):
            coverages.append([float(key1),cov[(key1,key2)][0]])
        else:
            coverages.append([float(key1),0.0])
    coverages=np.array(coverages)
    coverages.view(coverages.dtype.str+','+coverages.dtype.str).sort(order=['f0'], axis=0)

    all_coverages.append(coverages)
    print 'coverages',coverages
    ax2.plot(coverages[:,0],coverages[:,1],'-',label=f.replace('calc_std_settings','std'))
#ax2.semilogy(coverages[:,0],all_coverages[1][:,1]/all_coverages[0][:,1],'-',label=f)


ax1.set_xlabel('Voltage vs. RHE (V)')
ax2.set_xlabel('Voltage vs. RHE (V)')
#ax2.legend()
#ax1.legend()
ax1.set_ylabel('pH')
ax2.set_ylabel(r'$\theta_\mathrm{CO}$')
plt.show()
