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

voltage=-0.5 # vs. RHE
pH = 13
voltage-=0.0592*pH


colors=cycle(['orange','b']) #,'k','r','y'])
#colors=cycle(['orange','b','k','r','y'])
styles=cycle(['-','--']) #,'-.',':','o','d','x'])
folders=sys.argv[1:] #glob('test') #calc_std_settings*') #obj') #calc_std*')
tp=Transport()
oldlabel=''
fig=plt.figure(figsize=(9,7))
ax1=fig.add_subplot('221')
ax1.set_title('Surface Concentrations')
ax1.set_xlabel('Voltage vs. RHE (V)')
ax2=fig.add_subplot('223')
ax2.set_xlabel('Voltage vs. RHE (V)')
#ax2.set_title('Surface Concentrations (logy)')
ax3=fig.add_subplot('222')
ax3.set_xlabel('x (m)')
ax4=fig.add_subplot('224')
ax4.set_xlabel('x (m)')
all_coverages=[]

for ax in [ax1,ax2,ax3,ax4]:
    ax.set_ylabel('$c_i$ (mol/L)')

for f in folders:

    read_all(tp,f,only=['all_data','species','descriptors','comsol_outputs','comsol_outputs_data'])
    conc=[[]]
    for sp in tp.species:
        conc.append([])
    desc=[]
    xx=[]
    min_dist=np.inf
    for key in tp.all_data:
        for key2 in tp.all_data[key]:
            for isp,sp in enumerate(tp.species):
                conc[isp].append([float(key),tp.all_data[key][key2]['species'][sp]['concentration'][0]/1000.])
#                conc[isp].append([float(key),tp.all_data[key][key2]['system']['cout'][-1][isp]/1000.])
                if abs(float(key)-voltage)<min_dist and isp==0:
                    min_dist=abs(float(key)-voltage)
                    voltage_found=round(float(key),3)
                    conc_x=[[]]
                    for sp in tp.species:
                        conc_x.append([])
                if isp>0:
                    conc_x[isp]=tp.all_data[key][key2]['species'][sp]['concentration']/1000.
                desc.append(float(key))
    for isp,sp in enumerate(tp.species):
        data=np.array(conc[isp])
        data.view(data.dtype.str+','+data.dtype.str).sort(order=['f0'], axis=0)
        x=data[:,0]
        y=data[:,1]
        ax1.plot(x+0.0592*pH,y,label=sp)
        ax2.semilogy(x+0.0592*pH,y,label=sp)
        print np.shape(conc_x[isp])
        ax3.plot(conc_x[isp],label=sp)
        ax4.loglog(conc_x[isp],label=sp)
ax3.set_title('Concentrations at '+str(voltage_found+0.0592*pH)+' V vs. RHE')
#ax4.set_title('Concentrations at '+str(voltage_found-0.0592*pH)+' V vs. RHE (loglog)')
ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()
plt.show()
