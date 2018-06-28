import numpy as np
import sys
import re
import matplotlib.pyplot as plt
from itertools import cycle


colors=cycle(['C'+str(i) for i in range(10)])
fig=plt.figure(figsize=(7,6))
ax1=fig.add_subplot('211')
ax2=fig.add_subplot('212')
for arg in sys.argv[1:]:
    color=next(colors)
    accs=[]
    curs=[]
    phis=[]
    tmp=[]
    tmp2={}
    for line in open(arg+'/transport_id000.log','r'):
        phi=re.findall('Starting calculation for phiM = (-?\d+\.\d+)',line)
        if len(phi)>0:
            if len(tmp)>0:
                accs.append(tmp)
            if len(tmp2)>0:
                curs.append(tmp2)
            tmp=[]
            tmp2={}
            cphi=float(phi[0])
            phis.append(cphi)
        acc=re.findall('Current accuracy in current density = (.*) mA',line)
        if len(acc)>0:
            if acc[0] != 'inf': #float(acc[0])>0:
                tmp.append(float(acc[0]))
#        cur=re.findall('\| Current current density = (.*) ',line)
#        cur=re.findall('Current Density of ([a-zA-Z-^+0-9]+) = ([0-9.-e]+) mA\/cm\^2',line)
        cur=re.findall('Current Density of ([a-zA-Z-^+0-9]+) = (-?[0-9]+.?[0-9]+e?-?\d{0,2}) mA\/cm\^2',line)
        #Current Density of HCO3- = 0.0 mA/cm^2
        if len(cur)>0:
            cur=cur[0]
            if cur[1] != 'inf': #float(cur[0])>0:
                if cur[0] not in tmp2:
                    tmp2[cur[0]]=[]
                tmp2[cur[0]].append(float(cur[1]))
    accs.append(tmp)
    curs.append(tmp2)
    i=-1
    for acc in accs:
        i+=1
 #       if i>0:
 #           break
        ax1.semilogy(acc,'-o',color=color,label=str(phis[i]))
    ax1.legend()
    i=-1
    for cur in curs:
        i+=1
        print [key for key in cur]
#        if i>0:
#            break
        for key in cur:
            ax2.semilogy([abs(a) for a in cur[key]],'-o',label=key+','+str(phis[i]))
    ax2.legend()
ax1.set_ylabel('RMSD in Current Density (mV/cm^2)')
ax2.set_ylabel('Current Density (mV/cm^2)')
ax1.set_xlabel('Iteration')
ax2.set_xlabel('Iteration')
plt.show()
