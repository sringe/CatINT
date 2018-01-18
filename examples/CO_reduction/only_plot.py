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

SA=1

def plot_leis_high_surface_data():
    data=[]
    for line in open('data/COR_high_surface.txt'):
        if line.startswith('#'):
            header=line.split()
        else:
            data.append(line.split())
            nmol=len(line.split())
    data=np.array(data)
    for i in range(nmol):
        if header[i]=='etol':
            print data[:,i]
            plt.semilogy(data[:,0],[float(d)/SA for d in data[:,i]],'--o',color='red',label=header[i])

def plot_leis_data():
    data=np.loadtxt('data/CO2R_lei.txt')
    species=['CH4','CH2CH2','CH3CHOO-','EtOH','n-PrOH']
    for i in range(5):
        x=data[:,0]
        y=data[:,10+i]
        datap=[]
        for xx,yy in zip(x,y):
            if yy!=0.0:
                datap.append([xx,-yy])
        datap=np.array(datap)
        if species[i]=='CH4':
            color='b'
        else:
            color='orange'
        if species[i] not in ['CH4','EtOH']:
            continue
        plt.semilogy(datap[:,0],datap[:,1],'--o',color=color,label=species[i])

def plot_kanans_data():
    data=np.loadtxt('data/COR.txt')
    species=['H2','CH4','CH2OO-','CO','MeOH','CH2CH2','EtOH', 'GlyAl','AcAl','Acetate','EtGly','n-PrOH','AllylOH','PrAl','Acetone','OHAc']
    for i in range(len(species)):
        x=data[:,0]
        y=data[:,i+1]
        datap=[]
        for xx,yy in zip(x,y):
            if yy!=0.0:
                datap.append([xx,-yy])
        datap=np.array(datap)
        if species[i] not in ['CH4','EtOH']: #['CH2CH2','CH4','EtOH','n-PrOH','Acetate']:
            continue
        if species[i] in ['CH4','MeOH','CH2OO-']:
            color='b'
        elif species[i]=='H2':
            color='k'
        elif species[i]=='CO':
            color='y'
        else:
            color='orange'
        plt.semilogy(datap[:,0],datap[:,1],'-o',color=color) #,label=species[i])

rho_act=1.004495558139274e-05



def C1_rate(voltage, CO_cvg):
    Ga_CHO, Ga_CHOH = 1.11746219, 2.37467774
    kT, A = 0.02585199, 1.e13
    return rho_act*CO_cvg*A*np.exp(-max([Ga_CHO+0.5*voltage, Ga_CHOH+2*voltage])/kT)

def C2_rate(voltage, CO_cvg):
    Ga_OCCO, Ga_OCCOH = 0.578959276, 1.10495851
    kT, A = 0.02585199, 1.e13
    return rho_act*CO_cvg**2*A*np.exp(-max([Ga_OCCOH+0.5*voltage, Ga_OCCO])/kT)

def plot_xinyans_equation():
    voltage = np.linspace(-1.6, 0.4, 101)
    CO_cvg = 0.415139754
    plt.semilogy(voltage, [C1_rate(v, CO_cvg) for v in voltage], '-', color='b', label = 'C1_rate')
    plt.semilogy(voltage, [C2_rate(v, CO_cvg) for v in voltage], '-', color='orange', label = 'C2_rate')

colors=cycle(['orange','lightblue','r','blue','darkred','black']) #,'k','r','y'])
#colors=cycle(['orange','b','k','r','y'])
styles=cycle(['-','--']) #,'-.',':','o','d','x'])

#folders=glob('calc_std_settings*') #obj') #calc_std*')
folders=sys.argv[1:] #glob('test') #calc_std_settings*') #obj') #calc_std*')
tp=Transport()
oldlabel=''
for f in folders:
    print f
    pH=re.findall('pH(\d+)',f)
    if len(pH)>0:
        newlabel=re.findall('pH(\d+)',f)[0]
    else:
        newlabel='none'
    print newlabel
    if newlabel!=oldlabel:
        style=next(styles)
#    if newlabel=='6.8':
#        style=styles[0]
#    else: #if newlabel=='13.0':
#        style=styles[1]
    oldlabel=newlabel
#    style=next(styles)
#    p=Plot(transport=tp,init_from_file='calc_std_settings')
    read_all(tp,f,only=['electrode_reactions','species','comsol_outputs','comsol_outputs_data']) #system')
    for sp in tp.electrode_reactions:
        color=next(colors)
        jout=tp.electrode_reactions[sp]['electrode_current_density']
        data=[]
        for v1,v2 in jout:
            data.append([float(v1),jout[(v1,v2)]])
        data=np.array(data)
        if 'CD' in f:
            data[:,1]=data[:,1]/SA
#        data=np.sort(a.view(a.dtype.str+','+a.dtype.str), order=['f1'], axis=0).view(np.float32)
        data.view(data.dtype.str+','+data.dtype.str).sort(order=['f0'], axis=0)
#        print 'data',f,sp,data[0,1]
        plt.semilogy(data[:,0],data[:,1],style,color=color) #,label=f.replace('calc_std_settings','std')+', '+sp)
        #plt.semilogy(data[:,0],data[:,1],'o')
    symbols=['o','x','d','D']
    isp=-1
    for sp in [o[1] for o in tp.comsol_outputs if o[1] != 'coverage']:
        isp+=1
        symbol=symbols[isp]
        jout=tp.comsol_outputs_data[sp]
        data=[]
        for v1,v2 in jout:
            data.append([float(v1),jout[(v1,v2)][0]])
        data=np.array(data)
        data.view(data.dtype.str+','+data.dtype.str).sort(order=['f0'], axis=0)
        if sp in ['CHO','HCOH']:
            nel=6
            nc=1
        else:
            nel=8
            nc=2
#        plt.semilogy(data[:,0],data[:,1]*nel*unit_F/10,symbol,color=color,label=f.replace('calc_std_settings','std')+', '+sp)
plt.ylim([1e-8,5e2])
plt.xlim([-1.03,0.0])
plt.xlabel('Voltage vs RHE (V)')
plt.ylabel(r'$j$ (mA/cm$^2$)')
plot_leis_data()
#plot_leis_high_surface_data()
plot_kanans_data()
#plot_xinyans_equation()
plt.legend()
plt.show()
plt.close() #empty()
sys.exit()
for arg in sys.argv[1:]:
    p=Plot(transport=tp,init_from_file=arg)
    read_all(tp,f) #,only=['','species','system','all_data']) #electrode_reactions','species']) #system')
    p.plot(large_plots=['concentrations_electrolyte','concentrations_electrode'],small_plots=['potential'])
#    plt.show()
