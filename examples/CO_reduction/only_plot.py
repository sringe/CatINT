import numpy as np
import matplotlib.pyplot as plt
from catint.plot import Plot
from catint.transport import Transport
from glob import glob
from catint.io import read_all
from itertools import cycle
import sys

def plot_leis_data():
    data=np.loadtxt('data/CORR_lei.txt')
    species=['CH4','CH2CH2','CH3CHOO-','EtOH','n-PrOH']
    for i in range(5):
        print data[:,0],data[:,10+i]
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
        plt.semilogy(datap[:,0],datap[:,1],'-o',color=color,label=species[i])
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

colors=cycle(['orange','b']) #,'k','r','y'])
styles=cycle(['-','--','-.',':'])

#folders=glob('calc_std_settings*') #obj') #calc_std*')
folders=sys.argv[1:] #glob('test') #calc_std_settings*') #obj') #calc_std*')
tp=Transport()

for f in folders:
    style=next(styles)
#    p=Plot(transport=tp,init_from_file='calc_std_settings')
    read_all(tp,f,only=['electrode_reactions','species']) #system')
    for sp in tp.electrode_reactions:
        color=next(colors)
        jout=tp.electrode_reactions[sp]['electrode_current_density']
        data=[]
        for v1,v2 in jout:
            data.append([float(v1),jout[(v1,v2)]])
        data=np.array(data)
#        data=np.sort(a.view(a.dtype.str+','+a.dtype.str), order=['f1'], axis=0).view(np.float32)
        data.view(data.dtype.str+','+data.dtype.str).sort(order=['f0'], axis=0)
        plt.semilogy(data[:,0],data[:,1],linestyle=style,color=color,label=f.replace('calc_std_settings','std')+', '+sp)
        #plt.semilogy(data[:,0],data[:,1],'o')
plt.ylim([1e-8,5e2])
plt.xlim([-1.03,0.0])
plt.xlabel('Voltage vs RHE (V)')
plt.ylabel('j (mV/cm^2')
plot_leis_data()
#plot_xinyans_equation()
plt.legend()
plt.show()

for arg in sys.argv[1:]:
    p=Plot(transport=tp,init_from_file=arg)
#    read_all(tp,f,only=['electrode_reactions','species']) #system')
    p.plot(large_plots=['concentrations_electrolyte','concentrations_electrode'])
    plt.show()
