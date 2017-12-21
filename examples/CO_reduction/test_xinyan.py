#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

mathtext_prop = {'fontset' : 'custom',
                 'it' : 'serif:italic',
                 'sf' : 'Arial:bold',
                 'cal' : 'serif:italic:bold'}
rc('font', family='serif', serif='Arial')
rc('mathtext', **mathtext_prop)

def C1_rate(voltage, CO_cvg):
    Ga_CHO, Ga_CHOH = 1.11746219, 2.37467774
    kT, A = 0.02585199, 1.e13
    return CO_cvg*A*np.exp(-max([Ga_CHO+0.5*voltage, Ga_CHOH+2*voltage])/kT)

def CHO_rate(voltage, CO_cvg):
    Ga_CHO, Ga_CHOH = 1.11746219, 2.37467774
    kT, A = 0.02585199, 1.e13
    return CO_cvg*A*np.exp(-(Ga_CHO+0.5*voltage)/kT)

def HCOH_rate(voltage, CO_cvg):
    Ga_CHO, Ga_CHOH = 1.11746219, 2.37467774
    kT, A = 0.02585199, 1.e13
    return CO_cvg*A*np.exp(-(Ga_CHOH+2*voltage)/kT)

def OCCOH_rate(voltage, CO_cvg):
    Ga_OCCO, Ga_OCCOH = 0.578959276, 1.10495851
    kT, A = 0.02585199, 1.e13
    return CO_cvg**2*A*np.exp(-(Ga_OCCOH+0.5*voltage)/kT)

def OCCO_rate(voltage, CO_cvg):
    Ga_OCCO, Ga_OCCOH = 0.578959276, 1.10495851
    kT, A = 0.02585199, 1.e13
    return CO_cvg**2*A*np.exp(-(Ga_OCCO)/kT)

def C2_rate(voltage, CO_cvg):
    Ga_OCCO, Ga_OCCOH = 0.578959276, 1.10495851
    kT, A = 0.02585199, 1.e13
    return CO_cvg**2*A*np.exp(-max([Ga_OCCOH+0.5*voltage, Ga_OCCO])/kT)
    
voltage = np.linspace(-1.6, 0.4, 101)
CO_cvg = 0.415139754

fig = plt.figure(1, figsize=(7.5, 5))
ax = fig.add_subplot(1,1,1)
plt.plot(voltage, [CHO_rate(v, CO_cvg) for v in voltage], '-', color='k', label = 'CHO_rate')
plt.plot(voltage, [HCOH_rate(v, CO_cvg) for v in voltage], '-', color='r', label = 'HCOH_rate')
plt.plot(voltage, [OCCOH_rate(v, CO_cvg) for v in voltage], '-', color='g', label = 'OCCOH_rate')
plt.plot(voltage, [OCCO_rate(v, CO_cvg) for v in voltage], '-', color='y', label = 'OCCO_rate')
plt.plot(voltage, [C1_rate(v, CO_cvg) for v in voltage], 'o', color='k', label = 'C1_rate')
plt.plot(voltage, [C2_rate(v, CO_cvg) for v in voltage], 'o',  color='r', label = 'C2_rate')
plt.xlim((-1.6, 0.4))
plt.ylim((1e-5,1e3))
plt.xlabel(r'$Voltage$ [V]', fontsize=16)
plt.ylabel(r'$j\ [mA/cm^{2}]$', fontsize=16)
ax.set_yscale('log')
leg = plt.legend(loc='best', ncol=1, prop={'size':10}, fancybox=True, shadow=False, numpoints=1)
leg.get_frame().set_facecolor('none')
leg.get_frame().set_linewidth(0.0)
fig.savefig('effective_rate.pdf')
plt.show()
#plt.close()
