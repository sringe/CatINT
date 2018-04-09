#tools for plotting various properties from transport simulations
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

from catint.experimental import EXPDATA
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--folder', help='folder to be evaluated',nargs='+')
parser.add_argument('--prop',help='property to be plotted',nargs='+')
parser.add_argument('--desc',help='value of descriptor at which to plot property, nearest point is taken')
args = parser.parse_args()

args.desc=float(args.desc)

exp=EXPDATA()

cut=-2 #-0.78

def plot_leis_new_data(ax):
    #data=np.loadtxt('her_pcCu_lei.csv')
    i=0
    symbols=cycle(['x','1','o'])
    for f in ['her_pcCu_lei.csv','her_pcCu_lei_NF.csv','her_pcCu_lei_NF_norm.csv']:
        symbol=next(symbols)
        i+=1
        data=np.loadtxt(f)
        label='exp, '
        if 'NF' in f and 'norm' in f:
            color='C0'
            label+='x390, normalized'
        elif 'NF' in f:
            color='C3'
            label+='x390'
        else:
            color='k'
        ax.semilogy(data[:,0],10**(data[:,1]),symbol,color=color,label=label)
    return ax

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

colors=cycle(['C'+str(i) for i in range(10)])
linestyles=cycle(['-',':','--'])

folders=sys.argv[1:]
tp=Transport(only_plot=True)

l=len(args.prop)
#round to divisible of 2
l=int( 2 * round( l / 2. ))
n_col=l/2 #int(np.ceil(np.sqrt(l)))
if n_col==0:
    n_col=1
n_row=l/2 #int(np.ceil(np.sqrt(l)))
if n_row==0:
    n_row=1
fig=plt.figure(figsize=(5*n_row,3.5*n_col))
prop_inx={}
ax_list=[]
for i,p in enumerate(args.prop):
    ax_list.append(fig.add_subplot(str(n_col)+str(n_row)+str(i+1)))
    prop_inx[p]=i

def settings(ax,prop):
    label=p.replace('_',' ')
    if prop=='concentration':
        ylabel=r'$c_i$ (M)'
        xlabel=r'x ($\AA$)'
        label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
    elif prop=='electrode_current_density':
        xlabel='Voltage vs. RHE (V)'
        ylabel=r'i (mA/cm$^2$)'
    elif prop=='pH':
        ylabel='pH'
        xlabel=r'x ($\AA$)'
        label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
    elif prop=='surface_pH':
        xlabel='Voltage vs. RHE (V)'
        ylabel='Surface pH'
    else:
        xlabel=''
        ylabel=''
    ax.set_title(label)
    return xlabel,ylabel

def plot(prop,d_sel_inx,ls):
    """create a plot of the property of interest"""
    #reset colors
    a=''
    while a!='C9':
        a=next(colors)
    inx=prop_inx[prop]
    ax=ax_list[inx]
    xlabel,ylabel=settings(ax,prop)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    xmesh=tp.xmesh
    if prop in ['concentration']:
        #x: xmesh
        #species
        x=xmesh
        for sp in tp.species:
            color=next(colors)
            y=tp.alldata[d_sel_inx]['species'][sp][prop]
            y=[yy/1000. for yy in y]
            ax.semilogx(x,y,ls,color=color,label=sp)
    elif prop in ['electrode_current_density','electrode_flux']:
        #x: descriptors
        #species 
        color=next(colors)
        x=tp.descriptors['phiM']
        #plot vs. RHE
        x=[xx+0.0592*tp.system['bulk_pH'] for xx in x]
        for sp in tp.species:
            if sp not in tp.electrode_reactions:
                continue
            color=next(colors)
            y=[tp.alldata[i]['species'][sp][prop] for i in range(len(x))]
            print 'checking electrode_current density'
            print y
            ax.semilogy(x,y,ls,color=color,label=sp)
            if prop=='electrode_current_density':
                ax=plot_leis_new_data(ax)
    elif prop in ['pH']:
        #x: xmesh
        #system
        x=xmesh
        y=tp.alldata[d_sel_inx]['system'][prop]
        ax.semilogx(x,y,ls,color='k')
    else:
        #x: descriptors
        #system
        x=tp.descriptors['phiM']
        #plot vs. RHE
        x=[xx+0.0592*tp.system['bulk_pH'] for xx in x]
        y=[tp.alldata[i]['system'][prop] for i in range(len(x))]
        ax.plot(x,y,ls,color='k')

for iif,f in enumerate(args.folder):
    ls=next(linestyles)
    read_all(tp,f,only=['alldata','species','system','xmesh','descriptors','electrode_reactions'])
    min_d=np.inf
    for i,d in enumerate(tp.descriptors['phiM']):
        if abs(d-args.desc)<min_d:
            min_d=abs(d-args.desc)
            d_sel_inx=i
            d_sel=d
    for p in args.prop:
        print('Plotting {}'.format(p))
        plot(p,d_sel_inx,ls)
for ax in ax_list:
    ax.legend(prop={'size': 6})
plt.tight_layout()
plt.show()
