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
parser.add_argument('--desc',help='value of descriptor at which to plot property, nearest point is taken',nargs='+')
parser.add_argument('--norm',help='normalize electrode current density by RF',action='store_true')
args = parser.parse_args()

#args.desc=float(args.desc)

exp=EXPDATA()

cut=-2 #-0.78

def plot_leis_new_data(ax):
    #data=np.loadtxt('her_pcCu_lei.csv')
    i=0
    symbols=cycle(['x','1','o'])
    #for f in ['her_pcCu_lei.csv','her_pcCu_lei_NF.csv','her_pcCu_lei_NF_norm.csv']:
    ax.scatter(-0.33,0.125/390,marker='d',color='C3')
    ax.scatter(-0.33,0.05406/390,marker='D',color='C3')
    for f in ['her_pcCu_lei.csv','her_pcCu_lei_NF_norm.csv']:
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

c_list=['C'+str(i) for i in range(10)]
colors=cycle(c_list)
ls_list=['-',':','--']
linestyles=cycle(ls_list)
#m_list=['x','o','1','d','D','2']
m_list=['']
markerstyles=cycle(m_list)

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
fig=plt.figure(figsize=(4.6*n_row,3.2*n_col))
prop_inx={}
ax_list=[]
for i,p in enumerate(args.prop):
    ax_list.append(fig.add_subplot(str(n_col)+str(n_row)+str(i+1)))
    prop_inx[p]=i

def settings(ax,prop,d_sel):
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
    elif prop=='pKw':
        xlabel=r'x ($\AA$)'
        ylabel=r'pKw'
        label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
    elif prop=='efield':
        xlabel=r'x ($\AA$)'
        ylabel=r'E (V/m)'
        label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
    elif prop=='potential':
        xlabel=r'x ($\AA$)'
        ylabel=r'$\phi$ (V)'
        label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
    elif prop=='charge_density':
        xlabel=r'x ($\AA$)'
        ylabel=r'$\rho_{charge}$ (e/L)'
        label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
    else:
        xlabel=''
        ylabel=''
    ax.set_title(label)
    return xlabel,ylabel

def plot(prop):
    """create a plot of the property of interest"""
    inx=prop_inx[prop]
    ax=ax_list[inx]
    xmesh=tp.xmesh
    #if prop not in ['electrode_current_density','electrode_flux','surface_pH']:
    desc_iter=args.desc
    #else:
    #    desc_iter=[args.desc[0]]
    for desc in desc_iter:
        #reset colors
        a=''
        while a!=c_list[-1]:
            a=next(colors)
        m=next(markerstyles)
        d_sel=None
        if prop not in ['electrode_current_density','electrode_flux','surface_pH']:
            min_d=np.inf
            for i,d in enumerate(tp.descriptors['phiM']):
                if abs(d-float(desc))<min_d:
                    min_d=abs(d-float(desc))
                    d_sel_inx=i
                    d_sel=d
        xlabel,ylabel=settings(ax,prop,d_sel)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if prop in ['concentration','pKw']:
            #x: xmesh
            #species
            x=xmesh
            if prop=='pKw':
                color=next(colors)
                y=[a*b for a,b in zip(tp.alldata[d_sel_inx]['species']['H+']['concentration'],tp.alldata[d_sel_inx]['species']['OH-']['concentration'])]
                y=[-np.log10(yy/1000./1000.) for yy in y]
                ax.semilogx(x,y,ls+m,color=color,label='water diss')
            else:
                for sp in tp.species:
                    color=next(colors)
                    y=tp.alldata[d_sel_inx]['species'][sp][prop]
                    y=[yy/1000. for yy in y]
                    ax.semilogx(x,y,ls+m,color=color,label=sp)
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
                if args.norm and prop=='electrode_current_density':
                    RF=tp.system['RF']
                    y=[yy/RF for yy in y]
                ax.semilogy(x,y,ls+m,color=color,label=sp)
                if prop=='electrode_current_density':
                    ax=plot_leis_new_data(ax)
        elif prop in ['pH','potential','efield','charge_density','pKa']:
            #x: xmesh
            #system
            x=xmesh
            y=tp.alldata[d_sel_inx]['system'][prop]
            if prop=='charge_density':
                y=[yy/1000/unit_F for yy in y]
            #if prop=='pKa':
                #y=[-(14.5-8.49)/19.8*rho_c/unit_F/1000.+14.5 for rho_c in tp.alldata[d_sel_inx]['system']['charge_density']]
                
            ax.semilogx(x,y,ls+m,color='k')
        else:
            #x: descriptors
            #system
            x=tp.descriptors['phiM']
            #plot vs. RHE
            x=[xx+0.0592*tp.system['bulk_pH'] for xx in x]
            y=[tp.alldata[i]['system'][prop] for i in range(len(x))]
            ax.plot(x,y,ls+m,color='k')

for iif,f in enumerate(args.folder):
    print 'Working on folder ',f
    read_all(tp,f,only=['alldata','species','system','xmesh','descriptors','electrode_reactions'])
#    print tp.alldata[0]['system']['potential']
    ls=next(linestyles)
    a=''
    while a!=m_list[-1]:
        a=next(markerstyles)
    for p in args.prop:
        print('Plotting {}'.format(p))
        plot(p)
#for ax in ax_list:
#    ax.legend(prop={'size': 6})
plt.tight_layout()
plt.show()
