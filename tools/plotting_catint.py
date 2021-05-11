#tools for plotting various properties from transport simulations
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from catint.plot import Plot
from catint.transport import Transport
from glob import glob
from catint.catint_io import read_all
from itertools import cycle
import sys
from units import *
import re
from scipy.interpolate import UnivariateSpline

from catint.experimental import EXPDATA
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--file', help='folder to be evaluated',nargs='+')
parser.add_argument('--prop',help='property to be plotted',nargs='+')
#concentration
#potential
#field
#surface_concentration
parser.add_argument('--desc',help='value of descriptor at which to plot property, nearest point is taken',nargs='+')
parser.add_argument('--norm',help='normalize electrode current density by RF',action='store_true')
parser.add_argument('--scale',help='RHE or SHE')
parser.add_argument('--xfixed',help='x values at which to plot some quantities, called by "pH_at_x"',nargs='+')
args = parser.parse_args()

#args.desc=float(args.desc)

if args.xfixed is not None:
    args.xfixed=[float(xx) for xx in args.xfixed]

exp=EXPDATA()
color_species={}

if args.scale is None:
    args.scale='RHE'

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
ls_list=['-','--',':']
#ls_list=['-','--','--']
#ls_list=['--',':']
linestyles=cycle(ls_list)
#m_list=['x','o','1','d','D','2']
m_list=['']
markerstyles=cycle(m_list)

folders=sys.argv[1:]
tp=Transport(only_plot=True)

l=len(args.prop)
#round to divisible of 2
l=int( 2 * round( l / 2. ))
#l=int(l/2.)
n_row=int(np.ceil(np.sqrt(l)))
if n_row==0:
    n_row=1
n_col=l/n_row
if n_col==0:
    n_col=1
fig=plt.figure(figsize=(4.6*n_row,3.2*n_col))
prop_inx={}
ax_list=[]
for i,p in enumerate(args.prop):
    ax_list.append(fig.add_subplot(n_col,n_row,i+1)) #str(n_row)+str(i+1)))
    prop_inx[p]=i

def settings(ax,prop,d_sel):
    label=p.replace('_',' ')
    if prop=='concentration':
        ylabel=r'$c_i$ (M)'
        xlabel=r'x ($\AA$)'
        if args.scale=='RHE':
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
        else:
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel,3))
    elif prop=='electrode_current_density':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel=r'i (mA/cm$^2$)'
    elif prop=='surface concentration':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel=r'$c_i$ (x=0) (M)'
    elif prop=='pH':
        ylabel='pH'
        xlabel=r'x ($\AA$)'
        if args.scale=='RHE':
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
        else:
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel,3))
    elif prop=='surface_pH':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel='Surface pH'
    elif prop=='pH_at_x':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel='pH at x nm' #'+str(xfixed*1e9)+' nm'
    elif prop=='pKw':
        xlabel=r'x ($\AA$)'
        ylabel=r'pKw'
        if args.scale=='RHE':
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
        else:
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel,3))
    elif prop=='efield':
        xlabel=r'x ($\AA$)'
        ylabel=r'$E_x$ (V/\AA)'
        if args.scale=='RHE':
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
        else:
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel,3))
    elif prop=='potential':
        xlabel=r'x ($\AA$)'
        ylabel=r'$\phi$ (V)'
        if args.scale=='RHE':
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
        else:
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel,3))
    elif prop=='charge_density':
        xlabel=r'x ($\AA$)'
        ylabel=r'$\rho_{charge}$ (e/L)'
        if args.scale=='RHE':
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
        else:
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel,3))
    elif prop=='Stern_efield':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel=r'$E_{\mathrm{Stern},x}$ (V/\AA)'
    elif prop=='surface_potential':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel=r'phidg (V)'
    elif prop=='surface_charge_density':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel=r'$\sigma$ ($\mu$C/cm$^2$)'
    elif prop=='overpotential_dunwell':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel=r'$\eta_c$ ($V$)'
    elif prop=='Potential_drop':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel=r'$\phi^\mathrm{M}-\phi^\ddagger$ (V)'
    elif prop=='surface_efield':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel=r'$E^\ddagger$ (V)'
    elif prop=='Stern_epsilon_func':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel=r'$\varepsilon_\mathrm{S}$'
    elif prop=='activity':
        xlabel=r'x ($\AA$)'
        ylabel='$a_i = \gamma_i c_i$/(1 M)'
        if args.scale=='RHE':
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel+0.0592*tp.system['bulk_pH'],3))
        else:
            label+=r' at $\phi_M$ = {} V'.format(round(d_sel,3))
    elif prop=='surface_activity':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel='$a_i(x=0) = \gamma_i(x=0) c_i(x=0)$/(1 M)'
    elif prop=='surface_activity_coefficient':
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel='$\gamma_i(x=0)'
    elif 'concentration' in prop and 'at_x' in prop:
        xlabel='Voltage vs. '+args.scale+' (V)'
        ylabel=r'$c_\mathrm{'+prop.split('_')[2]+'}$ at x nm (M)'
    else:
        xlabel=''
        ylabel=''
    ax.set_title(label)
    return xlabel,ylabel

def plot(prop,color=None):
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
        #while a!=c_list[-1]:
        #    a=next(colors)
        #color=a
        if color is None:
            color=next(colors)
        else:
            color=color
        m=next(markerstyles)
        d_sel=None
        if prop not in ['electrode_current_density','electrode_flux','surface_pH','surface_concentration','reaction_rate','surface_activity','surface_activity_coefficient']:
            min_d=np.inf
            for i,d in enumerate(tp.descriptors['phiM']):
                if abs(d-float(desc))<min_d:
                    min_d=abs(d-float(desc))
                    d_sel_inx=i
                    d_sel=d
        if prop not in tp.alldata[0]['system'] and prop not in tp.alldata[0]['species'][[sp for sp in tp.alldata[0]['species']][0]] and (prop not in ['surface_charge_density','overpotential_dunwell','pH_at_x','capacitance','pKw','reaction_rate','Keq_buffer','surface_activity','surface_activity_coefficient','activity'] and not ('concentration' in prop and 'at_x' in prop)):
            return
        xlabel,ylabel=settings(ax,prop,d_sel)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        print('current prop',prop)
        if prop in ['concentration','pKw','activity','activity_coefficient']:
            #x: xmesh
            #species
            x=xmesh
            if prop=='pKw':
                if 'H+' not in tp.alldata[d_sel_inx]['species']:
                    continue
                if 'activity_coefficient' not in tp.alldata[d_sel_inx]['species']['OH-']:
                    continue
                y=[ca*ga*cb*gb for ca,ga,cb,gb in zip(tp.alldata[d_sel_inx]['species']['H+']['concentration'],tp.alldata[d_sel_inx]['species']['H+']['activity_coefficient'],tp.alldata[d_sel_inx]['species']['OH-']['concentration'],tp.alldata[d_sel_inx]['species']['OH-']['activity_coefficient'])]
                y=[-np.log10(yy/1000./1000.) for yy in y]
                ax.semilogx(x,y,ls+m,color=color,label='water diss (log)')
                ax2=ax.twiny()
                ax2.plot(x,y,ls+m,color='0.5',alpha=0.5,label='water diss')
            else:
                if prop in ['concentration','activity','activity_coefficient']:
                    func=ax.loglog #semilogy#plot #semilogx #loglog #semilogx #loglog
                    #func=ax.semilogy
                    #func=ax.plot
                else:
                    func=ax.semilogx
                for sp in tp.species:
                    if sp not in color_species:
                        ccolor=next(colors)
                        color_species[sp]=ccolor
                    else:
                        ccolor=color_species[sp]
                    print(d_sel_inx,'species',sp,prop)
                    if prop=='activity':
                        propread='concentration'
                    else:
                        propread=prop
                    y=tp.alldata[d_sel_inx]['species'][sp][propread]
                    if prop=='activity':
                        if 'activity_coefficient' not in tp.alldata[d_sel_inx]['species'][sp]:
                            continue
                        y=[a*b for a,b in zip(tp.alldata[d_sel_inx]['species'][sp]['concentration'],tp.alldata[d_sel_inx]['species'][sp]['activity_coefficient'])]
                    y=[yy/1000. for yy in y]
                    func(x,y,ls+m,color=ccolor,label=sp)
#                    ax.plot(x,y,ls+m,color=color,label=sp)
        elif prop in ['electrode_current_density','electrode_flux','surface_concentration','reaction_rate','Keq_buffer','surface_activity_coefficient','surface_activity']:
            #x: descriptors
            #species 
            x=tp.descriptors['phiM']
            #plot vs. RHE
            if args.scale=='RHE':
                x=[xx+0.0592*tp.system['bulk_pH'] for xx in x]
            else:
                x=[xx for xx in x]
            for sp in tp.species:
                if sp not in tp.electrode_reactions and prop not in ['surface_concentration','surface_activity_coefficient','surface_activity']:
                    continue
                if sp not in color_species:
                    ccolor=next(colors)
                    color_species[sp]=ccolor
                else:
                    ccolor=color_species[sp]
                if prop=='reaction_rate':
                    propread='electrode_current_density'
                elif prop=='surface_activity':
                    propread='surface_concentration'
                elif prop=='Keq_buffer':
                    propread='surface_concentration'
                else:
                    propread=prop
                if propread in tp.alldata[0]['species'][sp]:
                    y=[tp.alldata[i]['species'][sp][propread] for i in range(len(x))]
                else:
                    print(propread,'not found in all data, current keys = ',list(tp.alldata[0]['species'][sp].keys()))
                    print('skipping')
                    continue
                if args.norm and (prop=='electrode_current_density' or prop=='reaction_rate'):
                    RF=tp.system['RF']
                    y=[yy/RF for yy in y]
                    sp_diff='CO2' #HCO3- or CO2
                    #get index at certain x position
                    xfixed=0.#2e-10
                    d=np.inf
                    for ix,xx in enumerate(xmesh):
                        if d>abs(xx-xfixed):
                            d=abs(xx-xfixed)
                            inx=ix
                    #
                    #plot also transport limited current
                    ct=tp.alldata[i]['species'][sp_diff]['concentration'][:inx+3]
                    dc_dx=(ct[inx+1]-ct[inx])/(x[inx+1]-x[inx])
#                    dc_dt_2=(ct[2]-2*ct[1]+ct[0])/(x[1]-x[0])**2 #dc/dt
                    i=-dc_dx*tp.species[sp_diff]['diffusion'] #flux in mol/m^2/s
                    nel=2.
                    nprod=1.
                    j=i*nel*unit_F/nprod/10. #current density in mA/cm^2
                    print('slope',dc_dx,j)
                    print('limiting current = {}'.format(j))
                    ax.axhline(j,color='0.5')
                    func=ax.semilogy
                    if prop=='electrode_current_density':
                        func(x,y,ls+m,color=ccolor,label=sp)
                    elif prop=='reaction_rate':
                        ax.plot(tp.species['HCO3-']['bulk_concentration']/1000.,[yy for xx,yy in zip(x,y) if xx==-0.9][0],'o',color='k')
                elif prop=='surface_concentration':
                    y=[(yy)/1000. for yy in y]
                    func=ax.semilogy #ax.plot
                    func(x,y,ls+m,color=ccolor,label=sp)
                elif prop=='surface_activity_coefficient':
                    #y=[(yy)/1000. for yy in y]
                    func=ax.semilogy
                    func(x,y,ls+m,color=ccolor,label=sp)
                elif prop=='surface_activity':
                    try:
                        y=[tp.alldata[i]['species'][sp]['surface_activity_coefficient']*tp.alldata[i]['species'][sp]['surface_concentration'] for i in range(len(x))]
                        func=ax.semilogy
                        func(x,y,ls+m,color=ccolor,label=sp)
                    except:
                        pass
                if prop=='surface_concentration':
                    for xx,yy in zip(x,y):
                        print('sc',sp,xx,yy)
                #if prop=='electrode_current_density':
                #    ax=plot_leis_new_data(ax)
            if prop=='Keq_buffer':
                hco3=np.array([tp.alldata[i]['species']['HCO3-'][propread] for i in range(len(x))])
                co2=np.array([tp.alldata[i]['species']['CO2'][propread] for i in range(len(x))])
                co3=np.array([tp.alldata[i]['species']['CO32-'][propread] for i in range(len(x))])
                oh=np.array([tp.alldata[i]['species']['OH-'][propread] for i in range(len(x))])
                h=np.array([tp.alldata[i]['species']['H+'][propread] for i in range(len(x))])
#                y=[(yy)/1000. for yy in y]
                func=ax.semilogy
                func(x,hco3/(oh*co2),'-')
                ax.axhline(44400.0)
        elif prop in ['pH','potential','efield','charge_density','pKa']:
            #x: xmesh
            #system
            x=xmesh
            y=tp.alldata[d_sel_inx]['system'][prop]
            if prop=='charge_density':
                y=[yy/1000/unit_F for yy in y]
            elif prop=='efield':
                y=[yy*1e-10 for yy in y]
            #if prop=='pKa':
                #y=[-(14.5-8.49)/19.8*rho_c/unit_F/1000.+14.5 for rho_c in tp.alldata[d_sel_inx]['system']['charge_density']]
            func=ax.semilogx #plot
            func(x,y,a+m,color=color) #next(colors))
        else:
            #x: descriptors
            #system
            x=tp.descriptors['phiM']
            #plot vs. RHE
            x_she=x
            if args.scale=='RHE':
                x=[xx+0.0592*tp.system['bulk_pH'] for xx in x]
            else:
                x=[xx for xx in x]
            if prop=='surface_charge_density':
                #calculate surface charge density from surface potential
                if tp.system['charging_scheme']=='comsol':
                    y=[tp.system['Stern capacitance']*(x_she[i]-tp.alldata[i]['system']['surface_potential']-tp.system['phiPZC']) for i in range(len(x))]
                elif tp.system['charging_scheme']=='input':
                    #read the 'sigma_input' entry from catmap input file
                    if sigma_input_str.strip().split('=')[-1].startswith('['):
                        sigma_input=[float(sis) for sis in sigma_input_str.strip().split('=')[-1].strip('[').strip(']').split(',')]
                    p=np.poly1d(sigma_input)
                    y=[p(x_she[i]) for i in range(len(x_she))]
            elif prop=='capacitance':
                #plot the double layer capacitance
                #1) calculate surface charge density from surface potential
                y=[tp.system['Stern capacitance']*(x_she[i]-tp.alldata[i]['system']['surface_potential']-tp.system['phiPZC']) for i in range(len(x))]
                #2) spline sigma vs. v
                data=np.array(sorted([[xx,yy] for xx,yy in zip(x,y)]))
                x=data[:,0]
                y=data[:,1]
                print(x,y)
                try:
                    spl=UnivariateSpline(x,y,k=2,s=0)
                    #3) calculate double layer capacitance from derivative
                    spl_deriv=spl.derivative()
                    y=[spl_deriv(xx) for xx in x]
                except:
                    y=x
            elif 'concentration' in prop and 'at_x' in prop:
                sp=prop.split('_')[1]
                for xfixed in args.xfixed:
                    d=np.inf
                    for ix,xx in enumerate(xmesh):
                        if d>abs(xx-xfixed):
                            d=abs(xx-xfixed)
                            inx=ix
                    sp='HCO3-' #K+'
                    y=np.array([tp.alldata[i]['species'][sp]['concentration'][inx]/1000. for i in range(len(x))])
                    ax.plot(x,y,ls+m,label='at '+str(xfixed*1e9)+' nm')
            elif prop=='pH_at_x':
                for xfixed in args.xfixed:
                    d=np.inf
                    for ix,xx in enumerate(xmesh):
                        if d>abs(xx-xfixed):
                            d=abs(xx-xfixed)
                            inx=ix
                    y=np.array([tp.alldata[i]['system']['pH'][inx] for i in range(len(x))])
                    ax.plot(x,y,ls+m,label='at '+str(xfixed*1e9)+' nm')
            elif prop=='overpotential_dunwell':
                yco2_surf=np.array([tp.alldata[i]['species']['CO2']['surface_concentration'] for i in range(len(x))])
                yh_surf=1e-14/np.array([tp.alldata[i]['species']['OH-']['surface_concentration'] for i in range(len(x))])
                yco2_bulk=np.array([tp.species['CO2']['bulk_concentration'] for i in range(len(x))])
                yh_bulk=1e-14/np.array([tp.species['OH-']['bulk_concentration'] for i in range(len(x))])
                y=2.3*unit_R*unit_T/unit_F*np.log(yco2_bulk*yh_bulk/yco2_surf/yh_surf)
            else:
                y=[tp.alldata[i]['system'][prop] for i in range(len(x))]
            if prop in ['Stern_efield']:
                y=[yy*1e-10 for yy in y]
            if not prop=='pH_at_x':
                ax.plot(x,y,ls+m,color=color)
        ax.set_xlim(-1.25,-0.5)

for iif,f in enumerate(args.file):
    print('Working on folder ',f)
    color=next(colors)
    read_all(tp,f,only=['alldata','species','system','xmesh','descriptors','electrode_reactions'])
    #check how many real datapoints we actually have:
    n_conv=0
    for i in range(len(tp.alldata)):
        key0=[key for key in tp.alldata[i]['species']][0]
        if len(tp.alldata[i]['species'][key0])>0:
            n_conv+=1
    tp.descriptors['phiM']=tp.descriptors['phiM'][:n_conv]
    tp.alldata=tp.alldata[:n_conv]
    #read sigma_input, if charging_scheme is defined in catmap input file
    if tp.system['charging_scheme']=='input':
        for line in open(glob(f+'/catmap_input/desc*/catmap_CO2R.mkm')[0],'r'):
            if line.strip().startswith('sigma_input'):
                sigma_input_str=line
#    sys.exit()
#    print tp.alldata[0]['system']['potential']
    ls=next(linestyles)
    a=''
    while a!=m_list[-1]:
        a=next(markerstyles)
    for p in args.prop:
        print(('Plotting {}'.format(p)))
        plot(p,color=color)
for ax in ax_list:
    ax.legend(prop={'size': 6})
plt.tight_layout()
#plt.savefig('test.pdf')
plt.show()
