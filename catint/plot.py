import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle
from units import *
import sys
from matplotlib import colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
import os
from io import read_all

class Plot():

    def __init__(self,transport=None,logscale=False,init_from_file=None):
        if transport==None:
            print('No transport object provided for calculator. Stopping here.')
            sys.exit()
        else:
            self.tp=transport #transport object
        if init_from_file is not None:
            #load data from pickle
            read_all(self.tp,init_from_file)
        self.logscale=logscale

    def plot(self,large_plots=[],small_plots=[]):
        self.large_plots=large_plots
        self.small_plots=small_plots

        outdir='catint_plots'

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        if self.tp.descriptors is not None:
            data_desc={}
            #copy descriptor based data to single arrays for plotting
            desc_keys=[key for key in self.tp.descriptors]
            if any([a.startswith('desc_') for a in sum([self.large_plots,self.small_plots],[])]):
                #a plot of a quantity should be performed with the descriptor as x axis
                plot_names=['_'.join(a.split('_')[1:]) for a in sum([self.large_plots,self.small_plots],[]) if a.startswith('desc_')]
                for p in plot_names:
                    data_desc[p]=[]
            i1=0
            i2=0
            for value1 in self.tp.descriptors[desc_keys[0]]:
                i1+=1
                i2=0
                for value2 in self.tp.descriptors[desc_keys[1]]:
                    i2+=1
                    self.tp.potential=self.tp.all_data[str(value1)][str(value2)]['system']['potential'].copy()
                    self.tp.efield=self.tp.all_data[str(value1)][str(value2)]['system']['efield'].copy()
                    self.tp.current_density=self.tp.all_data[str(value1)][str(value2)]['system']['current_density'].copy()
                    self.tp.species=self.tp.all_data[str(value1)][str(value2)]['species'].copy()
                    self.tp.cout=self.tp.all_data[str(value1)][str(value2)]['system']['cout'].copy()
                    self.tp.electrode_flux=self.tp.all_data[str(value1)][str(value2)]['system']['electrode_flux'].copy()
                    if any([not a.startswith('desc_') for a in sum([self.large_plots,self.small_plots],[])]):
                        self.plot_single()
                        plt.savefig(outdir+'/results_'+str(i1)+str(i2)+'.pdf')
                        plt.savefig(outdir+'/results_'+str(i1)+str(i2)+'.png')
                        plt.close()
                    if any([a.startswith('desc_') for a in sum([self.large_plots,self.small_plots],[])]):
                        for p in plot_names:
                            if p=='current_density':
                                current_density=0
                                all_currents=[]
                                for sp in self.tp.electrode_reactions:
                                    nprod=len([a for a in self.tp.electrode_reactions[sp]['reaction'][1] if a==sp])
                                    isp=[i for i,sp2 in enumerate(self.tp.species) if sp2==sp][0]
                                    c_current_density=self.tp.electrode_flux[-1][isp*self.tp.nx]*self.tp.electrode_reactions[sp]['nel']*unit_F/nprod/10.
                                    current_density+=c_current_density
                                    all_currents.append(c_current_density)
                                data_desc[p].append(sum([[value1,value2,current_density],all_currents],[]))
            #convert to np.arrays
            for p in plot_names:
                data_desc[p]=np.array(data_desc[p])
            if any([a.startswith('desc_') for a in sum([self.large_plots,self.small_plots],[])]):
                for p in plot_names:
                    plt.figure()
                    plt.ylabel('i (mV/cm^2)')
                    plt.xlabel('Voltage vs RHE (V)') 
                    plt.semilogy(data_desc[p][:,0],data_desc[p][:,2],'-o',label='total')
                    for isp,sp in enumerate(self.tp.electrode_reactions):
                        plt.semilogy(data_desc[p][:,0],data_desc[p][:,3+isp],'-o',label=sp)
                    plt.legend()
                    plt.savefig(outdir+'/results_'+p+'.pdf')
                    plt.savefig(outdir+'/results_'+p+'.png')
        else:
            self.plot_single()
            plt.savefig(outdir+'/results.pdf')
            plt.savefig(outdir+'/results.png')
#        plt.show()

    def plot_single(self):

        cout=self.tp.cout
        electrode_flux=self.tp.electrode_flux

        colors=dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
        colorlist=cycle([colors[key] for key in colors])

        fig = plt.figure(figsize=(14,9)) #gcf()

        flatten = lambda l: [item for sublist in l for item in sublist]

        gs1 = gridspec.GridSpec(3,4)
        ax1 = fig.add_subplot(gs1[0:2,0:2])
        ax2 = fig.add_subplot(gs1[0:2,2:4])

        product_list=self.tp.educt_list+self.tp.product_list #['HCOO-','CO','H2','etol','propol','metol','C2H4','CH4','allyl']

        color_offset=0.3
#        print 'cout in plot', self.tp.all_data[str(desc_val[0])][str(desc_val[1])]['system']['cout'][0,9*self.tp.nx:10*self.tp.nx-1]
        ax3 = fig.add_subplot(gs1[8])
        ax4 = fig.add_subplot(gs1[9])
        ax5 = fig.add_subplot(gs1[10])
        ax6 = fig.add_subplot(gs1[11])

        ax_large_plots_used=[False,False] #ax1,ax2]
        ax_large_plots=[ax1,ax2]
        ax_small_plots_used=[False,False,False,False] #=[ax3,ax4,ax5,ax6]
        ax_small_plots=[ax3,ax4,ax5,ax6]

        #concentrations

        def get_current_ax(method):
            if method in self.large_plots:
                window='large'
            elif method in self.small_plots:
                window='small'
            if window=='large':
                for iw,w in enumerate(ax_large_plots_used):
                    if not w:
                        ax_large_plots_used[iw]=True
                        break
                return ax_large_plots[iw]
            elif window=='small':
                for iw,w in enumerate(ax_small_plots_used):
                    if not w:
                        ax_small_plots_used[iw]=True
                        break
                return ax_small_plots[iw]

        if any([a.startswith('concentrations') for a in sum([self.large_plots,self.small_plots],[])]) or 'pH' in sum([self.large_plots,self.small_plots],[]) or 'electrode_flux' in sum([self.large_plots,self.small_plots],[]):
            went_to_elec=False
            went_to_prod=False
            went_to_elflux=False
            for k,sp in enumerate(self.tp.species):
                if not any([a.startswith('concentrations') for a in sum([self.large_plots,self.small_plots],[])]) and sp not in ['OH-','H+'] and not 'electrode_flux' in sum([self.large_plots,self.small_plots],[]):
                    continue
                color=next(colorlist)
                if sp in ['OH-','H+']:
                    if not 'pH' in sum([self.large_plots,self.small_plots],[]):
                        continue
                    else:
                        ax_pH=get_current_ax('pH')
                        ax_pH.set_title('pH')
                        ax_pH.set_xlabel('x (m)')
                        ax_pH.set_ylabel('pH')
                if 'electrode_flux' in sum([self.large_plots,self.small_plots],[]):
                    if not went_to_elflux:
                        ax_elflux=get_current_ax('electrode_flux')
                        ax_elflux.set_title('Electrode Reaction Flux')
                        ax_elflux.set_xlabel('x (m)')
                        ax_elflux.set_ylabel('j (mol/m^2/s)')
                        went_to_elflux=True
                if sp in product_list:
                    if not 'concentrations_reaction' in sum([self.large_plots,self.small_plots],[]):
                        continue
                    else:
                        if not went_to_prod:
                            ax=get_current_ax('concentrations_reaction')
                            ax.set_title('Educt/Product Concentrations')
                            ax_p=ax
                        else:
                            ax=ax_p
                    went_to_prod=True
                else:
                    if not 'concentrations_electrolyte' in sum([self.large_plots,self.small_plots],[]):
                        continue
                    else:
                        if not went_to_elec:
                            ax=get_current_ax('concentrations_electrolyte')
                            ax.set_title('Electrolyte Concentrations')
                            ax_e=ax
                        else:
                            ax=ax_e
                    went_to_elec=True
                i=-1
                for c,j in zip(cout,electrode_flux):
                    i+=1
                    if i!=len(cout)-1:
                        lw=0.5
                        zorder=100
                    else:
                        lw=2.0
                        zorder=0
                    brightness=1.-(color_offset+(i*1.)/len(cout)*(1.-color_offset))
                    if k==0:
                        color=str(brightness)
                    elif k==1:
                        color=(brightness,1.,1.)

                    if sp in ['OH-','H+'] and 'pH' in sum([self.large_plots,self.small_plots],[]):
                        if i==len(cout)-1:
                            ax_pH.plot(self.tp.xmesh,14.+np.log10(c[k*self.tp.nx:(k+1)*self.tp.nx] /10**3),'-',color=color,linewidth=lw,zorder=zorder,label=r'$'+self.tp.species[sp]['symbol']+'$')
                                    #'z='+str(int(self.tp.charges[k]/unit_F)))
                        else:
                            ax_pH.plot(self.tp.xmesh,14.+np.log10(c[k*self.tp.nx:(k+1)*self.tp.nx] /10**3),'-',color=color,linewidth=lw,zorder=zorder)
                    if 'electrode_flux' in sum([self.large_plots,self.small_plots],[]):
                        if i==len(cout)-1:
                            ax_elflux.plot(self.tp.xmesh,j[k*self.tp.nx:(k+1)*self.tp.nx],'-',color=color,linewidth=lw,zorder=zorder,label=r'$'+self.tp.species[sp]['symbol']+'$')
                        else:
                            ax_elflux.plot(self.tp.xmesh,j[k*self.tp.nx:(k+1)*self.tp.nx],'-',color=color,linewidth=lw,zorder=zorder,label=r'$'+self.tp.species[sp]['symbol']+'$')
                    if any([a.startswith('concentrations') for a in sum([self.large_plots,self.small_plots],[])]):
                        if self.logscale:
                            func=ax.semilogy
                        else:
                            func=ax.plot
                        if i==len(cout)-1:
                            func(self.tp.xmesh,c[k*self.tp.nx:(k+1)*self.tp.nx] /10**3,'-',color=color,linewidth=lw,zorder=zorder,label=r'$'+self.tp.species[sp]['symbol']+'$')
                                    #'z='+str(int(self.tp.charges[k]/unit_F)))
                        else:
                            func(self.tp.xmesh,c[k*self.tp.nx:(k+1)*self.tp.nx] /10**3,'-',color=color,linewidth=lw,zorder=zorder)
                if any([a.startswith('concentrations') for a in sum([self.large_plots,self.small_plots],[])]):
                    ax.legend(ncol=2)
                    ax.set_xlabel('x (m)')
                    ax.set_ylabel('c (mol/L)')
                if 'electrode_flux' in sum([self.large_plots,self.small_plots],[]):
                    ax_elflux.legend(ncol=2)

        #c=self.integrate_FTCS(self.dt,self.dx)
        #for t in np.arange(0.0,1.,0.1):
        #    plt.plot(self.xmesh,c[int(t/self.dt),:],'-o',label=str(t))
        if 'efield' in sum([self.large_plots,self.small_plots],[]):
            ax=get_current_ax('efield')
            ax.plot(self.tp.xmesh,self.tp.efield/1e10,'-')
            ax.set_title('Electric field')
            ax.set_xlabel('x (m)')
            ax.set_ylabel('E (V/Ang)')
            ax.legend()
#        if 'vzeta' in self.tp.system:
#            ax2.plot(self.tp.xmesh,[-self.tp.gouy_chapman(x)[1]/1e10 for x in self.tp.xmesh],'-',color='r',linewidth=lw,label='Gouy-Chapman')

        if 'potential' in sum([self.large_plots,self.small_plots],[]):
            ax=get_current_ax('potential')
            ax.set_title('Potential')
            ax.plot(self.tp.xmesh,self.tp.potential,'-')
            #if 'vzeta' in self.tp.system:
            #    ax4.plot(self.tp.xmesh,[self.tp.gouy_chapman(x)[0] for x in self.tp.xmesh],'-',color='r',linewidth=lw,label='Gouy-Chapman')
            ax.set_ylabel('v (V)')
            ax.set_xlabel('x (m)')
            ax.legend()
#        ax4.plot(self.tp.xmesh[:-1],self.tp.external_charge[:-1]/10**3/unit_F, '-',label='n_ext')

        if 'total charge' in sum([self.large_plots,self.small_plots],[]):
            ax=get_current_ax('total charge')
            ax.plot(self.tp.xmesh,self.tp.total_charge/10**3/unit_F,'-',label='n_tot')
            ax.set_ylabel('n_ion (e*mol/L)')
            ax.legend()
            ax.set_title('Charge Density')

        if 'current_density' in sum([self.large_plots,self.small_plots],[]):
            ax=get_current_ax('current_density')
            ax.plot(self.tp.xmesh,self.tp.current_density/10,'-',label='J')
            ax.set_ylabel('J (mA/cm^2)')
            ax.legend()
            ax.set_title('Current Density')

        for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

        gs1.tight_layout(fig, rect=[0.0, 0, 1, 1], h_pad=1.0)

