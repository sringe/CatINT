import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle
from units import *
import sys
from matplotlib import colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick


class Plot():

    def __init__(self,transport=None):
        if transport==None:
            print('No transport object provided for calculator. Stopping here.')
            sys.exit()
        else:
            self.tp=transport #transport object

    def plot(self,cout=None):

        cout=self.tp.cout

        colors=dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
        colorlist=cycle([colors[key] for key in colors])

        fig = plt.figure(figsize=(14,9)) #gcf()

        flatten = lambda l: [item for sublist in l for item in sublist]

        gs1 = gridspec.GridSpec(3,4)
        ax1 = fig.add_subplot(gs1[0:2,0:2])
        ax2 = fig.add_subplot(gs1[0:2,2:4])

        product_list=['HCOO-','CO','H2','etol','propol','metol','C2H4','CH4','allyl']

        color_offset=0.3
        for k,sp in enumerate(self.tp.species):
            color=next(colorlist)
            if sp in product_list:
                ax=ax1
            else:
                ax=ax2
            i=-1
            for c in cout:
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
                if i==len(cout)-1:
                    ax.plot(self.tp.xmesh,c[k*self.tp.nx:(k+1)*self.tp.nx] /10**3,'-',color=color,linewidth=lw,zorder=zorder,label=r'$'+self.tp.species[sp]['symbol']+'$')
                            #'z='+str(int(self.tp.charges[k]/unit_F)))
                else:
                    ax.plot(self.tp.xmesh,c[k*self.tp.nx:(k+1)*self.tp.nx] /10**3,'-',color=color,linewidth=lw,zorder=zorder)
        for ax in [ax1,ax2]:
            ax.legend(ncol=2)
            ax.set_xlabel('x (m)')
            ax.set_ylabel('c (mol/L)')
        ax1.set_title('Product Concentrations')
        ax2.set_title('Electrolyte Concentrations')

        ax3 = fig.add_subplot(gs1[8])
        ax4 = fig.add_subplot(gs1[9])
        ax5 = fig.add_subplot(gs1[10])
        ax6 = fig.add_subplot(gs1[11])
        #c=self.integrate_FTCS(self.dt,self.dx)
        #for t in np.arange(0.0,1.,0.1):
        #    plt.plot(self.xmesh,c[int(t/self.dt),:],'-o',label=str(t))
        ax3.plot(self.tp.xmesh,self.tp.efield/1e10,'-')
        if 'vzeta' in self.tp.system:
            ax2.plot(self.tp.xmesh,[-self.tp.gouy_chapman(x)[1]/1e10 for x in self.tp.xmesh],'-',color='r',linewidth=lw,label='Gouy-Chapman')
        ax3.set_title('Electric field')
        ax3.set_xlabel('x (m)')
        ax3.set_ylabel('E (V/Ang)')
        ax3.legend()

        ax4.set_title('Potential')
        #self.potential=np.zeros([len(self.xmesh)])
        #integral=0.0
        #for i in range(len(self.xmesh)):
        #    integral-=self.efield[i]*self.dx
        #    self.potential[i]+=integral
        ax4.plot(self.tp.xmesh,self.tp.potential,'-')
        if 'vzeta' in self.tp.system:
            ax4.plot(self.tp.xmesh,[self.tp.gouy_chapman(x)[0] for x in self.tp.xmesh],'-',color='r',linewidth=lw,label='Gouy-Chapman')
        ax4.set_ylabel('v (V)')
        ax4.set_xlabel('x (m)')
        ax4.legend()
#        ax4.plot(self.tp.xmesh[:-1],self.tp.external_charge[:-1]/10**3/unit_F, '-',label='n_ext')
        ax5.plot(self.tp.xmesh,self.tp.total_charge/10**3/unit_F,'-',label='n_tot')
        ax5.set_ylabel('n_ion (e*mol/L)')
        ax5.legend()
        ax5.set_title('Charge Density')

        ax6.plot(self.tp.xmesh,self.tp.current_density/100**3,'-',label='J')
        ax6.set_ylabel('J (A/cm^2)')
        ax6.legend()
        ax6.set_title('Current Density')

        for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

        gs1.tight_layout(fig, rect=[0.0, 0, 1, 1], h_pad=1.0)

        plt.show()
