import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle
from units import *
import sys
from matplotlib import colors as mcolors

class Plot():

    def __init__(self,transport=None):
        if transport==None:
            print('No transport object provided for calculator. Stopping here.')
            sys.exit()
        else:
            self.tp=transport #transport object

        self.ax1=plt.subplot('311')
        self.ax2=plt.subplot('312')
        self.ax3=plt.subplot('313')

    def plot(self,cout):

#        colorlist=cycle(['b','k','])
        colors=dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
        colorlist=cycle([colors[key] for key in colors])

#        ax = plt.subplot(2, 2, 1)

        ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=2, rowspan=2,)
        ax2 = plt.subplot2grid((2, 4), (0, 2)) 
        ax3 = plt.subplot2grid((2, 4), (0, 3)) 
        ax4 = plt.subplot2grid((2, 4), (1, 2)) 
        ax5 = plt.subplot2grid((2, 4), (1, 3)) 

#        ax1=plt.subplot('221')
#        ax4=plt.subplot('222')
#        ax3=plt.subplot('223')
#        ax2=plt.subplot('224')
    
        color_offset=0.3
        for k,sp in enumerate(self.tp.species):
            color=next(colorlist)
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
                    ax1.plot(self.tp.xmesh,c[k*self.tp.nx:(k+1)*self.tp.nx] /10**3,'-',color=color,linewidth=lw,zorder=zorder,label=self.tp.species[sp]['name'])
                            #'z='+str(int(self.tp.charges[k]/unit_F)))
                else:
                    ax1.plot(self.tp.xmesh,c[k*self.tp.nx:(k+1)*self.tp.nx] /10**3,'-',color=color,linewidth=lw,zorder=zorder)
        ax1.legend(ncol=2)
        ax1.set_xlabel('x (m)')
        ax1.set_ylabel('c (mol/L)')
        ax1.set_title('Concentrations')
        #c=self.integrate_FTCS(self.dt,self.dx)
        #for t in np.arange(0.0,1.,0.1):
        #    plt.plot(self.xmesh,c[int(t/self.dt),:],'-o',label=str(t))
        ax2.plot(self.tp.xmesh,self.tp.efield/1e10,'-')
        if 'vzeta' in self.tp.system:
            ax2.plot(self.tp.xmesh,[-self.tp.gouy_chapman(x)[1]/1e10 for x in self.tp.xmesh],'-',color='r',linewidth=lw,label='Gouy-Chapman')
        ax2.set_title('Electric field')
        ax2.set_xlabel('x (m)')
        ax2.set_ylabel('E (V/Ang)')
        ax2.legend()

        ax3.set_title('Potential')
        #self.potential=np.zeros([len(self.xmesh)])
        #integral=0.0
        #for i in range(len(self.xmesh)):
        #    integral-=self.efield[i]*self.dx
        #    self.potential[i]+=integral
        ax3.plot(self.tp.xmesh,self.tp.potential,'-')
        if 'vzeta' in self.tp.system:
            ax3.plot(self.tp.xmesh,[self.tp.gouy_chapman(x)[0] for x in self.tp.xmesh],'-',color='r',linewidth=lw,label='Gouy-Chapman')
        ax3.set_ylabel('v (V)')
        ax3.set_xlabel('x (m)')
        ax3.legend()
#        ax4.plot(self.tp.xmesh[:-1],self.tp.external_charge[:-1]/10**3/unit_F, '-',label='n_ext')
        ax4.plot(self.tp.xmesh[:-1],self.tp.total_charge[:-1]/10**3/unit_F,'-',label='n_tot')
        ax4.set_ylabel('n_ion (e*mol/L)')
        ax4.legend()
        ax4.set_title('Charge Density')
        #plt.tight_layout()
        plt.show()
    #    self.integrate_pb()
