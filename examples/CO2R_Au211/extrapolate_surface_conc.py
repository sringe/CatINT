import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from scipy.optimize import curve_fit
from catint.comsol_reader import Reader
from glob import glob
import re
from copy import deepcopy

class extrapolate():
    def __init__(self,tp=None,extrapol_folder=None):
        self.extrapol_data={}
        self.extrapol_data['voltage_diff_drop']=[]
        for topdir in extrapol_folder:
            for folder in glob(topdir+'/comsol_results*'):
                #get the potential:
                for line in open(folder+'/pnp_transport.java','r'):
                    if 'phiM' in line and 'Metal Potential' in line:
                        print 'found!!!',line
                        sol=re.findall('phiM\",\"(-?[0-9.]{2,30}).*Metal Potential',line)
                        if len(sol)>0:
                            phiM=float(sol[0])
                        break
                if not os.path.exists(folder+'/concentrations.txt'):
                    continue
                comsol_reader=Reader(transport=tp,results_folder=folder,\
                    outputs=['concentrations','electrostatics','electrode_flux','rho_charge'],comsol_args=tp.comsol_args)
                comsol_reader.read_all()
                for sp in tp.species:
                    if sp not in self.extrapol_data:
                        self.extrapol_data[sp]=[]
                    self.extrapol_data[sp].append([phiM,tp.species[sp]['surface_concentration']])
                self.extrapol_data['voltage_diff_drop'].append([phiM,tp.system['potential'][0]])

        for key in self.extrapol_data:
            self.extrapol_data[key]=np.array(sorted(self.extrapol_data[key]))

        self.extrapol_func_list={}
        def exponential(x,*p):
            A,b,c,d=p #,d,e,f=p
            return A*np.tanh(b*(x-c))+d #*np.log(b*(x+c))+d+e*(x-f)
        def exponential_poly2(x,*p):
            A,b,c,d,e,f=p #,d,e,f=p
            return A*np.tanh(b*(x-c))+d+e*x+f*x**2 #*np.log(b*(x+c))+d+e*(x-f)
        def exponential_v2(x,p):
            A,b,c,d=p #,d,e,f=p
            return A*np.tanh(b*(x-c))+d #*np.log(b*(x+c))+d+e*(x-f)
        def exponential_poly2_v2(x,p):
            A,b,c,d,e,f=p #,d,e,f=p
            return A*np.tanh(b*(x-c))+d+e*x+f*x**2 #*np.log(b*(x+c))+d+e*(x-f)
        
        fig=plt.figure(figsize=(10,6))
        ax1=fig.add_subplot('131')
        ax2=fig.add_subplot('132')
        ax3=fig.add_subplot('133')
        
        colors=cycle(['C'+str(i) for i in range(10)])
        for ads in tp.species: #data:
            color=next(colors)
            x=self.extrapol_data[ads][:,0]
            y=self.extrapol_data[ads][:,1]
            y=np.log10(y)
            ax1.plot(x,y,'x',color=color)
            if ads == 'K+':
                p0=[0.1,-3.,-0.5,0.4] #,0.6,0.5,0.0,-0.2]
            elif ads == 'CO':
                p0=[3.,-3.,-0.5,-7.5]
            elif ads == 'OH-':
                p0=[3.,-3.,-0.5,-7.5]
            if ads in ['K+','CO','OH-']: #,'CO']:
                coeff, var_matrix = curve_fit(exponential, x,y, p0=p0)
                #coeff=p0
                if ads in ['CO','OH-']:
                    p0=list(coeff)+[0.1,0.5]
                    coeff, var_matrix = curve_fit(exponential_poly2, x,y, p0=p0)
                    self.extrapol_func_list[ads]=[coeff,exponential_poly2_v2]
                else:
                    self.extrapol_func_list[ads]=[coeff,exponential_v2]
            else:
                #if ads in ['OH-']:
                #    d=np.array([[xx,yy] for xx,yy in zip(x,y) if xx<-0.7])
                #    x=d[:,0]
                #    y=d[:,1]
                #    order=2
                #else:
                order=2
                z=np.polyfit(x,y,order)
                flambda=lambda x,c: np.poly1d(c)(x)
                self.extrapol_func_list[ads]=[z,flambda]
        #    coeff=p0
            xf=np.linspace(-5,5,1000)# -2.0,-0.4,1000)
            if ads=='OH-':
                ax3.plot(x,y-3.+14.,'x',color=color)
        #replace functions of coefficients and x with simple function of x
        colors=cycle(['C'+str(i) for i in range(10)])
        def make_closure(sp):
            return lambda x: self.extrapol_func_list[sp][1](x, self.extrapol_func_list[sp][0])
        self.extrapol_func={}
        for sp in tp.species:
            color=next(colors)
            self.extrapol_func[sp]=make_closure(sp)
            ax1.plot(xf,self.extrapol_func[sp](xf),'-',color=color,label=sp)
            if sp=='OH-':
                self.extrapol_func['surface_pH']=lambda x: self.extrapol_func['OH-'](x)-3.+14.
                ax3.plot(xf,self.extrapol_func['surface_pH'](xf),color=color)
        x=self.extrapol_data['voltage_diff_drop'][:,0]
        y=self.extrapol_data['voltage_diff_drop'][:,1]
        ax2.plot(x,y,'o')
        z=np.polyfit(x,y,2)
        p=np.poly1d(z)
        self.extrapol_func['voltage_diff_drop']=p
        xf=np.linspace(-2.0,0.0,1000)
        ax2.plot(xf,p(xf),'-')
        ax1.legend()
        ax1.set_ylabel('Surface Concentrations (mol/dm^3)')
        ax2.set_xlabel('phiM (V vs. SHE)')
        ax1.set_xlabel('phiM (V vs. SHE)')
        ax2.set_ylabel('Diffusional Potential Drop (V)')
        ax3.set_ylabel('Surface pH')
        ax3.set_xlabel('phiM (V vs. SHE)')
        #ax1.set_xlim([-1.4,-0.4])
        ax1.set_ylim([-17,6.])
        ax1.set_xlim([-2.0,1.0])
        plt.tight_layout()
    def plot(self):
        plt.show()
        sys.exit()
