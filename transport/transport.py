"""
transport class. defines the transport model.
---
Stefan Ringe (stefan.ringe.tum@gmail.com)
"""
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
#from ase import units
from units import *
from itertools import cycle
import sys
from copy import deepcopy
import collections

class Transport(object):

    def __init__(self, species=None,reactions=None,system=None):

        #go over input data and put in some defaults if none
        if not type(species)==dict:
            self.species={'species1':       {'name':r'K^+',
                                            'diffusion':1.96e-9,
                                            'bulk concentration':0.001*1000.},
                          'species2':       {'name':r'HCO_3^-',
                                            'diffusion':1.2e-9,
                                            'bulk concentration':0.001*1000.}}
        else:
            self.species=species
        self.nspecies=len(self.species)

        self.reactions=reactions
        if not type(system)==dict:
            self.system={'temperature': 298,
                         'epsilon':     80,
                         'vzeta':    -0.025}
        else:
            self.system=system

        self.eps = self.system['epsilon']*unit_eps0 #1.1e11 #*unit_eps0
        self.beta = 1./(self.system['temperature'] * unit_R)

        #determine charges and create arrays of charges and D's
        k=-1
        self.charges=[]
        self.D=[]
        for sp in self.species:
            k+=1
            string=self.species[sp]['name'].split('^')
            if len(string)>1:
                print string
                sign=string[-1][-1]
                if len(string[-1])>1:
                    no=int(string[-1][:-1])
                else:
                    no=1
                if sign=='-':
                    no*=(-1)
            else:
                no=0
            self.species[sp]['charge']=no
            #create a few shorter arrays
            self.charges.append(no*unit_F)
            self.D.append(self.species[sp]['diffusion'])
        self.charges=np.array(self.charges)
        self.D=np.array(self.D)
        self.mu = self.D * self.charges *self.beta  #ion mobilities according to Einstein relation

        #Debye-Hueckel screening length
        self.ionic_strength=0.0
        for isp,sp in enumerate(self.species):
            self.ionic_strength+=self.charges[isp]**2*self.species[sp]['bulk concentration']
        self.ionic_strength*=0.5
        self.debye_length = np.sqrt( self.eps/self.beta/2./self.ionic_strength ) #in m


        #THE MESH
        self.dx=self.debye_length/10.
        #self.xmax=80.e-6 #*1e-10
        self.xmax=10*self.debye_length
        self.xmesh=np.arange(0,self.xmax+self.dx,self.dx)
        self.nx=len(self.xmesh)

        #A FEW VARIABLES
        self.external_charge=np.zeros([len(self.xmesh)])
        #self.external_charge=self.gaussian(sigma=5e-6,z=1.0,mu=4.e-5,cmax=5e-5)+self.gaussian(sigma=5e-6,z=-1.0,mu=5.e-5,cmax=5e-5)
        self.count=1

        #BOUNDARY AND INITIAL CONDITIONS

        c_initial_general={}
        for isp,sp in enumerate(self.species):
            c_initial_general[str(isp)]=self.species[sp]['bulk concentration']

        self.set_initial_conditions(
                c_initial_general=c_initial_general)
                #c_initial_specific={'0':{'0':0.1},'1':{'0':0.1}})

        self.set_boundary_conditions(\
                flux_boundary={'0':{'l':0.0},'1':{'l':0.0}},         #in mol/s/cm^2
                dc_dt_boundary={'all':{'l':0.0}},     #in mol/l/s #give either left OR right boundary condition here
                #integration will start at the site where dc_dt is defined
                efield_boundary={'l':0.0})    #in V/Ang

    def set_initial_concentrations(self,func):
        """we can set the initial concentrations to a particular function"""
        if func=='Gouy-Chapman':
            if self.nspecies!=2:
                print('Gouy-Chapman limit only implemented for two species, cationic' 
                        'and anionic. Not applying initialization.')
                return
            function=self.gouy_chapman
        def tree():
            return collections.defaultdict(tree)

        c_initial_specific = tree() #collections.defaultdict(list)

        for k,sp in enumerate(self.species):
            for i in range(self.nx):
                c_initial_specific[str(k)][str(i)]=\
                    self.species[sp]['bulk concentration']*\
                    np.exp(-self.beta*self.charges[k]*function(self.xmesh[i]))
        self.set_initial_conditions(c_initial_specific=c_initial_specific)
        return

    def get_static_concentrations(self):
        """solves the PBE in order to get the static limit for the concentrations"""
        def d(z, x):
            zold=z
            z=[]
            #z[0] = phi'
            #z[1] = phi
            #first block defines z[0]' which is d^2/dx^2 phi as function of
            #second block defines z[1]' which is phi, so z[0]
            z.append(-1./self.eps*self.external_charge_func(x)) #(sum([charge*unit_F*0.1*10**3*np.exp(-self.beta*unit_F*charge*zold[1]) for charge in self.charges])+self.external_charge_func(x)))
            z.append(zold[0])
            return z

        v0=np.zeros([2]) #this is phi' and phi at x=0
        sol,output = integrate.odeint(d, v0, self.xmesh,full_output=True) #, args=(b, c))
        print np.shape(sol)
        print sol[:,0]
        #b_min = minimize(func, b0, args=(z0, m, k, g, a), constraints=cons)
        self.ax1.plot(self.xmesh,sol[:,0]/1e10,'-',label='E (V/Ang)')
        self.ax2.plot(self.xmesh,sol[:,1],'-',label='phi (V)')
        self.ax1.legend()
        self.ax2.legend()
        plt.show()
        sys.exit()
        return sol

       
    def gouy_chapman(self,x):
        term1 = 1.+np.tanh(self.system['vzeta']*self.beta*unit_F/4.)*\
            np.exp(-1./self.debye_length*x)
        term2 = 1.-np.tanh(self.system['vzeta']*self.beta*unit_F/4.)*\
            np.exp(-1./self.debye_length*x)
        return 2./(self.beta*unit_F)*np.log(term1/term2)


    def get_initial_conditions(self):
        return self.c0 

    def set_initial_conditions(self,\
        c_initial_general={},\
        c_initial_specific={}):
        """accepts the initial concentrations of all species in the c_general list (all in mol/L), 
        the c_specific dictionary allows to parse specific values of the concentrations at 
        specific grid points (1st key is species index, 2nd is xmesh index, 
        value is concentration)."""

        c0=np.zeros([self.nspecies*len(self.xmesh)])
        j=-1

        for k in range(self.nspecies):
            for i in range(self.nx):
                j+=1
                if str(k) in c_initial_general:
                    c0[j]=c_initial_general[str(k)]
                elif 'all' in c_initial_general:
                    c0[j]=c_initial_general['all']
                if (str(k) in c_initial_specific and str(i) in c_initial_specific[str(k)]):
                    c0[j]=c_initial_specific[str(k)][str(i)]
                if ('all' in c_initial_specific and str(i) in c_initial_specific['all']):
                    c0[j]=c_initial_specific['all'][str(i)]

        c0=np.array(c0)
        self.c0=c0


    def get_boundary_conditions(self):
        return self.c_bound,self.dc_dt_bound,self.efield_bound

    def set_boundary_conditions(self,\
        flux_boundary={},\
        dc_dt_boundary={},\
        efield_boundary={}):

        """allows to set boundary conditions in the concentrations themselves (c_boundary) 
        and the derivatives dc_dx_boundary and the efield. 
        species, direction, value"""

        if len(flux_boundary)>0:
            self.boundary_type='flux'
        elif len(dc_dt_boundary)>0:
            self.boundary_type='dc_dt'
        else:
            print('No boundary conditions defined, stopping here.')
            sys.exit()

        flux_bound=np.zeros([self.nspecies,2])
        dc_dt_bound=np.zeros([self.nspecies,2])
        efield_bound=np.array([None,None])

        def key_to_index(key):
            if key=='l':
                index=0
            elif key=='r':
                index=1
            return index

        def species_to_index(key):
            if key.isdigit():
                index=int(key)
            else:
                index=self.nspecies+1

        for key in efield_boundary:
            efield_bound[key_to_index(key)]=efield_boundary[key]

        def allocate_fluxes(array_in,array_out):
            for key1_raw in array_in:
                if key1_raw=='all':
                    keys1=map(str,range(self.nspecies))
                else:
                    keys1=[str(key1_raw)]
                for key1 in keys1:
                    for key2_raw in array_in[key1_raw]:
                        if key2_raw=='all':
                            keys2=map(str,range(2))
                        else:
                            keys2=[key_to_index(key2_raw)]
                        for key2 in keys2:
                            array_out[int(key1),int(key2)]=array_in[key1_raw][key2_raw]
            return array_out

        dc_dt_bound=allocate_fluxes(dc_dt_boundary,dc_dt_bound)
        flux_bound=allocate_fluxes(flux_boundary,flux_bound)

        for l in [dc_dt_bound,efield_bound,flux_bound]:
            l=np.array(l)

        self.flux_bound=flux_bound*100**2
        self.dc_dt_bound=dc_dt_bound*10**3
        self.efield_bound=[]
        for e in efield_bound:
            if e==None:
                self.efield_bound+=[None]
            else:
                self.efield_bound+=[e*1e10]
        self.efield_bound=np.array(self.efield_bound)

    def gaussian(self,sigma=0.01,z=1,mu=0.2,cmax=1):
        """define gaussian charge density. c is in molar"""
#        sigma=1./(np.sqrt(2.*np.pi))*1/cmax
        gaussian=np.zeros([len(self.xmesh)])
        for i in range(0,len(self.xmesh)):
            gaussian[i]=\
                    1./(sigma*np.sqrt(2.*np.pi))\
                    *np.exp(-0.5*((self.xmesh[i]-mu)/sigma)**2)
        gaussian=gaussian/max(gaussian)*cmax*10**3*z*unit_F
        return gaussian

    def external_charge_func(self,x):
        return self.gaussian_func(x,sigma=5e-6,z=1.0,mu=4.e-5,cmax=5e-5)+self.gaussian_func(x,sigma=5e-6,z=-1.0,mu=5.e-5,cmax=5e-5)

    def gaussian_func(self,x,sigma=0.01,z=1,mu=0.2,cmax=1):
        """define gaussian charge density. c is in molar"""
#        sigma=1./(np.sqrt(2.*np.pi))*1/cmax
        gaussian=\
             1./(sigma*np.sqrt(2.*np.pi))\
             *np.exp(-0.5*((x-mu)/sigma)**2)
        gaussian=gaussian*cmax*10**3*z*unit_F
        return gaussian

    def set_calculator(self,calc=None):
        """Attach calculator object."""
        self.calc=calc

#    def update():
#
#        self.rates = []
#        self.rates += self.migration_rates
#        self.rates += self.diffusion_rates
#        self.rates += self.convection_rates
#        self.rates += self.reaction_rates

    #def get_rates(self):
    #        for i in self.xmesh:
    #            rates[i] = self.D*self.z/unit.kB/self.T*concentraitions*
    #    return rates



#    def read_catmap():
#        """Purpose: read the catmap input file"""
         
