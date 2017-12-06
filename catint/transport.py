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
import logging
import os
class Transport(object):

    def __init__(self, species=None,reactions=None,system=None,pb_bound=None,nx=100,scf_bound=False,comsol_params={}):

        # set up logging to file - see previous section for more details
        logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='transport.log',
                    filemode='w')
        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        # set a format which is simpler for console use
        formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
        # tell the handler to use this format
        console.setFormatter(formatter)
        # add the handler to the root logger
        logging.getLogger('').addHandler(console)
        
        # Now, we can log to the root logger, or any other logger. First the root...
        logging.info('Starting Transport Calculation.')
        
        # Now, define a couple of other loggers which might represent areas in your
        # application:
        
        self.logger_db = logging.getLogger('transport.debug')
        self.logger = logging.getLogger('transport.info')
        
        #all the possible keys:
        species_keys=['bulk concentration', 'diffusion', 'name', 'symbol', 'zeff','flux','req']
        system_keys=['vzeta','temperature','pressure','water viscosity','electrolyte viscosity',\
                'epsilon','exclude species','migration','boundary thickness','educts',\
                'products']

        #go over input data and put in some defaults if none
        if species is None:
            self.species={'species1':       {'symbol':r'K^+',
                                            'name':'potassium',
                                            'diffusion':1.96e-9,
                                            'bulk concentration':0.001*1000.},
                          'species2':       {'symbol':r'HCO_3^-',
                                            'name':'bicarbonate',
                                            'diffusion':1.2e-9,
                                            'bulk concentration':0.001*1000.}}
        else:
            for sp in species:
                for key in species[sp]:
                    if key not in species_keys:
                        self.logger.error('No such key "'+key+'" in species list. Quitting here.')
                        sys.exit()
            self.species=species

        if pb_bound is None:
            pb_bound={
            'potential': {'wall':'zeta'},
            'gradient': {'bulk':0.0}}

        #self-consistent boundary condition
        self.scf_bound=scf_bound

        system_defaults={
                'epsilon':78.36,
                'Stern epsilon': 2.0, #dielectric permittivity in Stenr layer (for Robin BCs)
                'Stern capacitance':18e-2, #Stern layer capacitance for Robin BCs
                'temperature':298.14,
                'vzeta':0.0, #-0.0125,
                'pressure':1,
                'educts':[], #'CO2',\
                'products':[]}
        if system is None:
            self.system=system_defaults
        else:
            for key in system:
                if key not in system_keys:
                    self.logger.error('No such key "'+key+'" in system list. Quitting here.')
                    sys.exit()
            self.system=system

        for key in system_defaults:
           # ['epsilon','temperature','pressure','vzeta']:
            if key not in self.system:
                self.system[key]=system_defaults[key]

        #delete species that should not be considered here
        if 'exclude species' in self.system:
            for es in self.system['exclude species']:
                del self.species[es]

        self.nspecies=len(self.species)

        #defaults for system settings:
        self.reactions=reactions

        #add concentration = 0 for all species with no separately defined values
        for sp in self.species:
            if 'bulk concentration' not in self.species[sp]:
                self.species[sp]['bulk concentration']=0.0

            self.species[sp]['updated concentration']=self.species[sp]['bulk concentration']

        self.eps = self.system['epsilon']*unit_eps0 #1.1e11 #*unit_eps0
        self.beta = 1./(self.system['temperature'] * unit_R)

        #GET CHARGES AND NCATOMS FROM CHEMICAL SYMBOLS
        self.charges=self.symbol_reader(self.species)

        self.use_migration=True
        if 'migration' in self.system:
            if not self.system['migration']:
                self.use_migration=False

        self.use_convection=False

        self.use_reactions=False
        if self.reactions is not None:
            if any(['rates' in self.reactions[reaction] for reaction in self.reactions]):
                self.use_reactions=True
                self.logger.info('Found reactions.')
#        self.use_reactions=False
        #DIFFUSION CONSTANTS
        self.D=[]
        for sp in self.species:
            if 'diffusion' in self.species[sp]:
                self.D.append(self.species[sp]['diffusion'])
            else:
                self.D.append(0.0)
        self.D=np.array(self.D)

        if all([a in self.system for a in ['water viscosity','electrolyte viscosity']]):
            #rescale diffusion coefficients according to ionic strength (Stokes-Einstein):
            self.logger.info('Rescaling diffusion coefficients from water viscosity {} to electrolyte viscosity {}'.format(self.system['water viscosity'], self.system['electrolyte viscosity']))
            self.D=np.array([d*float(self.system['water viscosity'])/float(self.system['electrolyte viscosity']) for d in self.D])
            

        self.mu = self.D * self.charges *self.beta  #ion mobilities according to Einstein relation

        #Debye-Hueckel screening length
        self.ionic_strength=0.0
        for isp,sp in enumerate(self.species):
            self.ionic_strength+=self.charges[isp]**2*self.species[sp]['bulk concentration']
        self.ionic_strength*=0.5
        self.debye_length = np.sqrt( self.eps/self.beta/2./self.ionic_strength ) #in m

        if 'boundary thickness' in self.system:
            self.boundary_thickness=self.system['boundary thickness']

        #THE MESH
        self.nx=nx
        if 'boundary thickness' in self.system:
            self.xmax=self.boundary_thickness
            self.dx=self.xmax/(self.nx*1.)
        else:
            #use debye-hueckel length to get a reasonable estimate
            #WARNING: this can be way to small, depending on the system of interest
            nx_mod=max(1.,np.ceil(self.nx/10.))
            self.xmax=self.debye_length*nx_mod #*50. #*nx_mod
            self.dx=self.debye_length/nx_mod
        self.xmesh=np.arange(0,self.xmax+self.dx,self.dx)
        self.nx=len(self.xmesh)

#        if min(self.xmesh)<=0:
#            min_x=1e-15
#        else:
#            min_x=min(self.xmesh)
#        self.xmesh=np.logspace(np.log10(min_x),np.log10(self.xmax-self.dx),self.nx)
#        self.nx=len(self.xmesh)
#        self.nmax=max(self.xmesh)

        #A FEW VARIABLES
        self.external_charge=np.zeros([len(self.xmesh)])
        #self.external_charge=self.gaussian(sigma=5e-6,z=1.0,mu=4.e-5,cmax=5e-5)+self.gaussian(sigma=5e-6,z=-1.0,mu=5.e-5,cmax=5e-5)
        self.count=1

        #STEADY-STATE REACTION RATES/FLUXES
        for sp in self.species:
            if 'flux' not in self.species[sp]:
                self.species[sp]['flux']=0.0
            else:
                if type(self.species[sp]['flux'])==str:
                    self.logger.info('Flux of species '+sp+' is assumed to be an equation.')
                elif type(self.species[sp]['flux'])==dict:
                    self.logger.info('Flux of species '+sp+' will be evaluated from the fluxes of '+str(self.species[sp]['flux']['values']))
        #we need to calculate the flux of the educt as the sum of all the product rates
        #only do this, if the rate is not given as function
        if not any([type(self.species[sp]['flux'])==str for sp in self.species]):
            #functions work only for comsol and are implemented in comsol calculator
            self.evaluate_fluxes()


        #FLUX AND FARADAIC YIELD
        flux_bound={} 
        for isp,sp in enumerate(self.species):
            flux_bound[str(isp)]={}
        for isp,sp in enumerate(self.species):
            if 'flux' in self.species[sp]:
                pass
                #self.species[sp]['flux']*=-1
            else:
                self.species[sp]['flux']=0.0
            flux_bound[str(isp)]['l']=self.species[sp]['flux']

        #debugging fluxes
        self.logger_db.debug("DB: FLUXES")
        for isp, sp in enumerate(self.species): #['CO2','HCO3-','CO32-','OH-','H2','CO','CH4','C2H4','HCOO-','etol','propol','allyl','metol','acet','etgly','K']:
  #          for isp2,sp2 in enumerate(self.species):
  #              if sp==sp2:
  #                  break
            self.logger_db.debug('{} {}'.format(sp, flux_bound[str(isp)]['l']))

        #BOUNDARY AND INITIAL CONDITIONS
        self.set_boundary_and_initial_conditions(pb_bound,flux_bound)


        #INITIALIZE FIELD AND POTENTIAL
        self.efield=np.zeros([self.nx])
        self.potential=np.zeros([self.nx])
        self.total_charge=np.zeros([self.nx])

        #define additional parameters which should be used inside comsol
        self.comsol_params=comsol_params

    def evaluate_fluxes(self):
        self.logger.info('Evaluating fluxes of {} as sum of products/educts'.format([sp for sp in self.species if type(self.species[sp]['flux'])==dict]))
        tmp_rates=np.zeros([self.nspecies])
        for isp,sp in enumerate(self.species):
            if 'zeff' in self.species[sp]:
                if type(self.species[sp]['flux'])!=dict:
                    crate=self.species[sp]['flux']/(unit_F*self.species[sp]['zeff'])
                    self.species[sp]['flux']=crate
#                if sp not in sum([self.system['sum rates'][sp]['values'] for sp in self.system['sum rates']],[]):
#                    continue
                for isp2,sp2 in enumerate(self.species):
                    if type(self.species[sp2]['flux'])==dict:
                        if sp in self.species[sp2]['flux']['values']:
                            if self.species[sp2]['flux']['reqs']=='number_of_catoms':
                                req=self.number_of_catoms[isp]
                            elif self.species[sp2]['flux']['reqs']=='zeff':
                                req=self.species[sp]['zeff']
                            else:
                                req=self.species[sp2]['flux']['reqs'][isp]
                            if self.species[sp2]['flux']['kind']=='educt':
                                tmp_rates[isp2]-=crate*req
                            else:
                                tmp_rates[isp2]+=crate*req
#            elif sp=='HCO3-':
#                self.species[sp]['flux']=-self.species['H2']['flux']/(unit_F*self.species['H2']['zeff'])*2.
        for isp,sp in enumerate(self.species):
            if type(self.species[sp]['flux'])==dict:
                self.species[sp]['flux']=tmp_rates[isp]
#        for ieduct,educt in enumerate(self.system['educts']):
#            self.species[educt]['flux']=-tmp_educts[ieduct]
        if 'unknown' in self.species:
            #we needed this flux only to calculate the co2 and oh- fluxes
            self.species['unknown']['flux']=0.0

    def symbol_reader(self,species):
    #determine charges and create arrays of charges and D's
        k=-1
        charges=[]
        number_of_catoms=[]
        noc_requested=False
        products=[]
        educts=[]
        catom_requested_values=[]
        #first find species which is educt:
        for sp in self.species:
            if 'flux' in self.species[sp]:
                if type(self.species[sp]['flux'])==dict:
                    if self.species[sp]['flux']['kind']=='educt' and self.species[sp]['flux']['reqs']=='number_of_catoms':
                        educts.append(sp)
                        catom_requested_values+=self.species[sp]['flux']['values']
                    elif self.species[sp]['flux']['kind']=='product' and self.species[sp]['flux']['reqs']=='number_of_catoms':
                        products.append(sp)
                        catom_requested_values+=self.species[sp]['flux']['values']
                    if self.species[sp]['flux']['reqs']=='number_of_catoms':
                        noc_requested=True


        for sp in species:
            k+=1
            number_of_catoms.append(0.0)
            fstring=species[sp]['symbol']

            #first extract charges from array
            string=fstring.split('^')
            if not len(string)==1:
                string=string[-1]
                string=string.replace('{', '').replace('}','')
                if string[-1]=='-':
                    if len(string)>1:
                        no=-int(string[:-1])
                    else:
                        no=-1
                else:
                    if string[-1]=='+':
                        if len(string)>1:
                            no=int(string[:-1])
                        else:
                            no=1
                    else:
                        no=int(string)
            else:
                no=0
            species[sp]['charge']=no
            #create a few shorter arrays
            charges.append(no*unit_F)

            if sp in educts and len(educts)!=0:
                continue
            if sp not in catom_requested_values:
                continue
            #secondly extract number of c atoms
            check_number=False
            number=''
            for s in fstring:
                if check_number:
                    if s=='_':
                        continue
                    elif s.isdigit():
                        number+=s
                        continue
                    elif number!='':
                        number_of_catoms[k]-=1
                        number_of_catoms[k]+=int(number)
                        check_number=False
                    else:
                        check_number=False
                if s=='C':
                    number_of_catoms[k]+=1
                    check_number=True
                else:
                    check_number=False
                    number=''
        charges=np.array(charges)
        number_of_catoms=np.array(number_of_catoms)
        self.number_of_catoms=number_of_catoms
        if noc_requested:
            self.logger.info('Number of C-Atoms for Evaluation of '+str(educts)+str(products)+' fluxes:')
            for isp,sp in enumerate(self.species):
                if sp not in educts and sp in catom_requested_values:
                    self.logger.info('  {}: {}'.format(sp,self.number_of_catoms[isp]))
        #the automatic reader of the number of catoms is maybe not what we want
        #for the reaction equivalents:
        return charges

    def set_boundary_and_initial_conditions(self,pb_bound,flux_bound):

        c_initial_general={}
        for isp,sp in enumerate(self.species):
            c_initial_general[str(isp)]=self.species[sp]['bulk concentration']

        self.set_initial_conditions(
                c_initial_general=c_initial_general)
                #c_initial_specific={'0':{'0':0.1},'1':{'0':0.1}})

        if pb_bound is None:
            self.pb_bound={
                'potential':    {'wall':    self.tp.system['vzeta'],
                                'bulk':     None},
                'gradient':     {'wall':    None,
                                'bulk':     0.0},
                }
        else:
            self.pb_bound={} 
            for key1 in ['potential','gradient']:
                self.pb_bound[key1]={}
                if key1 in pb_bound:
                    for key2 in ['wall','bulk']:
                        if key2 in pb_bound[key1]:
                            if pb_bound[key1][key2]=='zeta':
                                value=self.system['vzeta']
                            else:
                                value=pb_bound[key1][key2]
                            self.pb_bound[key1][key2]=value
                        else:
                            self.pb_bound[key1][key2]=None
                else:
                    self.pb_bound[key1]['wall']=None
                    self.pb_bound[key1]['bulk']=None

        #if 'vzeta' in self.system:
        #    self.vzeta_init=self.system['vzeta']
        #else:
        self.vzeta_init=None

        self.set_boundary_conditions(\
                flux_boundary=flux_bound,         #in mol/s/cm^2
                dc_dt_boundary={'all':{'r':0.0}},     #in mol/l/s #give either left OR right boundary condition here
                #integration will start at the site where dc_dt is defined
                efield_boundary={'l':0.0})    #in V/Ang


    def set_initial_concentrations(self,func,vzeta=None):
        """we can set the initial concentrations to a particular function"""
        if func=='Gouy-Chapman':
            if vzeta is None:
                vzeta=self.system['vzeta']
            else:
                self.vzeta_init=vzeta
            if self.nspecies!=2:
                self.logger.error('Gouy-Chapman limit only implemented for two species, cationic'+\
                        'and anionic. Not applying initialization.')
                return
            function=self.gouy_chapman

        c_initial_specific={}
        for k,sp in enumerate(self.species):
            c_initial_specific[str(k)]={}
            for i in range(self.nx):
                c_initial_specific[str(k)][str(i)]=\
                    self.species[sp]['bulk concentration']*\
                    np.exp(-self.beta*self.charges[k]*function(self.xmesh[i],vzeta=vzeta)[0])
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
        self.logger_db.debug('Shape off solution vector = '.format(np.shape(sol)))
        self.logger_db.debug('Initial solution vector = '.format(sol[:,0]))
        #b_min = minimize(func, b0, args=(z0, m, k, g, a), constraints=cons)
        self.ax1.plot(self.xmesh,sol[:,0]/1e10,'-',label='E (V/Ang)')
        self.ax2.plot(self.xmesh,sol[:,1],'-',label='phi (V)')
        self.ax1.legend()
        self.ax2.legend()
        return sol

       
    def gouy_chapman(self,x,vzeta=None):
        if vzeta is None:
            vzeta=self.system['vzeta']
        def func(x):
            term1 = 1.+np.tanh(vzeta*self.beta*unit_F/4.)*\
                np.exp(-1./self.debye_length*x)
            term2 = 1.-np.tanh(vzeta*self.beta*unit_F/4.)*\
                np.exp(-1./self.debye_length*x)
            return 2./(self.beta*abs(self.charges[0]))*np.log(term1/term2)
        grad=(func(x+1e-10)-func(x-1e-10))/(2*1e-10)
        return func(x),grad

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
            self.logger.error('No boundary conditions defined, stopping here.')
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

        if not any([type(self.species[sp]['flux'])==str for sp in self.species]):
            dc_dt_bound=allocate_fluxes(dc_dt_boundary,dc_dt_bound)
            flux_bound=allocate_fluxes(flux_boundary,flux_bound)
            flux_bound=np.array(flux_bound)
            self.flux_bound=flux_bound #*100**2

        dc_dt_bound=np.array(dc_dt_bound)
        efield_bound=np.array(efield_bound)

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

