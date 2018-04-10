"""
transport class. defines the transport model.
---
Stefan Ringe (stefan.ringe.tum@gmail.com)
"""
import scipy.integrate as integrate
import numpy as np
#import matplotlib.pyplot as plt
#from ase import units
from units import *
from itertools import cycle
import sys
from copy import deepcopy
import collections
import logging
import os
import re
from io import sync_mpi,save_all #,MPIFileHandler
import subprocess
from glob import glob
from shutil import copy
import imp
use_mpi=False
try:
    imp.find_module('mpi4py')
    use_mpi=True
except ImportError:
    pass
use_mpi=False
if use_mpi:
    from mpi4py import MPI
use_mpi=False

class Transport(object):

    def __init__(self, species=None,electrode_reactions=None,electrolyte_reactions=None,\
            system=None,pb_bound=None,nx=100,\
            descriptors=None,model_name=None,\
            comsol_args={},catmap_args={},only_plot=False):
        """
        only_plot   only initialize transport without creating folders
        """

        if only_plot:
            return

        #MPI setup
        if use_mpi:
            self.mpi_comm = MPI.COMM_WORLD
            comm=self.mpi_comm
            self.mpi_rank = self.mpi_comm.Get_rank()
            rank=self.mpi_rank
            self.mpi_size = self.mpi_comm.Get_size()
            size=self.mpi_size
        else:
            self.mpi_rank=rank=0
            self.mpi_size=size=1

        ##############################################
        ###########FOLDERS AND FILES##################
        ##############################################

        #folder and file names
        if model_name is None:
            self.model_name='catint'
        else:
            self.model_name=model_name
        self.outputfoldername=self.model_name+'_results' #folder where all results will be saved with the self.save function
        self.inputfilename=sys.argv[0] #the input file

        if rank==0:
            if not os.path.exists(self.outputfoldername):
                os.makedirs(self.outputfoldername)
            else:
                existing_files=sorted([f for f in os.listdir('.') if re.search(self.outputfoldername+'_[0-9]+', f)])
                if len(existing_files)>0: #self.outputfoldername.split('_')[-1].isdigit():
                    number=int(existing_files[-1].split('_')[-1])+1
                else:
                    number=2
                self.outputfoldername='_'.join(sum([[self.outputfoldername],[str(number).zfill(4)]],[]))
                print 'making folder',rank,self.outputfoldername
                os.makedirs(self.outputfoldername)
            self.logfilename=self.outputfoldername+'/transport.log' # the log file
        else:
            self.logfilename=None
            self.outputfoldername=None

        if use_mpi:
            self.logfilename=sync_mpi(self.logfilename)
            self.outputfoldername=sync_mpi(self.outputfoldername)

        #copy input file for later reference
        if rank==0:
            copy(self.inputfilename,self.outputfoldername+'/'+self.inputfilename)

        #wait for all processors
        if use_mpi:
            comm.Barrier()

        logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='.'.join(self.logfilename.split('.')[:-1])+'_id'+str(rank).zfill(3)+'.log',
                    filemode='w')

        if rank==0:
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
        try:
            logging.info('Starting Transport Calculation. Current Version: {}'.format(subprocess.check_output(["git", "describe","--always"]).strip()))
        except:
            logging.info('Starting Transport Calculation. Current Version not available.')
        
        # Now, define a couple of other loggers which might represent areas in your
        # application:
        
        self.logger_db = logging.getLogger('transport.debug')
        self.logger = logging.getLogger('transport.info')

        ##############################################
        ###########DICTIONARY HANDLING################
        ##############################################
        
        #all the possible keys:
        species_keys=['bulk concentration', 'diffusion', 'name', 'symbol', 'flux','current density','flux-equation']
        system_keys=[
                'phiM',                     #V
                'Stern capacitance',        #mF/cm^2
                'bulk_pH',
                'phiPZC',                   #V
                'temperature',              #K
                'pressure',     
                'water viscosity',
                'electrolyte viscosity',
                'epsilon',                  #epsilon_0
                'migration',
                'electrode reactions',
                'electrolyte reactions',
                'boundary thickness',       #m
                'exclude species',
                'active site density',       #mol/m^2
                'current density',
                'RF',                       #roughness factor
                ]

        #go over input data and put in some defaults if none
        if species is None:
            self.species={'species1':       {'symbol':r'K^+',
                                            'name':'potassium',
                                            'diffusion':1.96e-9,
                                            'kind':'electrolyte',
                                            'bulk concentration':0.001*1000.},
                          'species2':       {'symbol':r'HCO_3^-',
                                            'name':'bicarbonate',
                                            'diffusion':1.2e-9,
                                            'kind':'electrolyte',
                                            'bulk concentration':0.001*1000.}}
        else:
            for sp in species:
                for key in species[sp]:
                    if key not in species_keys:
                        self.logger.error('No such key "'+key+'" in species list. Quitting here.')
                        sys.exit()
            self.species=species

        #convert to ordered dictionary to have consistent indices of species
        self.species=collections.OrderedDict(self.species)
    
        if pb_bound is None:
            pb_bound={
            'potential': {'wall':'phiM'},
            'gradient': {'bulk':0.0}}


        system_defaults={
                'epsilon':78.36,
                'Stern epsilon':2.0, #only important for mesh definition (no actual parameter in model)
                'Stern capacitance':18e-2, #Stern layer capacitance for Robin BCs
                'temperature':298.14,
                'phiM':0.0, #-0.0125, #potential vs SHE
                'phiPZC':0.0,       #PZC potential
                'electrode reactions': False,
                'electrolyte reactions': False,
                'exclude species': ['H2O','e-'],
                'pressure':1}
        if system is None:
            self.system=system_defaults
        else:
            for key in system:
                if key not in system_keys:
                    self.logger.error('No such key "'+key+'" in system list. Quitting here.')
                    sys.exit()
            self.system=system


        for key in system_defaults:
           # ['epsilon','temperature','pressure','phiM']:
            if key not in self.system:
                self.system[key]=system_defaults[key]

        if 'exclude species' not in self.system:
            self.system['exclude species']=['e-']
        elif 'e-' not in self.system['exclude species']:
            self.system['exclude species']+=['e-']
        if 'exclude species' not in self.system:
            self.system['exclude species']+=['H2O']
        elif 'H2O' not in self.system['exclude species']:
            self.system['exclude species']+=['H2O']
        #delete species which should not be considered for PNP dynamics
        for es in self.system['exclude species']:
            if es in self.species:
                del self.species[es]
        self.logger.info('Excluding {} from PNP transport. They will also not participate in reactions (activity = 1)'.format(self.system['exclude species']))

        self.nspecies=len(self.species)


        #add concentration = 0 for all species with no separately defined values
        for sp in self.species:
            if 'bulk concentration' not in self.species[sp]:
                self.species[sp]['bulk concentration']=0.0

        #get pH
        if 'bulk_pH' in self.system:
            self.logger.info('pH given in system list, updating H+ and OH- concentrations if species exist')
            if 'H+' in self.species:
                self.species['H+']['bulk concentration']=10**(-self.system['bulk_pH'])*1000.
            elif 'OH-' in self.species:
                self.species['OH-']['bulk concentration']=10**(-(14-self.system['bulk_pH']))*1000.
        else:
            if 'H+' in self.species:
                self.system['bulk_pH']=-np.log10(self.species['H+']['bulk concentration']/1000.)
            elif 'OH-' in self.species:
                self.system['bulk_pH']=14+np.log10(self.species['OH-']['bulk concentration']/1000.)
            else:
                #setting the pH to an arbitrary value, it is not relevant here
                self.system['bulk_pH']=7.0


        self.system['surface_pH']=self.system['bulk_pH']

        #initialize concentrations at electrode
        for sp in self.species:
            self.species[sp]['surface concentration']=self.species[sp]['bulk concentration']



        self.eps = self.system['epsilon']*unit_eps0 #1.1e11 #*unit_eps0
        self.beta = 1./(self.system['temperature'] * unit_R)

        #GET CHARGES AND NCATOMS FROM CHEMICAL SYMBOLS
        self.charges=self.symbol_reader(self.species)

        self.use_migration=True
        if 'migration' in self.system:
            if not self.system['migration']:
                self.use_migration=False

        self.use_convection=False

        ######################
        #REACTIONS
        ######################
        #if reactions are given as input, they will always be used!
        #use_electrolyte/electrode_reactions is not really needed 

        #Working on electrolyte reactions if requested
        self.electrolyte_reactions=electrolyte_reactions
        self.use_electrolyte_reactions=False
        if 'electrolyte reactions' in self.system:
            if not self.system['electrolyte reactions']:
                self.use_electrolyte_reactions=False
            elif self.system['electrolyte reactions']:
                self.use_electrolyte_reactions=True
        if self.electrolyte_reactions is not None:
            if any(['rates' in self.electrolyte_reactions[reaction] for reaction in self.electrolyte_reactions]):
                self.use_electrolyte_reactions=True
                self.logger.info('Found electrolyte reactions with specified rates. Switching electrolyte reactions on. Preparing...')
                self.electrolyte_reactions=self.initialize_reactions(self.electrolyte_reactions)
        if self.electrolyte_reactions is None and self.use_electrolyte_reactions:
            self.logger.error('Electrolyte reactions were requested by input, but no electrolyte reaction was defined. Define electrolyte reaction first.')
            sys.exit()

        if self.use_electrolyte_reactions:
            for el in self.electrolyte_reactions:
                if not 'rates' in self.electrolyte_reactions[el]:
                    self.logger.info('Reaction {} has no rates given. It will not be considered for PNP dynamics!'.format(el))
                for sp in sum(self.electrolyte_reactions[el]['reaction'],[]):
                    if sp not in self.species and sp not in self.system['exclude species']:
                        self.logger.error('Species {} has not been defined, but is used in the electrolyte reactions, define it first!'.format(sp))
                        self.logger.error('  This is the current species list:')
                        for sp in self.species:
                            self.logger.error('  {}'.format(sp))
                        sys.exit()

        #Working on electrode reactions if requested
        self.electrode_reactions=electrode_reactions
        self.use_electrode_reactions=False
        if 'electrode reactions' in self.system:
            if not self.system['electrode reactions']:
                self.use_electrode_reactions=False
            elif self.system['electrode reactions']:
                self.use_electrode_reactions=True

        if self.electrode_reactions is not None:
            self.use_electrode_reactions=True
            self.logger.info('Found electrode reactions. Switching electrode reactions on. Preparing...')
            self.electrode_reactions=self.initialize_reactions(self.electrode_reactions)
        elif self.electrode_reactions is None and self.use_electrode_reactions:
            self.logger.error('Electrode reactions were requested by input, but no electrode reaction was defined. Define electrode reaction first.')
            sys.exit()

        if self.use_electrode_reactions:
            for el in self.electrode_reactions:
                for sp in sum(self.electrode_reactions[el]['reaction'],[]):
                    print 'TEST',sp
                    if sp not in self.species and sp not in self.system['exclude species'] and not sp.startswith('*'):
                        self.logger.error('Species {} has not been defined, but is used in the electrode reactions, define it first!'.format(sp))
                        self.logger.error('  This is the current species list:')
                        for sp in self.species:
                            self.logger.error('  {}'.format(sp))
                        sys.exit()

        #sort different species into lists:
        self.product_list=[]
        self.educt_list=[]
        self.electrolyte_list=[]

        if self.use_electrode_reactions:
            for sp in self.electrode_reactions:
                self.product_list.append(sp)
                for rr in self.electrode_reactions[sp]['reaction'][0]:
                    if rr != sp and rr not in ['e-'] and rr not in self.system['exclude species'] and rr not in self.educt_list:
                        self.educt_list.append(rr)
                for rr in self.electrode_reactions[sp]['reaction'][1]:
                    if rr != sp and rr not in ['e-'] and rr not in self.system['exclude species'] and rr not in self.product_list:
                        self.product_list.append(rr)


        for sp in self.species:
            if sp not in self.product_list and sp not in self.educt_list:
                self.electrolyte_list.append(sp)

        if self.use_electrode_reactions:
            self.logger.info('Educts: {}'.format(self.educt_list))
            self.logger.info('Products: {}'.format(self.product_list))
        if self.use_electrolyte_reactions:
            self.logger.info('Electrolyte Components: {}'.format(self.electrolyte_list))

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

        self.xmesh_init=self.xmesh
        self.nx_init=self.nx
        self.xmax_init=self.xmax

        if self.debye_length>self.xmax/2.:
            self.logger.warning('Debye length is larger than 1/4th of the xmesh. Take care that the x discretization is not too coarse!.')
            self.logger.warning('Current xmesh: xmax={}, dx={} at debye_length={}'.format(self.xmax,self.dx,self.debye_length))

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


        #STATIONARY-STATE REACTION RATES/FLUXES
        self.initialize_fluxes()


        #FLUX AND FARADAIC YIELD
        flux_bound={} 
        for isp,sp in enumerate(self.species):
            flux_bound[str(isp)]={}
        for isp,sp in enumerate(self.species):
            #if 'flux' in self.species[sp]:
            #    pass
            #    #self.species[sp]['flux']*=-1
            #else:
            #    self.species[sp]['flux']=0.0
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

        self.initialize_descriptors(descriptors)

        if hasattr(self,'descriptors') and use_mpi:
            ntasks=np.prod(map(len,[self.descriptors[key] for key in self.descriptors]))
        else:
            ntasks=1
        if size!=ntasks and use_mpi:
            self.logger.error('# of CPUs is different from # of tasks. This is currently not supported.')
            sys.exit()

        #INITIALIZE PARAMETERS OF EXTERNAL SOFTWARE
        #-- CATMAP
        self.catmap_args=catmap_args
        #-- COMSOL
        self.initialize_comsol(comsol_args)

    def initialize_comsol(self,comsol_args):
        """
        Short description of keys:
            solver:             parametric or simple
            studies:            stat or time
            grid_factor:        factor that determines fine-ness of grid
            bin_path:           comsol executable
            global_variables:   variables defined on the entire geometry
            boundary_variables: variables defined on a boundary
            global_equations:   differential equations defined on entire geometry
            outputs:            outputs to be created
            nflux:              number of steps in increasing flux in parametric sweep
        """
        comsol_keys=['outputs','boundary_variables','global_variables','global_equations','parameter','bin_path','nflux',\
            'par_name','par_values','desc_method','model_type','solver','studies']

        #tp_dilute_species or porous_electrode
        if 'model_type' not in comsol_args:
            comsol_args['model_type']='tp_dilute_species'

        if 'studies' not in comsol_args:
            comsol_args['studies']=['stat']

        if 'solver' not in comsol_args:
            #parametric or simple
            comsol_args['solver']='parametric'

        for a in ['outputs','global_variables','boundary_variables','parameter']:
            if not a in comsol_args:
                comsol_args[a]={}
        if 'grid_factor' not in comsol_args['parameter']:
            comsol_args['parameter']['grid_factor']=['200','Fineness of Grid']

        #COMSOL Treatment of descriptors
        #if the descriptor is the potential, we can use it internally in COMSOL
        #as a parametric sweep. this kills to birds with a single shot, since
        #slowly turning on the potential is like turning on the flux, but at the 
        #same time evaluates the potential dependence of the current density, so
        #we do not need to recompile COMSOL for every potential value

        ### "desc_method" defines if/how to treat descriptors as parameter sweeps
        #internal: descriptor given here will be used for iterating internally in comsol
        #  if this selected descriptor is not the potential this could lead to convergence problems!
        #  internal-reinit: at each parameter set, the solutions are reinitialized (default)
        #  internal-cont: the solution of the previous parameter set is used to initialize the next
        #external: parameters are updated within this comsol.py routine and comsol 
        #  is recompiled and relaunched for each new parameter set
        ### "par_name" selects a parameter from self.descriptors to use as parameter sweep
        ### "par_values" lists the value of this parameter 

        if 'RF' not in self.system:
            if 'RF' in comsol_args['parameter']:
                print 'in RF'
                self.system['RF']=float(comsol_args['parameter']['RF'][0])
            else:
                comsol_args['parameter']['RF']=['1','Roughness Factor']
                self.system['RF']=1.0
        else:
            comsol_args['parameter']['RF']=[str(self.system['RF']),'Roughness Factor']

        if not 'nflux' in comsol_args:
            comsol_args['nflux']=32 #default to 32 steps for flux ramping
        for a in comsol_args:
            if a not in comsol_keys:
                self.logger.error('{} is not a standard key of COMSOL. Implement this first. Exiting here to be sure that this key is what you want'.format(a))
                sys.exit()

        if 'desc_method' not in comsol_args:
            #default is using flux sweep and external looping over descriptors
            comsol_args['desc_method']='external'
            if 'par_name' not in comsol_args:
                comsol_args['par_name']='flux_factor'
                comsol_args['par_values']=np.linspace(0,1,comsol_args['nflux'])
        if comsol_args['desc_method']=='external':
            if 'par_name' not in comsol_args:
                comsol_args['par_name']='flux_factor'
                comsol_args['par_values']=np.linspace(0,1,comsol_args['nflux'])

        if comsol_args['desc_method'].startswith('internal'):
            if self.use_catmap:
                self.logger.warning('CatMAP does not work together with passing descriptors as COMSOL'+
                    ' parameter sweep, changing desc_method to external')
                comsol_args['desc_method']='external'
                if 'par_name' not in comsol_args:
                    comsol_args['par_name']='flux_factor'
                    comsol_args['par_values']=np.linspace(0,1,comsol_args['nflux'])
        if comsol_args['desc_method'].startswith('internal'):
            if 'par_name' not in comsol_args:
                self.logger.error('Descriptor was requested to be used as parameter sweep inside COMSOL'+
                    ', however no descriptor was selected. Select one by setting the par_name key in'+
                    'the comsol_args dictionary')
                sys.exit()
            if comsol_args['par_name'] not in self.descriptors:
                self.logger.error('Selected descriptor for internal COMSOL ramping {} was not found'+
                    ' in global descriptor list.')
                sys.exit()
            if 'par_values' not in comsol_args:
                comsol_args['par_values']=self.descriptors[comsol_args['par_name']]
            if comsol_args['par_name']!='phiM':
                self.logger.warning('Selected Descriptor for COMSOL parameter sweep'+
                    'is not the potential phiM. This is not recommended, be sure that'+
                    'you do not get convergence issues')

        self.comsol_args=comsol_args


    def initialize_fluxes(self):

        if not self.use_electrode_reactions:
            for sp in self.species:
                self.species[sp]['flux']=0.0
            return

        self.use_catmap=False

        for sp in self.species:
            if 'flux' in self.species[sp]:
                if self.species[sp]['flux']=='catmap':
                    self.use_catmap=True
        if self.use_catmap:
            for sp in self.species:
                if not 'flux' in self.species[sp]:
                    self.species[sp]['flux']='catmap'
            self.logger.info('Found flux = catmap, all fluxes will be calculated by CatMAP.')
            return
                
        #some consistency check, either flux, current density or flux-equation should be given, NOT all
        for sp in self.species:
            count=0
            for key in self.species[sp]:
                if key in ['flux','current density','flux-equation']:
                    count+=1
            if count>1:
                self.logger.error('Flux of species {} has been defined by more than one method.'.format(sp))
                sys.exit()
                
        #another check, the current density method should be only selected if the species is a product:
        for sp in self.species:
            if 'current density' in self.species[sp]:
                if sp not in [prod.split('-')[0] for prod in self.electrode_reactions]:
                    self.logger.error('Flux of species {} has been given as current density but this species is not product.'.format(sp))

        #last check, check if more than one flux per equation has been defined which is not necessary!
        ers=self.electrode_reactions
        for er in ers:
            educts=ers[er]['reaction'][0]
            products=ers[er]['reaction'][1]
            count=0
            for ep in sum([educts,products],[]):
                if ep in self.system['exclude species']:
                    continue
                if '*' in ep:
                    continue
                if any([a in ['flux','current density','flux-equation'] for a in self.species[ep]]):
                    count+=1
            if count>1:
                self.logger.error('More than one flux has been defined for equation {}. Select one of the fluxes, the rest will be automatically calculated.'.format(ers[er]['reaction']))
                sys.exit()
            if count==0:
                self.logger.error('No flux defined in equation {}. Define one flux.'.format(ers[er]['reaction']))
                sys.exit()

        #first search fluxes. if any flux is given as equation, we create a comsol parameter first for
        #the fixed fluxes
        using_equations=False
        for sp in self.species:
            if 'flux-equation' in self.species[sp]:
                using_equations=True
                break

        #is this a reduction or oxidation?
        for er in ers:
            educts=ers[er]['reaction'][0]
            products=ers[er]['reaction'][1]
            if 'e-' in educts:
                e_in_educts=True
            elif 'e-' in products:
                e_in_educts=False
            else:
                self.logger.error('No electron found in the reactions.')
                sys.exit()

        #convert given current densities into fluxes
        if using_equations:
            #convert fixed fluxes to strings
            for sp in self.species:
                if 'flux' in self.species[sp]:
                    self.species[sp]['flux']=str(self.species[sp]['flux'])
                elif 'current density' in self.species[sp]:
                    if e_in_educts:
                        sign='(-1)'
                    else:
                        sign='1'
                    for reac in self.electrode_reactions:
                        if sp==reac.split('-')[0]:
                            nprod=len([a for a in self.electrode_reactions[reac]['reaction'][1] if a==sp])
                            self.species[sp]['flux']=sign+'*'+str(self.species[sp]['current density']/self.electrode_reactions[reac]['nel']/unit_F*nprod)
                elif 'flux-equation' in self.species[sp]:
                    self.species[sp]['flux']=self.species[sp]['flux-equation']
        else:
            for sp in self.species:
                if 'current density' in self.species[sp]:
                    if e_in_educts:
                        sign=-1
                    else:
                        sign=1
                    for reac in self.electrode_reactions:
                        if sp==reac.split('-')[0]:
                            nprod=len([a for a in self.electrode_reactions[reac]['reaction'][1] if a==sp])
                            self.species[sp]['flux']=sign*self.species[sp]['current density']/self.electrode_reactions[reac]['nel']/unit_F*nprod

        ers=self.electrode_reactions
        #first set the fluxes of all products
        missing_species=[]
        for er in ers:
            #check for undefined fluxes
            for reac in sum(ers[er]['reaction'],[]):
                if reac not in ers and reac not in ['e-'] and reac not in self.system['exclude species'] and reac not in missing_species and not reac.startswith('*'):
                    missing_species+=[reac]
            if len(missing_species)>0:
                self.logger.info('Calculating fluxes of {} as sum of other fluxes'.format(missing_species))
        ref_sp={}
        for er in ers:
            educts=ers[er]['reaction'][0]
            products=ers[er]['reaction'][1]
            #determine the species from which flux the remaining fluxes will be calculated
            for ep in sum([educts,products],[]):
                if ep in self.system['exclude species']:
                    continue
                if '*' in ep:
                    continue
                if any([a in ['flux','current density','flux-equation'] for a in self.species[ep]]):
                    ref_sp[er]=ep
                    break
        #calculate fluxes of remaining species by using the fluxes of products evaluated before
        for er in ers:
            educts=ers[er]['reaction'][0]
            products=ers[er]['reaction'][1]
            sp=ref_sp[er]
            #adjust fluxes of all missing species
            for missing in missing_species:
                if missing not in educts and missing not in products:
                    continue
                count_missing=max(educts.count(missing),products.count(missing))*1.
                count_ref=max(educts.count(sp),products.count(sp))*1.
                both_in_educts=(missing in educts)*(sp in educts)
                both_in_products=(missing in products)*(sp in products)
                if both_in_educts or both_in_products:
                    if using_equations:
                        factor='1'
                    else:
                        factor=1
                else:
                    if using_equations:
                        factor='(-1)'
                    else:
                        factor=(-1.)
                if not using_equations and 'flux' not in self.species[missing]:
                    self.species[missing]['flux']=0.0
                elif using_equations and 'flux' not in self.species[missing]:
                    self.species[missing]['flux']='0'
                if using_equations:
                    self.species[missing]['flux']+='+'+str(factor)+'*'+self.species[sp]['flux']+'*'+str(count_missing/count_ref)
                else:
                    self.species[missing]['flux']+=factor*self.species[sp]['flux']*count_missing/count_ref

        #finally set all remaining fluxes equal to zero
        for sp in self.species:
            if 'flux' not in self.species[sp]:
                if using_equations:
                    self.species[sp]['flux']='0.0'
                else:
                    self.species[sp]['flux']=0.0


    def initialize_reactions(self,reactions):
        """ replace the reaction string by a list"""
        for reaction in reactions:
            string=reactions[reaction]['reaction']
            if '<->' in string:
                si=sum([ab.split('->') for ab in string.split('<->')],[])
            else:
                si=string.split('->')
            new_si=[]
            for sii in si:
                new_si.append(sii.strip())
            si=new_si
            final_reactants=[]
            for isii,sii in enumerate(si):
                reactants=sii.split(' + ') #list of all reactants of current equation side
                reactants_out=[]
                for ir,reactant in enumerate(reactants):
                    nel=re.findall('([0-9]+)[ ]+e-',reactant)
                    if len(nel)>0:
                        #this is e-, have to count this
                        reactions[reaction]['nel']=int(nel[0])
                    N=re.findall('([0-9]+)[ ]+[*A-Za-z]+',reactant)
                    if len(N)<1:
                        reactants_out.append(reactant.strip())
#                        reactants_out.append([r.strip() for r in reactants if len(r)>0])
                        continue
                    else:
                        N=int(N[0])
#                        del reactants_out[ir]
                        for n in range(N):
                            reactants_out.append(reactant[len(str(N))+1:].strip())
#                        reactants_out=[r.strip() for r in reactants_out if len(r)>0]
                final_reactants.append(reactants_out)
            reactions[reaction]['reaction']=final_reactants
        return reactions


    def initialize_descriptors(self,descriptors):
        #if internal comsol iteration should be used the potential should be 
        #like 0,-0.1,-0.2,-0.3...
        #or 0,0.1,0.2,0.3...

        #list of descriptors over which to iterate
        if descriptors is not None:
            if any([type(descriptors[desc]) not in [list,np.array] for desc in descriptors]):
                self.logger.error('Descriptors must be given as list. Stopping here for safety')
                sys.exit()
            self.descriptors=descriptors
        else:
            #no descriptors given at input, add some here for convenience:
            self.descriptors={}
            self.logger.warning('No descriptor list given at input, performing single point calculation')
            self.descriptors['phiM']=[self.system['phiM']]
            self.descriptors['temperature']=[self.system['temperature']]
            return
#            self.descriptors
            #no descriptors specified. put one here just that the following routine work fine
#            self.descriptors={'voltage':[self.system['phiM']],'temperature':[self.system['temperature']]}

        desc_keys=[key for key in self.descriptors]
        if len(desc_keys)==1:
            self.logger.debug('Adding a dummy descriptor for convenience in the code')
            if 'temperature' not in desc_keys:
                self.descriptors['temperature']=[self.system['temperature']]
            else:
                self.descriptors['phiM']=[self.system['phiM']]


        desc_keys=[key for key in self.descriptors]
        if len(desc_keys)!=2:
            self.logger.error('Cannot use other than 2 descriptors')
            sys.exit()

        for desc in desc_keys:
            if desc not in self.system:
                self.logger.error(desc+' not found in system list, cannot evaluate other than system descriptors, yet')
                self.logger.error('Here is the current system list:\n{}'.format(self.system))

        #self.alldata_names
        #gives the values of the descriptors for the alldata array
        #suppose we have:
        #   self.descriptors['desc_1']=[1,2,3]
        #   self.descriptors['desc_2']=[10,11,12,13]
        #then
        #   self.alldata_name=[[1,10],[1,11],[1,12],...,[3,12],[3,13]]
        #self.alldata
        #contains all data on the same grid
        #each grid point contains
        self.alldata=[]
        self.alldata_names=[]
        i=-1
        for value1 in self.descriptors[desc_keys[0]]:
            for value2 in self.descriptors[desc_keys[1]]:
                i+=1
                self.alldata_names.append([value1,value2])
                self.alldata.append({'species':{},'system':{}})
                for sp in self.species:
                    self.alldata[i]['species'][sp]={}

#    def evaluate_fluxes(self):
#        self.logger.info('Evaluating fluxes of {} as sum of products/educts'.format([sp for sp in self.species if type(self.species[sp]['flux'])==dict]))
#        tmp_rates=np.zeros([self.nspecies])
#        for isp,sp in enumerate(self.species):
#            if 'zeff' in self.species[sp]:
#                if type(self.species[sp]['flux'])!=dict:
#                    crate=self.species[sp]['flux']/(unit_F*self.species[sp]['zeff'])
#                    self.species[sp]['flux']=crate
##                if sp not in sum([self.system['sum rates'][sp]['values'] for sp in self.system['sum rates']],[]):
##                    continue
#                for isp2,sp2 in enumerate(self.species):
#                    if type(self.species[sp2]['flux'])==dict:
#                        if sp in self.species[sp2]['flux']['values']:
#                            if self.species[sp2]['flux']['reqs']=='number_of_catoms':
#                                req=self.number_of_catoms[isp]
#                            elif self.species[sp2]['flux']['reqs']=='zeff':
#                                req=self.species[sp]['zeff']
#                            else:
#                                req=self.species[sp2]['flux']['reqs'][isp]
#                            if self.species[sp2]['flux']['kind']=='educt':
#                                tmp_rates[isp2]-=crate*req
#                            else:
#                                tmp_rates[isp2]+=crate*req
##            elif sp=='HCO3-':
##                self.species[sp]['flux']=-self.species['H2']['flux']/(unit_F*self.species['H2']['zeff'])*2.
#        for isp,sp in enumerate(self.species):
#            if type(self.species[sp]['flux'])==dict:
#                self.species[sp]['flux']=tmp_rates[isp]
##        for ieduct,educt in enumerate(self.system['educts']):
##            self.species[educt]['flux']=-tmp_educts[ieduct]
#        if 'unknown' in self.species:
#            #we needed this flux only to calculate the co2 and oh- fluxes
#            self.species['unknown']['flux']=0.0

    def update_alldata(self,desc_value1,desc_value2):
        """updates the alldata array with current system and species dictionaries for the descriptors [desc_value1,desc_value2]"""
        #get index
        index=self.alldata_names.index([desc_values[0],desc_values[1]])
        self.alldata[index]['system']=self.tp.system.copy()
        self.alldata[index]['species']=self.tp.species.copy()

    def symbol_reader(self,species):
    #determine charges and create arrays of charges and D's
        k=-1
        charges=[]

        for sp in species:
            k+=1
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

        charges=np.array(charges)
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
                'potential':    {'wall':    self.tp.system['phiM'],
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
                            if pb_bound[key1][key2]=='phiM':
                                value=self.system['phiM']
                            else:
                                value=pb_bound[key1][key2]
                            self.pb_bound[key1][key2]=value
                        else:
                            self.pb_bound[key1][key2]=None
                else:
                    self.pb_bound[key1]['wall']=None
                    self.pb_bound[key1]['bulk']=None

        #if 'phiM' in self.system:
        #    self.phiM_init=self.system['phiM']
        #else:
        self.phiM_init=None

        self.set_boundary_conditions(\
                flux_boundary=flux_bound,         #in mol/s/cm^2
                dc_dt_boundary={'all':{'r':0.0}},     #in mol/l/s #give either left OR right boundary condition here
                #integration will start at the site where dc_dt is defined
                efield_boundary={'l':0.0})    #in V/Ang


    def set_initial_concentrations(self,func,phiM=None):
        """we can set the initial concentrations to a particular function"""
        if func=='Gouy-Chapman':
            if phiM is None:
                phiM=self.system['phiM']
            else:
                self.phiM_init=phiM
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
                    np.exp(-self.beta*self.charges[k]*function(self.xmesh[i],phiM=phiM)[0])
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

       
    def gouy_chapman(self,x,phiM=None):
        if phiM is None:
            phiM=self.system['phiM']
        def func(x):
            term1 = 1.+np.tanh(phiM*self.beta*unit_F/4.)*\
                np.exp(-1./self.debye_length*x)
            term2 = 1.-np.tanh(phiM*self.beta*unit_F/4.)*\
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

    def save(self):
        if self.mpi_rank==0:
            self.logger.info('Saving all data into binary pickle files.')
            save_all(self)
