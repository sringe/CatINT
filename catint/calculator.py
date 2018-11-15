"""
calculator class
solves the transport equations
either by using finite difference techniques to solve the PNP equations
or by calling COMSOL
kinetics are implemented as fixed fluxes, rate equations or 
full microkinetic modeling via CatMAP
---
Stefan Ringe (stefan.ringe.tum@gmail.com)
"""
import scipy.integrate as integrate
from scipy.sparse import diags
from scipy.linalg import block_diag
from scipy.optimize import minimize
from scipy.integrate import ode
from scipy import interpolate
import numpy as np
#import matplotlib.pyplot as plt
#from ase import units
from units import *
import sys
from copy import deepcopy
import os
import math
from copy import deepcopy
import imp
from io import sync_mpi,reduce_dict_mpi
from comsol_wrapper import Comsol
from catmap_wrapper import CatMAP
from comsol_reader import Reader

#import mpi if available
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

#import odespy if available
try:
    imp.find_module('odespy')
    found = True
except ImportError:
    found = False
if found:
    import odespy

class Calculator():

    def __init__(self,transport=None,dt=None,tmax=None,ntout=1,calc=None,
            scale_pb_grid=None,tau_jacobi=1e-7,tau_scf=5e-5,mix_scf=0.5,mode='time-dependent'):
    
        self.mode=mode #calculation mode for comsol: time-dependent or stationary. the local
        #solvers here are all time-dependent
        self.tau_scf=tau_scf
        self.mix_scf=mix_scf

        if transport is None:
            print('No transport object provided for calculator. Stopping here.')
            sys.exit()
        else:
            self.tp=transport #transport object
        
        if calc is None:
            calc=self.tp.calc

        if self.tp.use_catmap:
            self.catmap=CatMAP(transport=self.tp)
        
        self.tp.ntout=ntout
#        if os.path.exists('results.txt'):
#            os.remove('results.txt')



        self.use_lax_friedrich=False
        string=calc.split('--')
        self.calc_method=None
        if len(string)>1:
            if string[-1]=='LF':
                self.use_lax_friedrich=True
            else:
                self.calc_method=string[-1]
            self.calc=string[0]
        else:
            self.calc=calc

        self.calc_list=['FTCS','Crank-Nicolson','odeint','vode','lsoda','dopri5','dop853','odeint','odespy','comsol']

        if self.calc not in self.calc_list:
            self.tp.logger.error('No calculator found with this name. Aborting.')
            sys.exit()
        
        #use a finer grid for PBE integration (factor of scale_pb_grid smaller than
        #normal xmesh
        self.scale_pb_grid=scale_pb_grid
        #accuracy of Jacoby iteration if potential needs to be outputted:
        self.tau_jacobi=tau_jacobi

        if dt is not None: # and hasattr(self.tp,'dt'):
            self.tp.dt=dt
        if tmax is not None: # and hasattr(self.tp,'tmax'):
            self.tp.tmax=tmax
       # if hasattr(self.tp,'dt') and hasattr(self.tp,'tmax'):
        if tmax is not None or dt is not None:
            self.tp.tmesh=np.arange(0,self.tp.tmax+self.tp.dt,self.tp.dt)
            self.tp.nt=len(self.tp.tmesh)
        else:
            self.tp.logger.warning('No time mesh given, defaulting to range(0,1,0.1)')
            self.tp.tmesh=np.arange(0,1,0.1)
            self.tp.nt=len(self.tp.tmesh)

        #save this initial mesh
        self.tp.tmesh_init=self.tp.tmesh
        self.tp.nt_init=self.tp.nt
        self.tp.dt_init=self.tp.dt
        self.tp.tmax_init=self.tp.tmax

        self.oldtime=np.inf
        
        #go over the times and decide which ones to output
        self.tp.itout=[]
        cc=0
        for it,t in enumerate(self.tp.tmesh):
            if it==self.tp.nt-1:
                cc+=1
                self.tp.itout.append(it)
            elif it>1 and it%int(self.tp.nt/float(self.tp.ntout)) == 0:
                cc+=1
                self.tp.itout.append(it)
        #set the output number correctly
        ntout=cc
        self.tp.ntout=ntout #

        self.initialize='bla'
        if self.calc == 'comsol':
            self.comsol=Comsol(transport=self.tp,mode=self.mode)
        return

    def get_rates(self,C):
        #get list of all species names:
        species_names=[sp for sp in self.tp.species]
        rates=np.zeros([self.tp.nspecies,self.tp.nx])
        for i in range(self.tp.nx):
            for r in self.tp.reactions:
                reaction=self.tp.reactions[r]
                if 'rates' not in reaction:
                    continue
                #rates of educts:
                for reactant in reaction['reactants'][0]:
                    if reactant not in self.tp.species:
                        continue
                    k=species_names.index(reactant)
                    rates[k,i]=0.0
                    prod=1.
                    for reactant2 in reaction['reactants'][0]:
                        if reactant2 not in self.tp.species:
                            continue
                        k2=species_names.index(reactant2)
                        prod*=C[k2,i]
                    rates[k,i]-=prod*reaction['rates'][0]
                    prod=1.
                    for reactant2 in reaction['reactants'][1]:
                        if reactant2 not in self.tp.species:
                            continue
                        k2=species_names.index(reactant2)
                        prod*=C[k2,i]
                    rates[k,i]+=prod*reaction['rates'][1]
                #rates of products
                for reactant in reaction['reactants'][1]:
                    if reactant not in self.tp.species:
                        continue
                    k=species_names.index(reactant)
                    rates[k,i]=0.0
                    prod=1.
                    for reactant2 in reaction['reactants'][0]:
                        if reactant2 not in self.tp.species:
                            continue
                        k2=species_names.index(reactant2)
                        prod*=C[k2,i]
                    rates[k,i]+=prod*reaction['rates'][0]
                    prod=1.
                    for reactant2 in reaction['reactants'][1]:
                        if reactant2 not in self.tp.species:
                            continue
                        k2=species_names.index(reactant2)
                        prod*=C[k2,i]
                    rates[k,i]-=prod*reaction['rates'][1]
        return rates

    def run(self):
        itask=0
        if self.tp.descriptors is not None and \
            (self.tp.comsol_args['desc_method'] == 'external' or self.tp.use_catmap):

            desc_keys=[key for key in self.tp.descriptors]
            i1=0
            i2=0
            for value1 in self.tp.descriptors[desc_keys[0]]:
                i1+=1
                i2=0
                for value2 in self.tp.descriptors[desc_keys[1]]:
                    i2+=1
                    itask+=1
                    #only proceed if this should be performed for current task
                    if itask%self.tp.mpi_size!=self.tp.mpi_rank:
                        continue

                    #update the iterating descriptor based property in system properties
                    self.tp.system[desc_keys[0]]=float(value1)
                    self.tp.system[desc_keys[1]]=float(value2)

                    self.tp.logger.info('| CI | -- | Starting calculation for {} = {} and {} = {} on CPU {} of {}'.format(desc_keys[0],value1,desc_keys[1],value2,self.tp.mpi_rank,self.tp.mpi_size))
                    label=str(i1).zfill(4)+'_'+str(i2).zfill(4)

                    if not self.tp.use_catmap:
                        self.run_single_step(label=label)
                    else:
                        self.run_scf_cycle(label=label)
                    
                    self.tp.save()
            #synchronize all_data over CPUs
#            self.tp.all_data=reduce_dict_mpi(self.tp.all_data)
        else:
        #internal treatment of parameter
            if not self.tp.use_catmap:
                self.run_single_step()
            else:
                self.tp.logger.error('| CI | -- | Catmap has been selected but current descriptor method is Comsol internal which does not work. Stopping here.')
                sys.exit()

            self.tp.save()

        if use_mpi:
            self.tp.comm.Barrier()

    def evaluate_accuracy(self,par,par_old,kind='relative'):
        #evaluate accuracy relative or absolute
        #automatic takes first the relative error, but if absolute error
        #is below treshhold also stops
        rmsd=0
        n=0
        par_sum=sum(par)
        #self.tp.logger.debug(' | Current current density = {} mA/cm^2'.format(par_sum))
        self.tp.logger.debug('| CI | -- | Flux = {} mol/m^2/s'.format(par_sum))
        par_old_sum=sum(par_old)
        errors=[]
        for p1,p2 in zip(par,par_old):
            if kind=='relative':
                if p1!=0:
                    errors.append(abs(p1-p2)/p1)
            else:
                errors.append(abs(p1-p2))
        error=max(errors)
#        rmsd=np.sqrt((par_sum-par_old_sum)**2)
#        for val1,val2 in zip(par,par_old):
#            n+=1
#            rmsd+=(val1-val2)**2
#        rmsd=np.sqrt(rmsd/(n*1.))
        return error

    def converged(self,old,new):
        cmrsd=0.0
        for isp in range(self.tp.nspecies):
            crmsd+=(old[isp*self.tp.nx]-new[isp*self.tp.nx])**2 #RMSD
        if np.sqrt(sum(crmsd)/len(crmsd))<self.tp.tau_scf:
            return True
        else:
            return False

    def run_scf_cycle(self,label=''):
        """run scf cycle to converge catmap-comsol iterations"""
        istep=0
        step_to_check=0
        scf_accuracy=np.inf
        self.tp.logger.info('| CI | -- | Starting iterative solution with CatMAP and COMSOL')
        self.tp.logger.info('| CI | -- |  using a current density accuracy cutoff of {} mV/cm^2 and a linear mixing parameter of {}'.format(self.tau_scf,self.mix_scf))

        #initialization
        if self.tp.system['init_folder'] is not None:
            self.tp.logger.info('| CI | -- | Initializing surface concentrations with {}'.format(self.tp.system['init_folder']))
            comsol_reader=Reader(transport=self.tp,results_folder=self.tp.system['init_folder'],\
                    outputs=['concentrations','electrostatics','electrode_flux','rho_charge'],comsol_args=self.tp.comsol_args)
            comsol_reader.read_all()
            for sp in self.tp.species:
                self.tp.logger.debug('| CI | -- | ci(x=0)_{} = {} M'.format(sp,self.tp.species[sp]['surface_concentration']/1000.))

        accuracies=[]
        desc_keys=self.tp.descriptors.keys()
        desc1_val=self.tp.system[desc_keys[0]]
        desc2_val=self.tp.system[desc_keys[1]]
        alldata_inx=self.tp.alldata_names.index([desc1_val,desc2_val])
        while scf_accuracy>self.tau_scf:
            istep+=1
            #evaluate how the accuracy during the last 50 steps, if it does not significantly decrease, reduce the mixing factor
            if istep-step_to_check>40: # and abs(accuracies[-2]-accuracies[-1])>1e-1:
                #still no convergence, try to decrease mixing factor
                self.mix_scf*=0.9
                self.tp.logger.info('| CI | -- | Accuracy is still < 1e-1, decreasing mixing factor in order to speed up the convergence')
                step_to_check=istep
            self.tp.logger.info('| CI | -- | Solving transport step {}. Current accuracy in current density = {} mA/cm^2'.format(istep,scf_accuracy))

            #linear mixing
            if istep>2:
                #mix fluxes
                for sp in self.tp.species:
#                    self.tp.species[sp]['flux']=self.mix_scf*self.tp.species[sp]['flux']+(1.-self.mix_scf)*fl_old[sp]
                    self.tp.species[sp]['surface_concentration']=self.mix_scf*self.tp.species[sp]['surface_concentration']+\
                            (1.-self.mix_scf)*sc_old[sp]
                    self.tp.alldata[alldata_inx]['species'][sp]['surface_concentration']=self.tp.species[sp]['surface_concentration']
            #update surface pH
            if 'H+' in self.tp.species:
                lss=self.tp.species['H+']['surface_concentration']
                self.tp.system['surface_pH']=-np.log10(float(lss)/1000.)
                self.tp.alldata[alldata_inx]['system']['surface_pH']=-np.log10(float(lss)/1000.)
            elif 'OH-' in self.tp.species:
                lss=self.tp.species['OH-']['surface_concentration']
                self.tp.system['surface_pH']=14+np.log10(float(lss)/1000.)
                self.tp.alldata[alldata_inx]['system']['surface_pH']=14+np.log10(float(lss)/1000.)

            #for linear mixing, safe all surface_concentrations
            sc_old={}
            for sp in self.tp.species:
                sc_old[sp]=self.tp.species[sp]['surface_concentration']

            #1) Microkinetic Model: CatMAP
            if istep==1:
                self.tp.logger.debug('| CI | -- | Electrode Fluxes:')
                for sp in self.tp.species:
                    self.tp.logger.debug('| CI | -- | {}: {}'.format(sp,self.tp.species[sp]['flux']))

            #run catmap. the fluxes will be updated automatically
            self.catmap.run()
            self.tp.logger.debug('| CI | -- | Electrode Fluxes:')
            for sp in self.tp.species:
                self.tp.logger.debug('| CI | -- | {}: {}'.format(sp,self.tp.species[sp]['flux']))

            #2) Transport: COMSOL

            self.tp.logger.debug('| CI | -- | Surface Concentrations:')
            for sp in self.tp.species:
                self.tp.logger.debug('| CI | -- | {}: {} mol/L'.format(sp,self.tp.species[sp]['surface_concentration']/1000.))

            #only_last updates the descriptor based dictionaries only for the last entry in self.tp.descriptors
            self.run_single_step(label=label)

            #3) Check Convergence
            current_density=[]
            for sp in self.tp.species: #self.tp.electrode_reactions:
#                current_density.append(self.tp.electrode_reactions[sp]['electrode_current_density'][(str(1.0),str(value2))])
#                current_density.append(self.tp.species[sp]['electrode_current_density'])
                if sp in self.tp.electrode_reactions:
                    nprod=len([a for a in self.tp.electrode_reactions[sp]['reaction'][1] if a==sp])
                    nel=self.tp.electrode_reactions[sp]['nel']
                else:
                    nprod=1
                    nel=1
                cd=self.tp.species[sp]['flux']*nel*unit_F/nprod/10.
                current_density.append(cd)
                self.tp.logger.debug('| CI | -- | Current Density of {} = {} mA/cm^2'.format(sp,cd))
            if istep>1:
                scf_accuracy=self.evaluate_accuracy(current_density,old_current_density)
            old_current_density=deepcopy(current_density)
            accuracies.append(scf_accuracy)
        if scf_accuracy<=self.tau_scf:
            self.tp.logger.info('| CI | -- | Iterative Solver converged in {} steps. Final accuracy in current density = {} mV/cm^2'.format(istep,scf_accuracy))

    def run_single_step(self,label=''):
        def nan_in_surface():
            if any([math.isnan(self.tp.species[sp]['surface_concentration']) for sp in self.tp.species]) or\
                math.isnan(self.tp.system['surface_pH']):
                return True
            else:
                return False
        #start with a single comsol calculation
        self.comsol.run(label=label)
        #check error
        error=self.comsol.check_error()
        #check for nan
        nan=nan_in_surface()
        was_nan=False
        #save the grid_factor and dflux before changing them in the loop
        grid_factor=float(self.tp.comsol_args['parameter']['grid_factor'][0])
        dramp=self.tp.comsol_args['solver_settings']['ramp']['dramp']
        dramp_step=dramp
        #if there were error's in the calculation or nan's in surface, run convergence loop
        nruns=0
        while nan or error:
            nruns+=1
            was_nan=True
            #rerun comsol with finer mesh
            if nan:
                self.tp.logger.warning('|    | CS | NaN appeared in surface_concentrations, rerunning COMSOL with slower ramping'+\
                    ' and finer resolution of grid')
            if error:
                self.tp.logger.warning('|    | CS | COMSOL encountered convergence problems')

            if nruns>25:
                self.tp.logger.error('| CI | -- | Restarting COMSOL with various settings and restarts did not help, stopping.')
                self.tp.logger.error('| CI | -- | Try to load the COMSOL file into the GUI and see how you can get convergence, maybe a denser grid '+
                        'can help, maybe a finer ramping of the load/non-linearity')
                sys.exit()
            elif nruns>20:
                self.tp.logger.warning('|    | CS | 5x restart did again not help, maximize all settings as a final try')
                self.tp.comsol_args['parameter']['grid_factor_bound'][0]=str(int(float(self.tp.comsol_args['parameter']['grid_factor_bound'][0])*1.5))
                self.tp.comsol_args['solver_settings']['ramp']['dramp']/=2.
                self.tp.comsol_args['par_values']='range(0,'+str(self.tp.comsol_args['solver_settings']['ramp']['dramp'])+',1)'
                self.tp.logger.info('|    | CS | Maximal discretization of boundary is now {}'.format(self.tp.comsol_args['parameter']['grid_factor_bound']))
                self.tp.logger.info('|    | CS | Load/Non-linearity ramping with an interval of {}'.format(self.tp.comsol_args['solver_settings']['ramp']['dramp']))
            elif nruns>15:
                self.tp.logger.warning('|    | CS | 5x restart with finer ramping did not help, finally try to increase both ramping and grid density')
                self.tp.comsol_args['parameter']['grid_factor_bound'][0]=str(int(float(self.tp.comsol_args['parameter']['grid_factor_bound'][0])*1.5))
                self.tp.logger.info('|    | CS | Maximal discretization of boundary is now {}'.format(self.tp.comsol_args['parameter']['grid_factor_bound']))
                self.tp.logger.info('|    | CS | Load/Non-linearity ramping with an interval of {}'.format(self.tp.comsol_args['solver_settings']['ramp']['dramp']))
            elif nruns>10:
                self.tp.logger.warning('|    | CS | 5x restart with finer grid did not help trying to reduce load/non-linearity ramping with initial grid')
                self.tp.comsol_args['solver_settings']['ramp']['dramp']/=2.
                self.tp.comsol_args['par_values']='range(0,'+str(self.tp.comsol_args['solver_settings']['ramp']['dramp'])+',1)'
                self.tp.comsol_args['parameter']['grid_factor_bound'][0]=grid_factor
                self.tp.logger.info('|    | CS | Load/Non-linearity ramping with an interval of {}'.format(self.tp.comsol_args['solver_settings']['ramp']['dramp']))
            elif nruns>5:
                self.tp.logger.warning('|    | CS | 5x restart did not help, trying to increase grid density')
                self.tp.comsol_args['parameter']['grid_factor_bound'][0]=str(int(float(self.tp.comsol_args['parameter']['grid_factor_bound'][0])*1.5))
                self.tp.logger.info('|    | CS | Maximal discretization of boundary is now {}'.format(self.tp.comsol_args['parameter']['grid_factor_bound']))

            #self.tp.comsol_args['dflux']/=2.
#           # self.tp.comsol_args['dflux']=int(self.tp.comsol_args['dflux'])
            #self.tp.logger.info('|    | CS | Flux descretization reduced to {}'.format(self.tp.comsol_args['dflux']))

            #if self.tp.comsol_args['dflux']/dflux_step == 2**3:
            #    dflux_step=self.tp.comsol_args['dflux']
            #    self.tp.logger.info('|    | CS | Discretization decrease of flux by factor of {} did not help, trying to decrease also minimal grid discretization'.format(2**3))
            #    #self.tp.comsol_args['parameter']['grid_factor_domain'][0]=str(int(float(self.tp.comsol_args['parameter']['grid_factor_domain'][0])*1.1))
            #    self.tp.comsol_args['parameter']['grid_factor_bound'][0]=str(int(float(self.tp.comsol_args['parameter']['grid_factor_bound'][0])*1.1))
            #    self.tp.logger.info('|    | CS | Maximal discretization of domain raised by {}x'.format(self.tp.comsol_args['parameter']['grid_factor_domain']))

            #if self.tp.comsol_args['dflux']<1e-5:
            #    self.tp.logger.error('|    | CS | Discretization of flux is smaller than 1e-5, stopping here, we will probably not get any convergence')
            #    sys.exit()

            #self.tp.comsol_args['par_values']='range(0,'+str(self.tp.comsol_args['dflux'])+',1)' #np.linspace(0,1,self.tp.comsol_args['nflux'])

#           # self.tp.logger.warning('|    | CS | Current x-axis resolution =  {} intervals.'.format(self.tp.comsol_args['nflux']))
            #self.tp.logger.warning('|    | CS | Current ramping with a flux decrement of  {}.'.format(self.tp.comsol_args['dflux']))

            self.comsol.run(label=label,restart=True)
            error=self.comsol.check_error()

            if nan:
                self.tp.logger.debug('Surface Concentrations:')
                for sp in self.tp.species:
                    self.tp.logger.debug('  {}: {} mol/L'.format(sp,self.tp.species[sp]['surface_concentration']/1000.))

            nan=nan_in_surface()

        #reset the grid_factor to original value:
        self.tp.comsol_args['solver_settings']['ramp']['dramp']=dramp
        self.tp.comsol_args['parameter']['grid_factor_domain'][0]=str(grid_factor)
