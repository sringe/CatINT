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

                    self.tp.logger.info('Starting calculation for {} = {} and {} = {} on CPU {} of {}'.format(desc_keys[0],value1,desc_keys[1],value2,self.tp.mpi_rank,self.tp.mpi_size))
                    label=str(i1).zfill(4)+'_'+str(i2).zfill(4)

                    if not self.tp.use_catmap:
                        self.run_single_step(label=label)
                    else:
                        self.run_scf_cycle(label=label)
            #synchronize all_data over CPUs
#            self.tp.all_data=reduce_dict_mpi(self.tp.all_data)
        else:
        #internal treatment of parameter
            if not self.tp.use_catmap:
                self.run_single_step()
            else:
                self.tp.logger.error('Catmap has been selected but current descriptor method is Comsol internal which does not work. Stopping here.')
                sys.exit()

        if use_mpi:
            self.tp.comm.Barrier()

    def evaluate_accuracy(self,par,par_old):
        rmsd=0
        n=0
        par_sum=sum(par)
        self.tp.logger.debug(' | Current current density = {} mV/cm^2'.format(par_sum))
        par_old_sum=sum(par_old)
        rmsd=np.sqrt((par_sum-par_old_sum)**2)
#        for val1,val2 in zip(par,par_old):
#            n+=1
#            rmsd+=(val1-val2)**2
#        rmsd=np.sqrt(rmsd/(n*1.))
        return rmsd

    def converged(self,old,new):
        cmrsd=0.0
        for isp in range(self.tp.nspecies):
            crmsd+=(old[isp*self.tp.nx]-new[isp*self.tp.nx])**2 #RMSD
        if np.sqrt(sum(crmsd)/len(crmsd))<self.tp.tau_scf:
            return True
        else:
            return False

    def run_single_step(self,label=''):
        """run a single transport calculation"""
#        if self.calc != 'comsol':
#            cout=self.integrate_pnp(self.tp.dx,self.tp.nx,self.tp.dt,\
#                len(self.tp.tmesh),self.tp.ntout,method=self.calc)
#            for i_sp,sp in enumerate(self.tp.species):
#                self.tp.species[sp]['concentration']=cout[-1,i_sp*self.tp.nx:(i_sp+1)*self.tp.nx]
#            self.tp.cout=cout
#        else:
        self.comsol.run(label=label)

    def run_scf_cycle(self,label=''):
        """run scf cycle to converge catmap-comsol iterations"""
        istep=0
        step_to_check=0
        scf_accuracy=np.inf
        self.tp.logger.info('Starting iterative solution with CatMAP and COMSOL')
        self.tp.logger.info('  using a current density accuracy cutoff of {} mV/cm^2 and a linear mixing parameter of {}'.format(self.tau_scf,self.mix_scf))
        accuracies=[]
        while scf_accuracy>self.tau_scf:
            istep+=1
            #evaluate how the accuracy during the last 50 steps, if it does not significantly decrease, reduce the mixing factor
            if istep-step_to_check>80 and abs(accuracy[-2]-accuracy[-1])>1e-1:
                #still no convergence, try to decrease mixing factor
                self.scf_mix*=0.9
                self.tp.logger.info(' | Accuracy is still < 1e-1, decreasing mixing factor in order to speed up the convergence')
                step_to_check=istep
            self.tp.logger.info(' | Solving transport step {}. Current accuracy in current density = {} mV/cm^2'.format(istep,scf_accuracy))

            #linear mixing
            if istep>2:
                #mix fluxes
                for sp in self.tp.species:
#                    self.tp.species[sp]['flux']=self.mix_scf*self.tp.species[sp]['flux']+(1.-self.mix_scf)*fl_old[sp]
                    self.tp.species[sp]['surface concentration']=self.mix_scf*self.tp.species[sp]['surface concentration']+\
                            (1.-self.mix_scf)*sc_old[sp]

            #for linear mixing, safe all surface concentrations
            sc_old={}
            for sp in self.tp.species:
                sc_old[sp]=self.tp.species[sp]['surface concentration']

            #1) Microkinetic Model: CatMAP
            if istep==1:
                self.tp.logger.debug('Electrode Fluxes:')
                for sp in self.tp.species:
                    self.tp.logger.debug('  {}: {}'.format(sp,self.tp.species[sp]['flux']))

            #run catmap. the fluxes will be updated automatically
            self.catmap.run()
            self.tp.logger.debug('Electrode Fluxes:')
            for sp in self.tp.species:
                self.tp.logger.debug('  {}: {}'.format(sp,self.tp.species[sp]['flux']))

            #2) Transport: COMSOL

            self.tp.logger.debug('Surface Concentrations:')
            for sp in self.tp.species:
                self.tp.logger.debug('  {}: {} mol/L'.format(sp,self.tp.species[sp]['surface concentration']/1000.))

            #only_last updates the descriptor based dictionaries only for the last entry in self.tp.descriptors
            self.run_single_step(label=label)

            was_nan=False
            grid_factor=float(self.tp.comsol_args['parameter']['grid_factor'][0])
            nflux=self.tp.comsol_args['nflux']
            nflux_step=0
            while any([math.isnan(self.tp.species[sp]['surface concentration']) for sp in self.tp.species]):
                was_nan=True
                #rerun comsol with finer mesh
                self.tp.logger.warning(' | CS | NaN appeared in surface concentrations, rerunning COMSOL with slower ramping'+\
                        ' and finer resolution of grid')
                self.tp.comsol_args['nflux']*=1.25
                self.tp.comsol_args['nflux']=int(self.tp.comsol_args['nflux'])
                self.tp.logger.info(' | CS | Number of flux ramping steps increased to {}'.format(self.tp.comsol_args['nflux']))

                if self.tp.comsol_args['nflux']/1000 != nflux_step:
                    nflux_step=self.tp.comsol_args['nflux']/1000
                    self.tp.logger.info(' | CS | Ramping increase of 1000 did not help, trying to decrease also minimal grid discretization')
                    self.tp.comsol_args['parameter']['grid_factor_domain'][0]=str(int(float(self.tp.comsol_args['parameter']['grid_factor_domain'][0])*1.1))
                    self.tp.logger.info(' | CS | Maximal discretization of domain raised by {}x'.format(self.tp.comsol_args['parameter']['grid_factor_domain']))

                if self.tp.comsol_args['nflux']>20000:
                    self.tp.logger.error(' | CS | ramping # nflux is larger than 20000, stopping here, we will probably not get any convergence')
                    sys.exit()

                self.tp.comsol_args['par_values']=np.linspace(0,1,self.tp.comsol_args['nflux'])

#                self.tp.logger.warning(' | CS | Current x-axis resolution =  {} intervals.'.format(self.tp.comsol_args['nflux']))
                self.tp.logger.warning(' | CS | Current ramping of flux =  {} intervals.'.format(self.tp.comsol_args['nflux']))

                self.run_single_step(label=label)

                self.tp.logger.debug('Surface Concentrations:')
                for sp in self.tp.species:
                    self.tp.logger.debug('  {}: {} mol/L'.format(sp,self.tp.species[sp]['surface concentration']/1000.))

            if was_nan:
                #reset the grid_factor to original value:
                self.tp.comsol_args['nflux']=nflux
                self.tp.comsol_args['parameter']['grid_factor_domain'][0]=str(grid_factor)

            #3) Check Convergence
            current_density=[]
            for sp in self.tp.electrode_reactions:
#                current_density.append(self.tp.electrode_reactions[sp]['electrode_current_density'][(str(1.0),str(value2))])
                current_density.append(self.tp.species[sp]['electrode_current_density'])
            if istep>1:
                scf_accuracy=self.evaluate_accuracy(current_density,old_current_density)
            old_current_density=deepcopy(current_density)
            accuracies.append(scf_accuracy)
