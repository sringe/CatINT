import re
from units import *
import numpy as np
import os
from subprocess import call
from io import replace_line,insert_line,mpi_make_dir
from glob import glob
import sys
sys.path.insert(0,'/scratch/users/sringe/software/catmap')
import pickle
from catmap.model import ReactionModel
from catmap import analyze
from string import Template
from shutil import copy
from tabulate import tabulate
def flatten_2d(output):
    "Helper function for flattening rate_control output"
    flat = []
    for x in output:
        flat+= x
    return flat

class Object(object):
    pass

class CatMAP():
    """This class modifies the concentrations in the catmap input file (defined by user) and runs catmap. The output is the read in and the boundary conditions adjusted properly"""

    def __init__(self,path=os.getcwd(),transport=None,model_name=None):
        #delta_desc is the descriptor delta with which the descriptor axes is resolved
        #in case of potential 0.05 is fine
        #max and min_desc are the  bounds of the descriptor region (here default for potential being the descriptor)



        self.use_interactions=False

        if transport is None:
            self.tp.logger.error('|    | CM | No transport object provided for catmap_wrapper. Stopping here.')
            sys.exit()
        else:
            self.tp=transport
        if 'n_inter' not in self.tp.catmap_args:
            self.n_inter=10 #number of steps to converge interactions
        else:
            self.n_inter=self.tp.catmap_args['n_inter']
        self.tp.path=path
        self.output_folder='results_catmap'
        if model_name is not None:
            self.model_name='catmap_'+model_name
        else:
            self.model_name='catmap_'+self.tp.model_name
        model_name=self.model_name
        self.catmap_model=model_name+'.mkm'
        self.output_base_folder=self.tp.outputfoldername+'/catmap_output'
        self.input_base_folder=self.tp.outputfoldername+'/catmap_input'
        #create output folder (only rank==0 can do this, others have to wait)
        mpi_make_dir(self.output_base_folder)
        #create input/running folder (only rank==0 can do this, others have to wait)
        mpi_make_dir(self.input_base_folder)

        if 'n_inter_max' in self.tp.catmap_args:
            self.n_inter_max=self.tp.catmap_args['n_inter_max']
        else:
            self.n_inter_max=150

        if 'desc_method' not in self.tp.catmap_args:
            self.tp.catmap_args['desc_method']='automatic'
            #possible settings:
            # - automatic
            # - single_point
            # - descriptor_range
            # - from_input

        #AUXILIARY DESCRIPTOR SWEEP

        #minimum descriptor value in auxiliary descriptor sweep
        if 'min_desc' in self.tp.catmap_args:
            self.min_desc=self.tp.catmap_args['min_desc']
        else:
            self.min_desc=-1.5
        #maximum descriptor value in auxiliary descriptor sweep
        if 'max_desc' in self.tp.catmap_args:
            self.max_desc=self.tp.catmap_args['max_desc']
        else:
            self.max_desc=0.2
        #maximum descriptor value in auxiliary descriptor sweep
        #relative to current descriptor
        if 'max_desc_delta' in self.tp.catmap_args:
            self.max_desc_delta=self.tp.catmap_args['max_desc_delta']
        else:
            self.max_desc_delta=None
        #minimum descriptor value in auxiliary descriptor sweep
        #relative to current descriptor
        if 'min_desc_delta' in self.tp.catmap_args:
            self.min_desc_delta=self.tp.catmap_args['min_desc_delta']
        else:
            self.min_desc_delta=None
        #discretization of auxiliary descriptor sweep
        if 'delta_desc' in self.tp.catmap_args:
            self.delta_desc=self.tp.catmap_args['delta_desc']
        else:
            self.delta_desc=0.05


    def run(self):
        desc_keys=[key for key in self.tp.descriptors]
        desc_val=[self.tp.system[key] for key in self.tp.descriptors]
        self.tp.logger.info('|    | CM | Starting CatMAP for {} = {}'.format(desc_keys[0],desc_val[0]))

        self.input_folder=self.input_base_folder+'/desc_'+str(desc_val[0])
        self.output_folder=self.output_base_folder+'/desc_'+str(desc_val[0])

        if not os.path.isdir(self.input_folder):
            os.makedirs(self.input_folder)
        if not os.path.isdir(self.output_folder):
            os.makedirs(self.output_folder)

        #go to output folder and run catmap there
        root=os.getcwd()
        os.chdir(self.input_folder)

        #setting up the mkm files
        mkm_template_file=root+'/'+self.model_name+'_template.mkm'
        if not os.path.exists(mkm_template_file):
            self.tp.logger.error('|    | CM | mkm file {} required for CatMAP does not exist.'.format(mkm_template_file))
            sys.exit()
        else:
            copy(mkm_template_file,self.model_name+'_template.mkm')
        copy(self.model_name+'_template.mkm',self.model_name+'.mkm')
#        mkm_template = Template(open(mkm_template_file).read())
#        mkm_text = mkm_template.substitute(pH_new = '7')
        mkm_file = self.model_name+'.mkm'
#        with open(mkm_file,'w') as f:
#            f.write(mkm_text)

        #update the mkm file with current concentrations
        self.update_input(desc_val)

        #energies file
        energies_file=root+'/'+self.model_name+'_energies.txt'
        if not os.path.exists(energies_file):
            self.tp.logger.error('|    | CM | energy data file {} required by CatMAP does not exist.'.format(energies_file))
            sys.exit()
        else:
            copy(energies_file,self.model_name+'_energies.txt')

        #setting up the model
        model = ReactionModel(setup_file = mkm_file, max_log_line_length=0) #we set the maximum line length to 0, because we work
        #with pickle files anyhow
        #output
        model.output_variables+=['consumption_rate','production_rate', 'free_energy', 'selectivity', 'interacting_energy','turnover_frequency','rate_control']
        def plot_fed(pressure_corr=False,coverage_corr=False,method=0):
            ma = analyze.MechanismAnalysis(model)
            ma.surface_colors = ['k','b','r','yellow','green','orange','cyan']
            ma.label_args['size'] = 14
            ma.energy_type = 'free_energy' #can also be free_energy/potential_energy
            ma.include_labels = True #way too messy with labels
            ma.pressure_correction = pressure_corr #assume all pressures are 1 bar (so that energies are the same as from DFT)
            ma.coverage_correction = coverage_corr #False
            min_a=np.inf
            desc_cm=None
            for a,b in model.descriptor_ranges:
                if min_a>abs(a-float(desc_val[0])):
                    min_a=abs(a-float(desc_val[0]))
                    desc_cm=a
            if min_a>1e-5 or desc_cm is None:
                self.tp.logger.warning('| CI | CM | Somehow, there was an error searching for the descriptor in the catmap'+\
                        'descriptor list, skipping the FED plot to not harm simulation')
                return
            if self.use_interactions:
                ma.energy_type = 'interacting_energy'
            if not pressure_corr and not coverage_corr:
                fig = ma.plot(save='FED.pdf',plot_variants=[desc_cm],method=method)
            elif pressure_corr and not coverage_corr:
                fig = ma.plot(save='FED_pressure_corrected.pdf',plot_variants=[desc_cm],method=method)
            elif pressure_corr and coverage_corr:
                fig = ma.plot(save='FED_pressure_and_cov_corrected.pdf',plot_variants=[desc_cm],method=method)
            elif not pressure_corr and coverage_corr:
                fig = ma.plot(save='FED_cov_corrected.pdf',plot_variants=[desc_cm],method=method)
        #plot_fed(True)
#        if not self.use_interactions:
#            plot_fed(False)
#        #    pass
#        else:
#            pass
        #sys.exit()
        #slowly ramp up the interactions if desired
        if type(self.n_inter)==int:
            if self.use_interactions:
                if self.n_inter==1:
                    inter=[self.interaction_strength]
                else:
                    inter=np.linspace(0,self.interaction_strength,self.n_inter)
            else:
                inter=[0]
        elif 'n_inter_min' in self.tp.catmap_args:
            if self.use_interactions:
                inter=np.linspace(0,self.interaction_strength,self.tp.catmap_args['n_inter_min'])
            else:
                inter=[0]
        else:
            if self.use_interactions:
                inter=[self.interaction_strength]
            else:
                inter=[0]
#        #redirect stdout to file
#        old = os.dup(1)
#        sys.stdout.flush()
#        os.close(1)
#        print os.getcwd()
#        os.open(self.input_folder+"/catmap_run.err", os.O_CREAT)
        desc_method=None
        jj=0
        while True:
            #for ii in inter:
            #    i=0
            #    if self.use_interactions:
            #        for line in open(self.catmap_model):
            #            i+=1
            #            if 'interaction_strength' in line:
            #                replace_line(mkm_file,i-1,'interaction_strength = '+str(ii))
            #        self.tp.logger.info('|    | CM | Running interaction_strength = {}'.format(ii))
            #    model = ReactionModel(setup_file = mkm_file, max_log_line_length=0)
            #    model.output_variables+=['consumption_rate','production_rate', 'free_energy', 'selectivity', 'interacting_energy','turnover_frequency']
            #    model.run()
            try:
                for ii in inter:
                    i=0
                    if self.use_interactions:
                        for line in open(self.catmap_model):
                            i+=1
                            if 'interaction_strength' in line:
                                replace_line(mkm_file,i-1,'interaction_strength = '+str(ii))
                        self.tp.logger.info('|    | CM | Running interaction_strength = {}'.format(ii))
                    model = ReactionModel(setup_file = mkm_file, max_log_line_length=0)
                    model.output_variables+=['consumption_rate','production_rate', 'free_energy', 'selectivity', 'interacting_energy','turnover_frequency','rate_control']
                    model.run()
            except:
                #SOME CRASH OF CATMAP
                # - check if due to high interaction strength or missing data points
                self.tp.logger.warning('|    | CM | CatMAP did not converge with interaction strength = {}'.format(ii))

                if desc_method is not None:
                    self.tp.logger.error('|    | CM | Running on discrete descriptor space did not help, CatMAP did not converge. Try to increase descriptor space manually.')

                desc_method=None
                
                if not self.use_interactions or (self.use_interactions and ii==0):
                    self.tp.logger.error('|    | CM | Unexpected end of CatMAP even though interactions_strength = 0, check error files for hints.')

                    if self.tp.catmap_args['desc_method']=='automatic':
                        self.tp.logger.info('|    | CM | Since desc_method = automatic, CatINT tries to solve this problem by')
                        self.tp.logger.info('|    | CM | running CatMAP on discrete descriptor space instead of single point.')
                        self.tp.logger.info('|    | CM | This will usually take much more time.')
                        desc_method='descriptor_range'

                        ######UPDATE INPUT########
                        #update the mkm file with current concentrations
                        self.update_input(desc_val,desc_method=desc_method)

                        #energies file
                        energies_file=root+'/'+self.model_name+'_energies.txt'
                        if not os.path.exists(energies_file):
                            self.tp.logger.error('|    | CM | energy data file {} required by CatMAP does not exist.'.format(energies_file))
                            sys.exit()
                        else:
                            copy(energies_file,self.model_name+'_energies.txt')

                        #setting up the model
                        model = ReactionModel(setup_file = mkm_file, max_log_line_length=0) #we set the maximum line length to 0, because we work
                        #with pickle files anyhow
                        #output
                        model.output_variables+=['consumption_rate','production_rate', 'free_energy', 'selectivity', 'interacting_energy','turnover_frequency']
                        ##########################
                    else:
                        if self.tp.catmap_args['desc_method']=='descriptor_range':
                            self.tp.logger.error('|    | CM | Descriptor method was already descriptor_range, but no convergence was achieved, try to increase max_bisections')
                            sys.exit()
                        elif self.tp.catmap_args['desc_method']=='single_point':
                            self.tp.logger.error('|    | CM | Try to rerun using desc_method = automatic settings.')
                            sys.exit()
                else:
                    if self.n_inter=='automatic':
                        jj+=1
                        if 10*jj>self.n_inter_max:
                            self.tp.logger.error('|    | CM | Adjusted interaction ramping range is larger than n_inter_max. Adjust n_inter_max in order to run finer range')
                            sys.exit()
                        self.tp.logger.warning('|    | CM | Adjusting interaction ramping to range 0 to 1 with {} steps'.format(10*jj))
                        inter=np.linspace(0,self.interaction_strength,10*jj)
                        pass
                    else:
                        break
            else:
                break

#        rate_constants = model.solver.get_rate_constants(rxn_parameters,coverages)
#        kfs, krs, dkfs, dkrs = model.rate_constants(rxn_parameters,coverages,
#        model._gas_energies,model._site_energies,
#        model.temperature,model.interaction_response_function,
#        model._mpfloat,model._matrix,model._math.exp)
        if self.use_interactions:
#            self.tp.logger.info(desc_val)
#            self.tp.logger.info(float(desc_val[0]))
            try:
                idx = [i for i in range(len(model.interacting_energy_map)) if abs(model.interacting_energy_map[i][0][0]-float(desc_val[0]))<1e-5][0]
                descrip, coverages = model.coverage_map[idx]
                rxn_parameters = model.scaler.get_rxn_parameters(descrip)
                self.tp.logger.info('|    | CM | --- Interaction energies ---')
                self.tp.logger.info('|    | CM |  - desc = {}'.format(desc_val[0]))
#                self.tp.logger.info(model.output_labels['interacting_energy'])
#                self.tp.logger.info(model.solver.get_interacting_energies(rxn_parameters))
                all_ads = model.adsorbate_names + model.transition_state_names
                N_ads = len(all_ads)
                energies = rxn_parameters[:N_ads]
                eps_vector = rxn_parameters[N_ads:]
                cvg = coverages + [0]*len(model.transition_state_names)
                #self.tp.logger.info('-- checking --')
                #self.tp.logger.info(model.interaction_function(cvg,energies,eps_vector,model.thermodynamics.adsorbate_interactions.interaction_response_function,False,False))
                #self.tp.logger.info('-- end checking --')
                self.tp.logger.info(tabulate(zip(model.output_labels['interacting_energy'],\
                        energies),\
                        headers=['species','energy']))
                #self.tp.logger.info(tabulate(zip(model.output_labels['interacting_energy'],\
                #        model.interaction_function(cvg,energies,eps_vector,\
                #        model.thermodynamics.adsorbate_interactions.interaction_response_function,False,False)[1]),\
                #        headers=['species','energy']))
                self.tp.logger.info('|    | CM | --- END OUTPUT ---')
                sys.exit()
            except:
                self.tp.logger.info('|    | CM | Error in reading interaction energies')
                pass
        try:
            plot_fed(False,False,method=2)
        except UserWarning:
            self.tp.logger.warning('|    | CM | Error in writing FED with no corrections')
        try:
            plot_fed(True,False,method=2)
        except UserWarning:
            self.tp.logger.warning('|    | CM | Error in writing FED with pressure correction')
        try:
            plot_fed(True,True,method=2)
        except UserWarning:
            self.tp.logger.warning('|    | CM | Error in writing FED with all corrections')
        try:
            plot_fed(False,True,method=2)
        except UserWarning:
            self.tp.logger.warning('|    | CM | Error in writing FED with coverage correction')
#        sys.stdout.flush()
#        os.close(1)
#        os.dup(old) # should dup to 1
#        os.close(old) # get rid of left overs

        #run!
#        stdout = sys.stdout
#        sys.stdout = open('std.log', 'w')
#        sys.stdout = stdout

        converged=False
        count=0
        max_iter=0
        for line in open(self.model_name+'.log','r'):
            if 'mapper_iteration_' in line:
                sol=int(re.findall('mapper_iteration_(\d+)',line)[0])
                max_iter=sol
            if 'status - 0 points do not have valid solution.' in line:
                converged=True
#        for line in open(self.model_name+'.mkm','r'):
#            if 'resolution' in line:
#                sol=re.findall('resolution[ ]*=[ ]*\[[ ]*(\d+)[ ]*,[\d\s]+\]',line)
#                number_of_jobs=int(sol[0])
#                break
        #go back to root folder
        os.chdir(root)

        if not converged:
            self.tp.logger.error('|    | CM | CatMAP did not converge')
        else:
            self.tp.logger.info('|    | CM | CatMAP finished successfully in {} steps'.format(max_iter))

        self.read_output(desc_val)
        #species_definitions['CO2_g'] = {'concentration':0.2}

    def update_input(self,desc_val,desc_method=None):

        #go over input file and check if single point or full descriptor range is used as method:
        found_resolution=False
        found_descriptor_range=False
        found_descriptors=False

        desc_method_in=None
        i=0
        for line in open(self.catmap_model):
            i+=1
            if line.lstrip().startswith('#'):
                continue
            if 'descriptor_range' in line:
                found_descriptor_range=True
            if 'resolution' in line:
                found_resolution=True
            if 'descriptors' in line:
                desc_method_in='single_point'
                found_descriptors=True
                idesc=i

        if desc_method is None:
            if self.tp.catmap_args['desc_method'] in ['single_point','descriptor_range']:
                desc_method=self.tp.catmap_args['desc_method']
            elif self.tp.catmap_args['desc_method'] == 'automatic':
                desc_method='single_point'
            elif self.tp.catmap_args['desc_method'] == 'from_input':
                if desc_method_in is not None:
                    desc_method=desc_method_in
                else:
                    desc_method='descriptor_range'
                self.tp.catmap_args['desc_method']=desc_method
                self.tp.logger.info('|    | CM | Descriptor method was determined from CatMAP input file to be {}'.format(desc_method))

        self.tp.logger.info('|    | CM | Running CatMAP using descriptor method = {}'.format(desc_method))

#        if self.method is None:
#            self.method='descriptor_range'
#        
#        #overwrite method by input parameter (preference)
#        if method is not None:
#            self.method=method

        if not found_resolution:
            insert_line(self.catmap_model,idesc,'resolution = [1,1]\n')
            found_resolution=True
        if desc_method=='single_point':
            if not found_descriptor_range:
                insert_line(self.catmap_model,idesc,'descriptor_range = [[0,0],[0,0]]')
                found_descriptor_range=True
        if not all([found_descriptor_range,found_resolution]):
                #or found_descriptors):
            self.tp.logger.error('|    | CM | Missing resolution definition or descriptor range')
            sys.exit()

        def convert(name_catint):
            if name_catint=='phiM':
                name_catmap='voltage'
            else:
                name_catmap=name_catint
            return name_catmap
        #1) work on input files
        #we need to replace the concentration lines in the file
        desc_keys=[key for key in self.tp.descriptors]
        desc_1=convert(desc_keys[0])
        desc_2=convert(desc_keys[1])

        if desc_method=='descriptor_range':
            desc_list=self.tp.descriptors[desc_keys[0]]
            if len(desc_list)<5 and desc_1 != 'voltage':
                self.tp.logger.warning('|    | CM | It could be that no pkl files are written out because the given descriptor axis is too coarse, consider a denser axis')

            if desc_1=='voltage':
                #use manual range here that contains current descriptor value
                #and starts
                val=desc_val[0]
                desc_list=[val]
                if self.max_desc_delta is not None:
                    max_desc=desc_val[0]+self.max_desc_delta
                else:
                    max_desc=max(self.max_desc,max(self.tp.descriptors['phiM']))
                while val<max_desc:
                    val+=self.delta_desc
                    desc_list.append(val)
                val=desc_val[0]
                if self.min_desc_delta is not None:
                    min_desc=desc_val[0]-self.min_desc_delta
                else:
                    min_desc=min(self.min_desc,min(self.tp.descriptors['phiM']))
                while val>min_desc:
                    val-=self.delta_desc
                    desc_list.append(val)
            desc_list=sorted(desc_list,key=float)

            min_desc=min(desc_list)
            max_desc=max(desc_list)
            n_desc=len(desc_list)
            delta_desc=desc_list[1]-desc_list[0]
            if delta_desc<1e-9:
                self.tp.logger.warning('|    | CM | The descriptor separation is smaller then 1e-9. The index assignment of the catmap data will fail, so I quit here')
                sys.exit()

        #check if interactions are desired, if yes read out the desired interaction strength
        if any([line.strip().startswith('adsorbate_interaction_model') for line in open(self.catmap_model)]): # or \
                #any(['max_coverage' in line for line in open(self.catmap_model)]):
            self.use_interactions=True

        if self.use_interactions:
            for line in open(self.catmap_model):
                sol=re.findall('interaction_strength[ ]*=[ ]*([0-9.]+)',line)
                if len(sol)>0:
                    self.interaction_strength=float(sol[0])

        i=0
        replaced_species=[]
        self.tp.logger.debug('|    | CM | in | Diffusion drop passed to CatMAP is {} V'.format(self.tp.system['potential'][0]))
        for sp in self.tp.species:
            self.tp.logger.info('|    | CM | in | c(x=0) of {} = {} M'.format(sp,self.tp.species[sp]['surface_concentration']/1000.))
        for line in open(self.catmap_model):
            i+=1
            if line.lstrip().startswith('#'):
                continue
            if line.strip().startswith('input_file'):
                replace_line(self.catmap_model,i-1,'input_file = \''+self.model_name+'_energies.txt\'')
                continue
            if line.strip().startswith('data_file'):
                replace_line(self.catmap_model,i-1,'data_file = \''+self.model_name+'.pkl\'')
                continue
            if line.strip().startswith('bulk_ph'):
                replace_line(self.catmap_model,i-1,'bulk_ph = '+str(self.tp.system['bulk_pH']))
                continue
            if line.strip().startswith('Upzc'):
                replace_line(self.catmap_model,i-1,'Upzc = '+str(self.tp.system['phiPZC']))
                continue
            if line.strip().startswith('field'):
                #replace_line(self.catmap_model,i-1,'field = '+str(self.tp.system['efield'][0]*1e-10*self.tp.system['epsilon']/self.tp.system['Stern epsilon']))
                replace_line(self.catmap_model,i-1,'field = '+str(self.tp.system['Stern_efield']*1e-10))
            if self.tp.system['charging_scheme'] == 'comsol' and line.strip().startswith('sigma'):
                #replace_line(self.catmap_model,i-1,'field = '+str(self.tp.system['efield'][0]*1e-10*self.tp.system['epsilon']/self.tp.system['Stern epsilon']))
                replace_line(self.catmap_model,i-1,'sigma_input = '+str(((self.tp.system['phiM']-self.tp.system['phiPZC'])-self.tp.system['surface_potential'])*self.tp.system['Stern capacitance']))
            if line.strip().startswith('voltage_diff_drop') and self.tp.system['potential drop']=='Stern':
                self.tp.logger.info('|    | CM | in | Running catmap with voltage drop = {}'.format(self.tp.system['phiM']-self.tp.system['potential'][0]))
                replace_line(self.catmap_model,i-1,'voltage_diff_drop = '+str(self.tp.system['potential'][0])) #potential drop']))
            for sp in self.tp.species:
                sp_cm=self.species_to_catmap(sp)
                sol=re.findall('species_definitions\[\''+sp_cm+'_g\'\].*{\'pressure\':.*}',line)
                if len(sol)>0:
                    #shortly check if concentrations are very negative, stop if they are:
                    if self.tp.species[sp]['surface_concentration']<-1.:
                        self.tp.logger.warning('|    | CM | Surface concentration of {} is more negative than 1e-3 mol/L, stopping to be safe.'.format(sp))
                        sys.exit()
                    if sp in ['OH-','H+']:
                        activity=1 #this is irrelevant, since it will be controlled by the pH
                    elif not self.tp.system['use_activities']:
                        activity=self.tp.species[sp]['surface_concentration']/1000.
                    else:
                        #activity=self.tp.species[sp]['surface_concentration']/self.tp.system['reference_gas_concentration'] #/(self.tp.species[sp]['Henry constant']*self.tp.system['pressure']) #self.tp.species[sp]['bulk_concentration']
                        activity=self.tp.species[sp]['surface_concentration']/self.tp.species[sp]['Henry constant'] #*self.tp.system['pressure']) #self.tp.species[sp]['bulk_concentration']
                        #need to multiply with activity coefficient
                        activity*=self.tp.species[sp]['surface_activity_coefficient']
                    self.tp.logger.debug('|    | CM | a_{}(x=0) = {}'.format(sp,activity))
                    replace_line(self.catmap_model,i-1,"species_definitions['"+sp_cm+"_g'] = {'pressure':"+str(max(0.,activity))+"}")
                    replaced_species.append(sp)
            #set pressure for water
            #originally, this should be 0.035 bar, because this corresponds to the partial pressure of water in the gas phase
            #   at equilibrium, we have mu_l=mu_g=mu_g^0+RTln(P/RT), where P is the partial pressure of water. since the reference pressure
            #   is 1 bar, we have delta mu_l = delta mu_g^0 + RTln(P/Pref), with Pref=1bar, so effectively we just need to multiply all rates
            #   with the vapor pressure of water
            #problem is only, there is a bug in catmap, in the process when setting the OH- energy equal to H2O + H+ where the H2O energy
            #   does not include the pressure change! this means we are better off directly correcting the energy of water and set the activity
            #   to 1 here
            sp_cm='H2O'
            sol=re.findall('species_definitions\[\''+sp_cm+'_g\'\].*{\'pressure\':.*}',line)
            if len(sol)>0:
                replace_line(self.catmap_model,i-1,"species_definitions['"+sp_cm+"_g'] = {'pressure':1.0}")
            if desc_method=='descriptor_range':
                if 'descriptor_range' in line:
                    replace_line(self.catmap_model,i-1,'descriptor_ranges = [['+str(min_desc)+','+str(max_desc)+'],['+str(desc_val[1])+','+str(desc_val[1])+']]')
                if line.startswith('descriptors'):
                    #this line should not exist if descriptor_range is used
                    replace_line(self.catmap_model,i-1,'')
                if line.strip().startswith('resolution'):
                    replace_line(self.catmap_model,i-1,'resolution = ['+str(n_desc)+', 1]') #descriptor_names= [\''+desc_1+'\', \''+desc_2+'\']')
            elif desc_method=='single_point':
                if 'descriptor_range' in line:
                    replace_line(self.catmap_model,i-1,'descriptor_ranges = [['+str(desc_val[0])+','+str(desc_val[0])+'],['+str(desc_val[1])+','+str(desc_val[1])+']]')
                if line.strip().startswith('resolution'):
                    replace_line(self.catmap_model,i-1,'resolution = [1,1]') #descriptor_names= [\''+desc_1+'\', \''+desc_2+'\']')
                if 'descriptors' in line:
                    replace_line(self.catmap_model,i-1,'descriptors = ['+str(desc_val[0])+','+str(desc_val[1])+']')
            if 'descriptor_names' in line:
                replace_line(self.catmap_model,i-1,'descriptor_names= [\''+desc_1+'\', \''+desc_2+'\']')
            sol=re.findall('pH[ ]*=[ ]*\d',line)
            if len(sol)>0:
                replace_line(self.catmap_model,i-1,'pH = '+str(self.tp.system['surface_pH'])+'')
                self.tp.logger.debug('|    | CM | surface pH = {}'.format(self.tp.system['surface_pH']))
#                replace_line(self.catmap_model,i-1,'pH = '+str(max(self.tp.system['pH']))+'')
        for sp in self.tp.species:
            if sp not in replaced_species and sp not in self.tp.system['exclude species'] and sp not in self.tp.electrolyte_list:
                self.tp.logger.warning('|    | CM | Pressure of species {} not in catmap mkm file, give an arbitrary pressure of all species needed which will then be replaced by CatINT'.format(sp))
                sys.exit()
    #SETTINGS
    def convert_TOF(self,A): # Given a list, convert all the TOF to j(mA/cm2) using 0.161*TOF(According to Heine's ORR paper)
        B = [-0.161*rate for rate in A]
        return B

    def read_output(self,desc_val):
        root=os.getcwd()
        os.chdir(self.input_folder)
        log_file = self.model_name+'.log'
        model = ReactionModel(setup_file = log_file)
        pickle_file = self.model_name+'.pkl'
        data = self.get_data(pickle_file,model)
        self.tp.catmap_data=data
        os.chdir(root)
        for sp in self.tp.species:
            self.tp.species[sp]['flux']=0.0
        #get reaction mechanism
        mechanisms = model.rxn_mechanisms.values()
        mechanisms_names = model.rxn_mechanisms.keys()
        for sp in self.tp.species: #prod in self.tp.electrode_reactions:
            ###############
            #1) the TOF's save them as fluxes for the individual species
            ###############
            tof=None
            name=self.species_to_catmap(sp)
            name+='_g'
            idx=None

            #TOF
            #"tof" is the signed rate of conversion/active site/s
            if name in data.turnover_frequency_names:
                idx=data.turnover_frequency_names.index(name)
                tof=data.turnover_frequency[:,idx]
            elif sp not in self.tp.system['exclude species'] and sp not in self.tp.electrolyte_list:
                self.tp.logger.warning('|    | CM | No CatMAP TOF data was found for species {}. Check your species definition names. CatINT uses the CatINT equation names to map to the CatMAP names!'.format(sp))
            if tof is None:
                continue
            data_ref=np.column_stack((data.voltage, tof))
            voltages=data_ref[np.argsort(data_ref[:, 0])][:,0]
            if sp in self.tp.electrode_reactions:
                nprod=len([a for a in self.tp.electrode_reactions[sp]['reaction'][1] if a==sp])
                nel=self.tp.electrode_reactions[sp]['nel']
            else:
                nprod=1
                nel=1

            rates=data_ref[np.argsort(data_ref[:, 0])][:,1]*self.tp.system['active site density'] #rate in [mol/m^2/s]
            current_densities=rates*nel*unit_F/nprod #current density in [A/m^2]
            current_densities/=10. #convert A/m^2 to mA/cm^2
#            currents=self.convert_TOF(data_ref[np.argsort(data_ref[:, 0])][:,1])
            pol_file=self.output_folder+'/j_'+name.split('_')[0]+'.tsv'
            np.savetxt(pol_file, np.array([voltages,current_densities]).T)
            iv=-1
            for v in voltages:
                iv+=1
                if abs(v-float(desc_val[0]))<1e-10:
                    index=iv
                    break
            self.tp.species[sp]['flux']=rates[index]
            self.tp.species[sp]['electrode_current_density']=current_densities[index]

            #update also alldata array:
            alldata_inx=self.tp.alldata_names.index([desc_val[0],desc_val[1]])
            self.tp.alldata[alldata_inx]['species'][sp]['electrode_current_density']=self.tp.species[sp]['electrode_current_density']
            self.tp.alldata[alldata_inx]['species'][sp]['electrode_flux']=self.tp.species[sp]['flux']

            ###############
            #2) rate control
            ###############
            if name in data.rate_control_names[0]:
                idx=data.rate_control_names[0].index(name)
                rate_control=data.rate_control[:,idx]
            elif sp not in self.tp.system['exclude species'] and sp not in self.tp.electrolyte_list:
                self.tp.logger.warning('|    | CM | No CatMAP rate control data was found for species {}. Check your species definition names. CatINT uses the CatINT equation names to map to the CatMAP names!'.format(sp))
            #loop over influencing states
            for name2 in data.rate_control_names[1]:
                data_rc=np.column_stack((data.voltage,rate_control[:,data.rate_control_names[1].index(name2)]))
                data_rc=data_rc[np.argsort(data_rc[:,0])][:,1]
                rc_file=self.output_folder+'/rc_'+name.split('_')[0]+'_'+name2.split('_')[0]+'.tsv'
                np.savetxt(rc_file, np.array([voltages,data_rc]).T)

            ##############
            #3) the Coverages
            ##############
    
            name=name.split('_')[0]
            cov_names=[cov.split('_')[0] for cov in data.coverage_names]
            if name in cov_names:
                idx=cov_names.index(name)
                data_ref=np.column_stack((data.voltage, data.coverage[:,idx]))
                voltages=data_ref[np.argsort(data_ref[:, 0])][:,0]
                coverages=data_ref[np.argsort(data_ref[:, 0])][:,1]
                cov_file=self.output_folder+'/cov_'+name.split('_')[0]+'.tsv'
                np.savetxt(cov_file,np.array([voltages,coverages]).T)
                self.tp.species[sp]['coverage']=coverages[index]
        #print all intermediate coverages also
        for name in cov_names:
            if name in cov_names and name not in self.tp.species:
                idx=cov_names.index(name)
                data_ref=np.column_stack((data.voltage, data.coverage[:,idx]))
                voltages=data_ref[np.argsort(data_ref[:, 0])][:,0]
                coverages=data_ref[np.argsort(data_ref[:, 0])][:,1]
                cov_file=self.output_folder+'/cov_'+name.split('_')[0]+'.tsv'
                np.savetxt(cov_file,np.array([voltages,coverages]).T)
        ###############
        #4) current densities associated with elementary steps
        ###############
        rate=None
        idx=None
        labels=model.rxn_expressions_names
        for idx,name in enumerate(data.rate_names):
            #"tof" is the signed rate of conversion/active site/s
            rate=data.rate[:,idx]
            data_ref=np.column_stack((data.voltage, rate))
            voltages=data_ref[np.argsort(data_ref[:, 0])][:,0]
            #count ele_g
            #if sp in self.tp.electrode_reactions:
            #    nprod=len([a for a in self.tp.electrode_reactions[sp]['reaction'][1] if a==sp])
            #    nel=self.tp.electrode_reactions[sp]['nel']
            #else:
            #    nprod=1
            #    nel=1
            sol=re.findall('(\d+) ele_g',model.rxn_expressions[idx])
            if len(sol)>0:
                nel=int(sol[0])
            else:
                nel=0
            rates=data_ref[np.argsort(data_ref[:, 0])][:,1]*self.tp.system['active site density']
            current_densities=rates*nel*unit_F/10.
            pol_file=self.output_folder+'/jelem_'+labels[idx]+'.tsv'
            #first write the reaction in the first line
            with open(pol_file,'w') as of:
#                of.write('{} \n'.format(name))
#                print zip(voltages,current_densities)
                for v,c in zip(voltages,current_densities):
                    of.write('{} {}\n'.format(v,c))
#            np.savetxt(pol_file, np.array([voltages,current_densities]).T)


    def get_data(self,pickle_file,model):
        a = pickle.load(open(pickle_file))
        data = Object()
        #COVERAGES
        data.coverage_names = model.output_labels['coverage']
        coverage_map = np.array(a['coverage_map'])
        data.voltage = []
        scaler_array = coverage_map[:,0]
        for s in scaler_array:
            data.voltage.append(s[0])
        coverage_mpf = coverage_map[:,1]
        data.coverage = np.zeros((len(coverage_mpf),len(data.coverage_names)))    
        for i in range(0,len(coverage_mpf)):
            for j in range(0,len(coverage_mpf[i])):
                float_rate = float(coverage_mpf[i][j])
                data.coverage[i][j]=float_rate
        #PRODUCT NAMES
        data.prod_names = model.output_labels['production_rate']
        data.cons_names = model.output_labels['consumption_rate']
        data.turnover_frequency_names = model.output_labels['turnover_frequency']
        data.rate_control_names=model.output_labels['rate_control']

        production_rate_map = np.array(a['production_rate_map'])
        consumption_rate_map = np.array(a['consumption_rate_map'])
        turnover_frequency_map = np.array(a['turnover_frequency_map'])
        rate_control_map=np.array(a['rate_control_map'])

        production_rate_mpf = production_rate_map[:,1]
        consumption_rate_mpf = consumption_rate_map[:,1]
        turnover_frequency_mpf = turnover_frequency_map[:,1]
        rate_control_mpf = rate_control_map[:,1]

        data.production_rate = np.zeros((len(production_rate_mpf),len(data.prod_names)))
        data.consumption_rate = np.zeros((len(consumption_rate_mpf),len(data.cons_names)))
        data.turnover_frequency = np.zeros((len(turnover_frequency_mpf),len(data.turnover_frequency_names)))
        data.rate_control = np.zeros((len(rate_control_mpf),len(data.rate_control_names[0]),len(data.rate_control_names[1])))

        data.voltage = np.zeros((len(production_rate_mpf),1))
        for i in range(0,len(production_rate_mpf)):
            data.voltage[i][0] = production_rate_map[:,0][i][0]
            for j in range(0,len(data.prod_names)):
                float_rate = float(production_rate_mpf[i][j])
                data.production_rate[i][j]=float_rate
        for i in range(0,len(consumption_rate_mpf)):
            for j in range(0,len(data.cons_names)):
                float_rate = float(consumption_rate_mpf[i][j])
                data.consumption_rate[i][j]=float_rate
#                float_rate = float(turnover_frequency_mpf[i][j])
#                data.turnover_frequency[i][j]=float_rate
        for i in range(0,len(turnover_frequency_mpf)):
            for j in range(0,len(data.turnover_frequency_names)):
                float_rate = float(turnover_frequency_mpf[i][j])
                data.turnover_frequency[i][j]=float_rate
        #rate control is
            #2nd index: products
            #3rd index: controlling states
        for i in range(0,len(rate_control_mpf)):
            for j in range(0,len(data.rate_control_names[0])):
                for k in range(0,len(data.rate_control_names[1])):
                    rate_control = float(rate_control_mpf[i][j][k])
                    data.rate_control[i][j][k]=rate_control

        #RATES
        data.rate_names = model.output_labels['rate']
        rate_map = np.array(a['rate_map'])
        rate_mpf = rate_map[:,1]
        data.rate = np.zeros((len(rate_mpf),len(data.rate_names)))
        for i in range(0,len(rate_mpf)):
            for j in range(0,len(rate_mpf[i])):
                float_rate = float(rate_mpf[i][j])
                data.rate[i][j]=float_rate
        return data

    def format_output(self,output_variable):
        """ taken from https://github.com/SUNCAT-Center/catmap/wiki/Accessing-and-reformatting-output Reads in CATMAP log file and reformates the output """
        from catmap.model import ReactionModel
#        output_variable = 'rate' #sys.argv[1]
        logfile = self.model_name+'.log'
        if len(logfile) > 1:
            raise InputError('Ambiguous logfile. Ensure that only one file ends with .log')
        try:
            model = ReactionModel(setup_file=logfile[0])
        except:
            self.tp.logger.error('|    | CM | Error in reading CatMAP log file {}'.format(logfile[0]))
            sys.exit()
        
        if output_variable == 'rate_control':
            dim = 2
        else:
            dim = 1
        
        labels = model.output_labels[output_variable]

        #flatten rate_control labels
        if output_variable == 'rate_control':
            flat_labels = []
            for i in labels[0]:
                for j in labels[1]:
                    flat_labels.append('d'+i+'/d'+j)
            labels = flat_labels
        
        #flatten elementary-step specific labels
        if output_variable in ['rate','rate_constant','forward_rate_constant','reverse_rate_constant']:
            str_labels = []
            for label in labels:
                states = ['+'.join(s) for s in label]
                if len(states) == 2:
                    new_label = '<->'.join(states)
                else:
                    new_label = states[0]+'<->'+states[1]+'->'+states[2]
                str_labels.append(new_label)
            labels = str_labels
        
        table = '\t'.join(list(['descriptor-'+d for d in model.descriptor_names])+list(labels))+'\n'
        
        for pt, output in getattr(model,output_variable+'_map'):
            if dim == 2:
                output = flatten_2d(output)
            table += '\t'.join([str(float(i)) for i in pt+output])+'\n'
        
        f = open(output_variable+'_table.txt','w')
        f.write(table)
        f.close()

    def species_to_catmap(self,sp):
        species=self.tp.species[sp]['symbol']
        if 'catmap_symbol' in self.tp.species[sp]:
            species=self.tp.species[sp]['catmap_symbol']
        species=species.replace('^','')
        species=species.replace('_','')
        if species=='H+':
            species='H'
        elif species=='OH-':
            species='OH'
        elif species=='HCOO-':
            species='HCOOH'
        return species

