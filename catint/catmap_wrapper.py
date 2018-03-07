import re
from units import *
import numpy as np
import os
from subprocess import call
from io import replace_line,insert_line,mpi_make_dir
from glob import glob
import sys
import pickle
from catmap.model import ReactionModel
from catmap import analyze
from string import Template
from shutil import copy

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

    def __init__(self,path=os.getcwd(),transport=None,model_name=None,delta_desc=0.05,min_desc=-1.5,max_desc=0.2,n_inter=1):
        #delta_desc is the descriptor delta with which the descriptor axes is resolved
        #in case of potential 0.05 is fine
        #max and min_desc are the  bounds of the descriptor region (here default for potential being the descriptor)
        self.delta_desc=delta_desc
        self.max_desc=max_desc
        self.min_desc=min_desc

        self.n_inter_max=100


        self.use_interactions=False

        if transport is None:
            self.tp.logger.error('No transport object provided for calculator. Stopping here.')
            sys.exit()
        else:
            self.tp=transport
        if not 'n_inter' in self.tp.catmap_args:
            self.n_inter=n_inter #number of steps to converge interactions
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

        self.method=None

    def run(self,desc_val):
        desc_keys=[key for key in self.tp.descriptors]
        self.tp.logger.info(' | CM | Starting CatMAP for {} = {}'.format(desc_keys[0],desc_val[0]))

        self.input_folder=self.input_base_folder+'/desc_'+str(desc_val[0])
        self.output_folder=self.output_base_folder+'/desc_'+str(desc_val[0])

        if not os.path.isdir(self.input_folder):
            os.makedirs(self.input_folder)
        if not os.path.isdir(self.output_folder):
            os.makedirs(self.output_folder)

        #go to output folder and run catmap there
        root=os.getcwd()
        os.chdir(root+'/'+self.input_folder)

        #setting up the mkm files
        mkm_template_file=root+'/'+self.model_name+'_template.mkm'
        if not os.path.exists(mkm_template_file):
            self.tp.logger.error('mkm file {} required for CatMAP does not exist.'.format(mkm_template_file))
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
            self.tp.logger.error('energy data file {} required by CatMAP does not exist.'.format(energies_file))
            sys.exit()
        else:
            copy(energies_file,self.model_name+'_energies.txt')

        #setting up the model
        model = ReactionModel(setup_file = mkm_file, max_log_line_length=0) #we set the maximum line length to 0, because we work
        #with pickle files anyhow
        #output
        model.output_variables+=['consumption_rate','production_rate', 'free_energy', 'selectivity', 'interacting_energy','turnover_frequency']
        def plot_fed(corr):
            ma = analyze.MechanismAnalysis(model)
            ma.surface_colors = ['k','b','r','yellow','green','orange','cyan']
            ma.label_args['size'] = 14
            ma.energy_type = 'free_energy' #can also be free_energy/potential_energy
            ma.include_labels = True #way too messy with labels
            ma.pressure_correction = corr #assume all pressures are 1 bar (so that energies are the same as from DFT)
            ma.coverage_correction = False
            ma.include_labels = True
            if not corr:
                fig = ma.plot(save='FED.pdf',plot_variants=[float(desc_val[0])])
            else:
                fig = ma.plot(save='FED_pressure_corrected.pdf',plot_variants=[float(desc_val[0])])
#        plot_fed(True)
        if not self.use_interactions:
            plot_fed(False)
        else:
            pass
           # idx = [i for i in range(len(model.interacting_energy_map)) if model.interacting_energy_map[i][0][0] == -0.88][0]
           # descrip, coverages = model.coverage_map[idx]
           # rxn_parameters = model.scaler.get_rxn_parameters(descrip)
           # rate_constants = model.solver.get_rate_constants(rxn_parameters,coverages)
           # kfs, krs, dkfs, dkrs = model.rate_constants(rxn_parameters,coverages,
           # model._gas_energies,model._site_energies,
           # model.temperature,model.interaction_response_function,
           # model._mpfloat,model._matrix,model._math.exp)
           # model.solver.get_interacting_energies(rxn_parameters)
           # all_ads = model.adsorbate_names + model.transition_state_names
           # N_ads = len(all_ads)
           # energies = rxn_parameters[:N_ads]
           # eps_vector = rxn_parameters[N_ads:]
           # cvg = coverages + [0]*len(model.transition_state_names)
           # model.interaction_function(cvg,energies,eps_vector,model.thermodynamics.adsorbate_interactions.interaction_response_function,False,False)
        
        #slowly ramp up the interactions if desired
        if self.n_inter.isdigit():
            if self.use_interactions:
                inter=np.linspace(0,self.interaction_strength,self.n_inter)
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

        jj=0
        while True:
            try:
                for ii in inter:
                    i=0
                    if self.use_interactions:
                        for line in open(self.catmap_model):
                            i+=1
                            if 'interaction_strength' in line:
                                replace_line(mkm_file,i-1,'interaction_strength = '+str(ii))
                        self.tp.logger.info(' | CM | Running interaction_strength = {}'.format(ii))
                    model = ReactionModel(setup_file = mkm_file, max_log_line_length=0)
                    model.output_variables+=['consumption_rate','production_rate', 'free_energy', 'selectivity', 'interacting_energy','turnover_frequency']
                    model.run()
            except:
                self.tp.logger.warning('CatMAP did not converge with interaction strength = {}'.format(ii))
                if self.n_inter=='automatic':
                    jj+=1
                    if 10*jj>self.n_inter_max:
                        self.tp.logger.error('Adjusted interaction ramping range is larger than n_inter_max. Adjust n_inter_max in order to run finer range')
                        sys.exit()
                    self.tp.logger.warning('Adjusting interaction ramping to range 0 to 1 with {} steps'.format(10*jj))
                    inter=np.linspace(0,self.interaction_strength,10*jj)
                    pass
                else:
                    break
            else:
                break

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
            self.tp.logger.error('CatMAP did not converge')
        else:
            self.tp.logger.info(' | CM | CatMAP finished successfully in {} steps'.format(max_iter))

        self.read_output(desc_val)
        #species_definitions['CO2_g'] = {'concentration':0.2}

    def update_input(self,desc_val,method='single_point'):

        #go over input file and check if single point or full descriptor range is used as method:
        found_resolution=False
        found_descriptor_range=False
        found_descriptors=False
        self.method=method

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
                self.method='single_point'
                found_descriptors=True
                idesc=i
        if self.method=='single_point':
            if not found_resolution:
                insert_line(self.catmap_model,idesc,'resolution = [1,1]\n')
                found_resolution=True
            if not found_descriptor_range:
                insert_line(self.catmap_model,idesc,'descriptor_range = [[0,0],[0,0]]')
                found_descriptor_range=True
        if not all([found_descriptor_range,found_resolution]):
                #or found_descriptors):
            self.tp.logger.error('Missing resolution definition or descriptor range')
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

        if self.method=='descriptor_range':
            desc_list=self.tp.descriptors[desc_keys[0]]
            if len(desc_list)<5 and desc_1 != 'voltage':
                self.tp.logger.warning('It could be that no pkl files are written out because the given descriptor axis is too coarse, consider a denser axis')
            if desc_1=='voltage':
                #use manual range here that contains current descriptor value
                val=desc_val[0]
                desc_list=[val]
                while val<max(self.max_desc,max(self.tp.descriptors['phiM'])):
                    val+=self.delta_desc
                    desc_list.append(val)
                val=desc_val[0]
                while val>min(self.min_desc,min(self.tp.descriptors['phiM'])):
                    val-=self.delta_desc
                    desc_list.append(val)
            desc_list=sorted(desc_list,key=float)
            min_desc=min(desc_list)
            max_desc=max(desc_list)
            n_desc=len(desc_list)
            delta_desc=desc_list[1]-desc_list[0]
            if delta_desc<1e-9:
                self.tp.logger.warning('The descriptor separation is smaller then 1e-9. The index assignment of the catmap data will fail, so I quit here')
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
                replace_line(self.catmap_model,i-1,'bulk_ph = '+str(self.tp.system['pH']))
                continue
            for sp in self.tp.species:
                sp_cm=self.species_to_catmap(sp)
                sol=re.findall('species_definitions\[\''+sp_cm+'_g\'\].*{\'pressure\':.*}',line)
                if len(sol)>0:
                    #shortly check if concentrations are very negative, stop if they are:
                    if self.tp.species[sp]['surface concentration']<-1.:
                        self.tp.logger.warning('Surface concentration of {} is more negative than 1e-3 mol/L, stopping to be safe.'.format(sp))
                        sys.exit()
                    replace_line(self.catmap_model,i-1,"species_definitions['"+sp_cm+"_g'] = {'pressure':"+str(max(0.,self.tp.species[sp]['surface concentration']/1000.))+"}")
            if self.method=='descriptor_range':
                if 'descriptor_range' in line:
                    replace_line(self.catmap_model,i-1,'descriptor_ranges = [['+str(min_desc)+','+str(max_desc)+'],['+str(desc_val[1])+','+str(desc_val[1])+']]')
                if line.strip().startswith('resolution'):
                    replace_line(self.catmap_model,i-1,'resolution = ['+str(n_desc)+', 1]') #descriptor_names= [\''+desc_1+'\', \''+desc_2+'\']')
            elif self.method=='single_point':
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
                replace_line(self.catmap_model,i-1,'pH = '+str(self.tp.system['surface pH'])+'')

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
        os.chdir(root)
        for sp in self.tp.species:
            self.tp.species[sp]['flux']=0.0
        for sp in self.tp.species: #prod in self.tp.electrode_reactions:
            ###############
            #1) the TOF's save them as fluxes for the individual species
            ###############
            tof=None
            name=self.species_to_catmap(sp)
            name+='_g'
            idx=None
            #"tof" is the signed rate of conversion/active site/s
            if name in data.turnover_frequency_names:
                idx=data.turnover_frequency_names.index(name)
                tof=data.turnover_frequency[:,idx]
            elif sp not in self.tp.system['exclude species'] and sp not in self.tp.electrolyte_list:
                self.tp.logger.warning('No CatMAP TOF data was found for species {}. Check your species definition names. CatINT uses the CatINT equation names to map to the CatMAP names!'.format(sp))
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
            rates=data_ref[np.argsort(data_ref[:, 0])][:,1]*self.tp.system['active site density']*nel*unit_F/nprod/10.
#            currents=self.convert_TOF(data_ref[np.argsort(data_ref[:, 0])][:,1])
            pol_file=self.output_folder+'/j_'+name.split('_')[0]+'.tsv'
            np.savetxt(pol_file, np.array([voltages,rates]).T)
            iv=-1
            for v in voltages:
                iv+=1
                if abs(v-float(desc_val[0]))<1e-10:
                    index=iv
                    break
            self.tp.species[sp]['flux']=rates[index]
            ##############
            #2) the Coverages
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


    def get_data(self,pickle_file,model):
        a = pickle.load(open(pickle_file))
        data = Object()
        #COVERAGES
        data.coverage_names = model.output_labels['coverage']
        print 'data.coverage_names names',data.coverage_names
        print 'interacting names',model.output_labels['interacting_energy']
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
        production_rate_map = np.array(a['production_rate_map'])
        consumption_rate_map = np.array(a['consumption_rate_map'])
        turnover_frequency_map = np.array(a['turnover_frequency_map'])
        production_rate_mpf = production_rate_map[:,1]
        consumption_rate_mpf = consumption_rate_map[:,1]
        turnover_frequency_mpf = turnover_frequency_map[:,1]
        data.production_rate = np.zeros((len(production_rate_mpf),len(data.prod_names)))
        data.consumption_rate = np.zeros((len(consumption_rate_mpf),len(data.cons_names)))
        data.turnover_frequency = np.zeros((len(turnover_frequency_mpf),len(data.turnover_frequency_names)))
        data.voltage = np.zeros((len(production_rate_mpf),1))
        for i in range(0,len(production_rate_mpf)):
            data.voltage[i][0] = production_rate_map[:,0][i][0]
            for j in range(0,len(data.prod_names)):
                float_rate = float(production_rate_mpf[i][j])
                data.production_rate[i][j]=float_rate
        for i in range(0,len(production_rate_mpf)):
            for j in range(0,len(data.prod_names)):
                float_rate = float(consumption_rate_mpf[i][j])
                data.consumption_rate[i][j]=float_rate
#                float_rate = float(turnover_frequency_mpf[i][j])
#                data.turnover_frequency[i][j]=float_rate
        for i in range(0,len(turnover_frequency_mpf)):
            for j in range(0,len(data.turnover_frequency_names)):
                float_rate = float(turnover_frequency_mpf[i][j])
                data.turnover_frequency[i][j]=float_rate
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
            self.tp.logger.error('Error in reading CatMAP log file {}'.format(logfile[0]))
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
        species=species.replace('^','')
        species=species.replace('_','')
        if species=='H+':
            species='H'
        elif species=='OH-':
            species='OH'
        return species

