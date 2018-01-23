from units import *
import numpy as np
import os
from subprocess import call
from io import replace_line
from glob import glob
import sys
import pickle
from catmap.model import ReactionModel
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

    def __init__(self,path=os.getcwd(),transport=None,model_name=None):
        if transport is None:
            self.tp.logger.error('No transport object provided for calculator. Stopping here.')
            sys.exit()
        else:
            self.tp=transport
        self.tp.path=path
        self.results_folder='results_catmap'
        self.catmap_model=model_name+'.mkm'
        if model_name is not None:
            self.model_name=model_name
        else:
            self.model_name='catmap'
        self.results_folder=self.tp.outputfoldername+'/catmap_results'

    def run(self):
        #update the mkm file with current concentrations
        self.update_input()

        #go to output folder and run catmap there
        root=os.getcwd()
        os.chdir(root+'/'+self.tp.outputfoldername)

        #setting up the mkm files
        mkm_template_file=root+'/'+self.model_name+'_template.mkm'
        if not os.path.exists(mkm_template_file):
            self.tp.logger.error('mkm file {} required for CatMAP does not exist.'.format(mkm_template_file))
            sys.exit()
        else:
            copy(mkm_template_file,self.model_name+'_template.mkm')
        mkm_template = Template(open(mkm_template_file).read())
        mkm_text = mkm_template.substitute(pH_new = '7')
        mkm_file = self.model_name+'.mkm'
        with open(mkm_file,'w') as f:
            f.write(mkm_text)

        #energies file
        energies_file=root+'/'+self.model_name+'_energies.txt'
        if not os.path.exists(energies_file):
            self.tp.logger.error('energy data file {} required by CatMAP does not exist.'.format(energies_file))
            sys.exit()
        else:
            copy(energies_file,self.model_name+'_energies.txt')

        #setting up the model
        model = ReactionModel(setup_file = mkm_file)
        #output
        model.output_variables+=['production_rate', 'free_energy', 'selectivity', 'interacting_energy']
        #run!
        model.run()
    
        #go back to root folder
        os.chdir(root)
        self.read_output()
        #species_definitions['CO2_g'] = {'concentration':0.2}

    def update_input(self):
        #1) work on input files
        #we need to replace the concentration lines in the file
        i=0
        for line in open(self.catmap_model):
            i+=1
            for sp in self.tp.species:
                species=self.tp.species[sp]['symbol']
                if all([a in line for a in ['species_definitions',species,'concentration']]):
                    replace_line(self.catmap_model,i,"species_definitions['"+species+"_g'] = {'concentration':"+self.tp.species[sp]['updated concentration']+"}")
        #2) create results folder
        if not os.path.isdir(self.results_folder):
            os.makedirs(self.results_folder)

    #SETTINGS
    def convert_TOF(self,A): # Given a list, convert all the TOF to j(mA/cm2) using 0.161*TOF(According to Heine's ORR paper)
        B = [-0.161*rate for rate in A]
        return B

    def read_output(self):
        log_file = self.results_folder+'/'+self.model_name+'.log'
        model = ReactionModel(setup_file = log_file)
        pickle_file = self.results_folder+'/'+self.model_name+'.pkl'
        data = self.get_data(pickle_file,model)
        idx=data.prod_names.index('H2_g')
        for prod in sum([self.tp.electrode_reactions,['H']],[]):
            idx=data.prod_names.index(prod+'_g')
            data_ref=np.column_stack((data.voltage, data.production_rate[:,idx]))
            voltages=data_ref[np.argsort(data_ref[:, 0])][:,0]
            currents=self.convert_TOF(data_ref[np.argsort(data_ref[:, 0])][:,1])
            pol_file=self.results_folder+'/j_'+prod+'.tsv'
            np.savetxt(pol_file, np.array([voltages,currents]).T)
            #set the rate of the reactions to result
            self.tp.electrode_reactions[prod]['rates']=rates

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
        production_rate_map = np.array(a['production_rate_map'])
        production_rate_mpf = production_rate_map[:,1]
        data.production_rate = np.zeros((len(production_rate_mpf),len(data.prod_names)))
        data.voltage = np.zeros((len(production_rate_mpf),1))
        for i in range(0,len(production_rate_mpf)):
            data.voltage[i][0] = production_rate_map[:,0][i][0]
            for j in range(0,len(data.prod_names)):
                float_rate = float(production_rate_mpf[i][j])
                data.production_rate[i][j]=float_rate
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

