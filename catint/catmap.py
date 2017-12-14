from units import *
import numpy as np
import os
from subprocess import call
from io import replace_line
from glob import glob
import sys
from catmap.model import ReactionModel

def flatten_2d(output):
    "Helper function for flattening rate_control output"
    flat = []
    for x in output:
        flat+= x
    return flat

class Catmap():
    """This class modifies the concentrations in the catmap input file (defined by user) and runs catmap. The output is the read in and the boundary conditions adjusted properly"""

    def __init__(self,path=os.getcwd(),transport=None):
        if transport is None:
            self.tp.logger.error('No transport object provided for calculator. Stopping here.')
            sys.exit()
        else:
            self.tp=transport
        self.tp.path=path
        self.results_folder='results_catmap'
        self.catmap_model='catmap.mkm'

    def run(self):
        self.update_input()
        call('python catmap_job.py',shell=True)
        self.read_output()
        #species_definitions['CO2_g'] = {'concentration':0.2}

    def update_input(self):
        #we need to replace the concentration lines in the file
        i=0
        for line in open(self.catmap_model):
            i+=1
            for sp in self.tp.species:
                species=self.tp.species[tp][symbol]
                if all([a in line for a in ['species_definitions',species,'concentration']]):
                    replace_line(self.catmap_model,i,"species_definitions['"+species+"_g'] = {'concentration':"+self.tp.species[sp]['updated concentration']+"}")

    def read_output(self,output_variable):
        with open(output_variable+'_table.txt','r') as infile:
            for line in infile:


    def format_output(self,output_variable):
        """ taken from https://github.com/SUNCAT-Center/catmap/wiki/Accessing-and-reformatting-output Reads in CATMAP log file and reformates the output """
#        output_variable = 'rate' #sys.argv[1]
        logfile = glob('*.log')
        if len(logfile) > 1:
            raise InputError('Ambiguous logfile. Ensure that only one file ends with .log')
        model = ReactionModel(setup_file=logfile[0])
        
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

