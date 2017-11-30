from units import *
import numpy as np
import os
from subprocess import call
from io import replace_line

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

    def read_output(self):
        
        self.flux_bound=...
