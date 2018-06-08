#Comsol class of transport model
from units import *
import numpy as np
import os
from subprocess import call
import re 
from shutil import copyfile as copy
from copy import deepcopy
from comsol_model import Model
from comsol_reader import Reader
from shutil import rmtree

class Comsol():
    """This class does all operations need to write input files for comsol and read output"""
    def __init__(self,path=os.getcwd(),transport=None,mode='time-dependent'):
        ##IMPORT TP
        if transport is None:
            print('No transport object provided for COMSOL wrapper. Stopping here.')
            sys.exit()
        else:
            self.tp=transport #transport object
        ##COMSOL EXE
        if 'bin_path' in self.tp.comsol_args:
            exe_path=self.tp.comsol_args['bin_path']
        elif os.path.exists('/share/PI/suncat/COMSOL/comsol53a/multiphysics/bin/comsol'):
            exe_path='/share/PI/suncat/COMSOL/comsol53a/multiphysics/bin/comsol'
        else:
            exe_path='/Applications/COMSOL53a/Multiphysics/bin/comsol'
        if not os.path.exists(exe_path):
            self.tp.logger.error(' | CS | Binary path {} of COMSOL does not exist'.format(exe_path))
        self.tp.path=path
        self.exe=exe_path
        root=os.getcwd()
        ##RESULTS FOLDER
        self.results_folder_base=root+'/'+self.tp.outputfoldername+'/comsol_results_id'+str(self.tp.mpi_rank).zfill(3)
        self.results_folder=self.results_folder_base
        self.outputs=['concentrations','electrostatics','electrode_flux','rho_charge']
        for out in self.tp.comsol_args['outputs']:
            self.outputs.append(out)
        self.mode=mode
#        #additionally specified comsol output
#        if not hasattr(self.tp,'comsol_outputs_data'):
#            self.tp.comsol_outputs_data={}

    def run(self,label=''):
        #only_last: if True, update only the data in the global arrays and dictionaries
        #corresponding to the last parameter in the parameter list

        java_name='pnp_transport.java'
        class_name=java_name.split('.')[0]+'.class'
        model_name=java_name.split('.')[0]+'.mph'

        self.results_folder=self.results_folder_base+'_'+label

        #create results folder
        if not os.path.isdir(self.results_folder):
            os.makedirs(self.results_folder)
        else:
            self.tp.logger.info(' | CS | Removing old folder and creating new one')
            rmtree(self.results_folder)
            os.makedirs(self.results_folder)

        root=os.getcwd()

        ##WRITE INPUT
        self.write_input(file_name=java_name,\
            model_name=model_name)

        os.chdir(self.results_folder)

        ##COMPILE
        self.tp.logger.info(' | CS | Compiling COMSOL.')
        self.tp.logger.info(self.exe+" compile "+os.getcwd()+'/'+java_name+' | tee comsol_compile.out')
        call(self.exe+" compile "+os.getcwd()+'/'+java_name+' | tee comsol_compile.out',shell=True)
        for line in open(os.getcwd()+'/'+java_name):
            self.tp.logger.debug(line+'\n')
        #call(self.exe+" compile "+'/'.join([root,file_name])+' | tee '+self.results_folder+'/comsol_compile.out',shell=True)
        self.tp.logger.info(' | CS | Compiling {}'.format(os.getcwd()+'/'+java_name))

        ##RUN
        self.tp.logger.info(' | CS | Starting COMSOL.')
        self.tp.logger.debug(self.exe+" batch -inputfile "+os.getcwd()+'/'+class_name+' | tee comsol_run.out')
        #call(self.exe+" batch -inputfile "+'/'.join([root,'.'.join(file_name.split('.')[:-1])+".class"])+' | tee '+self.results_folder+'/comsol_run.out',shell=True)
        call(self.exe+" batch -inputfile "+os.getcwd()+'/'+class_name+' | tee comsol_run.out',shell=True)
        #self.tp.logger.info('Running {}'.format('/'.join([root,'.'.join(file_name.split('.')[:-1])+".class"])))
        self.tp.logger.info(' | CS | Running {}'.format(os.getcwd()+'/'+class_name))

        os.chdir(root)

        ##READING RESULTS
        self.tp.logger.info(' | CS | Reading COMSOL output.')
        self.read_output()

    def write_input(self,file_name='',model_name='',model_comments=''):
        model=Model(file_name,model_name,model_comments,outputs=self.outputs,\
            transport=self.tp,comsol_args=self.tp.comsol_args,\
            results_folder=self.results_folder)
        #build model from components
        model.build_all()
        #write to file
        model.write()
        root=os.getcwd()
        copy(root+'/'+file_name,self.results_folder+'/'+file_name)

    def read_output(self):
        reader=Reader(outputs=self.outputs,transport=self.tp,\
            comsol_args=self.tp.comsol_args,\
            results_folder=self.results_folder)
        reader.read_all()
