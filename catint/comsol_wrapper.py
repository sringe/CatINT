#Comsol class of transport model
from .units import *
import numpy as np
import os
from subprocess import call
import re 
from shutil import copyfile as copy
from copy import deepcopy
from .comsol_model import Model
from .comsol_reader import Reader
from shutil import rmtree
from subprocess import check_output,Popen,PIPE
import datetime
import time
import sys

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

        if not hasattr(self,'exe'):
            if 'bin_path' in self.tp.comsol_args:
                exe_path=self.tp.comsol_args['bin_path']
                self.version=self.tp.comsol_args['bin_version']
            else:
                #try to find executable
                try:
                    s=check_output('locate Multiphysics/bin/comsol',shell=True).decode(sys.stdout.encoding).split('\n')
                    sr=[ss for ss in s if ss.endswith('comsol')]
                    if len(sr)==1:
                        exe_path=sr[0]
                    elif len(sr)==0:
                        raise LookupError('Could not find any comsol binary, when locating Multiphysics/bin/comsol')
                    elif len(sr)>1:
                        raise LookupError('|    | CS | More than one binary found for COMSOL, specify via comsol_args[\'bin_path\']')
                    self.version=re.findall('COMSOL([0-9a-zA-Z.]+)',exe_path)[0]
                except:
                    try:
                        s=check_output('locate multiphysics/bin/comsol',shell=True).decode(sys.stdout.encoding).split('\n')
                        sr=[ss for ss in s if ss.endswith('comsol')]
                        if len(sr)==1:
                            exe_path=sr[0]
                        elif len(sr)==0:
                            self.tp.logger.error('|    | CS | Could not find binary for COMSOL, neither by locating Multiphysics/bin/comsol nor by locating multiphysics/bin/comsol, specify via comsol_args[\'bin_path\']')
                            sys.exit()
                        elif len(sr)>1:
                            self.tp.logger.error('|    | CS | More than one binary found for COMSOL, specify via comsol_args[\'bin_path\']')
                            sys.exit()
                        sol=re.findall('([0-9]{1,2}.?[0-9]{1,2})',exe_path)
                        if len(sol)>0:
                            v=sol[0]
                            if '.' not in v:
                                v=v[0]+'.'+v[1:]
                            if v+'a' in exe_path:
                                v+='a'
                            elif v+'b' in exe_path:
                                v+='b'
                            elif v+'c' in exe_path:
                                v+='c'
                            self.version=v
                        else:
                            self.version=None
                    except:
                        self.tp.logger.error('|    | CS | Could not find binary for COMSOL, possibly locate did not work or try recreating your locate database. specify via comsol_args[\'bin_path\']')
            self.exe=exe_path
            #get COMSOL version, /Applications/COMSOL53a/Multiphysics/bin/comsol
            self.tp.logger.info('|    | CS | Found COMSOL at {}'.format(self.exe))
            self.tp.logger.info('|    | CS | Running COMSOL v{}'.format(self.version))
            if self.version is not None:
                self.tp.logger.info('|    | CS | Check this version number carefully, because some features/settings in COMSOL might have changed and thus COMSOL might crash or this can lead to convergence issues if convergence settings were changed. CatINT has been tested to run succesfully with COMSOL v5.3a and v5.5')
            elif self.version is None:
                self.tp.logger.error('|    | CS | The COMSOl version could not be identified from the executable folder name. Provide it as input to the COMSOl variables manually via comsol_args[\'bin_version\']')
                sys.exit()
            #convert version to a float
            self.tp.comsol_args['bin_version']=float(self.version.replace('a','1').replace('b','2').replace('c','3'))

#        elif os.path.exists('/share/PI/suncat/COMSOL/comsol53a/multiphysics/bin/comsol'):
#            exe_path='/share/PI/suncat/COMSOL/comsol53a/multiphysics/bin/comsol'
#        else:
#            exe_path='/Applications/COMSOL53a/Multiphysics/bin/comsol'
        if not os.path.exists(exe_path):
            self.tp.logger.error('|    | CS | Binary path {} of COMSOL does not exist'.format(exe_path))
        self.tp.path=path
        root=os.getcwd()
        ##RESULTS FOLDER
        self.results_folder_base=self.tp.outputfoldername+'/comsol_results_id'+str(self.tp.mpi_rank).zfill(3)
        self.results_folder=self.results_folder_base
        self.outputs=['concentrations','electrostatics','electrode_flux','rho_charge']
        for out in self.tp.comsol_args['outputs']:
            self.outputs.append(out)
        self.mode=mode
#        #additionally specified comsol output
#        if not hasattr(self.tp,'comsol_outputs_data'):
#            self.tp.comsol_outputs_data={}

    def run(self,label='',restart=False):
        #only_last: if True, update only the data in the global arrays and dictionaries
        #corresponding to the last parameter in the parameter list

        java_name='pnp_transport.java'
        class_name=java_name.split('.')[0]+'.class'
        model_name=java_name.split('.')[0]+'.mph'

        self.results_folder=self.results_folder_base+'_'+label

        if not restart:
        #create results folder
            if not os.path.isdir(self.results_folder):
                os.makedirs(self.results_folder)
            else:
                self.tp.logger.info('|    | CS | Removing old folder and creating new one')
##                ts= time.time()
##                st = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M%S')
#                os.rename(self.results_folder+'/'+java_name,self.results_folder+'/../'+'comsol_backup_'+st+'_'+java_name)
                rmtree(self.results_folder)
                os.makedirs(self.results_folder)

        root=os.getcwd()

        ##WRITE INPUT
        if not restart:
            self.write_input(file_name=java_name,\
                model_name=model_name)

            #copy input file
            ts= time.time()
            st = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M%S')
            copy(java_name,self.results_folder+'/../'+'comsol_backup_'+st+'_'+java_name)


        os.chdir(self.results_folder)

        if not restart:
            ##COMPILE
            self.tp.logger.info('|    | CS | Compiling COMSOL.')
            self.tp.logger.info('|    | CS | Compiling {}'.format(os.getcwd()+'/'+java_name))
            self.tp.logger.info(self.exe+" compile "+os.getcwd()+'/'+java_name+' | tee comsol_compile.out')
            call(self.exe+" compile "+os.getcwd()+'/'+java_name+' | tee comsol_compile.out',shell=True)
#        for line in open(os.getcwd()+'/'+java_name):
#            self.tp.logger.debug(line+'\n')
        #call(self.exe+" compile "+'/'.join([root,file_name])+' | tee '+self.results_folder+'/comsol_compile.out',shell=True)

        ##RUN
        if restart:
            self.tp.logger.info('|    | CS | Restarting COMSOL (starting from previous calculation)')
        else:
            self.tp.logger.info('|    | CS | Starting COMSOL.')


        #recover and continue restart the run, if it failed before starting from the current result
        call(self.exe+" batch -recoverydir $HOME/.comsol/v"+self.version+"/recoveries -inputfile "+os.getcwd()+'/'+class_name+' | tee '+self.results_folder+'/comsol_run.out',shell=True)

        os.chdir(root)

        error=self.check_error()
        if error:
            self.tp.logger.warning('|    | CS | Error appeared in COMSOL calculation.')
        else:
            ##READING RESULTS
            self.tp.logger.info('|    | CS | Reading COMSOL output.')
            self.read_output()

    def check_error(self):
        error=False
        file_name=self.results_folder+'/comsol_run.out'
        for line in open(file_name.replace('\\\\','\\').replace('\\','/'),'r'):
#            if 'Not all parameter steps returned.' in line:
            if 'Failed to find a solution for all parameters' in line or 'Failed to find a solution for the initial parameter' in line:
                error=True
        return error

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
