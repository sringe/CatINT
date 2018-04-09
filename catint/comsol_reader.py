#Reads all the COMSOL output
import re
import numpy as np
from units import *

class Reader():
    def __init__(self,transport=None,outputs=[],comsol_args={},\
            results_folder=''):

        ##IMPORT TP
        if transport is None:
            print('No transport object provided for COMSOL Reader. Stopping here.')
            sys.exit()
        else:
            self.tp=transport

        self.outputs=outputs
        self.results_folder=results_folder
        self.comsol_args=comsol_args


    def read_all(self):
        for output in self.outputs:
            self.tp.logger.info(' | CS | Reading output {}.txt'.format(output))
            if output in ['concentrations','electrode_flux']:
                species_var=True
            else:
                species_var=False
            self.read_single_file(results_folder=self.results_folder,output=output,species_var=species_var)

    def cs_to_ci(self,var_name):
        """maps COMSOL argument names to CATINT argument names"""
        if var_name=='cp':
            return 'concentration'
        elif var_name=='j':
            return 'electrode_flux'
        elif var_name=='phi':
            return 'potential'
        elif var_name=='Ex':
            return 'efield'
        else:
            return var_name

    def read_single_file(self,results_folder='',output='',species_var=True):
        """
        Read a single COMSOL output file. species_var specifies if the variable in this file is a species variable (important because numbers need to be replaced by species name.
        results_folder:     folder where outputs files are located
        output:             name of output file + .txt
        species_var:        is the output variable in the output file a species dependent variable?
        """
        start_reading=False
        initialized=False
        jj=-1
        desc_keys=[key for key in self.tp.descriptors]
        if output in self.comsol_args['outputs']:
            toutput=output[1]
        else:
            toutput=output

        #first read of file to determine if this is a variable tabulated on domain or boundary
        for line in open(results_folder+'/'+toutput+'.txt', 'r'):
            if line.startswith('% Nodes'):
                nodes=int(re.findall('Nodes:\s+(\d+)',line)[0])
                if nodes==2:
                    geo='point'
                else:
                    geo='domain'
                break

        if geo != 'point':
            xmesh=[]
        for line in open(results_folder+'/'+toutput+'.txt', 'r'):
            if initialized:
                jj+=1 #counter for relevant lines
                ls=line.split()
                x=float(ls[0])
                if geo!='point':
                    xmesh.append(x)
                for i,lss in enumerate(ls[1:]):
                    if geo=='point' and x!=0:
                        #we need to read only the first line, otherwise break
                        break
                    var_name=variable_names[i]
                    par_name=par_names[i]
                    par_value=float(par_values[i])
                    update_last=False
                    if self.comsol_args['par_values'][-1]==par_value:
                        update_last=True
                    if self.comsol_args['par_name']=='flux':
                        alldata_inx=0
                    else:
                        alldata_inx=(i)/len(set(variable_names))
                    if self.comsol_args['par_name']=='flux' and not update_last:
                        #if the par name is the flux, this means we run comsol with a single descriptor set and we are only interested in the last value
                        continue
                    if species_var:
                        var=re.findall('([a-zA-Z]+)',var_name)[0]
                        ind_sp=int(re.findall('[a-zA-Z]+(\d+)',var_name)[0])
                        sp=[sp2 for j,sp2 in enumerate(self.tp.species) if j+1==ind_sp][0]
                        var_name=self.cs_to_ci(var)
                        if jj==0 and geo=='domain':
                            self.tp.species[sp][var_name]=[]
#                            if sp not in self.tp.alldata[alldata_inx]['species']:
#                                self.tp.alldata[alldata_inx]['species'][sp]={}
                            self.tp.alldata[alldata_inx]['species'][sp][var_name]=[]
                            if var=='cp':
                                self.tp.system['pH']=[]
                                self.tp.alldata[alldata_inx]['system']['pH']=[]
                        #variables derived from concentrations:
                        if var=='cp':
                            if 'H+' in self.tp.species:
                                if sp=='H+':
                                    if update_last:
                                        self.tp.system['pH'].append(-np.log10(float(lss)/1000.))
                                    self.tp.alldata[alldata_inx]['system']['pH'].append(-np.log10(float(lss)/1000.))
                            elif 'OH-' in self.tp.species:
                                if sp=='OH-':
                                    if update_last:
                                        self.tp.system['pH'].append(14+np.log10(float(lss)/1000.))
                                    self.tp.alldata[alldata_inx]['system']['pH'].append(14+np.log10(float(lss)/1000.))

                        if x==0.0 and var=='cp':
                            if update_last:
                                self.tp.species[sp]['surface_concentration']=float(lss)
                            self.tp.alldata[alldata_inx]['species'][sp]['surface_concentration']=float(lss)
                            if 'H+' in self.tp.species:
                                if sp=='H+':
                                    if update_last:
                                        self.tp.system['surface_pH']=-np.log10(float(lss)/1000.)
                                    self.tp.alldata[alldata_inx]['system']['surface_pH']=-np.log10(float(lss)/1000.)
                            elif 'OH-' in self.tp.species:
                                if sp=='OH-':
                                    if update_last:
                                        self.tp.system['surface_pH']=14+np.log10(float(lss)/1000.)
                                    self.tp.alldata[alldata_inx]['system']['surface_pH']=14+np.log10(float(lss)/1000.)
                        if update_last:
                            if geo=='domain':
                                self.tp.species[sp][var_name].append(float(lss))
                            else:
                                self.tp.species[sp][var_name]=float(lss)
                        if geo=='domain':
                            self.tp.alldata[alldata_inx]['species'][sp][var_name].append(float(lss))
                        else:
                            self.tp.alldata[alldata_inx]['species'][sp][var_name]=float(lss)
                        #calculate electrode current density from electrode flux
                        if var=='j' and sp in self.tp.electrode_reactions:
                            nprod=len([a for a in self.tp.electrode_reactions[sp]['reaction'][1] if a==sp])
                            if update_last:
                                self.tp.species[sp]['electrode_current_density']=float(lss)*self.tp.electrode_reactions[sp]['nel']*unit_F/nprod/10.
                            self.tp.alldata[alldata_inx]['species'][sp]['electrode_current_density']=float(lss)*self.tp.electrode_reactions[sp]['nel']*unit_F/nprod/10.
                    else:
                        var_name=self.cs_to_ci(var_name)
                        if jj==0 and geo=='domain':
                            self.tp.system[var_name]=[]
                            self.tp.alldata[alldata_inx]['system'][var_name]=[]
                        if update_last:
                            if geo=='domain':
                                self.tp.system[var_name].append(float(lss))
                            else:
                                self.tp.system[var_name]=float(lss)
                        if geo=='domain':
                            self.tp.alldata[alldata_inx]['system'][var_name].append(float(lss))
                        else:
                            self.tp.alldata[alldata_inx]['system'][var_name]=float(lss)
                continue

            if start_reading:
#                variable_names=re.findall('([a-zA-Z0-9]+)\s+\(.*\)\s+@\s+[a-zA-Z]+\s?=\s?-?\d+\.?\d?',line)
#                variable_names=re.findall('([a-zA-Z]+\d+)',line)
                if species_var:
                    variable_names=re.findall('([a-zA-Z]+[0-9]+)',line)
                else:
                    variable_names=re.findall('([a-zA-Z]+)\s*\(',line)
#                par_names=re.findall('[a-zA-Z0-9]+\s+\(.*\)\s+@\s+([a-zA-Z]+)\s?=\s?-?\d+\.?\d?',line)
                par_names=re.findall('@\s?([a-zA-Z0-9]+)',line)
                #par_values=re.findall('[a-zA-Z0-9]+\s+\(.*\)\s+@\s+[a-zA-Z]+\s?=\s?(-?\d+\.?\d?)',line)
                par_values=re.findall('[a-zA-Z0-9]+\s?=\s?(-?\d+.?\d*)',line)
                initialized=True
            if not line.startswith('% Description'):
                continue
            else:
                start_reading=True
        if geo!='point':
            self.tp.xmesh=np.array(xmesh)
            self.tp.xmax=max(xmesh)
            self.tp.nx=len(xmesh)
            self.tp.dx=self.tp.xmax/self.tp.nx

