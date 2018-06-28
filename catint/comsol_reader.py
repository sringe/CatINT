#Reads all the COMSOL output
import re
import numpy as np
from units import *
from scipy.optimize import fsolve,basinhopping

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
            desc_keys=[key for key in self.tp.descriptors]
            if output in self.comsol_args['outputs']:
                toutput=output[0]
            else:
                toutput=output
            self.tp.logger.info('|    | CS | Reading output {}.txt'.format(toutput))
            if output in ['concentrations','electrode_flux']:
                species_var=True
            else:
                species_var=False
            self.read_single_file(results_folder=self.results_folder,output=toutput,species_var=species_var)

    def cs_to_ci(self,var_name):
        """maps COMSOL argument names to CATINT argument names"""
        if var_name=='cp':
            return 'concentration'
        elif var_name=='j':
            return 'electrode_flux'
        elif var_name=='phi':
            return 'potential'
        elif var_name=='es.Ex':
            return 'efield'
        elif var_name=='rho_charge':
            return 'charge_density'
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
        toutput=output

        #####some functions, if nonlinear permittivity function in Stern layer is used
        #we have to solve a possibly non-linear equation system
        if self.tp.system['Stern epsilon']=='Booth':
            def coth(x):
                return 1.+2./(np.exp(2*x)-1)
            def eps(E):
                #parameter for water at zero field:
                n=1.33
                beta=1.41e-8
                #eps0=self.tp.system['Stern_epsilon'] #self.tp.system['epsilon']
                eps0=self.tp.system['epsilon']
                if E>=1e7:
                    return n**2+(eps0-n**2)*3./(beta*E)*(coth(beta*E)-1./(beta*E))
                else:
                    return eps0
        def equations(p):
            E=p
            #Eout,=args
            eq=(E-Eout*self.tp.system['epsilon']/eps(abs(E)))**2
            return eq
        ####end

        jj=-1
        #first read of file to determine if this is a variable tabulated on domain or boundary
        for line in open(results_folder+'/'+toutput+'.txt', 'r'):
            if line.startswith('% Nodes'):
                nodes=int(re.findall('Nodes:\s+(\d+)',line)[0])
                if nodes==2:
                    geo='point'
                else:
                    geo='domain'
                break
        if self.comsol_args['par_name']=='flux_factor':
            #get the current values of the two descriptors
            desc_keys=self.tp.descriptors.keys()
            desc1_val=self.tp.system[desc_keys[0]]
            desc2_val=self.tp.system[desc_keys[1]]
            #get the index:
            alldata_inx=self.tp.alldata_names.index([desc1_val,desc2_val])

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
                    if not self.comsol_args['solver_settings']['solver_sequence']=='tds_elstat':
                        par_name=par_names[i]
                        par_value=float(par_values[i])
                    update_last=False
                    if self.comsol_args['par_name']!='flux_factor':
                        alldata_inx=(i)/len(set(variable_names))
            
                    #print par_values_list
                    #print variable_names
                    if (i+1)>len(par_values_list)*len(set(variable_names))-len(set(variable_names)) or\
                        self.comsol_args['solver_settings']['solver_sequence']=='tds_elstat':
                        update_last=True
                    if self.comsol_args['par_name']=='flux_factor' and not update_last:
                        #if the par name is the flux, this means we run comsol with a single descriptor set and we are only interested in the last value
                        continue
                    if species_var:
                        var=re.findall('([a-zA-Z]+)',var_name)[0]
                        ind_sp=int(re.findall('[a-zA-Z]+(\d+)',var_name)[0])
                        sp=[sp2 for j,sp2 in enumerate(self.tp.species) if j+1==ind_sp][0]
                        var_name=self.cs_to_ci(var)
                        if var=='j' and self.tp.use_catmap:
                            #corresponding variables will be updated by catmap instead
                            continue
                        update_pH=None
                        if var=='cp':
                            if 'H+' in self.tp.species:
                                if sp=='H+':
                                    update_pH='H+'
                            elif 'OH-' in self.tp.species:
                                if sp=='OH-':
                                    update_pH='OH-'
                        if jj==0 and geo=='domain':
                            self.tp.species[sp][var_name]=[]
#                            if sp not in self.tp.alldata[alldata_inx]['species']:
#                                self.tp.alldata[alldata_inx]['species'][sp]={}
                            self.tp.alldata[alldata_inx]['species'][sp][var_name]=[]
                            if update_pH is not None:
                                self.tp.system['pH']=[]
                                self.tp.alldata[alldata_inx]['system']['pH']=[]
                        #variables derived from concentrations:
                        if var=='cp':
                            if update_pH=='H+':
                                if update_last:
                                    self.tp.system['pH'].append(-np.log10(float(lss)/1000.))
                                self.tp.alldata[alldata_inx]['system']['pH'].append(-np.log10(float(lss)/1000.))
                            elif update_pH=='OH-':
                                if update_last:
                                    self.tp.system['pH'].append(14+np.log10(float(lss)/1000.))
                                self.tp.alldata[alldata_inx]['system']['pH'].append(14+np.log10(float(lss)/1000.))

                        if x==0.0 and var=='cp':
                            if update_last:
                                self.tp.species[sp]['surface_concentration']=float(lss)
                            self.tp.alldata[alldata_inx]['species'][sp]['surface_concentration']=float(lss)
                            #our reference for the pH is the SHE at 1M concentrations, so the activity is defined as
                            #local concentration of protons devided by 1M
                            if 'H+' in self.tp.species:
                                if sp=='H+':
                                    if float(lss)<0:
                                        self.tp.logger.warning('|    | CS | COMSOL returned negative proton concentrations, pH cannot be evaluated.')
                                        self.tp.logger.warning('|    | CS | Do not update pH here to enable futher calculation.')
                                    else:
                                        if update_last:
                                            self.tp.system['surface_pH']=-np.log10(float(lss)/1000.) #+self.tp.system['bulk_pH']
                                        self.tp.alldata[alldata_inx]['system']['surface_pH']=-np.log10(float(lss)/1000.)
                            elif 'OH-' in self.tp.species:
                                if sp=='OH-':
                                    if float(lss)<0:
                                        self.tp.logger.warning('|    | CS | COMSOL returned negative hydroxide concentrations, pH cannot be evaluated.')
                                        self.tp.logger.warning('|    | CS | Do not update pH here to enable futher calculation.')
                                    else:
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
                        if var=='j' and sp in self.tp.electrode_reactions and not self.tp.use_catmap:
                            #(in case we use catmap, this will be responsible for evaluating the current densities and rates)
                            nprod=len([a for a in self.tp.electrode_reactions[sp]['reaction'][1] if a==sp])
                            if update_last:
                                self.tp.species[sp]['electrode_current_density']=float(lss)*self.tp.electrode_reactions[sp]['nel']*unit_F/nprod/10.
                            self.tp.alldata[alldata_inx]['species'][sp]['electrode_current_density']=float(lss)*self.tp.electrode_reactions[sp]['nel']*unit_F/nprod/10.
                        elif var=='j' and sp in self.tp.electrode_reactions:
                            self.tp.alldata[alldata_inx]['species'][sp]['electrode_current_density']=self.tp.species[sp]['electrode_current_density']
                    else:
                        var_name=self.cs_to_ci(var_name)
                        if jj==0 and geo=='domain':
                            self.tp.system[var_name]=[]
                            self.tp.alldata[alldata_inx]['system'][var_name]=[]
                        if update_last:
                            if geo=='domain':
                                self.tp.system[var_name].append(float(lss))
                                if x==0.0 and var_name=='efield':
                                    self.tp.system['surface_efield']=float(lss)
                                    if type(self.tp.system['Stern epsilon'])!=str:
                                        self.tp.system['Stern_efield']=float(lss)*self.tp.system['epsilon']/self.tp.system['Stern epsilon']
                                        self.tp.system['Stern_epsilon_func']=self.tp.system['Stern epsilon']
                                    else:
                                        #self.tp.system['Stern_efield'],=fsolve(equations,(float(lss),),args=(float(lss),))
                                        E0=float(lss)*self.tp.system['epsilon']/10. #self.tp.system['Stern epsilon']
#                                        self.tp.system['Stern_efield'],=fsolve(equations,(float(lss),),args=(E0,))
                                        minimizer_kwargs = {"method": "BFGS"}
                                        Eout=float(lss)
                                        self.tp.system['Stern_efield'],=basinhopping(equations, E0, minimizer_kwargs=minimizer_kwargs,
                                                niter=200).x
                                        self.tp.system['Stern_epsilon_func']=eps(abs(self.tp.system['Stern_efield']))
                                elif x==0.0 and var_name=='potential':
                                    self.tp.system['surface_potential']=float(lss)
                            else:
                                self.tp.system[var_name]=float(lss)
                        if geo=='domain':
                            self.tp.alldata[alldata_inx]['system'][var_name].append(float(lss))
                            if x==0.0 and var_name=='efield':
                                self.tp.alldata[alldata_inx]['system']['surface_efield']=float(lss)
                                if type(self.tp.system['Stern epsilon'])!=str:
                                    self.tp.alldata[alldata_inx]['system']['Stern_efield']=float(lss)*self.tp.system['epsilon']/self.tp.system['Stern epsilon']
                                    self.tp.alldata[alldata_inx]['system']['Stern_epsilon_func']=self.tp.system['Stern epsilon']
                                else:
                                    E0=float(lss)*self.tp.system['epsilon']/10. #self.tp.system['Stern epsilon']
                                    #self.tp.alldata[alldata_inx]['system']['Stern_efield'],=fsolve(equations,(float(lss),),args=(E0,))
                                    minimizer_kwargs = {"method": "BFGS"}
                                    Eout=float(lss)
                                    self.tp.alldata[alldata_inx]['system']['Stern_efield'],=basinhopping(equations, E0, minimizer_kwargs=minimizer_kwargs,
                                            niter=200).x
                                    self.tp.alldata[alldata_inx]['system']['Stern_epsilon_func']=eps(abs(self.tp.system['Stern_efield']))
                            elif x==0.0 and var_name=='potential':
                                self.tp.alldata[alldata_inx]['system']['surface_potential']=float(lss)
                        else:
                            self.tp.alldata[alldata_inx]['system'][var_name]=float(lss)
                continue

            if start_reading:
                if self.comsol_args['solver_settings']['solver_sequence']=='tds_elstat':
                    #we do not use the flux_factor in the final step
                    a=re.findall('([a-zA-Z.\_-]{1,50}\d{0,4})(?:\s)(\(.*?\))?\s*',line)
                    a=a[1:]
                    par_names='flux_factor'
                    par_values=[1.0]
                    par_values_list=[1.0]
                else:
                    a=re.findall('\s+([a-zA-Z.\_-]{1,50}\d{0,4})(?:\s)(\(.*?\))?\s*@',line)
                    par_names=re.findall('@\s?([a-zA-Z0-9]+)',line)
                    par_values=re.findall('[a-zA-Z0-9]+\s?=\s?(-?\d+.?\d*)',line)
                    par_values_list=sorted(list(set([float(v) for v in set(par_values)])))
                variable_names,variable_units=map(list, zip(*a))
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
