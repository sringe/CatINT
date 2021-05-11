#COMSOL Model class, providing functions analogeously to the java class functions
import sys
import re
import numpy as np

class Model():

    def __init__(self,file_name='pnp_transport.java',model_name='pnp_transport.mph',model_comments='',comsol_args={},transport=None,outputs=[],\
            results_folder=''):
        """
        outputs:    outputs to be created
        comsol_args:comsol arguments. by default taken from input file, but the Model class can 
            be also called with a new dictionary of comsol_args if wished, e.g.:
            - studies:    stat or time
            - par_name:   parameter used for parametric sweep
            - par_values: values of parameter used for parametric sweep
        transport:  the CatINT transport instance with all system and species dictionaries
        """

        ##IMPORT TP
        if transport is None:
            print('No transport object provided for COMSOL model. Stopping here.')
            sys.exit()
        else:
            self.tp=transport

        self.outputs=outputs

        self.comsol_args=comsol_args
        self.file_name=file_name
        self.par_name=comsol_args['par_name']
        self.par_values=comsol_args['par_values']
        self.par_method=comsol_args['par_method']
        self.studies=comsol_args['studies']
        self.results_folder=results_folder
        self.version=comsol_args['bin_version']

        self.js=''
        self.js+='import com.comsol.model.*;\n'
        self.js+='import com.comsol.model.util.*;\n'
        self.js+='/*\n'
        self.js+=' *Transport Model as created by Transport Framework\n'
        self.js+=' */\n'
        self.js+='public class '+str('.'.join(file_name.split('.')[:-1]))+' {\n'
        self.js+='  public static Model run() {\n'
        self.js+='    Model model = ModelUtil.create("Model");\n'
        self.js+='    model.modelPath("'+self.tp.path+'");\n'
        self.js+='    model.label("'+model_name+'");\n'
        self.js+='    model.comments("'+model_comments+'");\n'

    def build_all(self):

        #add all classes need for model_type
        if self.comsol_args['model_type']=='tp_dilute_species':
            #first make a list of all classes
            classes=[]
            print('self.tp',self.tp)
            classes.append(self.param(self.tp,comsol_args=self.comsol_args))
            classes.append(self.components())
            classes.append(self.variables(self.tp,index=1,geo='b1',\
                comsol_args=self.comsol_args))
            classes.append(self.variables(self.tp,index=2,geo='d1',\
                comsol_args=self.comsol_args))
            classes.append(self.cpl(index=1,geo='b1'))
            classes.append(self.cpl(index=2,geo='d1'))
            classes.append(self.physics(self.tp,methods=['es','tds','ge'],\
                comsol_args=self.comsol_args))
            classes.append(self.multiphysics(method='tds-es'))
            classes.append(self.probe())
            classes.append(self.mesh())
            classes.append(self.std(self.tp,studies=self.studies,par_name=self.par_name,\
                par_values=self.par_values,par_method=self.par_method,comsol_args=self.comsol_args))
            classes.append(self.output(self.tp,self.outputs,self.results_folder,self.comsol_args))
            #now get all the strings from the different classes
            for c in classes:
                self.js+=c.get()

        #finalize model:
        self.js+='    return model;\n'
        self.js+='  }\n'
        self.js+='\n'
        self.js+='  public static Model run2(Model model) {\n'
        self.js+='  return model;\n'
        self.js+='  }\n'
        self.js+='\n'
        self.js+='  public static void main(String[] args) {\n'
        self.js+='    Model model = run();\n'
        self.js+='    run2(model);\n'
        self.js+='  }\n'
        self.js+='\n'
        self.js+='}\n'

    def get(self):
        return self.js

    def write(self):
        with open(self.file_name,'w') as inp:
            inp.write(self.js)

    ########################################################################
    ########################################################################
    ##################C O M S O L #####  C L A S S E S######################
    ########################################################################
    ########################################################################

    class output():
        def __init__(self,tp=None,outputs=[],results_folder='',comsol_args={}):
            if tp is None:
                return
            else:
                self.tp=tp
            self.comsol_args=comsol_args
            self.results_folder=results_folder
            self.outputs=outputs
            self.s=''
            self.s+='/*\n'
            self.s+=' *OUTPUTS\n'
            self.s+=' */\n'
            label_append=''
            if 'concentrations' in outputs:
                self.export_data('cp','mol/m^3', 'Concentrations'+label_append, self.results_folder+'/concentrations.txt',export_count=1)
            if 'electrostatics' in outputs:
                self.export_data(['phi','es.Ex'],['V','V/m'],'Potential, Field'+label_append, self.results_folder+'/electrostatics.txt',export_count=2)
            if 'electrode_flux' in outputs:
                self.export_data('j','mol/m^2/s','Electrode flux'+label_append, self.results_folder+'/electrode_flux.txt',export_count=3,geo='b1')
            if 'rho_charge' in outputs:
                self.export_data(['rho_charge'],['e*mol/m^3'],'Local Charge Density'+label_append, self.results_folder+'/rho_charge.txt',export_count=4)
            i=0
            for out in outputs:
                if out in self.comsol_args['outputs']:
                    i+=1
                    self.export_data([out[0]],[out[1]],out[0]+label_append, self.results_folder+'/'+out[0]+'.txt',i+4)

        def export_data(self,var_name='cp',unit_name='mol/m^3',label='Concentrations',file_name='results/concentrations.txt',export_count=1,geo='d1',dset_count=1):
            """exports a quantity of interest. var_name/unit_name/label can be either a single string, in which case
            it is assumed that var_name should be exported for each individual species. or it can be a list, over which
            to iterate the output"""

            self.s+='    model.result().export().create("data'+str(export_count)+'", "Data");\n'
            c_str=""
            unit_str=""
            name_str=""
            if type(var_name)==str:
                for sp,i in zip(self.tp.species,list(range(len(self.tp.species)))):
                    c_str+="\""+var_name+str(i+1)+"\""
                    unit_str+="\""+unit_name+"\""
                    name_str+="\""+label+" "+self.tp.species[sp]['name']+"\""
                    if i != len(self.tp.species)-1:
                        c_str+=", "
                        unit_str+=", "
                        name_str+=", "
            elif type(var_name)==list:
                i=-1
                for v,u in zip(var_name,unit_name):
                    i+=1
                    c_str+="\""+v+"\""
                    unit_str+="\""+u+"\""
                    name_str+="\""+v+"\""
                    if i != len(var_name)-1:
                        c_str+=", "
                        unit_str+=", "
                        name_str+=", "
            if self.comsol_args['solver_settings']['solver_sequence']=='tds_elstat':
                #assign corret dataset to this export
                self.s+='    model.result().export("data'+str(export_count)+'").set("data", "dset5");\n'
            self.s+='    model.result().export("data'+str(export_count)+'").set("expr", new String[]{'+c_str+'});\n'
            self.s+='    model.result().export("data'+str(export_count)+'").set("unit", new String[]{'+unit_str+'});\n'
            self.s+='    model.result().export("data'+str(export_count)+'").set("descr", new String[]{'+name_str+'});\n'
            self.s+='    model.result().export("data'+str(export_count)+'").set("filename", "'+file_name+'");\n'
            if geo.startswith('b'):
                self.s+='    model.result().export("data'+str(export_count)+'").set("level", "point");\n'
            self.s+='    model.result().export("data'+str(export_count)+'").run();\n'

        def get(self):
            return self.s



    class std():
        """Studies"""
        def __init__(self,tp=None,studies=[],par_name='',par_values=[],par_method='internal',comsol_args={}):
            self.comsol_args=comsol_args
            self.s=''
            if tp is None:
                return
            else:
                self.tp=tp
            self.s+='/*\n'
            self.s+=' *STUDIES & SOLVERS\n'
            self.s+=' */\n'
            self.std_index=0
            self.sol_index=0
            std_sol_index={}

            if self.comsol_args['solver_settings']['solver_sequence'] is not None:
                #a particular solver sequence has been requested
                #in this case all the studies input will be neglected
                if self.comsol_args['solver_settings']['solver_sequence']=='tds_elstat':
                    studies=self.tds_elstat(par_values=par_values,par_unit='')
            else:

                #1st define all studies
                for i,study in enumerate(studies):
                    #setup study and assign required solver
                    getattr(self,study)(i+1)
        
                #2nd define all solvers required for studies
                for study in studies:
                    #set up specific solver for study
                    self.sol(study)

                    #add parametric sweep if requested
                    if self.comsol_args['solver']=='parametric':
                        self.add_feature(study,'param',\
                            par_name=par_name,par_values=par_values,\
                            par_unit='',par_method=par_method)

                    self.sol_rm_feature('fcDef')

                #run all solver's with assigned studies
                for i,study in enumerate(studies):
                    #run solver
                    self.run(i+1)

        def tds_elstat(self,**kwargs):
            study='stat'
            if kwargs is not None:
                locals().update(kwargs)
                if 'initstep' not in kwargs:
                    initstep=0.01
                if 'minstep' not in kwargs:
                    minstep=1e-6
                if 'maxiter' not in kwargs:
                    if study=='stat':
                        maxiter=50
                    elif study=='time':
                        maxiter=8
                if 'damp' not in kwargs:
                    damp=0.9
                if 'stol' not in kwargs:
                    stol=1e-3
                par_values=kwargs['par_values']
                par_unit=kwargs['par_unit']
            else:
                initstep=0.01
                minstep=1e-6
                if study=='stat':
                    maxiter=50
                elif study=='time':
                    maxiter=8
                damp=0.9
                stol=1e-4
            #create a list of internal descriptors for which COMSOL will do the calculation
            if self.comsol_args['desc_method'].startswith('internal'):
                par_values_str=''
                count=0
                for p in par_values:
                    count+=1
                    par_values_str+=str(p)+' '
                    if count%5==0:
                        par_values_str+='\"\n      +\"'
                par_values_str=par_values_str[:-1]+'\"'
            else:
                par_values_str=par_values
            def ffc(number):
                number=str(number)
                number.replace('e','E')
                return number
            self.s+='    model.study().create("std1");\n'
            self.s+='    model.study("std1").create("stat2", "Stationary");\n'
            self.s+='    model.study("std1").create("stat", "Stationary");\n'
            self.s+='    model.study("std1").feature("stat2")\n'
            self.s+='         .set("activate", new String[]{"es", "off", "tds", "on", "frame:spatial1", "on"});\n'
            self.s+='    model.study("std1").feature("stat2").set("activateCoupling", new String[]{"pc1", "off", "scdc1", "off"});\n'
            self.s+='    model.study("std1").feature("stat")\n'
            self.s+='         .set("activate", new String[]{"es", "on", "tds", "off", "frame:spatial1", "on"});\n'
            self.s+='    model.study("std1").feature("stat").set("activateCoupling", new String[]{"pc1", "on", "scdc1", "off"});\n'
            self.s+='    model.study().create("std2");\n'
            self.s+='    model.study("std2").create("stat", "Stationary");\n'

            self.s+='    model.sol().create("sol1");\n'
            self.s+='    model.sol("sol1").study("std1");\n'
            self.s+='    model.sol("sol1").attach("std1");\n'
            self.s+='    model.sol("sol1").create("st1", "StudyStep");\n'
            self.s+='    model.sol("sol1").create("v1", "Variables");\n'
            self.s+='    model.sol("sol1").create("s1", "Stationary");\n'
            self.s+='    model.sol("sol1").create("su1", "StoreSolution");\n'
            self.s+='    model.sol("sol1").create("st2", "StudyStep");\n'
            self.s+='    model.sol("sol1").create("v2", "Variables");\n'
            self.s+='    model.sol("sol1").create("s2", "Stationary");\n'
            self.s+='    model.sol("sol1").feature("s1").create("p1", "Parametric");\n'
            self.s+='    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");\n'
            self.s+='    model.sol("sol1").feature("s1").create("d1", "Direct");\n'
            self.s+='    model.sol("sol1").feature("s1").feature().remove("fcDef");\n'
            self.s+='    model.sol("sol1").feature("s2").create("p1", "Parametric");\n'
            self.s+='    model.sol("sol1").feature("s2").create("fc1", "FullyCoupled");\n'
            self.s+='    model.sol("sol1").feature("s2").feature().remove("fcDef");\n'
            self.s+='    model.sol().create("sol3");\n'
            self.s+='    model.sol("sol3").study("std2");\n'
            self.s+='    model.sol("sol3").attach("std2");\n'
            self.s+='    model.sol("sol3").create("st1", "StudyStep");\n'
            self.s+='    model.sol("sol3").create("v1", "Variables");\n'
            self.s+='    model.sol("sol3").create("s1", "Stationary");\n'
            self.s+='    model.sol("sol3").feature("s1").create("fc1", "FullyCoupled");\n'
            self.s+='    model.sol("sol3").feature("s1").create("d1", "Direct");\n'
            self.s+='    model.sol("sol3").feature("s1").feature().remove("fcDef");\n'

            self.s+='    model.study("std1").feature("stat2").label("Stationary TDS");\n'
            self.s+='    model.study("std1").feature("stat2").set("useparam", true);\n'
            self.s+='    model.study("std1").feature("stat2").set("pname", new String[]{"flux_factor"});\n' #par_name
            self.s+='    model.study("std1").feature("stat2").set("plistarr", new String[]{"'+par_values_str+'"});\n'
            self.s+='    model.study("std1").feature("stat2").set("punit", new String[]{"'+par_unit+'"});\n'
            self.s+='    model.study("std1").feature("stat").label("Stationary ES");\n'
            self.s+='    model.study("std1").feature("stat").set("useparam", true);\n'
            self.s+='    model.study("std1").feature("stat").set("pname", new String[]{"PZC_factor"});\n'
            self.s+='    model.study("std1").feature("stat").set("plistarr", new String[]{"'+par_values_str+'"});\n'
            self.s+='    model.study("std1").feature("stat").set("punit", new String[]{""});\n'
            self.s+='    model.study("std2").feature("stat").set("useinitsol", true);\n'
            self.s+='    model.study("std2").feature("stat").set("initmethod", "sol");\n'
            self.s+='    model.study("std2").feature("stat").set("initstudy", "std1");\n'
            self.s+='    model.study("std2").feature("stat").set("solnum", "auto");\n'

            #study step 1
            self.s+='    model.sol("sol1").attach("std1");\n'
            self.s+='    model.sol("sol1").feature("v1").set("clistctrl", new String[]{"p1"});\n'
            self.s+='    model.sol("sol1").feature("v1").set("cname", new String[]{"flux_factor"});\n'
            self.s+='    model.sol("sol1").feature("v1").set("clist", new String[]{"'+par_values_str+'"});\n'

            #solution 1: pure tds problem, tune up flux_factor
            self.s+='    model.sol("sol1").feature("s1").set("probesel", "none");\n'
            self.s+='    model.sol("sol1").feature("s1").feature("p1").set("pname", new String[]{"flux_factor"});\n'
            self.s+='    model.sol("sol1").feature("s1").feature("p1").set("plistarr", new String[]{"'+par_values_str+'"});\n'
            self.s+='    model.sol("sol1").feature("s1").feature("p1").set("punit", new String[]{""});\n'
            self.s+='    model.sol("sol1").feature("s1").feature("fc1").set("initstep", '+ffc(initstep)+');\n'
            self.s+='    model.sol("sol1").feature("s1").feature("fc1").set("minstep", '+ffc(minstep)+');\n'
            self.s+='    model.sol("sol1").feature("s1").feature("fc1").set("maxiter", '+ffc(maxiter)+');\n'
            self.s+='    model.sol("sol1").feature("s1").feature("d1").set("linsolver", "pardiso");\n'

            #study step 2
            self.s+='    model.sol("sol1").feature("st2").set("studystep", "stat");\n'
            #initialize with solution of std1
            self.s+='    model.sol("sol1").feature("v2").set("initmethod", "sol");\n'
            self.s+='    model.sol("sol1").feature("v2").set("initsol", "sol1");\n'
            self.s+='    model.sol("sol1").feature("v2").set("solnum", "auto");\n'
            self.s+='    model.sol("sol1").feature("v2").set("notsolmethod", "sol");\n'
            self.s+='    model.sol("sol1").feature("v2").set("notsol", "sol1");\n'
            self.s+='    model.sol("sol1").feature("v2").set("notsolnum", "auto");\n'
            self.s+='    model.sol("sol1").feature("v2").set("clistctrl", new String[]{"p1"});\n'
            self.s+='    model.sol("sol1").feature("v2").set("cname", new String[]{"PZC_factor"});\n'
            self.s+='    model.sol("sol1").feature("v2").set("clist", new String[]{"'+par_values_str+'"});\n'
            self.s+='    model.sol("sol1").feature("s2").set("probesel", "none");\n'
            self.s+='    model.sol("sol1").feature("s2").feature("p1").set("pname", new String[]{"PZC_factor"});\n'
            self.s+='    model.sol("sol1").feature("s2").feature("p1").set("plistarr", new String[]{"'+par_values_str+'"});\n'
            self.s+='    model.sol("sol1").feature("s2").feature("p1").set("punit", new String[]{""});\n'
            self.run(1)
            self.s+='    model.sol("sol3").attach("std2");\n'
            self.s+='    model.sol("sol3").feature("v1").set("initmethod", "sol");\n'
            self.s+='    model.sol("sol3").feature("v1").set("initsol", "sol1");\n'
            self.s+='    model.sol("sol3").feature("v1").set("solnum", "auto");\n'
            self.s+='    model.sol("sol3").feature("s1").set("stol", "'+ffc(stol)+'");\n' #relative tolerance (default 1e-3)
            self.s+='    model.sol("sol3").feature("s1").feature("fc1").set("initstep", '+ffc(initstep)+');\n'
            self.s+='    model.sol("sol3").feature("s1").feature("fc1").set("minstep", '+ffc(minstep)+');\n'
            self.s+='    model.sol("sol3").feature("s1").feature("fc1").set("maxiter", '+ffc(maxiter)+');\n'
            self.s+='    model.sol("sol3").feature("s1").feature("d1").set("nliniterrefine", true);\n'
            self.s+='    model.result().dataset().create("dset5", "Solution");\n'
            self.s+='    model.result().dataset("dset5").set("solution", "sol3");\n'
            self.run(3)

        def stat(self,sol_index):
            self.std_index+=1
            j=self.std_index
            self.s+='    model.study().create("std'+str(j)+'");\n'
            self.s+='    model.study("std'+str(j)+'").create("stat", "Stationary");\n'
            #now we can activate here certain physics that we are interested in
#            if self.comsol_args['model_type']=='tp_dilute_species':
#                self.s+='    model.study("std'+str(j)+'").feature("stat").activate("es", true);\n'
#                self.s+='    model.study("std'+str(j)+'").feature("stat").activate("tds", true);\n'
#                self.s+='    model.study("std'+str(j)+'").feature("stat").activate("pc1", true);\n'
#                self.s+='    model.study("std'+str(j)+'").feature("stat").activate("scdc1", true);\n'

        def time(self,sol_index):
            self.std_index+=1
            j=self.std_index
            self.s+='    model.study().create("std'+str(j)+'");\n'
            self.s+='    model.study("std'+str(j)+'").create("time", "Transient");\n'
            self.s+='    model.study("std'+str(i)+'").feature("time").set("tlist", "range('+str(min(self.tp.tmesh_init))+','+str(self.tp.dt_init)+','+str(self.tp.tmax_init)+')");\n'

        def add_feature(self,study,feature,**kwargs):
            getattr(self,'add_'+feature)(study,**kwargs)

        def add_param(self,study,**kwargs):
            if kwargs is not None:
                par_name=kwargs['par_name']
                par_values=kwargs['par_values']
                par_unit=kwargs['par_unit']
                par_method=kwargs['par_method'] #either internal (auxiliary) or external parameter sweep can be used in comsol

            j=self.std_index
            i=self.sol_index

            #create a list of internal descriptors for which COMSOL will do the calculation
            if self.comsol_args['desc_method'].startswith('internal'):
                par_values_str=''
                count=0
                for p in par_values:
                    count+=1
                    par_values_str+=str(p)+' '
                    if count%5==0:
                        par_values_str+='\"\n      +\"'
                par_values_str=par_values_str[:-1]+'\"'
            else:
                par_values_str=par_values

            self.s+='    model.sol("sol'+str(i)+'").feature("s1").create("p1", "Parametric");\n'
            self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("p1").active(true);\n'
            self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("p1").set("pname", new String[]{\"'+par_name+'\"});\n'
            self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("p1").set("plistarr", new String[]{\"'+par_values_str+'\"});\n'
            self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("p1").set("punit", new String[]{\"'+par_unit+'\"});\n'
            self.s+='    model.sol("sol'+str(i)+'").feature("s1").set("probesel", "none");\n'
            if par_method == 'external':

                self.s+='    model.study("std'+str(j)+'").create("param", "Parametric");\n'
                self.s+='    model.study("std'+str(j)+'").feature("param").set("pname", new String[]{\"'+par_name+'\"});\n'
                self.s+='    model.study("std'+str(j)+'").feature("param").set("plistarr", new String[]{\"'+par_values_str+'\"});\n'
                self.s+='    model.study("std'+str(j)+'").feature("param").set("punit", new String[]{\"'+par_unit+'\"});\n'
                #set "param" parametric sweep instance as controling element for sweep
                self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("p1").set("control", "param");\n'
                #if comsol crashes during sweep, store empty solutions or stop
                self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("p1").set("ponerror", "empty");\n'
                if 'internal' in self.comsol_args['desc_method']:
                    if self.tp.comsol_args['desc_method'].split('-')[1]=='reinit':
                        self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("p1").set("pcontinuationmode", "no");\n'
                    else:
                        self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("p1").set("preusesol", "yes");\n'
                else:
                    self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("p1").set("pcontinuationmode", "no");\n'
            elif par_method == 'internal':
                #let parameter sweep to be controlled by study
                self.s+='    model.sol("sol'+str(i)+'").feature("s1").set("control", "'+study+'");\n'
                self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("p1").set("control", "'+study+'");\n'
                #pDef is the parametric sweep inside the solver (control = user) in the above line
#                self.s+='    model.sol("sol'+str(i)+'").feature("s1").create("pDef", "Parametric");\n'
#                self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("pDef").set("pname", new String[]{\"'+par_name+'\"});\n'
#                self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("pDef").set("plistarr", new String[]{\"'+par_values_str+'\"});\n'
#                self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("pDef").set("punit", new String[]{\"'+par_unit+'\"});\n'
#                self.s+='    model.sol("sol'+str(i)+'").feature("s1").feature("pDef").set("uselsqdata", false);\n'
                self.s+='    model.study("std'+str(j)+'").feature("'+study+'").set("useparam", true);\n'
                self.s+='    model.study("std'+str(j)+'").feature("'+study+'").set("pname", new String[]{\"'+par_name+'\"});\n'
                self.s+='    model.study("std'+str(j)+'").feature("'+study+'").set("plistarr", new String[]{\"'+par_values_str+'\"});\n'
                self.s+='    model.study("std'+str(j)+'").feature("'+study+'").set("punit", new String[]{\"'+par_unit+'\"});\n'

            #general settings
            #variables needed for parametric solver
            self.s+='    model.sol("sol'+str(i)+'").feature("v1").set("clistctrl", new String[]{"p1"});\n'
            self.s+='    model.sol("sol'+str(i)+'").feature("v1").set("cname", new String[]{\"'+par_name+'\"});\n'
            self.s+='    model.sol("sol'+str(i)+'").feature("v1").set("clist", new String[]{\"'+par_values_str+'\"});\n'


        def run(self,sol_index):
            self.s+='    model.sol("sol'+str(sol_index)+'").runAll();\n'

        def sol_rm_feature(self,feature):
            self.s+='    model.sol("sol'+str(self.sol_index)+'").feature("s1").feature().remove("'+feature+'");\n'

        def sol(self,study,**kwargs):
            if kwargs is not None:
                locals().update(kwargs)
                if 'initstep' not in kwargs:
                    initstep=0.01
                if 'minstep' not in kwargs:
                    minstep=1e-6
                if 'maxiter' not in kwargs:
                    if study=='stat':
                        maxiter=50
                    elif study=='time':
                        maxiter=8
                if 'damp' not in kwargs:
                    damp=0.9
                if 'stol' not in kwargs:
                    stol=1e-3
            else:
                initstep=0.01
                minstep=1e-6
                if study=='stat':
                    maxiter=50
                elif study=='time':
                    maxiter=8
                damp=0.9
                stol=1e-4
            self.sol_index+=1
            j=self.sol_index
            i=self.std_index
            def ffc(number):
                number=str(number)
                number.replace('e','E')
                return number
            #create solver
            self.s+='    model.sol().create("sol'+str(j)+'");\n'
            self.s+='    model.sol("sol'+str(j)+'").study("std'+str(i)+'");\n'
            self.s+='    model.sol("sol'+str(j)+'").attach("std'+str(i)+'");\n'
            #setup solver
            self.s+='    model.sol("sol'+str(j)+'").create("st1", "StudyStep");\n'
            self.s+='    model.sol("sol'+str(j)+'").create("v1", "Variables");\n'
            if study=='stat':
                self.s+='    model.sol("sol'+str(j)+'").create("s1", "Stationary");\n'
                self.s+='    model.sol("sol'+str(j)+'").feature("s1").create("d1", "Direct");\n'
                ##fc1
                self.s+='    model.sol("sol'+str(j)+'").feature("s1").create("fc1", "FullyCoupled");\n'
                self.s+='    model.sol("sol'+str(j)+'").feature("s1").feature("fc1").set("initstep", '+ffc(initstep)+');\n'
                self.s+='    model.sol("sol'+str(j)+'").feature("s1").feature("fc1").set("minstep", '+ffc(minstep)+');\n'
                self.s+='    model.sol("sol'+str(j)+'").feature("s1").feature("fc1").set("maxiter", '+ffc(maxiter)+');\n'
                self.s+='    model.sol("sol'+str(j)+'").feature("s1").set("stol", "'+ffc(stol)+'");\n' #relative tolerance (default 1e-3)
                if 'direct' in self.comsol_args['solver_settings'] and \
                        'nliniterrefine' in self.comsol_args['solver_settings']['direct'] and\
                        self.comsol_args['solver_settings']['direct']['nliniterrefine']:
                    self.s+='    model.sol("sol'+str(j)+'").feature("s1").feature("d1").set("nliniterrefine", true);\n'
                ##end fc1
            elif study=='time':
                self.s+='    model.sol("sol'+str(j)+'").feature("v1").set("clist", new String[]{"'+ffc(min(self.tp.tmesh_init))+','+ffc(self.tp.dt_init)+', '+ffc(self.tp.tmax_init)+'"});\n'
                ##t1 -- Time
                self.s+='    model.sol("sol'+str(j)+'").create("t1", "Time");\n'
                self.s+='    model.sol("sol'+str(j)+'").feature("t1").create("fc1", "FullyCoupled");\n'
                self.s+='    model.sol("sol'+str(j)+'").feature("t1").create("d1", "Direct");\n'
                self.s+='    model.sol("sol'+str(j)+'").feature("t1").feature().remove("fcDef");\n'
                self.s+='    model.sol("sol'+str(j)+'").feature("t1").set("tlist", "range('+ffc(min(self.tp.tmesh_init))+','+ffc(self.tp.dt_init)+','+ffc(self.tp.tmax_init)+')");\n'
                self.s+='    model.sol("sol'+str(j)+'").feature("t1").set("maxorder", 2);\n'
                  ##fc1
                self.s+='    model.sol("sol'+str(j)+'").feature("t1").feature("fc1").set("maxiter", '+ffc(maxiter)+');\n'
                self.s+='    model.sol("sol'+str(j)+'").feature("t1").feature("fc1").set("damp", '+ffc(damp)+');\n'
                self.s+='    model.sol("sol'+str(j)+'").feature("t1").feature("fc1").set("jtech", "once");\n'
                  ##end fc1
                ##end t1 -- Time

        def get(self):
            return self.s

    class probe():
        """Defining probes"""
        #TODO: work on this more
        def __init__(self):
            self.s=''
            self.s+='/*\n'
            self.s+=' *PROBES\n'
            self.s+=' */\n'
            #create domain point probe group pdom1
            group='pdom1'
            self.s+='    model.component("comp1").probe().create("'+group+'", "DomainPoint");\n'

            self.add_probe(group='pdom1',index=2,expr='phiM',unit='V',descr='Electric Potential',coord=0)

        def add_probe(self,group='',index=1,expr='',unit='',descr='',active=True,snap_closest=True,coord=0):
            self.s+='    model.component("comp1").probe("'+group+'").create("ppb'+str(index)+'", "PointExpr");\n'
            self.s+='    model.component("comp1").probe("pdom1").setIndex("coords1", "'+str(coord)+'", 0);\n'
            if snap_closest:
                self.s+='    model.component("comp1").probe("'+group+'").set("bndsnap1", true);\n'
            if not active:
                self.s+='    model.component("comp1").probe("'+group+'").feature("ppb'+str(index)+'").active(false);\n'
            #create table
            self.s+='    model.result().table().create("tbl'+str(index)+'", "Table");\n'
            #assign table
            self.s+='    model.component("comp1").probe("'+group+'").feature("ppb'+str(index)+'").set("table", "tbl'+str(index)+'");\n'
            #label probe table
            self.s+='    model.result().table("tbl'+str(index)+'").label("Probe Table '+str(index)+'");\n'
            #assign window
            self.s+='    model.component("comp1").probe("'+group+'").feature("ppb'+str(index)+'").set("window", "window'+str(index)+'");\n'
            #content
            self.s+='    model.component("comp1").probe("'+group+'").feature("ppb'+str(index)+'").set("expr", "'+expr+'");\n'
            self.s+='    model.component("comp1").probe("'+group+'").feature("ppb'+str(index)+'").set("unit", "'+unit+'");\n'
            self.s+='    model.component("comp1").probe("'+group+'").feature("ppb'+str(index)+'").set("descr", "'+descr+'");\n'

        def get(self):
            return self.s

    class mesh():
        """Defining the comsol mesh"""
        def __init__(self):
            self.s=''
            self.s+='/*\n'
            self.s+=' *MESH\n'
            self.s+=' */\n'
            self.s+='    model.component("comp1").mesh("mesh1").create("edg1", "Edge");\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").create("size1", "Size");\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").create("size2", "Size");\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").selection().geom("geom1", 0);\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").selection().set(new int[]{1, 2});\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hauto", 1);\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("custom", "on");\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hmax", "L_cell/grid_factor_domain");\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hmaxactive", true);\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hnarrow", 5000);\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hnarrowactive", false);\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("custom", "on");\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("hmax", "lambdaD/grid_factor_bound");\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("hmaxactive", true);\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("hmin", 3.26E-12);\n'
            self.s+='    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("hminactive", false);\n'
            self.s+='    model.component("comp1").mesh("mesh1").run();\n'
        def get(self):
            return self.s

    class multiphysics():
        """MULTIPHYSICS: tds-es COUPLING (potential and space charge density)"""
        def __init__(self,method="tds-es"):
            self.s=''
            self.s+='/*\n'
            self.s+=' *MULTIPHYSICS\n'
            self.s+=' */\n'
            if method=="tds-es":
                self.s+='    model.component("comp1").multiphysics().create("pc1", "PotentialCoupling", 1);\n'
                self.s+='    model.component("comp1").multiphysics("pc1").selection().all();\n'
                self.s+='    model.component("comp1").multiphysics().create("scdc1", "SpaceChargeDensityCoupling", 1);\n'
                self.s+='    model.component("comp1").multiphysics("scdc1").selection().all();\n'
                self.s+='    model.component("comp1").physics("es").feature("sfcd1").set("rhoqs", "rho_s");\n'
        def get(self):
            return self.s

    class physics():
        def __init__(self,tp=None,methods=[],comsol_args={}):
            if tp is None:
                return
            else:
                self.tp=tp
            self.comsol_args=comsol_args
            self.version=comsol_args['bin_version']
            self.s=''
            self.s+='/*\n'
            self.s+=' *PHYSICS\n'
            self.s+=' */\n'
            for method in methods:
                getattr(self,method)()
        def ge(self):
            if 'global_equations' not in self.comsol_args:
                return self.s
            else:
                self.s+='/*\n'
                self.s+=' *GLOBAL EQUATIONS\n'
                self.s+=' */\n'

                self.s+='    model.component("comp1").physics().create("ge", "GlobalEquations", "geom1");\n'
                for ge_name in self.comsol_args['global_equations']:
                    for func in self.comsol_args['global_equations'][ge_name]:
                        if type(self.comsol_args['global_equations'][ge_name][func])==str:
                            self.s+='    model.component("comp1").physics("ge").feature("'+ge_name+'").'+func+'(\"'+self.comsol_args['global_equations'][ge_name][func]+'\");\n'
                        else:
                            if 'DependentVariableQuantity' not in self.comsol_args['global_equations'][ge_name][func]:
                                self.s+='    model.component("comp1").physics("ge").feature("'+ge_name+'").'+func+'(\"DependentVariableQuantity\", \"none\");\n'
                                if 'CustomDependentVariableUnit' not in self.comsol_args['global_equations'][ge_name][func]:
                                    self.tp.logger.warning('No unit given for dependent variable in global COMSOL equation, could lead to error')
                            if 'SourceTermQuantity' not in self.comsol_args['global_equations'][ge_name][func]:
                                self.s+='    model.component("comp1").physics("ge").feature("'+ge_name+'").'+func+'(\"SourceTermQuantity\", \"none\");\n'
                                if 'SourceTermVariableUnit' not in self.comsol_args['global_equations'][ge_name][func]:
                                    self.tp.logger.warning('No unit given for source term in global COMSOL equation, could lead to error')
                            for key in self.comsol_args['global_equations'][ge_name][func]:
                                self.s+='    model.component("comp1").physics("ge").feature("'+ge_name+'").'+func+'(\"'+key+'\", \"'+self.comsol_args['global_equations'][ge_name][func][key]+'\");\n'

        def es(self):
            """Electrostatics"""
            self.s+='    model.component("comp1").physics().create("es", "Electrostatics", "geom1");\n'
            self.s+='    model.component("comp1").physics("es").field("electricpotential").field("phi");\n'
            self.s+='    model.component("comp1").physics("es").create("sfcd1", "SurfaceChargeDensity", 0);\n'
            self.s+='    model.component("comp1").physics("es").feature("sfcd1").selection().set(new int[]{1});\n'
            self.s+='    model.component("comp1").physics("es").create("pot1", "ElectricPotential", 0);\n'
            self.s+='    model.component("comp1").physics("es").feature("pot1").selection().set(new int[]{2});\n'
            #only for dirichlet:
            #self.s+='    model.component("comp1").physics("es").create("pot2", "ElectricPotential", 0);\n'
            #self.s+='    model.component("comp1").physics("es").feature("pot2").selection().set(new int[]{1});\n'
            self.s+='    model.component("comp1").physics("es").create("df1", "DisplacementField", 0);\n'
            self.s+='    model.component("comp1").physics("es").feature("df1").selection().set(new int[]{1});\n'
            #floating requires one of: AC/DC Module, MEMS Module, Plasma Module, Acoustics Module, Structural Mechanics Module, Semiconductor Module
#            self.s+='    model.component("comp1").physics("es").create("fp1", "FloatingPotential", 0);\n'
            self.s+='    model.component("comp1").physics("es").feature("ccn1").set("epsilonr_mat", "userdef");\n'
            self.s+='    model.component("comp1").physics("es").feature("ccn1").set("epsilonr", new String[][]{{"eps_r"}, {"0"}, {"0"}, {"0"}, {"eps_r"}, {"0"}, {"0"}, {"0"}, {"eps_r"}});\n'
                #for dirichlet BC's
                #self.s+='    model.component("comp1").physics("es").feature("pot2").set("V0", "phiM");\n'
                #self.s+='    model.component("comp1").physics("es").feature("pot2").active(false);\n'
                #end dirichlet
            self.s+='    model.component("comp1").physics("es").feature("df1").active(false);\n'
            if self.version>5.31:
#                self.s+='    model.component("comp1").physics("es").prop("PortSweepSettings").set("PortParamName", "PortName");\n'
                self.s+='    model.component("comp1").physics("es").feature("ccn1").set("minput_numberdensity", 0);\n'
                self.s+='    model.component("comp1").physics("es").feature("ccn1").set("minput_temperature_src", "userdef");\n'
        def tds(self):
            """Transport of Dilute Species"""
            self.s+='    model.component("comp1").physics().create("tds", "DilutedSpecies", "geom1");\n'
            self.s+='    model.component("comp1").physics("tds").field("concentration").field("cp1");\n'
            conc_str=""
            for i in range(len(self.tp.species)):
                conc_str+="\"cp"+str(i+1)+"\""
                if i != len(self.tp.species)-1:
                    conc_str+=", "
            self.s+='    model.component("comp1").physics("tds").field("concentration").component(new String[]{'+conc_str+'});\n'
            #if self.version>5.31:
            #    self.s+='    model.component("comp1").physics("tds").create("fl1", "FluxBoundary", 0);\n'
            self.s+='    model.component("comp1").physics("tds").create("fl1", "Fluxes", 0);\n'
            self.s+='    model.component("comp1").physics("tds").feature("fl1").selection().set(1);\n'
            self.s+='    model.component("comp1").physics("tds").create("conc1", "Concentration", 0);\n'
            self.s+='    model.component("comp1").physics("tds").feature("conc1").selection().set(2);\n'


            #### REACTIONS

            if self.tp.use_electrolyte_reactions:
                self.s+='    model.component("comp1").physics("tds").create("reac1", "Reactions", 1);\n'
                self.s+='    model.component("comp1").physics("tds").feature("reac1").selection().all();\n'
            self.s+='    model.component("comp1").physics("tds").prop("ShapeProperty").set("order_concentration", 2);\n'
            if self.tp.use_convection or self.tp.use_mpb:
                self.s+='    model.component("comp1").physics("tds").prop("TransportMechanism").set("Convection", true);\n'
            else:
                self.s+='    model.component("comp1").physics("tds").prop("TransportMechanism").set("Convection", false);\n'
            if self.tp.use_migration:
                self.s+='    model.component("comp1").physics("tds").prop("TransportMechanism").set("Migration", true);\n'
            else:
                self.s+='    model.component("comp1").physics("tds").prop("TransportMechanism").set("Migration", false);\n'
            if self.version>5.31:
                self.s+='    model.component("comp1").physics("tds").prop("TransportMechanism").set("MassTransferInPorousMedia", false);\n'

            #add diffusion coefficients to diffusion model
            self.s+='/*\n'
            self.s+=' *-> ADD diffusion coefficients\n'
            self.s+=' */\n'
            charge_str=""
            for i in range(len(self.tp.species)):
                charge_str+="{\"Z"+str(i+1)+"\"}"
                if i != len(self.tp.species)-1:
                    charge_str+=", "

            self.s+='    model.component("comp1").physics("tds").feature("cdm1").set("z", new String[][]{'+charge_str+'});\n'
            self.s+='    model.component("comp1").physics("tds").feature("cdm1").set("minput_temperature", "T");\n'
            if self.version>5.31:
                self.s+='    model.component("comp1").physics("tds").feature("cdm1").set("minput_temperature_src", "userdef");\n'

            for i in range(len(self.tp.species)):
                self.s+='    model.component("comp1").physics("tds").feature("cdm1").set("D_cp'+str(i+1)+'", new String[][]{{"D'+str(i+1)+'"}, {"0"}, {"0"}, {"0"}, {"D'+str(i+1)+'"}, {"0"}, {"0"}, {"0"}, {"D'+str(i+1)+'"}});\n'
    
            self.s+='/*\n'
            self.s+=' *-> ADD initial concentrations\n'
            self.s+=' */\n'
            #add initial concentrations to diffusion model
            ci_str=""
            for i in range(len(self.tp.species)):
                ci_str+="{\"ci"+str(i+1)+"\"}"
                if i != len(self.tp.species)-1:
                    ci_str+=", "
            self.s+='    model.component("comp1").physics("tds").feature("init1").set("initc", new String[][]{'+ci_str+'});\n'

            self.s+='/*\n'
            self.s+=' *-> ADD BOUNDARIES\n'
            self.s+=' */\n'
            #boundaries
            b_str=""
            for i in range(len(self.tp.species)):
                b_str+="{1}"
                if i != len(self.tp.species)-1:
                    b_str+=", "
            self.s+='    model.component("comp1").physics("tds").feature("fl1").set("species", new int[][]{'+b_str+'});\n'

            self.s+='/*\n'
            self.s+=' *-> ADD FLUXES\n'
            self.s+=' */\n'
            #fluxes
            f_str=""
            for i in range(len(self.tp.species)):
                f_str+="{\"j"+str(i+1)+"\"}"
                if i != len(self.tp.species)-1:
                    f_str+=", "

            #if self.version>5.31:
            #    self.s+='    model.component("comp1").physics("tds").feature("fl1").set("J0", new String[][]{'+f_str+'});\n'
            #else:
            self.s+='    model.component("comp1").physics("tds").feature("fl1").set("N0", new String[][]{'+f_str+'});\n'
            self.s+='    model.component("comp1").physics("tds").feature("conc1").set("species", new int[][]{'+b_str+'});\n'
            self.s+='    model.component("comp1").physics("tds").feature("conc1").set("c0", new String[][]{'+ci_str+'});\n'
            if self.version>5.31:
                self.s+='    model.component("comp1").physics("tds").feature("nflx1").set("IncludeConvection", true);\n'

            #global constraint
#            self.s+='    model.component("comp1").physics("tds").create("gconstr1", "GlobalConstraint", -1);\n'
#            self.s+='    model.component("comp1").physics("tds").feature("gconstr1").set("constraintExpression", "intop2(cm)-(cref*L)");\n'
#            self.s+='    model.component("comp1").physics("tds").feature("gconstr1").active(false);\n'

            def get_rates():
                #get list of all species names:
                species_names=[sp for sp in self.tp.species]
                rates=[""]*self.tp.nspecies
                jj=0
                for r in self.tp.electrolyte_reactions:
                    reaction=self.tp.electrolyte_reactions[r]
                    if 'rates' not in reaction:
                        continue
                    else:
                        jj+=1
                    #rates of educts:
                    for reactant in reaction['reaction'][0]:
                        if reactant not in self.tp.species:
                            #prod+="conc_std"
                            continue
                           # continue
                        else:
                            k=species_names.index(reactant)
                    #    rates[k]=""
                        prod=""
                        for reactant2 in reaction['reaction'][0]:
                            if reactant2 not in self.tp.species:
                                prod+="conc_std*"
                                #continue
                            else:
                                k2=species_names.index(reactant2)
                                prod+="ap"+str(k2+1)+"*"
                        rates[k]+="-"+prod+"k"+str(jj)+"f"
                        prod=""
                        for reactant2 in reaction['reaction'][1]:
                            if reactant2 not in self.tp.species:
                                prod+="conc_std*"
                                #pass
                                #continue
                            else:
                                k2=species_names.index(reactant2)
                                prod+="ap"+str(k2+1)+"*"
                        rates[k]+="+"+prod+"k"+str(jj)+"r"
                    #rates of products
                    for reactant in reaction['reaction'][1]:
                        if reactant not in self.tp.species:
                            #prod+="conc_std"
                            continue
                            #continue
                        else:
                            k=species_names.index(reactant)
                     #   rates[k]=""
                        prod=""
                        for reactant2 in reaction['reaction'][0]:
                            if reactant2 not in self.tp.species:
                                prod+="conc_std*"
                                #continue
                            else:
                                k2=species_names.index(reactant2)
                                prod+="ap"+str(k2+1)+"*"
                        rates[k]+="+"+prod+"k"+str(jj)+"f"
                        prod=""
                        for reactant2 in reaction['reaction'][1]:
                            if reactant2 not in self.tp.species:
                                prod+="conc_std*"
                                #continue
                            else:
                                k2=species_names.index(reactant2)
                                prod+="ap"+str(k2+1)+"*"
                        rates[k]+="-"+prod+"k"+str(jj)+"r"

                rates_new=[]
                i=0
                for rate in rates:
                    i+=1
                    if len(rate)>0:
                        rates_new.append([i,rate])
                    else:
                        rates_new.append([i,'0.0'])
                return rates_new
            if self.tp.use_electrolyte_reactions:
                self.s+='/*\n'
                self.s+=' *REACTIONS\n'
                self.s+=' */\n'
                rate_str=get_rates()
                species_names=[sp for sp in self.tp.species]
                for index,rate in rate_str:
                    if 'reactions' in self.comsol_args['solver_settings']['ramp']['names']:
                        self.s+='    model.component("comp1").physics("tds").feature("reac1").set("R_cp'+str(index)+'", "flux_factor*('+rate+')");\n'
                    else:
                        self.s+='    model.component("comp1").physics("tds").feature("reac1").set("R_cp'+str(index)+'", "'+rate+'");\n'
                        
                self.s+='    model.component("comp1").physics("tds").feature("reac1").active(true);\n'
            if self.tp.use_convection or self.tp.use_mpb:
                self.s+='    model.component("comp1").physics("tds").prop("AdvancedSettings").set("ConvectiveTerm", "cons");\n'
            if self.tp.use_mpb:
                #we have to define a fake convective term defined separately for each species:
                species_names=[sp for sp in self.tp.species]
                for sp in self.tp.species:
                    index=species_names.index(sp)+1
                    self.s+='    model.component("comp1").physics("tds").feature("cdm1").featureInfo("info").set("tds.cflux_cp'+str(index)+'x", new String[]{"cp'+str(index)+'*tds.u'+str(index)+'"});\n'
                for sp in self.tp.species:
                    index=species_names.index(sp)+1
                    #if self.version>5.31:
                    #    self.s+='    model.component("comp1").physics("tds").feature("cdm1").featureInfo("info").set("tds.Res_cp'+str(index)+'", new String[]{"d(-cp'+str(index)+'*tds.z_cp'+str(index)+'*tds.um_cp'+str(index)+'xx*F_const*d(tds.V,x),x)+tds.u'+str(index)+'*cp'+str(index)+'x-tds.R_cp'+str(index)+'"});\n'
                    #else:
                    self.s+='    model.component("comp1").physics("tds").feature("cdm1").featureInfo("info").set("tds.Res_cp'+str(index)+'", new String[]{"-tds.D_cp'+str(index)+'xx*cp'+str(index)+'xx+d(cp'+str(index)+'*(tds.u'+str(index)+'-tds.z_cp'+str(index)+'*tds.um_cp'+str(index)+'xx*F_const*d(tds.V,x)),x)-tds.R_cp'+str(index)+'"});\n'
                ii=0
                for sp in self.tp.species:
                    ii+=1
                    index=species_names.index(sp)+1
                    #$13 ???
                    #if self.version>5.31:
                    #    self.s+='    model.component("comp1").physics("tds").feature("cdm1").featureInfo("info").set("root.comp1.tds.cdm1.weak$'+str(len(species_names)*2+ii)+'", new String[]{"-tds.u'+str(index)+'*cp'+str(index)+'x*test(cp'+str(index)+')*(isScalingSystemDomain==0)", "2"});\n'
                    #else:
                    ####
                    #!!!!!!!!!!!!!!!!!!!!!
                    #VERY IMPORTANT. In the new version (tested 5.5) the indices of the weak expressions changed! With incorrect indices COMSOL crashes!
                    if self.version>5.31:
                        self.s+='    model.component("comp1").physics("tds").feature("cdm1").featureInfo("info").set("root.comp1.tds.cdm1.weak$'+str(len(species_names)*2+ii*2-1)+'", new String[]{"cp'+str(index)+'*tds.u'+str(index)+'*test(cp'+str(index)+'x)*(isScalingSystemDomain==0)", "4"});\n'
                    else:
                        self.s+='    model.component("comp1").physics("tds").feature("cdm1").featureInfo("info").set("root.comp1.tds.cdm1.weak$'+str(len(species_names)*2+ii)+'", new String[]{"cp'+str(index)+'*tds.u'+str(index)+'*test(cp'+str(index)+'x)*(isScalingSystemDomain==0)", "4"});\n'
                    #!!!!!!!!!!!!!!!!!!!!!
                    ####
            u=''
            if self.tp.use_convection:
                u+=str(self.tp.system['flow rate'])
#                    u+="-N_A_const*a^3*(D"+inx+"*cp"+inx+"x)/(1-phi_zero)"
#                if self.tp.system['MPB']['species']=='all':
#                    #this model is for z:z electrolyte
#                    if len(u)>0:
#                        u+='+'
#                    #get index of positive and negative ion:
#                    inx_cat=str([i+1 for i,sp in enumerate(self.tp.species) if self.tp.species[sp]['charge']>0][0])
#                    inx_an=str([i+1 for i,sp in enumerate(self.tp.species) if self.tp.species[sp]['charge']<0][0])
#                    #MPB model due to Kilic, https://www.ncbi.nlm.nih.gov/pubmed/17358343
#                    u+="-N_A_const*a^3*(D"+inx_cat+"*cp"+inx_cat+"x+D"+inx_an+"*cp"+inx_an+"x)/(1-a^3*cp"+inx_cat+"*N_A_const-a^3*cp"+inx_an+"*N_A_const)"
#                else:
#                    #this model corrects for finite size of only the selected species (should be the counter-ions)
#                    inx=str([i+1 for i,sp in enumerate(self.tp.species) if sp == self.tp.system['MPB']['species']][0])
#                    u+="-N_A_const*a^3*(D"+inx+"*cp"+inx+"x)/(1-a^3*cp"+inx+"*N_A_const)"
            if len(u)>0:
                self.s+='    model.component("comp1").physics("tds").feature("cdm1").set("u", new String[][]{{"'+u+'"}, {"0"}, {"0"}});\n'

        def get(self):
            return self.s


    class cpl():
        """defines integrations"""
        def __init__(self,index=1,geo='d1'):
            self.s=''
            self.s+='/*\n'
            self.s+=' *CPL\n'
            self.s+=' */\n'
            self.s=''
            self.s+='    model.component("comp1").cpl().create("intop'+str(index)+'", "Integration");\n'
            if geo.startswith('b'):
                self.s+='    model.component("comp1").cpl("intop'+str(index)+'").selection().geom("geom1",0);\n'
            self.s+='    model.component("comp1").cpl("intop'+str(index)+'").selection().set(new int[]{'+geo[-1]+'});\n'
        def get(self):
            return self.s

    class variables():
        def __init__(self,tp=None,index=1,geo='d1',comsol_args={}):

            self.comsol_args=comsol_args
            self.s=''

            if tp is None:
                return
            else:
                self.tp=tp

            self.index=index

            self.s+='/*\n'
            self.s+=' *VARIABLES '+str(index)+' -- '+geo+'\n'
            self.s+=' */\n'

            self.s+='    model.component("comp1").variable().create("var'+str(index)+'");\n'
            if geo.startswith('b'):
                #additional electrode boundary variables from input
                self.s+='    model.component("comp1").variable("var'+str(index)+'").selection().geom("geom1", 0);\n'
                self.s+='    model.component("comp1").variable("var'+str(index)+'").selection().set(new int[]{'+geo[-1]+'});\n'
                for pa in self.comsol_args['boundary_variables']:
                    pa_val=self.comsol_args['boundary_variables'][pa][0]
                    pa_des=self.comsol_args['boundary_variables'][pa][1]
                    matches=re.findall('\[\[(.*?)\]\]',pa_val,re.DOTALL)
                    for match in matches:
                        pa_val=pa_val.replace('[['+match+']]',str([ii+1 for ii,spp in enumerate(self.tp.species) if spp==match][0]))
                    self.set(pa, pa_val, pa_des)
            else:
                #additional global variables from input
#                self.s+='    model.component("comp1").variable("var'+str(index)+'").selection().set(new int[]{'+geo[-1]+'});\n'
                for pa in self.comsol_args['global_variables']:
                    pa_val=self.comsol_args['global_variables'][pa][0]
                    pa_des=self.comsol_args['global_variables'][pa][1]
                    matches=re.findall('\[\[(.*?)\]\]',pa_val,re.DOTALL)
                    for match in matches:
                        pa_val=pa_val.replace('[['+match+']]',str([ii+1 for ii,spp in enumerate(self.tp.species) if spp==match][0]))
                    self.set(pa, pa_val, pa_des)

            ##SOME default variables
            if geo.startswith('d'):
                self.set("rho_s", "((phiM-phiPZC)-phi)*CS", "Surface charge density")

                #### FLUXES

                for i,sp in enumerate(self.tp.species):
#                    if self.tp.system['kinetics']=='rate-equation':
                    if type(self.tp.species[sp]['flux'])==str:
                        mod_str=''
                        #flux is given as equation, have to replace concentrations with correct number here
                        string=self.tp.species[sp]['flux']
                        string.strip() #remove white spaces
                        string.replace('\n','')
                        string.replace('\t','')
                        iss=0
                        matches=re.findall('\[\[(.*?)\]\]',string,re.DOTALL)
                        for match in matches:
                            string=string.replace('[['+match+']]','cp'+str([ii+1 for ii,spp in enumerate(self.tp.species) if spp==match][0]))
                        self.tp.species[sp]['flux']=string
                        self.set("j{}".format(i+1), 'RF*flux_factor*('+self.tp.species[sp]['flux']+')', "{} flux".format(self.tp.species[sp]['name']))
            if geo.startswith('d'):#else:
                # - local charge density
                ldc_str='F_const*('
                for i in range(len(self.tp.species)):
                    if i!=0:
                        ldc_str+='+'
                    ldc_str+='cp'+str(i+1)+'*tds.z_cp'+str(i+1)
                ldc_str+=')'
                self.set("rho_charge",ldc_str,"Local Charge Density")
                # - conductivity
                cond_str='F_const^2*('
                for i in range(len(self.tp.species)):
                    if i!=0:
                        cond_str+='+'
                    cond_str+='tds.z_cp'+str(i+1)+'^2*cp'+str(i+1)+'*tds.um_cp'+str(i+1)+'xx'
                cond_str+=')'
                self.set("rho_c",cond_str,"Electrolyte conductivity")
                # - electrolyte current density
                cur_str='F_const*('
                for i in range(len(self.tp.species)):
                    if i!=0:
                        cur_str+='+'
                    cur_str+='tds.tflux_cp'+str(i+1)+'x*tds.z_cp'+str(i+1)
                cur_str+=')'
                dflux_charge_str='F_const*('
                for i in range(len(self.tp.species)):
                    if i!=0:
                        dflux_charge_str+='+'
                    dflux_charge_str+='tds.z_cp'+str(i+1)+'*tds.dflux_cp'+str(i+1)+'x'
                dflux_charge_str+=')'
                self.set("i_el",cur_str,"Electrolyte Current Density")
                self.set("delta_phi_iRx","-i_el/rho_c","Ohmic Loss Derivative")
                self.set("delta_phi_iR","intop2(delta_phi_iRx*(x<=dest(x)))","Spatially Dependent Ohmic Loss")
                self.set("dflux_times_charge",dflux_charge_str,"Diffusion Flux times charge")
                self.set("delta_phi_diffx","dflux_times_charge/rho_c","Diffusion Flux Derivative")
                self.set("delta_phi_diff","intop2(delta_phi_diffx*(x<=dest(x)))","Spatially Dependent Diffusion Loss")
                self.set("delta_phi_diff_inf","comp1.at0("+str(self.tp.system['boundary thickness'])+",delta_phi_diff)","Total Diffusion drop over the whole cell")
                self.set("delta_phi_iR_inf","comp1.at0("+str(self.tp.system['boundary thickness'])+",delta_phi_iR)","Total iR drop over the whole cell")
                self.set("delta_phi_inf","comp1.at0("+str(self.tp.system['boundary thickness'])+",phi)-comp1.at0(0,phi)","Total potential drop over the whole cell")
                self.set("delta_phi_inf_min_iR","delta_phi_inf-delta_phi_iR_inf","Total potential drop over the whole cell minux iR drop")
                if self.tp.use_mpb:
                    #volume fraction phi_0 = NA*sum_i a_i c_i
                    phi_zero='N_A_const*('
                    phi_zero_grad='N_A_const*('
                    species_names=[sp for sp in self.tp.species]
                    for sp in self.tp.species:
                        inx=str(species_names.index(sp)+1)
                        if 'MPB_radius' in self.tp.species[sp]:
                            self.set("a"+inx,str(self.tp.species[sp]['MPB_radius'])+"[m]","Ion size in MPB Model")
                            phi_zero+="+a"+inx+"^3*cp"+inx
                            phi_zero_grad+="+a"+inx+"^3*cp"+inx+"x"
                        if True: #'MPB_radius' in self.tp.species[sp]: #abs(self.tp.species[sp]['charge'])>0:
                            self.set('tds.u'+str(inx),'-D'+str(inx)+'*phi_zero_grad/(1-phi_zero)','Ion size contribution to ion flux')
                        else:
                            self.set('tds.u'+str(inx),'0.0','Ion size contribution to ion flux')
                        self.set('ap'+str(inx),'cp'+str(inx)+'*gamma','Activity of species {}'.format(sp))
                    phi_zero+=')'
                    phi_zero_grad+=')'
                    self.set("phi_zero",phi_zero,'Volume fraction of space occupied by ions in MPB Model')
                    self.set("phi_zero_grad",phi_zero_grad,'Gradient of the Volume fraction of space occupied by ions in MPB Model')
                    self.set('gamma','1./(1.-phi_zero)','Activitity Coefficient')
                else:
                    self.set('gamma','1','Activity Coefficient')
                ##REACTION RATES
                i=0
                if self.tp.use_electrolyte_reactions: 
                    for reaction in self.tp.electrolyte_reactions:
                        if 'rates' in self.tp.electrolyte_reactions[reaction]:
                            i+=1
                            educts=''
                            products=''
                            unit='1/s'
                            unit_f=unit+'*m^3/mol'*(len(self.tp.electrolyte_reactions[reaction]['reaction'][0])-1)
                            unit_r=unit+'*m^3/mol'*(len(self.tp.electrolyte_reactions[reaction]['reaction'][1])-1)
                            for reactant in self.tp.electrolyte_reactions[reaction]['reaction'][0]:
                                educts+=reactant
                                if reactant!=self.tp.electrolyte_reactions[reaction]['reaction'][0][-1]:
                                    educts+=' + '
                            for reactant in self.tp.electrolyte_reactions[reaction]['reaction'][1]:
                                products+=reactant
                                if reactant!=self.tp.electrolyte_reactions[reaction]['reaction'][1][-1]:
                                    products+=' + '
                            self.set('k'+str(i)+'f', "{} [{}]".format(self.tp.electrolyte_reactions[reaction]['rates'][0],unit_f), "rate constant: {}".format(educts+' -> '+products))
                            self.set('k'+str(i)+'r', "{} [{}]".format(self.tp.electrolyte_reactions[reaction]['rates'][1],unit_r), "rate constant: {}".format(products+' -> '+educts))

        def set(self,par_name,par_exp,par_desc):
            self.s+='    model.component("comp1").variable("var'+str(self.index)+'").set(\"{}\",\"{}\",\"{}\");\n'.format(par_name,\
                par_exp,par_desc)

        def get(self):
            return self.s

    class param():
        def __init__(self,tp=None,comsol_args={}):

            self.comsol_args=comsol_args
            self.s=''

            if tp is None:
                return
            else:
                self.tp=tp

            self.s+='/*\n'
            self.s+=' *PARAMETER\n'
            self.s+=' */\n'

            ##1) parameters from input
            for pa in self.comsol_args['parameter']:
                pa_val=self.comsol_args['parameter'][pa][0]
                pa_des=self.comsol_args['parameter'][pa][1]
                self.set(pa,str(pa_val),pa_des)
    
            ##2) default parameters for model_types
            param_set=['diffusion','charges','fluxes',\
                        'init_conc',\
                        'geo_bl']
            if self.comsol_args['model_type']=='tp_dilute_species':
                pass
            elif self.tp.model_type=='porous_electrode':
                param_set+=['porous_params']
    
            if 'diffusion' in param_set:
                for i,sp in enumerate(self.tp.species):
                    self.set('D'+str(i+1),str(self.tp.species[sp]['diffusion'])+'[m^2/s]',self.tp.species[sp]['name']+' Diffusion coefficient')
            if 'charges' in param_set:
                for i,sp in enumerate(self.tp.species):
                    self.set('Z'+str(i+1),self.tp.species[sp]['charge'],self.tp.species[sp]['name']+' ionic charge')
    
            if 'fluxes' in param_set:
                for i,sp in enumerate(self.tp.species):
                    if type(self.tp.species[sp]['flux'])!=str:
                        #define fluxes here as parameters
                        self.set("j{}".format(i+1), "RF*flux_factor*{}[mol/m^2/s]".format(self.tp.species[sp]['flux']),"{} flux".format(self.tp.species[sp]['name']))
    
            if 'init_conc' in param_set:
                for i,sp in enumerate(self.tp.species):
                    self.set("ci{}".format(i+1), "{}[mol/m^3]".format(self.tp.species[sp]['bulk_concentration']), "{} bulk_concentrations".format(self.tp.species[sp]['name']))
    
            if 'geo_bl' in param_set:
                self.set("L_cell", self.tp.xmax, "Cell length")
    
            ##standard parameter for all model_type:
            self.set("T", "{}[K]".format(self.tp.system['temperature']), "Temperature")
            self.set("RT", "R_const*T", "Molar gas constant * Temperature")
            self.set("phiM", str(self.tp.system['phiM'])+'[V]', "Metal Potential")
            if 'PZC' in self.comsol_args['solver_settings']['ramp']['names']:
                #metal potential should work just fine
                PZC_init=self.tp.system['phiM']
                if self.comsol_args['solver_settings']['solver_sequence']=='tds_elstat':
                    #self.set("phiPZC", '('+str(self.tp.system('phiPZC')-1.)+'+PZC_factor)[V]', "Metal PZC Potential")
                    self.set("phiPZC", '('+str(PZC_init)+'+PZC_factor*('+str(self.tp.system['phiPZC'])+'-'+str(PZC_init)+'))[V]', "Metal PZC Potential")
                else:
                    #self.set("phiPZC", '(-0.8+flux_factor)[V]', "Metal PZC Potential")
                    self.set("phiPZC", '('+str(PZC_init)+'+flux_factor*('+str(self.tp.system['phiPZC'])+'-'+str(PZC_init)+'))[V]', "Metal PZC Potential")
            else:
                self.set("phiPZC", str(self.tp.system['phiPZC'])+'[V]', "Metal PZC Potential")
            self.set("lambdaD", str(self.tp.debye_length)+'[m]', "Debye length")
            if 'CS' in self.comsol_args['solver_settings']['ramp']['names']:
                #2e-5 F/cm^2 (20 mikroF/cm^2) usually works fine, so ramp it up from there to the value that was requested
                CS_init=2e-5
                if self.comsol_args['solver_settings']['solver_sequence']=='tds_elstat':
                    self.set("CS", '('+str(CS_init)+'+PZC_factor*'+str(self.tp.system['Stern capacitance']/1e6-CS_init)+')[F/cm^2]', "Stern layer capacitance")
                else:
                    self.set("CS", '('+str(CS_init)+'+flux_factor*'+str(self.tp.system['Stern capacitance']/1e6-CS_init)+')[F/cm^2]', "Stern layer capacitance")
            else:
                self.set("CS", str(self.tp.system['Stern capacitance']/1e6)+'[F/cm^2]', "Stern layer capacitance")
            self.set("eps_r", self.tp.system['epsilon'], "relative permittivity")
#            self.set("epsS", 'epsilon0_const*'+str(self.tp.system['Stern epsilon']), "Stern layer effective permittivity")
#            self.set("lambdaS", "epsS/CS", "Stern layer thickness")
            self.set("conc_std", "1 [mol/m^3]", "Standard concentration (1mol/l)")
    
            self.set("flux_factor", "1", "factor scaling the flux")
            self.set("PZC_factor", "1", "factor scaling the flux")


        def set(self,par_name,par_exp,par_desc):
            self.s+='    model.param().set(\"{}\",\"{}\",\"{}\");\n'.format(par_name,\
                par_exp,par_desc)

        def get(self):
            return self.s
    
    class components():
        def __init__(self):
            #COMPONENTS: CREATE THE MODEL
            self.s=''
            self.s+='/*\n'
            self.s+=' *COMPONENTS & GEOMETRY\n'
            self.s+=' */\n'
            self.s+='    model.component().create("comp1", true);\n'
            self.s+='    model.component("comp1").geom().create("geom1", 1);\n'
            self.s+='    model.component("comp1").mesh().create("mesh1");\n'
            self.s+='    model.component("comp1").geom("geom1").create("i1", "Interval");\n'
            self.s+='    model.component("comp1").geom("geom1").feature("i1").set("p2", "L_cell");\n'
            self.s+='    model.component("comp1").geom("geom1").run();\n'
        def get(self):
            return self.s
    
