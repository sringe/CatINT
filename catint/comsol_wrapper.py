#Comsol class of transport model
from units import *
import numpy as np
import os
from subprocess import call
import re 
from shutil import copyfile as copy
from copy import deepcopy

class Comsol():
    """This class does all operations need to write input files for comsol and read output"""
    def __init__(self,path=os.getcwd(),transport=None,exe_path='/Applications/COMSOL53a/Multiphysics/bin/comsol',mode='time-dependent'): 
        if transport is None:
            self.tp.logger.error('No transport object provided for calculator. Stopping here.')
            sys.exit()
        else:
            self.tp=transport #transport object
        self.tp.path=path
        self.exe=exe_path
        self.results_folder_base=self.tp.outputfoldername+'/comsol_results'
        self.results_folder=self.results_folder_base
        self.outputs=['concentrations','electrostatics','current_density','electrode_flux']
        for out in self.tp.comsol_args['outputs']:
            self.outputs.append(out)
        self.mode=mode
        #additionally specified comsol output
        if not hasattr(self.tp,'comsol_outputs_data'):
            self.tp.comsol_outputs_data={}
        

    def run(self,label='',desc_val=[],
            studies=None,only_last=False):
        #only_last: if True, update only the data in the global arrays and dictionaries
        #corresponding to the last parameter in the parameter list
        self.results_folder=self.results_folder_base+'_'+label
        desc_keys=[key for key in self.tp.descriptors]
        if studies is None:
            studies={}
#            studies['time-dependent']={'None':{'None':None}}
            if self.tp.desc_method=='external':
                studies[self.mode]={'None':{'None':None}}
            else:
                #so far we can only have a single descriptor for comsol, choose the first one here:
                #idesc=-1
                #for desc in self.tp.descriptors:
                #    idesc+=1
                #    if len(self.tp.descriptors[desc])>1:
                #        if idesc>0:
                #            self.tp.logger.error('Please chose first descriptor as parametric COMSOL sweep')
                #            sys.exit()
                #        desc_choice=desc
                #        desc_val=self.tp.descriptors[desc_choice]
                #        break
                desc_choice=desc_keys[0]
                c_desc_val=self.tp.descriptors[desc_choice] #[desc_val[0]]
                studies[self.mode]={'parametric':{desc_choice:c_desc_val}}
        for study in studies:
            if study=='time-dependent' and 'parametric' in studies[study]:
                self.tp.logger.error('Internal parametric sweep in COMSOL does currently only work with stationary solver')
                sys.exit()
#            studies['stationary']={'None':{'None':None}}
#            studies['time-dependent']['None']=None
#            studies['static']={'parametric':{'epsilon':[0.01,0.1]}}
        self.studies=studies
        for study in studies:
            for method in studies[study]:
                self.tp.logger.info('Will perform {} calculation with {} setting COMSOL'.format(study,method))
                for pars in studies[study][method]:
                    self.tp.logger.info('Parameter set passed to COMSOL is {}'.format(studies[study][method]))
        if 'cont' in self.tp.desc_method:
            self.tp.logger.info('Solutions will be initialized with solutions of previous parameter set')
        else:
            self.tp.logger.info('Solutions will be reinitialized for each parameter set')
#        self.write_parameter_file()
        self.write_input()
        self.tp.logger.info('Compiling COMSOL.')
        call(self.exe+" compile "+'/'.join([os.getcwd(),self.inp_file_name]),shell=True)
        self.tp.logger.info('Starting COMSOL.')
        call(self.exe+" batch -inputfile "+'/'.join([os.getcwd(),'.'.join(self.inp_file_name.split('.')[:-1])+".class"]),shell=True)
        #~/software/transport/examples/diffuse_double_layer_with_charge_transfer_nonstatic_2.java
        self.tp.logger.info('Reading COMSOL output.')
        self.read_output(desc_val,only_last=only_last)

    def write_parameter_file(self,file_name='comsol_parameters.txt'):
        variable_names=[['D'+str(i+1) for i in range(len(self.tp.D))],\
                        ['Z'+str(i+1) for i in range(len(self.tp.D))],\
                        ['k'+str(i+1) for i in range(len(self.tp.D))],\
                        ['ci'+str(i+1) for i in range(len(self.tp.D))]]
        units=[[str(d)+'[m^2/s]' for d in self.tp.D],\
                [str(int(c/unit_F)) for c in self.tp.charges],\
                [str(-self.tp.flux_bound[isp][0])+'[m/s]' for isp in range(self.tp.nspecies)],\
                [str(self.tp.species[sp]['bulk concentration'])+'[mol/m^3]' for sp in self.tp.species]]
        labels=[[n+' Diffusion coefficient' for n in self.tp.species],\
                [n+' charge' for n in self.tp.species],\
                [n+' flux' for n in self.tp.species],\
                [n+' bulk concentrations' for n in self.tp.species]]
        if self.tp.pb_bound['potential']['wall'] is not None:
            potential_r=self.tp.pb_bound['potential']['wall']
        else:
            potential_r=0.0
        lines=[
            'L_cell '+str(max(self.tp.xmesh))+'[m] Cell length',
#            lambdaD/epsilon Cell length',
            'epsilon lambdaD/L_cell Dimensionless Debye length scale',
            'T '+str(self.tp.system['temperature'])+'[K] Temperature',
            'RT R_const*T Molar gas constant * Temperature',
            'delta lambdaS/lambdaD Dimensionless Stern layer thickness'
            'cM 1[mol/m^3] Metal reference concentration',
            'alphac 0.5 Cathodic charge transfer coefficient',
            'alphaa 1-alphac Anodic charge transfer coefficient',
            'V 0.5 Potential',
            #'jr 0.01 Dimensionless anodic reaction current_density',
            #'kc 0.01 Dimensionless cathodic rate coefficient',
            #'J 0.9 Dimensionless cell current_density',
            #'Kc kc*4*Dp/L Cathodic rate constant',
            #'Ka jr*4*Dp*cref/(L*cM) Anodic rate constant',
            #'id 4*Z*F_const*Dp*cref/L Nernst limiting current_density',
            #'icell J*id Cell current_density',
            'lambdaD '+str(self.tp.debye_length)+'[m] Debye length',
            'CS '+str(self.tp.system['Stern capacitance'])+'*1e-6[F/cm^2] Stern layer capacitance',
            'delta lambdaS/lambdaD Dimensionless Stern layer thickness',
            'eps_r '+str(self.tp.system['epsilon'])+' relative permittivity',
            'epsS epsilon0_const*eps_r Stern layer effective permittivity',
            'lambdaS epsS/CS Stern layer thickness'
            'V '+str(potential_r)+'[V] Electrode potential']
        if os.path.exists(file_name):
            os.remove(file_name)
        with open(file_name,'a') as outfile:
            for v,u,l in zip(variable_names, units, labels):
                for i in range(len(v)):
                    outfile.write('{} {} {}\n'.format(v[i],u[i],l[i]))
            for line in lines:
                outfile.write(line+'\n')

    def write_input(self,file_name='pnp_transport.java',model_name='pnp_transport.mph',model_comments=''):
        self.inp_file_name=file_name
        with open(file_name,'w') as inp:

            inp.write('import com.comsol.model.*;\n')
            inp.write('import com.comsol.model.util.*;\n')

            inp.write('/*\n')
            inp.write(' *Transport Model as created by Transport Framework\n')
            inp.write(' */\n')
            inp.write('public class '+str('.'.join(file_name.split('.')[:-1]))+' {\n')
            inp.write('  public static Model run() {\n')
            inp.write('    Model model = ModelUtil.create("Model");\n')
            inp.write('    model.modelPath("'+self.tp.path+'");\n')
            inp.write('    model.label("'+model_name+'");\n')
            inp.write('    model.comments("'+model_comments+'");\n')

            #PARAMETER

            #additional parameters from input
            for pa in self.tp.comsol_args['parameter']:
                pa_val=self.tp.comsol_args['parameter'][pa][0]
                pa_des=self.tp.comsol_args['parameter'][pa][1]
                inp.write('    model.param().set("'+pa+'", "'+str(pa_val)+'",   "'+pa_des+'");\n')

            #diffusion coefficients
            for i,sp in enumerate(self.tp.species):
                inp.write('    model.param().set("D{}", "{}[m^2/s]", "{} Diffusion coefficient");\n'.format(\
                        i+1, self.tp.species[sp]['diffusion'],self.tp.species[sp]['name']))
            #charges
            for i,sp in enumerate(self.tp.species):
                inp.write('    model.param().set("Z{}", "{}", "{} ionic charge");\n'.format(\
                        i+1, self.tp.species[sp]['charge'],self.tp.species[sp]['name']))

            for i,sp in enumerate(self.tp.species):
                if type(self.tp.species[sp]['flux'])!=str:
                    #define fluxes here as parameters
                    inp.write('    model.param().set("j{}", "{}[mol/m^2/s]", "{} flux");\n'.format(\
                        i+1, self.tp.species[sp]['flux'],self.tp.species[sp]['name']))


            #initial concentrations
            for i,sp in enumerate(self.tp.species):
                inp.write('    model.param().set("ci{}", "{}[mol/m^3]", "{} bulk concentrations");\n'.format(\
                        i+1, self.tp.species[sp]['bulk concentration'],self.tp.species[sp]['name']))
            #some system parameters
            inp.write('    model.param().set("L_cell", "{}", "Cell length");\n'.format(self.tp.xmax))
            inp.write('    model.param().set("epsilon", "lambdaD/L_cell", "Dimensionless Debye length scale");\n')
            inp.write('    model.param().set("T", "{}[K]", "Temperature");\n'.format(self.tp.system['temperature']))
            inp.write('    model.param().set("RT", "R_const*T", "Molar gas constant * Temperature");\n')
            inp.write('    model.param().set("delta", "lambdaS/lambdaD", "Dimensionless Stern layer thickness");\n')
#            inp.write('    model.param().set("cM", "1[mol/m^3]", " Metal reference concentration");\n')
            inp.write('    model.param().set("alphac", "0.5", "Cathodic charge transfer coefficient");\n')
            inp.write('    model.param().set("alphaa", "1-alphac", "Anodic charge transfer coefficient");\n')
            if self.tp.pb_bound['potential']['wall'] is not None:
                potential_r=self.tp.pb_bound['potential']['wall']
            else:
                potential_r=0.0
            inp.write('    model.param().set("phiM", "'+str(self.tp.system['phiM'])+' [V]", "Metal Potential");\n')
            inp.write('    model.param().set("phiPZC", "'+str(self.tp.system['phiPZC'])+' [V]", "Metal PZC Potential");\n')
            inp.write('    model.param().set("lambdaD", "'+str(self.tp.debye_length)+'[m]", "Debye length");\n')
            inp.write('    model.param().set("CS", "'+str(self.tp.system['Stern capacitance']/1e6)+'[F/cm^2]", "Stern layer capacitance");\n')
            inp.write('    model.param().set("eps_r", "'+str(self.tp.system['epsilon'])+'", "relative permittivity");\n')
            inp.write('    model.param().set("epsS", "epsilon0_const*'+str(self.tp.system['Stern epsilon'])+'", "Stern layer effective permittivity");\n')
            inp.write('    model.param().set("lambdaS", "epsS/CS", "Stern layer thickness");\n')
            inp.write('    model.param().set("conc_std", "1 [mol/m^3]", "Standard concentration (1mol/l)");\n')
#            inp.write('    model.param().set("V", "0.2[V]", "Electrode potential");\n')


            #reaction rates
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
                        inp.write('    model.param().set("{}", "{} [{}]", "rate constant: {}");\n'.format(\
                                'k'+str(i)+'f',self.tp.electrolyte_reactions[reaction]['rates'][0],unit_f,educts+' -> '+products))
                        inp.write('    model.param().set("{}", "{} [{}]", "rate constant: {}");\n'.format(\
                                'k'+str(i)+'r',self.tp.electrolyte_reactions[reaction]['rates'][1],unit_r,products+' -> '+educts))
            inp.write('    model.param().set("flux_factor", "1", "factor scaling the flux");\n')
            inp.write('    model.param().set("grid_factor", "'+str(self.tp.nx_init-1)+'", "minimal number of x mesh points (boundary and domain)");\n')

            inp.write('/*\n')
            inp.write(' *MODEL DEFINITION\n')
            inp.write(' */\n')
            #COMPONENTS: CREATE THE MODEL
            inp.write('    model.component().create("comp1", true);\n')
            inp.write('    model.component("comp1").geom().create("geom1", 1);\n')

            inp.write('    model.result().table().create("tbl1", "Table");\n')

            inp.write('    model.component("comp1").mesh().create("mesh1");\n')

            inp.write('    model.component("comp1").geom("geom1").create("i1", "Interval");\n')
            inp.write('    model.component("comp1").geom("geom1").feature("i1").set("p2", "L_cell");\n')
            inp.write('    model.component("comp1").geom("geom1").run();\n')

            inp.write('/*\n')
            inp.write(' *VARIABLES\n')
            inp.write(' */\n')
            #VARIABLES
            inp.write('    model.component("comp1").variable().create("var1");\n')
            
#            inp.write('    model.component("comp1").variable("var1").set("deltaphi", "phiM-phi", "Metal - reaction plane potential difference");\n')
            #inp.write('    model.component("comp1").variable("var1").set("rho_s", "epsS*deltaphi/lambdaS", "Surface charge density");\n')
            inp.write('    model.component("comp1").variable("var1").set("rho_s", "((phiM-phiPZC)-phi)*CS", "Surface charge density");\n')
            #inp.write('    model.component("comp1").variable("var1").set("rho_s", "((phiM-phiPZC))*CS", "Surface charge density");\n')
            if self.tp.pb_bound['potential']['bulk'] is not None:
                potential_r=self.tp.pb_bound['potential']['bulk']
            else:
                potential_r=0.0
#            inp.write('    model.component("comp1").variable("var1").set("phiM", "'+self.tp.system['phiM']+'[V]", "Metal phase potential (cell voltage)");\n')
            inp.write('    model.component("comp1").variable("var1").selection().geom("geom1", 0);\n')
            inp.write('    model.component("comp1").variable("var1").selection().set(new int[]{1});\n')
            inp.write('    model.component("comp1").variable().create("var2");\n')
            #additional parameters from input
            for pa in self.tp.comsol_args['variables']:
                pa_val=self.tp.comsol_args['variables'][pa][0]
                pa_des=self.tp.comsol_args['variables'][pa][1]
                matches=re.findall('\[\[(.*?)\]\]',pa_val,re.DOTALL)
                for match in matches:
                    pa_val=pa_val.replace('[['+match+']]','cp'+str([ii+1 for ii,spp in enumerate(self.tp.species) if spp==match][0]))


                inp.write('    model.component("comp1").variable("var2").set("'+pa+'", "'+str(pa_val)+'",   "'+pa_des+'");\n')
#            inp.write('    model.component("comp1").variable("var2").set("phiM", "V", "Metal phase potential (ground)");\n')
#            inp.write('    model.component("comp1").variable("var2").selection().set(new int[]{1});\n')
#            inp.write('    model.component("comp1").variable().create("var3");\n')
            #inp.write('    model.component("comp1").variable("var3").selection().geom("geom1", 0);\n')
            #inp.write('    model.component("comp1").variable("var3").selection().set(new int[]{2});\n')


            #rates
            for i,sp in enumerate(self.tp.species):
#                if self.tp.system['kinetics']=='rate-equation':
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
                    inp.write('    model.component("comp1").variable("var2").set("j{}", "{}", "{} flux");\n'.format(\
                        i+1, self.tp.species[sp]['flux'],self.tp.species[sp]['name']))
                elif i==0:
                    inp.write('    model.component("comp1").variable("var2").selection().geom("geom1", 0);\n')
                    inp.write('    model.component("comp1").variable("var2").selection().set(new int[]{1});\n')

            inp.write('/*\n')
            inp.write(' *INTEGRATION\n')
            inp.write(' */\n')
            #INTEGRATION
            inp.write('    model.component("comp1").cpl().create("intop1", "Integration");\n')
            inp.write('    model.component("comp1").cpl().create("intop2", "Integration");\n')
            inp.write('    model.component("comp1").cpl("intop1").selection().geom("geom1", 0);\n')
            inp.write('    model.component("comp1").cpl("intop1").selection().set(new int[]{2});\n')
            inp.write('    model.component("comp1").cpl("intop2").selection().set(new int[]{1});\n')

            #ELECTROSTATICS
            inp.write('/*\n')
            inp.write(' *ELECTROSTATICS\n')
            inp.write(' */\n')
            inp.write('    model.component("comp1").physics().create("es", "Electrostatics", "geom1");\n')
            inp.write('    model.component("comp1").physics("es").field("electricpotential").field("phi");\n')
            inp.write('    model.component("comp1").physics("es").create("sfcd1", "SurfaceChargeDensity", 0);\n')
            inp.write('    model.component("comp1").physics("es").feature("sfcd1").selection().set(new int[]{1});\n')
            inp.write('    model.component("comp1").physics("es").create("pot1", "ElectricPotential", 0);\n')
            inp.write('    model.component("comp1").physics("es").feature("pot1").selection().set(new int[]{2});\n')
            #only for dirichlet:
            #inp.write('    model.component("comp1").physics("es").create("pot2", "ElectricPotential", 0);\n')
            #inp.write('    model.component("comp1").physics("es").feature("pot2").selection().set(new int[]{1});\n')
            inp.write('    model.component("comp1").physics("es").create("df1", "DisplacementField", 0);\n')
            inp.write('    model.component("comp1").physics("es").feature("df1").selection().set(new int[]{1});\n')
            #floating requires one of: AC/DC Module, MEMS Module, Plasma Module, Acoustics Module, Structural Mechanics Module, Semiconductor Module
#            inp.write('    model.component("comp1").physics("es").create("fp1", "FloatingPotential", 0);\n')
            inp.write('    model.component("comp1").physics().create("tds", "DilutedSpecies", "geom1");\n')

            inp.write('    model.component("comp1").physics("tds").field("concentration").field("cp1");\n')
            conc_str=""
            for i in range(len(self.tp.species)):
                conc_str+="\"cp"+str(i+1)+"\""
                if i != len(self.tp.species)-1:
                    conc_str+=", "
            inp.write('    model.component("comp1").physics("tds").field("concentration").component(new String[]{'+conc_str+'});\n')
            inp.write('    model.component("comp1").physics("tds").create("fl1", "Fluxes", 0);\n')
            inp.write('    model.component("comp1").physics("tds").feature("fl1").selection().set(new int[]{1});\n')
            inp.write('    model.component("comp1").physics("tds").create("gconstr1", "GlobalConstraint", -1);\n')
            inp.write('    model.component("comp1").physics("tds").create("conc1", "Concentration", 0);\n')
            inp.write('    model.component("comp1").physics("tds").feature("conc1").selection().set(new int[]{2});\n')

            if self.tp.use_electrolyte_reactions:
                inp.write('    model.component("comp1").physics("tds").create("reac1", "Reactions", 1);\n')
                inp.write('    model.component("comp1").physics("tds").feature("reac1").selection().all();\n')
            inp.write('    model.component("comp1").physics().create("ge", "GlobalEquations", "geom1");\n')

            inp.write('    model.component("comp1").multiphysics().create("pc1", "PotentialCoupling", 1);\n')
            inp.write('    model.component("comp1").multiphysics("pc1").selection().all();\n')
            inp.write('    model.component("comp1").multiphysics().create("scdc1", "SpaceChargeDensityCoupling", 1);\n')
            inp.write('    model.component("comp1").multiphysics("scdc1").selection().all();\n')

            inp.write('    model.component("comp1").mesh("mesh1").create("edg1", "Edge");\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").create("size1", "Size");\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").create("size2", "Size");\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").selection().geom("geom1", 0);\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").selection().set(new int[]{1, 2});\n')

            inp.write('    model.component("comp1").probe().create("pdom1", "DomainPoint");\n')
            inp.write('    model.component("comp1").probe("pdom1").create("ppb2", "PointExpr");\n')
            inp.write('    model.component("comp1").probe("pdom1").create("ppb3", "PointExpr");\n')

            inp.write('    model.result().table("tbl1").label("Probe Table 1");\n')

#            inp.write('    model.component("comp1").variable("var3").active(false);\n')

            inp.write('    model.component("comp1").view("view1").axis().set("xmin", 0);\n')

            inp.write('    model.component("comp1").physics("es").feature("ccn1").set("epsilonr", new String[][]{{"eps_r"}, {"0"}, {"0"}, {"0"}, {"eps_r"}, {"0"}, {"0"}, {"0"}, {"eps_r"}});\n')
            inp.write('    model.component("comp1").physics("es").feature("sfcd1").set("rhoqs", "rho_s");\n')
                #for dirichlet BC's
                #inp.write('    model.component("comp1").physics("es").feature("pot2").set("V0", "phiM");\n')
                #inp.write('    model.component("comp1").physics("es").feature("pot2").active(false);\n')
                #end dirichlet
            inp.write('    model.component("comp1").physics("es").feature("df1").active(false);\n')
#                inp.write('    model.component("comp1").physics("es").feature("fp1").active(false);\n')
            inp.write('    model.component("comp1").physics("tds").prop("ShapeProperty").set("order_concentration", 2);\n')
            if self.tp.use_convection:
                inp.write('    model.component("comp1").physics("tds").prop("TransportMechanism").set("Convection", true);\n')
            else:
                inp.write('    model.component("comp1").physics("tds").prop("TransportMechanism").set("Convection", false);\n')
            if self.tp.use_migration:
                inp.write('    model.component("comp1").physics("tds").prop("TransportMechanism").set("Migration", true);\n')
            else:
                inp.write('    model.component("comp1").physics("tds").prop("TransportMechanism").set("Migration", false);\n')

            #add diffusion coefficients to diffusion model
            inp.write('/*\n')
            inp.write(' *-> ADD diffusion coefficients\n')
            inp.write(' */\n')
            charge_str=""
            for i in range(len(self.tp.species)):
                charge_str+="{\"Z"+str(i+1)+"\"}"
                if i != len(self.tp.species)-1:
                    charge_str+=", "

            inp.write('    model.component("comp1").physics("tds").feature("cdm1").set("z", new String[][]{'+charge_str+'});\n')
            inp.write('    model.component("comp1").physics("tds").feature("cdm1").set("minput_temperature", "T");\n')

            for i in range(len(self.tp.species)):
                inp.write('    model.component("comp1").physics("tds").feature("cdm1").set("D_cp'+str(i+1)+'", new String[][]{{"D'+str(i+1)+'"}, {"0"}, {"0"}, {"0"}, {"D'+str(i+1)+'"}, {"0"}, {"0"}, {"0"}, {"D'+str(i+1)+'"}});\n')
    
            inp.write('/*\n')
            inp.write(' *-> ADD initial concentrations\n')
            inp.write(' */\n')
            #add initial concentrations to diffusion model
            ci_str=""
            for i in range(len(self.tp.species)):
                ci_str+="{\"ci"+str(i+1)+"\"}"
                if i != len(self.tp.species)-1:
                    ci_str+=", "
            inp.write('    model.component("comp1").physics("tds").feature("init1").set("initc", new String[][]{'+ci_str+'});\n')

            inp.write('/*\n')
            inp.write(' *-> ADD BOUNDARIES\n')
            inp.write(' */\n')
            #boundaries
            b_str=""
            for i in range(len(self.tp.species)):
                b_str+="{1}"
                if i != len(self.tp.species)-1:
                    b_str+=", "
            inp.write('    model.component("comp1").physics("tds").feature("fl1").set("species", new int[][]{'+b_str+'});\n')

            inp.write('/*\n')
            inp.write(' *-> ADD FLUXES\n')
            inp.write(' */\n')
            #fluxes
            f_str=""
            for i in range(len(self.tp.species)):
                f_str+="{\"j"+str(i+1)+"*flux_factor\"}"
                if i != len(self.tp.species)-1:
                    f_str+=", "

            inp.write('    model.component("comp1").physics("tds").feature("fl1").set("N0", new String[][]{'+f_str+'});\n')
            inp.write('    model.component("comp1").physics("tds").feature("gconstr1").set("constraintExpression", "intop2(cm)-(cref*L)");\n')
            inp.write('    model.component("comp1").physics("tds").feature("gconstr1").active(false);\n')
            inp.write('    model.component("comp1").physics("tds").feature("conc1").set("species", new int[][]{'+b_str+'});\n')
            inp.write('    model.component("comp1").physics("tds").feature("conc1").set("c0", new String[][]{'+ci_str+'});\n')

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
                                prod+="cp"+str(k2+1)+"*"
                        rates[k]+="-"+prod+"k"+str(jj)+"f"
                        prod=""
                        for reactant2 in reaction['reaction'][1]:
                            if reactant2 not in self.tp.species:
                                prod+="conc_std*"
                                #pass
                                #continue
                            else:
                                k2=species_names.index(reactant2)
                                prod+="cp"+str(k2+1)+"*"
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
                                prod+="cp"+str(k2+1)+"*"
                        rates[k]+="+"+prod+"k"+str(jj)+"f"
                        prod=""
                        for reactant2 in reaction['reaction'][1]:
                            if reactant2 not in self.tp.species:
                                prod+="conc_std*"
                                #continue
                            else:
                                k2=species_names.index(reactant2)
                                prod+="cp"+str(k2+1)+"*"
                        rates[k]+="-"+prod+"k"+str(jj)+"r"

                rates_new=[]
                i=0
                for rate in rates:
                    i+=1
                    if len(rate)>0:
                        rates_new.append([i,rate])
                return rates_new

            inp.write('/*\n')
            inp.write(' *REACTIONS\n')
            inp.write(' */\n')
            #reactions
            if self.tp.use_electrolyte_reactions:
                rate_str=get_rates()
                for index,rate in rate_str:
                    inp.write('    model.component("comp1").physics("tds").feature("reac1").set("R_cp'+str(index)+'", "'+rate+'");\n')
                        
                inp.write('    model.component("comp1").physics("tds").feature("reac1").active(true);\n')

            inp.write('    model.component("comp1").physics("ge").active(false);\n')
            inp.write('    model.component("comp1").physics("ge").feature("ge1").set("name", "icell");\n')

            inp.write('/*\n')
            inp.write(' *GLOBAL EQUATIONS\n')
            inp.write(' */\n')
            #global equations
            inp.write('    model.component("comp1").physics("ge").feature("ge1").set("equation", "intop1(iloc)-icell");\n')
            inp.write('    model.component("comp1").physics("ge").feature("ge1").set("description", "Cell Voltage");\n')
            inp.write('    model.component("comp1").physics("ge").feature("ge1").set("DependentVariableQuantity", "electricpotential");\n')
            inp.write('    model.component("comp1").physics("ge").feature("ge1").set("SourceTermQuantity", "currentdensity");\n')
            inp.write('    model.component("comp1").physics("ge").feature("ge1")\n')
            inp.write('         .label("Solve for potential that results in given current_density icell");\n')

            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hauto", 1);\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("custom", "on");\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hmax", "L_cell/grid_factor");\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hmaxactive", true);\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hnarrow", 5000);\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size1").set("hnarrowactive", false);\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("custom", "on");\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("hmax", "lambdaD/grid_factor");\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("hmaxactive", true);\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("hmin", 3.26E-12);\n')
            inp.write('    model.component("comp1").mesh("mesh1").feature("edg1").feature("size2").set("hminactive", false);\n')
            inp.write('    model.component("comp1").mesh("mesh1").run();\n')

            inp.write('/*\n')
            inp.write(' *PROBES\n')
            inp.write(' */\n')

            #probes
            inp.write('    model.component("comp1").probe("pdom1").set("bndsnap1", true);\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb1").active(false);\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb1").set("table", "tbl1");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb1").set("window", "window1");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb2").set("expr", "es.Jdx");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb2").set("unit", "A/m^2");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb2").set("descr", "Displacement current_density, x component");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb2").set("table", "tbl1");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb2").set("window", "window1");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb3").active(false);\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb3").set("expr", "es.Jdx-r*F_const*charge");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb3").set("unit", "A/m^2");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb3").set("descr", "es.Jdx-r*F_const*charge");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb3").set("table", "tbl1");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb3").set("window", "window1");\n')

            inp.write('    model.component("comp1").physics("es").feature("ccn1").set("epsilonr_mat", "userdef");\n')

            inp.write('/*\n')
            inp.write(' *STUDIES\n')
            inp.write(' */\n')

            #solutions
            j=0
            ip=0
            for i,study in enumerate(self.studies):
                for method in self.studies[study]:
                    j+=1
                    stype=[key for key in self.studies[study][method]][0]
                    if study=='stationary':
                        inp.write('    model.study().create("std'+str(i+1)+'");\n')
                        inp.write('    model.study("std'+str(i+1)+'").create("stat", "Stationary");\n')
                    else:
                        inp.write('    model.study().create("std'+str(i+1)+'");\n')
                        inp.write('    model.study("std'+str(i+1)+'").create("time", "Transient");\n')
                    #create a new solver
                    inp.write('    model.sol().create("sol'+str(j)+'");\n')
                    inp.write('    model.sol("sol'+str(j)+'").study("std'+str(i+1)+'");\n')
                    inp.write('    model.sol("sol'+str(j)+'").attach("std'+str(i+1)+'");\n')
                    inp.write('    model.sol("sol'+str(j)+'").create("st1", "StudyStep");\n')
                    inp.write('    model.sol("sol'+str(j)+'").create("v1", "Variables");\n')
                    if study=='stationary':
                        inp.write('    model.sol("sol'+str(j)+'").create("s1", "Stationary");\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("s1").create("fc1", "FullyCoupled");\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("s1").create("d1", "Direct");\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("s1").feature().remove("fcDef");\n')
                    elif study=='time-dependent':
                        inp.write('    model.sol("sol'+str(j)+'").create("t1", "Time");\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("t1").create("fc1", "FullyCoupled");\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("t1").create("d1", "Direct");\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("t1").feature().remove("fcDef");\n')
                    ##not sure why there are these lines, just pasting them from the comsol java file:
                    #inp.write('    model.result().dataset().create("dset2", "Solution");\n')
                    #inp.write('    model.result().dataset().create("cpt1", "CutPoint1D");\n')
                    #inp.write('    model.result().dataset("dset2").set("probetag", "pdom1");\n')
                    #inp.write('    model.result().dataset("cpt1").set("probetag", "pdom1");\n')
                    #inp.write('    model.result().dataset("cpt1").set("data", "dset2");\n')
                    #inp.write('    model.result().numerical().create("pev1", "EvalPoint")\n');
                    #inp.write('    model.result().numerical().create("pev2", "EvalPoint");\n')
                    #inp.write('    model.result().numerical().create("pev3", "EvalPoint");\n')
                    #inp.write('    model.result().numerical("pev1").set("probetag", "pdom1/ppb1");\n')
                    #inp.write('    model.result().numerical("pev2").set("probetag", "pdom1/ppb2");\n')
                    #inp.write('    model.result().numerical("pev3").set("probetag", "pdom1/ppb3");\n')
                    #inp.write('    model.component("comp1").probe("pdom1").genResult(null);\n')
                    #inp.write('    model.result().dataset("dset2").label("Probe Solution 2");\n')
                    #inp.write('    model.result().numerical("pev3").set("unit", new String[]{""});\n')
                    #inp.write('    model.result().remove("pg1");\n')
                    ##end unknown lines
                    if method=='parametric':
                        ip+=1

                        par_name=''
                        par_unit=''
                        par_val=''
#                        for istype,stype in enumerate(self.studies[study][method]):
                        par_name+='\"'+stype+'\"'
                        par_unit+='\"\"'
                        par_val+='\"'
                        for val in self.studies[study][method][stype]:
                            par_val+=' '+str(val)+' '
                        par_val+='\"'
                        #if istype!=len(self.studies[study][method])-1:
                        #    par_name+=', '
                        #    par_val+=', '
                        #    par_unit+=', '

                        #setup parametric solver
                        #add parametric solver to study 1
                        inp.write('    model.study("std'+str(i+1)+'").create("param", "Parametric");\n')
                        inp.write('    model.study("std'+str(i+1)+'").feature("param").set("pname", new String[]{'+par_name+'});\n')
                        inp.write('    model.study("std'+str(i+1)+'").feature("param").set("plistarr", new String[]{'+par_val+'});\n')
                        inp.write('    model.study("std'+str(i+1)+'").feature("param").set("punit", new String[]{'+par_unit+'});\n')

                        #variables needed for parametric solver
                        inp.write('    model.sol("sol'+str(j)+'").feature("v1").set("clistctrl", new String[]{"p1"});\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("v1").set("cname", new String[]{'+par_name+'});\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("v1").set("clist", new String[]{'+par_val+'});\n')

                        #stationary solver config. select the param feature to control the parametric solver and 
                        #if needed enable solver continuation (usage of old solution)
                        inp.write('    model.sol("sol'+str(j)+'").feature("s1").create("p1", "Parametric");\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("s1").set("probesel", "none");\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("s1").feature("p1").set("control", "param");\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("s1").feature("p1").set("pname", new String[]{'+par_name+'});\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("s1").feature("p1").set("plistarr", new String[]{'+par_val+'});\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("s1").feature("p1").set("punit", new String[]{'+par_unit+'});\n')
                        inp.write('    model.sol("sol1").feature("s1").feature("p1").set("ponerror", "empty");\n')

                        if 'internal' in self.tp.desc_method:
                            if self.tp.desc_method.split('-')[1]=='reinit':
                                inp.write('    model.sol("sol1").feature("s1").feature("p1").set("pcontinuationmode", "no");\n')
                            else:
                                inp.write('    model.sol("sol1").feature("s1").feature("p1").set("preusesol", "yes");\n')
                        else:
                            inp.write('    model.sol("sol1").feature("s1").feature("p1").set("pcontinuationmode", "no");\n')

                    if study=='time-dependent':
                        inp.write('    model.study("std'+str(i+1)+'").feature("time").set("tlist", "range('+str(min(self.tp.tmesh_init))+','+str(self.tp.dt_init)+','+str(self.tp.tmax_init)+')");\n')

            #now set numeric parameters and run!
            j=0
            ip=0
            for i,study in enumerate(self.studies):
                for method in self.studies[study]:
                    j+=1
                    if study=='stationary':
                        inp.write('    model.sol("sol'+str(j)+'").attach("std'+str(i+1)+'");\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("s1").feature("fc1").set("initstep", 0.01);\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("s1").feature("fc1").set("minstep", 1.0E-6);\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("s1").feature("fc1").set("maxiter", 50);\n')
                        inp.write('    model.sol("sol'+str(j)+'").runAll();\n')
                    elif study=='time-dependent':
                        inp.write('    model.sol("sol'+str(j)+'").attach("std'+str(i+1)+'");\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("v1").set("clist", new String[]{"'+str(min(self.tp.tmesh_init))+','+str(self.tp.dt_init)+', '+str(self.tp.tmax_init)+'"});\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("t1").set("tlist", "range('+str(min(self.tp.tmesh_init))+','+str(self.tp.dt_init)+','+str(self.tp.tmax_init)+')");\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("t1").set("maxorder", 2);\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("t1").feature("fc1").set("maxiter", 8);\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("t1").feature("fc1").set("damp", 0.9);\n')
                        inp.write('    model.sol("sol'+str(j)+'").feature("t1").feature("fc1").set("jtech", "once");\n')
                        inp.write('    model.sol("sol'+str(j)+'").runAll();\n')

            inp.write('/*\n')
            inp.write(' *OUTPUTS\n')
            inp.write(' */\n')

            if not os.path.isdir(self.results_folder):
                os.makedirs(self.results_folder)


            root=os.getcwd()


            #finally output results
            def export_data(var_name='cp',unit_name='mol/m^3',label='Concentrations',file_name='results/concentrations.txt',export_count=1):
                """exports a quantity of interest. var_name/unit_name/label can be either a single string, in which case
                it is assumed that var_name should be exported for each individual species. or it can be a list, over which
                to iterate the output"""

                inp.write('    model.result().export().create("data'+str(export_count)+'", "Data");\n')
                c_str=""
                unit_str=""
                name_str=""
                if type(var_name)==str:
                    for sp,i in zip(self.tp.species,range(len(self.tp.species))):
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
                inp.write('    model.result().export("data'+str(export_count)+'").set("expr", new String[]{'+c_str+'});\n')
                inp.write('    model.result().export("data'+str(export_count)+'").set("unit", new String[]{'+unit_str+'});\n')
                inp.write('    model.result().export("data'+str(export_count)+'").set("descr", new String[]{'+name_str+'});\n')
                inp.write('    model.result().export("data'+str(export_count)+'").set("filename", "'+root+'/'+file_name+'");\n')
                inp.write('    model.result().export("data'+str(export_count)+'").run();\n')

            label_append=''
            if 'concentrations' in self.outputs:
                #concentrations
                export_data('cp','mol/m^3', 'Concentrations'+label_append, self.results_folder+'/concentrations.txt',1)
            if 'electrostatics' in self.outputs:
                #electrostatics
                export_data(['phi','es.Ex'],['V','V/m'],'Potential, Field'+label_append, self.results_folder+'/electrostatics.txt',2)
            if 'current_density' in self.outputs:
                export_data(['es.Jx'],['A/m^2'],'Current density'+label_append, self.results_folder+'/current_density.txt',3)
            if 'electrode_flux' in self.outputs:
                export_data('j','mol/m^2/s','Electrode flux'+label_append, self.results_folder+'/electrode_flux.txt',4)
            i=0
            for out in self.outputs:
                if out in self.tp.comsol_args['outputs']:
                    i+=1
                    export_data([out[0]],[out[1]],out[0]+label_append, self.results_folder+'/'+out[1]+'.txt',i+4)

            inp.write('    return model;\n')
            inp.write('  }\n')
            inp.write('\n')

            inp.write('  public static Model run2(Model model) {\n')
            inp.write('  return model;\n')
            inp.write('  }\n')

            inp.write('\n')
            inp.write('  public static void main(String[] args) {\n')
            inp.write('    Model model = run();\n')
            inp.write('    run2(model);\n')
            inp.write('  }\n')
            inp.write('\n')
            inp.write('}\n')

        #make a copy of the input file into the results folder
        copy(root+'/'+file_name,self.results_folder+'/'+file_name)

    def read_output(self,desc_val=[],only_last=False):
        """reads the output files of COMSOL written to the results folder"""
        #modify the Description line of the output with proper concentrations:
        if self.tp.descriptors is not None:
            label_append=' || '
            keys=[key for key in self.tp.descriptors]
            vals=[]
            print self.tp.descriptors
            for key in keys:
                vals.append(self.tp.descriptors[key][0])
            label_append+=keys[0]+'='+str(vals[0])+', '
            label_append+=keys[1]+'='+str(vals[1])
        else:
            label_append=''
#            label_append=''
        #if internal descriptor iteration create descriptor value list:
        if 'internal' in self.tp.desc_method:
            int_desc_list=[self.studies['stationary']['parametric'][d1] for d1 in self.studies['stationary']['parametric']][0]
            int_desc=[d1 for d1 in self.studies['stationary']['parametric']][0]
            #get the first value of the descriptor which is not used as parametric sweep
            int_desc_non=[self.tp.descriptors[desc][0] for desc in self.tp.descriptors if desc!=int_desc][0]
        else:
            int_desc_list=[]
            int_desc=None
        self.tp.potential=[]
        self.tp.efield=[]
        self.tp.current_density=[]
        if 'internal' in self.tp.desc_method:
            for id1,d1 in enumerate(self.tp.all_data):
                for id2,d2 in enumerate(self.tp.all_data[d1]):
                    self.tp.all_data[str(d1)][str(d2)]['system']['potential']=[]
                    self.tp.all_data[str(d1)][str(d2)]['system']['efield']=[]
                    self.tp.all_data[str(d1)][str(d2)]['system']['current_density']=[]
                    for iout,oout in enumerate(self.tp.comsol_args['outputs']):
                        out=oout[1]
                        if out not in self.tp.comsol_outputs_data:
                            self.tp.comsol_outputs_data[out]={(str(d1),str(d2)):[]}
                        else:
                            self.tp.comsol_outputs_data[out][(str(d1),str(d2))]=[]
        comsol_outputs_data=[]
        for i in range(len(self.tp.comsol_args['outputs'])):
            comsol_outputs_data.append([])
        for ioutput,output in enumerate(self.outputs):
            #add some more info to the description line
            if output in self.tp.comsol_args['outputs']:
                toutput=output[1]
            else:
                toutput=output
            file_lines=''
            with open(self.results_folder+'/'+toutput+'.txt', 'r') as f:
                for x in f.readlines():
                    if 'Description' in x:
                        file_lines+=''.join([x.strip(), label_append, '\n'])
                    else:
                        file_lines+=x
            with open(self.results_folder+'/'+toutput+'.txt', 'w') as f:
                f.writelines(file_lines)
            # end add info
            if ioutput==0:
                #first get the new xmesh:
                self.tp.xmesh=[]
                i=0

                for line in open(self.results_folder+'/'+toutput+'.txt','r'):
                    i+=1
                    if i>8:
                        ls=line.split()
                        for j,lss in enumerate(ls):
                            if j==0:
                                self.tp.xmesh.append(float(lss))
                self.tp.xmax=max(self.tp.xmesh)
                self.tp.nx=len(self.tp.xmesh)
                self.tp.dx=self.tp.xmax/self.tp.nx
                cout=np.zeros([self.tp.ntout,self.tp.nspecies*self.tp.nx])
                cout_tmp=np.zeros([self.tp.nt+1,self.tp.nspecies*self.tp.nx]) ##check here, why do we put nt+1?? if we put nt, we need to correct the stationary solver below
                electrode_flux=np.zeros_like(cout)
                electrode_flux_tmp=np.zeros_like(cout_tmp)
            if 'internal' in self.tp.desc_method and output=='concentrations':
                for desc in int_desc_list:
                    self.tp.all_data[str(desc)][str(int_desc_non)]['system']['cout']=np.zeros_like(cout_tmp) #([self.tp.nt,self.tp.nspecies*self.tp.nx])
            if 'internal' in self.tp.desc_method and output=='electrode_flux':
                for desc in int_desc_list:
                    self.tp.all_data[str(desc)][str(int_desc_non)]['system']['electrode_flux']=np.zeros_like(cout_tmp)
            #now read in all results
            i=0
            self.tp.logger.info('Reading COMSOL output from'+self.results_folder+'/'+toutput+'.txt')

            for line in open(self.results_folder+'/'+toutput+'.txt','r'):
                i+=1
                if i>8 and [key for key in self.studies][0]=='stationary':
                    if output in ['concentrations','electrode_flux']:
                        ls=line.split()
                        for j,lss in enumerate(ls):
                            if j>0:
                                i_sp=(j-1)%(self.tp.nspecies)
                                i_de=(j-1-i_sp)/(self.tp.nspecies)
                                if 'internal' in self.tp.desc_method:
                                    if output=='concentrations':
                                        self.tp.all_data[str(int_desc_list[i_de])][str(int_desc_non)]['system']['cout'][-2,i_sp*self.tp.nx+i-9]=float(lss)
                                    elif output=='electrode_flux':
                                        self.tp.all_data[str(int_desc_list[i_de])][str(int_desc_non)]['system']['electrode_flux'][-2,i_sp*self.tp.nx+i-9]=float(lss)
                                if output=='concentrations':
                                    cout_tmp[-2,i_sp*self.tp.nx+i-9]=float(lss)
                                elif output=='electrode_flux':
                                    electrode_flux_tmp[-2,i_sp*self.tp.nx+i-9]=float(lss)
                    elif output in ['electrostatics','current_density']:
                        ls=line.split()
                        for j,lss in enumerate(ls):
                            if j>0:
                                if output=='electrostatics':
                                    i_sp=(j-1)%2
                                    i_de=(j-1-i_sp)/2
                                    if 'internal' in self.tp.desc_method:
                                        if i_sp==0:
                                            self.tp.all_data[str(int_desc_list[i_de])][str(int_desc_non)]['system']['potential'].append(float(lss))
                                        elif i_sp==1:
                                            self.tp.all_data[str(int_desc_list[i_de])][str(int_desc_non)]['system']['efield'].append(float(lss))
                                    if i_sp==0:
                                        self.tp.potential.append(float(lss))
                                    elif i_sp==1:
                                        self.tp.efield.append(float(lss))
                                elif output=='current_density':
                                    i_sp=(j-1)%1
                                    i_de=(j-1-i_sp)
                                    if 'internal' in self.tp.desc_method:
                                        if i_sp==0:
                                            self.tp.all_data[str(int_desc_list[i_de])][str(int_desc_non)]['system']['current_density'].append(float(lss))
                                    if i_sp==0:
                                        self.tp.current_density.append(float(lss))
                    elif output in self.tp.comsol_args['outputs']:
                        ls=line.split()
                        for j,lss in enumerate(ls):
                            if j>0:
                                i_de=j-1
                                if 'internal' in self.tp.desc_method:
                                    self.tp.comsol_outputs_data[output[1]][(str(int_desc_list[i_de]),str(int_desc_non))].append(float(lss))
                                comsol_outputs_data[self.tp.comsol_args['outputs'].index(output)].append(float(lss))
                elif i>8 and [key for key in self.studies][0]=='time-dependent':
                    if output in ['concentrations','electrode_flux']:
                        #read a selection of time steps here
                        ls=line.split()
                        i_t=np.zeros([self.tp.nspecies],dtype=int)-1
                        for j,lss in enumerate(ls):
                            if j>0:
                                i_sp=(j-1)%self.tp.nspecies #index of species
                                i_t[i_sp]+=1 #index of current time step of species
                                n=i_t[i_sp]
#                                if n % int(self.tp.nt/float(self.tp.ntout)) == 0: # or n==self.tp.nt-1:
                                if n in self.tp.itout:
                                    if output=='concentrations':
                                        cout_tmp[n,i_sp*self.tp.nx+i-9]=float(lss)
                                    elif output=='electrode_flux':
                                        electrode_flux_tmp[n,i_sp*self.tp.nx+i-9]=float(lss)

                    elif output in ['electrostatics','current_density']:
                        #read only the last time step here
                        ls=line.split()
                        i_t=np.zeros([2],dtype=int)-1
                        for j,lss in enumerate(ls):
                            if j>0:
                                if output=='electrostatics':
                                    i_sp=(j-1)%2
                                else:
                                    i_sp=(j-1)%1
                                i_t[i_sp]+=1
                                n=i_t[i_sp]
                                if n==self.tp.nt-1: #n % int(self.tp.nt/float(self.tp.ntout)) == 0 or n==self.tp.nt-1:
                                    if output=='electrostatics':
                                        if i_sp==0:
                                            self.tp.potential.append(float(lss))
                                        elif i_sp==1:
                                            self.tp.efield.append(float(lss))
                                    elif output=='current_density': 
                                        if i_sp==0:
                                            self.tp.current_density.append(float(lss))
                    elif output in self.tp.comsol_args['outputs']:
                        ls=line.split()
                        i_t=np.zeros([2],dtype=int)-1
                        for j,lss in enumerate(ls):
                            if j>0:
                                i_sp=(j-1)%1
                                i_t[i_sp]+=1
                                n=i_t[i_sp]
                                if n==self.tp.nt-1:
                                    if i_sp==0:
                                        comsol_outputs_data[self.tp.comsol_args['outputs'].index(output)].append(float(lss))
        #compress cout_tmp to cout (save only relevant time steps)
        if 'internal' in self.tp.desc_method:
            for desc in int_desc_list:
                if only_last and desc!=int_desc_list[-1]:
                    continue
                tmp_ar=deepcopy(self.tp.all_data[str(desc)][str(int_desc_non)]['system']['cout'])
                self.tp.all_data[str(desc)][str(int_desc_non)]['system']['cout']=np.zeros([self.tp.ntout,self.tp.nspecies*self.tp.nx])
                nn=-1
                for n,cc in enumerate(tmp_ar):
                    if n in self.tp.itout:
                        nn+=1
                        self.tp.all_data[str(desc)][str(int_desc_non)]['system']['cout'][nn,:]=cc
                
                self.tp.all_data[str(desc)][str(int_desc_non)]['system']['cout']=np.array(self.tp.all_data[str(desc)][str(int_desc_non)]['system']['cout'])

                tmp_ar=deepcopy(self.tp.all_data[str(desc)][str(int_desc_non)]['system']['electrode_flux'])
                self.tp.all_data[str(desc)][str(int_desc_non)]['system']['electrode_flux']=np.zeros([self.tp.ntout,self.tp.nspecies*self.tp.nx])
                nn=-1
                for n,cc in enumerate(tmp_ar):
                    if n in self.tp.itout:
                        nn+=1
                        self.tp.all_data[str(desc)][str(int_desc_non)]['system']['electrode_flux'][nn,:]=cc

                self.tp.all_data[str(desc)][str(int_desc_non)]['system']['electrode_flux']=np.array(self.tp.all_data[str(desc)][str(int_desc_non)]['system']['electrode_flux'])
                
                
        nn=-1
        for n,cc in enumerate(cout_tmp):
#            if n % int(self.tp.nt/float(self.tp.ntout))==0:
            if n in self.tp.itout:
                nn+=1
                cout[nn,:]=cc
        nn=-1
        for n,cc in enumerate(electrode_flux_tmp):
            if n in self.tp.itout:
                nn+=1
                electrode_flux[nn,:]=cc


        self.tp.electrode_flux=electrode_flux

        self.tp.logger.info('Wrote out concentrations at '+str(nn+1)+' steps.')

        if 'internal' in self.tp.desc_method:
            for desc in int_desc_list:
                if only_last and desc!=int_desc_list[-1]:
                    continue
                for par in ['potential','efield','current_density']:
                    self.tp.all_data[str(desc)][str(int_desc_non)]['system'][par]=np.array(self.tp.all_data[str(desc)][str(int_desc_non)]['system'][par])
                for oout in self.tp.comsol_args['outputs']:
                    out=oout[1]
                    self.tp.comsol_outputs_data[out][(str(desc),str(int_desc_non))]=np.array(self.tp.comsol_outputs_data[out][(str(desc),str(int_desc_non))])

        self.tp.potential=np.array(self.tp.potential)
        self.tp.efield=np.array(self.tp.efield)
        self.tp.current_density=np.array(self.tp.current_density)
        comsol_outputs_data=np.array(comsol_outputs_data)
        if self.tp.descriptors is not None and not 'internal' in self.tp.desc_method:
            #save all of them in descriptor based dictionary:
            self.tp.all_data[str(desc_val[0])][str(desc_val[1])]['system']['potential']=self.tp.potential
            self.tp.all_data[str(desc_val[0])][str(desc_val[1])]['system']['efield']=self.tp.efield
            self.tp.all_data[str(desc_val[0])][str(desc_val[1])]['system']['current_density']=self.tp.current_density

        cout=np.array(cout)
        if not 'internal' in self.tp.desc_method:
            self.tp.all_data[str(desc_val[0])][str(desc_val[1])]['system']['cout']=cout
            self.tp.all_data[str(desc_val[0])][str(desc_val[1])]['system']['electrode_flux']=electrode_flux
            for i_sp,sp in enumerate(self.tp.species):
                self.tp.species[sp]['concentration']=cout[-1,i_sp*self.tp.nx:(i_sp+1)*self.tp.nx]
                if self.tp.descriptors is not None:
                    self.tp.all_data[str(desc_val[0])][str(desc_val[1])]['species'][sp]['concentration']=self.tp.species[sp]['concentration']
            self.tp.cout=cout
        else:
            for desc in int_desc_list:
                if only_last and desc!=int_desc_list[-1]:
                    continue
                for i_sp,sp in enumerate(self.tp.species):
                    self.tp.species[sp]['concentration']=self.tp.all_data[str(desc)][str(int_desc_non)]['system']['cout'][-1,i_sp*self.tp.nx:(i_sp+1)*self.tp.nx]
                    self.tp.all_data[str(desc)][str(int_desc_non)]['species'][sp]['concentration']=self.tp.species[sp]['concentration']

        
        #update total charge density
        self.tp.total_charge=np.zeros([self.tp.nx])
        for i_sp,sp in enumerate(self.tp.species):
            for ix,xx in enumerate(self.tp.xmesh):
                self.tp.total_charge[ix]+=self.tp.charges[i_sp]*self.tp.species[sp]['concentration'][ix]

        #update surface concentrations
        for i_sp,sp in enumerate(self.tp.species):
            self.tp.species[sp]['surface concentration']=self.tp.species[sp]['concentration'][0]

        #update surface pH
        if 'H+' in self.tp.species:
            self.tp.system['surface pH']=-np.log10(self.tp.species['H+']['surface concentration']/1000.)
        elif 'OH-' in self.tp.species:
            self.tp.system['surface pH']=14+np.log10(self.tp.species['OH-']['surface concentration']/1000.)

        #evaluate and save properties for each descriptor (so far the electrode current density and the comsol outputs - self.tp.comsol_args['outputs']):
        data=[]
        i1=0
        i2=0

        for reac in self.tp.electrode_reactions:
            if 'electrode_current_density' not in self.tp.electrode_reactions[reac]:
                self.tp.electrode_reactions[reac]['electrode_current_density']={}
        if not 'internal' in self.tp.desc_method:
            i2+=1
            for sp in self.tp.electrode_reactions:
                nprod=len([a for a in self.tp.electrode_reactions[sp]['reaction'][1] if a==sp])
                isp=[i for i,sp2 in enumerate(self.tp.species) if sp2==sp][0] #index in species array
                
                c_current_density=self.tp.electrode_flux[-1][isp*self.tp.nx]*self.tp.electrode_reactions[sp]['nel']*unit_F/nprod/10.
                self.tp.electrode_reactions[sp]['electrode_current_density'][(str(desc_val[0]),str(desc_val[1]))]=c_current_density
            
            for iout,oout in enumerate(self.tp.comsol_args['outputs']):
                out=oout[1]
                if out not in self.tp.comsol_outputs_data:
                    self.tp.comsol_outputs_data[out]={(str(desc_val[0]),str(desc_val[1])):comsol_outputs_data[iout]}
                else:
                    self.tp.comsol_outputs_data[out][(str(desc_val[0]),str(desc_val[1]))]=comsol_outputs_data[iout]
        else:
            i2+=1
            for desc in int_desc_list:
                if only_last and desc!=int_desc_list[-1]:
                    continue
                for sp in self.tp.electrode_reactions:
                    nprod=len([a for a in self.tp.electrode_reactions[sp]['reaction'][1] if a==sp])
                    isp=[i for i,sp2 in enumerate(self.tp.species) if sp2==sp][0] #index in species array
                    c_current_density=self.tp.all_data[str(desc)][str(int_desc_non)]['system']['electrode_flux'][-1][isp*self.tp.nx]*self.tp.electrode_reactions[sp]['nel']*unit_F/nprod/10.
                    self.tp.electrode_reactions[sp]['electrode_current_density'][(str(desc),str(int_desc_non))]=c_current_density

