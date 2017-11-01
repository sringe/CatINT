#Comsol class of transport model
from units import *
import numpy as np
import os

class Comsol():
    """This class does all operations need to write input files for comsol and read output"""
    def __init__(self,path='.',transport=None,method=[\
            ['static'],['time-dependent','parametric','epsilon']]):
        if transport is None:
            self.tp.logger.error('No transport object provided for calculator. Stopping here.')
            sys.exit()
        else:
            self.tp=transport #transport object
        self.tp.path=path
        self.method=method #static or time-dependent

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

        lines=[
            'L_cell '+str(max(self.tp.xmesh))+' Cell length',
#            lambdaD/epsilon Cell length',
            'epsilon lambdaD/L_cell Dimensionless Debye length scale',
            'T '+str(self.tp.system['temperature'])+'[K] Temperature',
            'RT R_const*T Molar gas constant * Temperature',
            'delta lambdaS/lambdaD Dimensionless Stern layer thickness'
            'cM 1[mol/m^3] Metal reference concentration',
            'alphac 0.5 Cathodic charge transfer coefficient',
            'alphaa 1-alphac Anodic charge transfer coefficient',
            'V 0.5 Potential',
            #'jr 0.01 Dimensionless anodic reaction current density',
            #'kc 0.01 Dimensionless cathodic rate coefficient',
            #'J 0.9 Dimensionless cell current density',
            #'Kc kc*4*Dp/L Cathodic rate constant',
            #'Ka jr*4*Dp*cref/(L*cM) Anodic rate constant',
            #'id 4*Z*F_const*Dp*cref/L Nernst limiting current density',
            #'icell J*id Cell current density',
            'lambdaD '+str(self.tp.debye_length)+'[m] Debye length',
            'CS 18*1e-6[F/cm^2] Stern layer capacitance',
            'delta lambdaS/lambdaD Dimensionless Stern layer thickness',
            'eps_r '+str(self.tp.system['epsilon'])+' relative permittivity',
            'epsS epsilon0_const*2 Stern layer effective permittivity',
            'lambdaS epsS/CS Stern layer thickness'
            'V 0.2[V] Electrode potential']

        if os.path.exists(file_name):
            os.remove(file_name)
        with open(file_name,'a') as outfile:
            for v,u,l in zip(variable_names, units, labels):
                for i in range(len(v)):
                    outfile.write('{} {} {}\n'.format(v[i],u[i],l[i]))
            for line in lines:
                outfile.write(line+'\n')

    def write_input(self,file_name='model.java',model_name='transport',model_comments=''):
        with open(file_name,'w') as inp:
            inp.write('/*\n')
            inp.write(' *Transport Model as created by Transport Framework\n')
            inp.write(' */\n')
            inp.write('public class transport_model {\n')
            inp.write('  public static Model run()\n')
            inp.write('    Model model = ModelUtil.create("Model");\n')
            inp.write('    model.modelPath("'+self.tp.path+'");\n')
            inp.write('    model.label("'+model_name+'");\n')
            inp.write('    model.comments("'+model_comments+'")\n')

            #PARAMETER
            #diffusion coefficients
            for i,sp in enumerate(self.tp.species):
                inp.write('    model.param().set("D{}", "{}[m^2/s]", "{} Diffusion coefficient");\n'.format(\
                        i+1, self.tp.species[sp]['diffusion'],self.tp.species[sp]['name']))
            #charges
            for i,sp in enumerate(self.tp.species):
                inp.write('    model.param().set("Z{}", "{}", "{} ionic charge");\n'.format(\
                        i+1, self.tp.species[sp]['charge'],self.tp.species[sp]['name']))
            #rates
            for i,sp in enumerate(self.tp.species):
                inp.write('    model.param().set("k{}", "{}[mol/m^2/s]", "{} flux");\n'.format(\
                        i+1, self.tp.species[sp]['flux'],self.tp.species[sp]['name']))
            #initial concentrations
            for i,sp in enumerate(self.tp.species):
                inp.write('    model.param().set("ci{}", "{}[mol/m^3]", "{} bulk concentrations");\n'.format(\
                        i+1, self.tp.species[sp]['bulk concentration'],self.tp.species[sp]['name']))
            #some system parameters
            inp.write('    model.param().set("L_cell", "{}", "Cell length")\n'.format(self.tp.xmax))
            inp.write('    model.param().set("epsilon", "lambdaD/L_cell", "Dimensionless Debye length scale");\n')
            inp.write('    model.param().set("T", "''[K]", "Temperature");\n'.format(self.tp.system['temperature']))
            inp.write('    model.param().set("RT", "R_const*T", "Molar gas constant * Temperature");\n')
            inp.write('    model.param().set("delta", "lambdaS/lambdaD", "Dimensionless Stern layer thicknesscM 1[mol/m^3] Metal reference concentration");\n')
            inp.write('    model.param().set("alphac", "0.5", "Cathodic charge transfer coefficient");\n')
            inp.write('    model.param().set("alphaa", "1-alphac", "Anodic charge transfer coefficient");\n')
            inp.write('    model.param().set("V", "'+str(self.tp.system['vzeta'])+' [V]", "Potential");\n')
            inp.write('    model.param().set("lambdaD", "'+str(self.tp.debye_length)+'[m]", "Debye length");\n')
            inp.write('    model.param().set("CS", "'+str(self.tp.system['Stern capacitance']/100**2)+'[F/cm^2]", "Stern layer capacitance");\n')
            inp.write('    model.param().set("eps_r", "'+str(self.tp.system['epsilon'])+'", "relative permittivity");\n')
            inp.write('    model.param().set("epsS", "epsilon0_const*'+str(self.tp.system['Stern epsilon'])+'", "Stern layer effective permittivity");\n')
            inp.write('    model.param().set("lambdaS", "epsS/CS", "Stern layer thicknessV 0.2[V] Electrode potential");\n')

            #reaction rates
            i=0
            for reaction in self.tp.reactions:
                if 'rates' in self.tp.reactions[reaction]:
                    i+=1
                    educts=''
                    products=''
                    unit='1/s'
                    unit_f=unit+'*m^3/mol'*(len(self.tp.reactions[reaction]['reactants'][0])-1)
                    unit_r=unit+'*m^3/mol'*(len(self.tp.reactions[reaction]['reactants'][1])-1)
                    for reactant in self.tp.reactions[reaction]['reactants'][0]:
                        educts+=reactant
                        if reactant!=self.tp.reactions[reaction]['reactants'][0][-1]:
                            educts+=' + '
                    for reactant in self.tp.reactions[reaction]['reactants'][1]:
                        products+=reactant
                        if reactant!=self.tp.reactions[reaction]['reactants'][1][-1]:
                            products+=' + '
                    inp.write('    model.param().set("{}", "{} [{}]", "rate constant: {}");\n'.format(\
                            'k'+str(i)+'f',self.tp.reactions[reaction]['rates'][0],unit_f,educts+' -> '+products))
                    inp.write('    model.param().set("{}", "{} [{}]", "rate constant: {}");\n'.format(\
                            'k'+str(i)+'r',self.tp.reactions[reaction]['rates'][1],unit_r,products+' -> '+educts))
            inp.write('    model.param().set("flux_factor", "1", "factor scaling the flux");\n')
            inp.write('    model.param().set("grid_factor", "'+str(self.tp.xmax/self.tp.debye_length)+'", "minimal number of x mesh points (boundary and domain)");\n')

            #COMPONENTS: CREATE THE MODEL
            inp.write('    model.component().create("comp1", true);\n')
            inp.write('    model.component("comp1").geom().create("geom1", 1);\n')

            inp.write('    model.result().table().create("tbl1", "Table");\n')

            inp.write('    model.component("comp1").mesh().create("mesh1");\n')

            inp.write('    model.component("comp1").geom("geom1").create("i1", "Interval");\n')
            inp.write('    model.component("comp1").geom("geom1").feature("i1").set("p2", "L_cell");\n')
            inp.write('    model.component("comp1").geom("geom1").run();\n')

            #VARIABLES
            inp.write('    model.component("comp1").variable().create("var1");\n')
            
            inp.write('    model.component("comp1").variable("var1").set("deltaphi", "phiM-phi", "Metal - reaction plane potential difference");\n')
            inp.write('    model.component("comp1").variable("var1").set("rho_s", "epsS*deltaphi/lambdaS", "Surface charge density");\n')
            inp.write('    model.component("comp1").variable("var1").selection().geom("geom1", 0);\n')
            inp.write('    model.component("comp1").variable("var1").selection().set(new int[]{1});\n')
            inp.write('    model.component("comp1").variable().create("var2");\n')
            inp.write('    model.component("comp1").variable("var2").set("phiM", "V", "Metal phase potential (ground)");\n')
            inp.write('    model.component("comp1").variable("var2").selection().geom("geom1", 0);\n')
            inp.write('    model.component("comp1").variable("var2").selection().set(new int[]{1});\n')
            inp.write('    model.component("comp1").variable().create("var3");\n')
            if self.tp.pb_bound['potential']['bulk'] is not None:
                potential_r=self.tp.pb_bound['potential']['bulk']
            else:
                potential_r=0.0
            inp.write('    model.component("comp1").variable("var3").set("phiM", "'+str(potential_r)+' [V]", "Metal phase potential (cell voltage)");\n')
            inp.write('    model.component("comp1").variable("var3").selection().geom("geom1", 0);\n')
            inp.write('    model.component("comp1").variable("var3").selection().set(new int[]{2});\n')

            #INTEGRATION
            inp.write('    model.component("comp1").cpl().create("intop1", "Integration");\n')
            inp.write('    model.component("comp1").cpl().create("intop2", "Integration");\n')
            inp.write('    model.component("comp1").cpl("intop1").selection().geom("geom1", 0);\n')
            inp.write('    model.component("comp1").cpl("intop1").selection().set(new int[]{2});\n')
            inp.write('    model.component("comp1").cpl("intop2").selection().set(new int[]{1});\n')

            #ELECTROSTATICS
            inp.write('    model.component("comp1").physics().create("es", "Electrostatics", "geom1");\n')
            inp.write('    model.component("comp1").physics("es").field("electricpotential").field("phi");\n')
            inp.write('    model.component("comp1").physics("es").create("sfcd1", "SurfaceChargeDensity", 0);\n')
            inp.write('    model.component("comp1").physics("es").feature("sfcd1").selection().set(new int[]{1});\n')
            inp.write('    model.component("comp1").physics("es").create("pot1", "ElectricPotential", 0);\n')
            inp.write('    model.component("comp1").physics("es").feature("pot1").selection().set(new int[]{2});\n')
            inp.write('    model.component("comp1").physics("es").create("pot2", "ElectricPotential", 0);\n')
            inp.write('    model.component("comp1").physics("es").feature("pot2").selection().set(new int[]{1});\n')
            inp.write('    model.component("comp1").physics("es").create("df1", "DisplacementField", 0);\n')
            inp.write('    model.component("comp1").physics("es").feature("df1").selection().set(new int[]{1});\n')
            inp.write('    model.component("comp1").physics("es").create("fp1", "FloatingPotential", 0);\n')
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

            inp.write('    model.component("comp1").variable("var3").active(false);\n')

            inp.write('    model.component("comp1").view("view1").axis().set("xmin", 0);\n')

            inp.write('    model.component("comp1").physics("es").feature("ccn1")\n')
            inp.write('         .set("epsilonr", new String[][]{{"eps_r"}, {"0"}, {"0"}, {"0"}, {"eps_r"}, {"0"}, {"0"}, {"0"}, {"eps_r"}});\n')
            inp.write('    model.component("comp1").physics("es").feature("sfcd1").set("rhoqs", "rho_s");\n')
            inp.write('    model.component("comp1").physics("es").feature("pot2").set("V0", "V");\n')
            inp.write('    model.component("comp1").physics("es").feature("pot2").active(false);\n')
            inp.write('    model.component("comp1").physics("es").feature("df1").active(false);\n')
            inp.write('    model.component("comp1").physics("es").feature("fp1").active(false);\n')
            inp.write('    model.component("comp1").physics("tds").prop("ShapeProperty").set("order_concentration", 2);\n')
            inp.write('    model.component("comp1").physics("tds").prop("TransportMechanism").set("Convection", false);\n')
            inp.write('    model.component("comp1").physics("tds").prop("TransportMechanism").set("Migration", true);\n')

            #add diffusion coefficients to diffusion model
            charge_str=""
            for i in range(len(self.tp.species)):
                charge_str+="{\"Z"+str(i+1)+"\"}"
                if i != len(self.tp.species)-1:
                    charge_str+=", "
            inp.write('    model.component("comp1").physics("tds").feature("cdm1").set("z", new String[][]{'+charge_str+'});\n')
            inp.write('    model.component("comp1").physics("tds").feature("cdm1").set("minput_temperature", "T");\n')
            for i in range(len(self.tp.species)):
                inp.write('    model.component("comp1").physics("tds").feature("cdm1").set("D_cp'+str(i+1)+'", new String[][]{{"D'+str(i+1)+'"}, {"0"}, {"0"}, {"0"}, {"D'+str(i+1)+'"}, {"0"}, {"0"}, {"0"}, {"D'+str(i+1)+'"}});\n')

            #add initial concentrations to diffusion model
            ci_str=""
            for i in range(len(self.tp.species)):
                ci_str+="\"ci"+str(i+1)+"\""
                if i != len(self.tp.species)-1:
                    ci_str+=", "
            inp.write('    model.component("comp1").physics("tds").field("init1").set("initc", new String{}[]{'+ci_str+'});\n')

            #boundaries
            b_str=""
            for i in range(len(self.tp.species)):
                b_str+="{1}"
                if i != len(self.tp.species)-1:
                    b_str+=", "
            inp.write('    model.component("comp1").physics("tds").feature("fl1").set("species", new int[][]{'+b_str+'});\n')

            #fluxes
            f_str=""
            for i in range(len(self.tp.species)):
                f_str+="{\"k"+str(i)+"*flux_factor\"}"
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
                for r in self.tp.reactions:
                    reaction=self.tp.reactions[r]
                    if 'rates' not in reaction:
                        continue
                    else:
                        jj+=1
                    #rates of educts:
                    for reactant in reaction['reactants'][0]:
                        if reactant not in self.tp.species:
                            continue
                        k=species_names.index(reactant)
                    #    rates[k]=""
                        prod=""
                        for reactant2 in reaction['reactants'][0]:
                            if reactant2 not in self.tp.species:
                                continue
                            k2=species_names.index(reactant2)
                            prod+="cp"+str(k2+1)+"*"
                        rates[k]+="-"+prod+"k"+str(jj)+"f"
                        prod=""
                        for reactant2 in reaction['reactants'][1]:
                            if reactant2 not in self.tp.species:
                                continue
                            k2=species_names.index(reactant2)
                            prod+="cp"+str(k2+1)+"*"
                        rates[k]+="+"+prod+"k"+str(jj)+"r"
                    #rates of products
                    for reactant in reaction['reactants'][1]:
                        if reactant not in self.tp.species:
                            continue
                        k=species_names.index(reactant)
                     #   rates[k]=""
                        prod=""
                        for reactant2 in reaction['reactants'][0]:
                            if reactant2 not in self.tp.species:
                                continue
                            k2=species_names.index(reactant2)
                            prod+="cp"+str(k2+1)+"*"
                        rates[k]+="+"+prod+"k"+str(jj)+"f"
                        prod=""
                        for reactant2 in reaction['reactants'][1]:
                            if reactant2 not in self.tp.species:
                                continue
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

            #reactions
            rate_str=get_rates()
            for index,rate in rate_str:
                inp.write('    model.component("comp1").physics("tds").feature("reac1").set("R_cp'+str(index)+'", "'+rate+'");\n')
                        
            inp.write('    model.component("comp1").physics("tds").feature("reac1").active(false);\n')
            inp.write('    model.component("comp1").physics("ge").active(false);\n')
            inp.write('    model.component("comp1").physics("ge").feature("ge1").set("name", "icell");\n')

            #global equations
            inp.write('    model.component("comp1").physics("ge").feature("ge1").set("equation", "intop1(iloc)-icell");\n')
            inp.write('    model.component("comp1").physics("ge").feature("ge1").set("description", "Cell Voltage");\n')
            inp.write('    model.component("comp1").physics("ge").feature("ge1").set("DependentVariableQuantity", "electricpotential");\n')
            inp.write('    model.component("comp1").physics("ge").feature("ge1").set("SourceTermQuantity", "currentdensity");\n')
            inp.write('    model.component("comp1").physics("ge").feature("ge1")\n')
            inp.write('         .label("Solve for potential that results in given current density icell");\n')

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

            #probes
            inp.write('    model.component("comp1").probe("pdom1").set("bndsnap1", true);\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb1").active(false);\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb1").set("table", "tbl1");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb1").set("window", "window1");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb2").set("expr", "es.Jdx");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb2").set("unit", "A/m^2");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb2").set("descr", "Displacement current density, x component");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb2").set("table", "tbl1");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb2").set("window", "window1");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb3").active(false);\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb3").set("expr", "es.Jdx-r*F_const*charge");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb3").set("unit", "A/m^2");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb3").set("descr", "es.Jdx-r*F_const*charge");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb3").set("table", "tbl1");\n')
            inp.write('    model.component("comp1").probe("pdom1").feature("ppb3").set("window", "window1");\n')

            inp.write('    model.component("comp1").physics("es").feature("ccn1").set("epsilonr_mat", "userdef");\n')

            #studies
            if 'static' in sum(self.method,[]):
                i=self.method.index('static')+1
                inp.write('    model.study().create("std'+str(i)+'");\n')
                inp.write('    model.study("std'+str(i)+'").create("stat", "Stationary");\n')
            if 'time-dependent' in sum(self.method,[]):
                i=self.method.index('time-dependent')+1
                inp.write('    model.study().create("std'+str(i)+'");\n')
                inp.write('    model.study("std'+str(i)+'").create("param", "Parametric");\n')
                inp.write('    model.study("std'+str(i)+'").create("time", "Transient");\n')

            #solutions
            i=0
            j=0
            for i in range(len(self.method)):
                method=self.method[i][0]
                i+=1
                if method=='static':
                    j+=1
                    inp.write('    model.sol().create("sol'+str(i)+'");\n')
                    inp.write('    model.sol("sol'+str(j)+'").study("std'+str(i)+'");\n')
                    inp.write('    model.sol("sol'+str(j)+'").attach("std'+str(i)+'");\n')
                    inp.write('    model.sol("sol'+str(j)+'").create("st'+str(i)+'", "StudyStep");\n')
                    inp.write('    model.sol("sol'+str(j)+'").create("v'+str(i)+'", "Variables");\n')
                    inp.write('    model.sol("sol'+str(j)+'").create("s'+str(i)+'", "Stationary");\n')
                    inp.write('    model.sol("sol'+str(j)+'").feature("s'+str(i)+'").create("fc1", "FullyCoupled");\n')
                    inp.write('    model.sol("sol'+str(j)+'").feature("s'+str(i)+'").create("d1", "Direct");\n')
                    inp.write('    model.sol("sol'+str(j)+'").feature("s'+str(i)+'").feature().remove("fcDef");\n')
                elif method=='time-dependent':
                    inp.write('    model.sol().create("sol'+str(i)+'");\n')
                    inp.write('    model.sol("sol'+str(i)+'").study("std2");\n')
                    inp.write('    model.sol("sol'+str(i)+'").attach("std2");\n')
                    inp.write('    model.sol("sol'+str(i)+'").create("st1", "StudyStep");\n')
                    inp.write('    model.sol("sol'+str(i)+'").create("v1", "Variables");\n')
                    inp.write('    model.sol("sol'+str(i)+'").create("t1", "Time");\n')
                    inp.write('    model.sol("sol'+str(i)+'").feature("t1").create("fc1", "FullyCoupled");\n')
                    inp.write('    model.sol("sol'+str(i)+'").feature("t1").create("d1", "Direct");\n')
                    inp.write('    model.sol("sol'+str(i)+'").feature("t1").feature().remove("fcDef");\n')
                    #inp.write('    model.sol().create("sol5");\n')
                    #inp.write('    model.sol("sol5").study("std2");\n')
                    #inp.write('    model.sol("sol5").label("Parametric Solutions 3");\n')

                    #inp.write('    model.batch().create("p1", "Parametric");\n')
                    #inp.write('    model.batch("p1").create("so1", "Solutionseq");\n')
                    #inp.write('    model.batch("p1").study("std2");\n')

            
