from catint.transport import Transport
from catint.calculator import Calculator 
from catint.plot import Plot
import numpy as np
import sys
from units import *
from read_data import read_data

pH_i=13.0
nobuffer=True #False #True #False #True #False #True #False #True #False #True #False #True #False #True 

educt='CO' #CO2 or CO

nx=200 #200
nphi=520 #260 #130

SA=1
use_elreac=True
if nobuffer:
    use_elreac=False
###########################################################################
#REACTIONS
###########################################################################
#reactants:             (first line are the educts, second the products)
#constant:              (dimensionless (mol/m^3))
#rates:                 (forward and backward rates)

#all constants & rates at room temperature
#Millero1997: http://www-naweb.iaea.org/napc/ih/documents/global_cycle/vol%20I/cht_i_09.pdf
electrolyte_reactions=\
    {
    ##############################################################################################################
    'buffer-acid':      {   'reaction':            'CO2 + H2O <-> HCO3- + H+', 
                            #PURE WATER, Emerson
                            'constant':             0.000445,                               #Gupta:  0.000444 
                             'rates':               [3.7e-2,8.3333]},                       #
                            #salinity S=35, Schulz2006
#                             'constant':             0.00138951310,
#                             'rates':               [3.71e-2,26.7]},
    ##############################################################################################################
    'buffer-acid2':     {   'reaction':             'HCO3- <-> CO32- + H+',
                            #PURE WATER, Millero1997
                            #'constant':              4.79e-8,
                            'constant':             3.5317025629468759e-07,              #https://www.iaea.org/ocean-acidification/act7/Guide%20best%20practices%20low%20res.pdf
                            'rates':                [59.44,1.68304093e8]},                  #assuming Schulz2006 for hin-reactio !!!!!!!NOT SALINITY CORRECTED!!!!!!!
##                            'rates':                [59.44,12.409e8]},                  #assuming Schulz2006 for hin-reaction !!!!!!!NOT SALINITY CORRECTED!!!!!!!
                            #salinity S=35, Schulz2006
#                            'constant':             1.1888e-06,                            #Emerson: 1.0715e-6
#                            'rates':                [59.44,5e7]},                          
    ##############################################################################################################
    'buffer-base':      {   'reaction':            'CO2 + OH- <-> HCO3-',
                            #PURE WATER, Emerson
                            'constant':             43750.0,                                #Gupta:  44400.0 
                            'rates':                [7.0,16e-5]},                           #Gupta:  [5.93,13.4e-5]
                            #salinity S=35, Schulz2006
#                             'constant':             22966.014418, #Schulz2006, m^3/mol
#                             'rates':                [2.23,9.71e-5]}, #Schulz2006
    ##############################################################################################################
    'buffer-base2':     {   'reaction':            'HCO3- + OH- <-> CO32- + H2O', 
                            #PURE WATER ???? Gupta
                            'constant':              4.66,
                            'rates':                [1.0e5,21459.2274]},
                            #salinity S=35, Schulz2006
#                             'constant':             19.60784, #Schulz2006, m^3/mol
#                             'rates':                [6e6,306000]}, #Schulz2006
    ##############################################################################################################
    'self-dissociation of water':            {   'reaction':             'H2O <-> OH- + H+',
                            'constant':             1e-8, #(mol/m^3)^2
                            'rates':                [2.4e-5*1000.,2.4e-5/1e-14/1000.]} #Singh
    ##############################################################################################################
    }

electrode_reactions={
    'C1':   {'reaction': 'CO + 5 H2O + 6 e- -> C1 + 6 OH-'}, #methane
    'C2':   {'reaction': '2 CO + 7 H2O + 8 e- -> C2 + 8 OH-'}, #ethanol
#    'H2':   {   'reaction':            '2 H2O + 2 e- -> H2 + 2 OH-'}}
    #'H2':           {   'reaction':             '2 HCO3- + 2 e- -> H2 + 2 CO32-'},
    }
#reactions=\
#    {
#    'buffe':           {   'reactants':            cp[['CO2','H2O'],['H2CO3']],
#                            'constant':             2.63e-3},                               #KH
#    'buffer-acid':      {   'reactants':            cp[['CO2','H2O'],['HCO3-','H+']],
#                            'constant':             (4.44e-7)*1000.0},                      #K1a
#    'buffer-base':      {   'reactants':            cp[['CO2','OH-'],['HCO3-']],
#                            'constant':             (4.44e7)/1000.0,                        #K1b
#                            'rates':                [(5.93e3)/1000.0,(5.93e3)/(4.44e7)]},   #"k1f, k1r"
#    'buffer-base2':     {   'reactants':            cp[['HCO3-','OH-'],['CO32-','H2O']],
#                            'constant':             (4.66e3)/1000.0,
#                            'rates':                [(1.0e8)/1000.0,(1.0e8)/(4.66e3)]},     #"k2f,k2r"
#    'buffer2':          {   'reactants':            cp[['CO2','CO32-','H2O'],['HCO3-','HCO3-']],
#                            'constant':             9.52e3}                                 #K3
#    }
###########################################################################


###########################################################################
#THERMODYNAMIC VARIABLES
###########################################################################
system=\
    {
    'temperature':  298,     #K
    'pressure':     1.,      #atm
    'water viscosity':  8.90e-004, #Pa*s at 25C
    #calculate the electrolyte viscosity. This will be used to rescale the diffusion coefficients
    #according to Einstein-Stokes relation: D_in_electrolyte = D_in_water * mu0/mu
    'epsilon': 78.36,
#    'exclude species': ['CO32-','HCO3-'], #exclude this species from PNP equations
    'migration': True,
    'electrode reactions': True,
    'electrolyte reactions': use_elreac, #False,
    'phiPZC': -0.07, #+unit_R*298.14/unit_F*pH_i*np.log(10.), #value at SHE: https://www.sciencedirect.com/science/article/pii/S002207280300799X
    'Stern capacitance': 20 #std: 20
    }
###########################################################################

###########################################################################
#READ DATA FILE
###########################################################################

data_fluxes,boundary_thickness,viscosity,bic_i=read_data()

###########################################################################
#INITIAL CONCENTRATIONS
###########################################################################
#set up the initial concentrationss from this constants:
CO2_i = 0.03419*system['pressure']*1000. #initial CO2(aq) bulk concentrations at t=0 and Pressure P in [mol/m3] units
#                        #from Henry constant (29.41 atm/M
CO_i = 9.5e-4*system['pressure']*1000.
#CO32m_i = ((2*bic_i+electrolyte_reactions['buffer2']['constant']*CO2_i)-\
#            (np.sqrt((2*bic_i+electrolyte_reactions['buffer2']['constant']*CO2_i)**2\
#            -4.0*(bic_i)**2)))/2  #initial (CO3)2- bulk concentrations at t=0 [mol/m3]

## Initial composition of the bulk electrolyte at t=0
#HCO3m_i = bic_i-CO32m_i #initial HCO3- bulk concentrations at t=0 [mol/m3]
#K_i = bic_i #initial K+ bulk concentrations at t=0 [mol/m3]
#OHm_i = HCO3m_i/electrolyte_reactions['buffer-base']['constant']/CO2_i #initial OH- bulk concentrations at t=0 [mol/m3]
#pH_i = 14+np.log10(OHm_i/1000.0) #initial pH (in log. arg must be conc in M)

##1) option: initialize with CO2_i and OHm_i
OHm_i=10**(pH_i-14.)*1000.0
HCO3m_i=electrolyte_reactions['buffer-base']['constant']*CO2_i*OHm_i
CO32m_i=electrolyte_reactions['buffer-base2']['constant']*HCO3m_i*OHm_i
#print 'HCO3m_i OHm_i CO2_i CO32m_i'
#print 'HCO3m_i',HCO3m_i, OHm_i, CO2_i, CO32m_i
##2) option: initialize with HCO3m_i and OHm_i #!currently used!!
#print 'CO2 before',CO2_i
#OHm_i=10**(pH_i-14.)*1000.0
#HCO3m_i=0.1*1000.
#CO32m_i=electrolyte_reactions['buffer-base2']['constant']*HCO3m_i*OHm_i
#CO2_i=HCO3m_i/OHm_i/electrolyte_reactions['buffer-base']['constant']
#print 'HCO3m_i',HCO3m_i, OHm_i, CO2_i, CO32m_i
#3) option: initialize with CO2_i and HCO3m_i
#HCO3m_i=0.5*1000
#OHm_i=HCO3m_i/CO2_i/electrolyte_reactions['buffer-base']['constant']
pH_i=14+np.log10(OHm_i/1000.)
print 'pH',pH_i


#sys.exit()

Hm_i=10**(-pH_i)*1000.0

#HCO3m_i=OHm_i*CO2_i*electrolyte_reactions['buffer-base']['constant']
#CO32m_i = HCO3m_i*OHm_i*electrolyte_reactions['buffer-base2']['constant'] #electrolyte_reactions['buffer-base2']['constant']**2*OHm_i**2*electrolyte_reactions['buffer2']['constant']*CO2_i
#bic_i = np.sqrt(electrolyte_reactions['buffer2']['constant']*CO2_i*CO32m_i)

#OHm_i=10**(pH_i-14.)*1000.0
#HCO3m_i=0.1*1000. #set to 0.1M
#CO2_i=HCO3m_i/OHm_i/electrolyte_reactions['buffer-base']['constant']
#CO32m_i = electrolyte_reactions['buffer-base2']['constant']*HCO3m_i*OHm_i
#CO32m_i = electrolyte_reactions['buffer-base2']['constant']**2*OHm_i**2*electrolyte_reactions['buffer2']['constant']*CO2_i
#bic_i=np.sqrt(electrolyte_reactions['buffer2']['constant']*CO2_i*CO32m_i)

#HCO3m_i=OHm_i*electrolyte_reactions['buffer-base']['constant']*CO2_i
#CO32m_i = electrolyte_reactions['buffer-base2']['constant']**2*OHm_i**2*electrolyte_reactions['buffer2']['constant']*CO2_i
#bic_i = np.sqrt(electrolyte_reactions['buffer2']['constant']*CO2_i*CO32m_i)

if nobuffer:
    K_i = OHm_i #+0.1/1000.
else:
    K_i = HCO3m_i+CO32m_i*2+OHm_i-Hm_i

#print 'new bic_i=',bic_i
#sys.exit()

#CO32m_i=
###########################################################################

###########################################################################
#SPECIES DATA
###########################################################################
#(the charges are determined automatically from the species name)
#diffusion          [m/s]       (infinite dilution in water at 25C)
#bulk concentrations [mol/m^3]
#zeff number of electrons need to derive product from CO2
#req number of educt molecules needed to derive product
#Henry              [atm/M] (from http://butane.chem.uiuc.edu/pshapley/GenChem1/L23/web-L23.pdf)
species=\
    {
    'K':                {   'symbol':               'K^+',
                            'name':                 'potassium',
                            'diffusion':            1.957e-9,
                            'bulk concentration':   K_i},
    'CO2':              {   'symbol':               'CO_2',
                            'name':                 'carbon dioxide',
                            'diffusion':            1.91e-9,
                            'bulk concentration':   CO2_i},
    'OH-':              {   'symbol':               'OH^-',
                            'name':                 'hydroxyl',
                            'diffusion':            5.273e-9,
                            'bulk concentration':   OHm_i},
    'H+':               {   'symbol':               'H^+',
                            'name':                 'hydronium',
                            'diffusion':            9.311e-9,   #CRC handbook, IONIC CONDUCTIVITY AND DIFFUSION AT INFINITE DILUTION
                            'bulk concentration':   Hm_i},
#    'H2':               {   'symbol':               'H_2',
#                            'name':                 'hydrogen',
#                            'diffusion':            4.50e-009},
    'CO':               {   'symbol':               'CO',
                            'name':                 'carbon monoxide',
                            'diffusion':            2.03e-9,
                            'bulk concentration':   CO_i},
    'C1':               {   'symbol':               'CH_4',
                            'name':                 'C1-methane',
                            'diffusion':            1.49e-9},
    'C2':               {   'symbol':               'C_2H_5OH',
                            'name':                 'C2-ethanol',
                            'diffusion':            0.84e-9},
    }

if not nobuffer:
    species['CO32-']={   'symbol':               'CO_3^{2-}',
                            'name':                 'carboxylate',
                            'diffusion':            9.23e-010,
                            'bulk concentration':   CO32m_i}
    species['HCO3-']={   'symbol':               'HCO_3^-',
                            'name':                 'bicarbonate',
                            'diffusion':            1.185e-009,
                            'bulk concentration':   HCO3m_i}

#else:
#    species['Cl-']={    'symbol':               'Cl^-',
#                        'name':                 'chlorine',
#                        'diffusion':            1.185e-009,
#                        'bulk concentration':   0.1/1000.}
    
###########################################################################

###########################################################################
#Additional COMSOL variable definitions
###########################################################################

comsol_args={}
comsol_args['parameter']={}
comsol_args['parameter']['A']=['1.e13[1/s]','Exponential prefactor']
#comsol_args['parameter']['A2']=['1.e13[1/s]','Exponential prefactor']
comsol_args['parameter']['Ga_CHO']=[str(1.11746219*unit_F)+'[J/mol]','CHO Activation Energy'] #[str(1.11746219*unit_F)+'[J/mol]','CHO Activation Energy']
comsol_args['parameter']['Ga_CHOH']=[str(2.37467774*unit_F)+'[J/mol]','CHOH Activation Energy'] #[str(2.37467774*unit_F)+'[J/mol]','CHOH Activation Energy']
comsol_args['parameter']['Ga_OCCO']=[str(0.578959276*unit_F)+'[J/mol]','OCCO Activation Energy'] #[str(0.578959276*unit_F)+'[J/mol]','OCCO Activation Energy']
comsol_args['parameter']['Ga_OCCOH']=[str(1.10495851*unit_F)+'[J/mol]','OCCOH Activation Energy'] #[str(1.10495851*unit_F)+'[J/mol]','OCCOH Activation Energy']
comsol_args['parameter']['DG_CHOH']=[str(0.0*unit_F)+'[J/mol]', 'CHOH free energy barrier']
comsol_args['parameter']['DG_OCCOH']=[str(0.0*unit_F)+'[J/mol]','OCCOH free energy barrier']
#comsol_args['parameter']['Ga_H']=[str(0.9*unit_F)+'[J/mol]','H Activation Energy']
#comsol_args['parameter']['eVToJmol']=[str(eVTokcal*1000*calToJ)+'[J/eV/mol]','eV to J/mol Conversion factor']
comsol_args['parameter']['alpha_CHO']=['0.5','Butler-Volmer Parameter'] #std 0.5
comsol_args['parameter']['alpha_CHOH']=['0.8','Butler-Volmer Parameter'] #std 2.0
comsol_args['parameter']['alpha_OCCOH']=['0.9','Butler-Volmer Parameter'] #std 0.5
comsol_args['parameter']['alpha_OCCO']=['0.5','Butler-Volmer Parameter'] #std 0.5
comsol_args['parameter']['n_CHO']=['1','Butler-Volmer Parameter'] #std 0.5
#comsol_args['parameter']['alpha_H']=['0.5','Butler-Volmer Parameter'] #std 0.5
comsol_args['parameter']['n_CHOH']=['1','Butler-Volmer Parameter'] #std 2.0
comsol_args['parameter']['n_OCCOH']=['0','Butler-Volmer Parameter'] #std 0.5
comsol_args['parameter']['n_OCCO']=['0','Butler-Volmer Parameter'] #std 0.5
#comsol_args['parameter']['n_H']=['0','Butler-Volmer Parameter'] #std 0.5
comsol_args['parameter']['e0']=['1[C]','electronic charge']
comsol_args['parameter']['erho']=['80.3e-6[C/cm^2]','surface density of active sites x elementary charge']
comsol_args['parameter']['Lmol']=['1[l/mol]','conversion factor']
#active site density. singh paper: 7.04e-6[mol/m^2]
#here: 3x3 Cu211 cell as example. area=6.363x7.794*1e-20, active sites=3 (step top sites, 1.004495558139274e-05), 9 (all top sites, 3.013486674417822e-05)
comsol_args['parameter']['rho_act']=['1.004495558139274e-05[mol/m^2]','Density of Active Sites'] #from Singh paper: 7.04e-6
comsol_args['parameter']['Ga_CO_ads']=[str(-0.3*unit_F)+'[J/mol]','Adsorption barrier for CO on Cu211']
comsol_args['parameter']['Ga_CO2_ads']=[str(-0.859729294*unit_F)+'[J/mol]','Adsorption barrier for CO2 on Cu211'] #https://smartech.gatech.edu/bitstream/handle/1853/43652/fergusson_alexander_i_201205_mast.pdf
#comsol_args['parameter']['G_H_ads']=['??','Adsorption Barrier for H on Cu211. this is an electrochemical adsorption process!!!']
comsol_args['parameter']['Kads_CO']=['exp(-Ga_CO_ads/RT)','Equilibrium constant for CO adsorption']
comsol_args['parameter']['Kads_CO2']=['exp(-Ga_CO2_ads/RT)','Equilibrium constant for CO adsorption']
comsol_args['parameter']['max_coverage']=['0.44','Maximal coverage with which the Langmuir isotherm will be scaled'] #0.44
#comsol_args['parameter']['Kads_H']=['exp(-G_H_ads/RT)','Equilibrium constant for HCO3- -> *H']
#comsol_args['parameter']['max_coverage_H']=['1.','Maximal coverage with which the Langmuir isotherm will be scaled'] #0.44
comsol_args['parameter']['SA']=[SA,'Surface Area Enhancement Factor']
###########################################################################
#RATE EQUATIONS/FLUXES
###########################################################################



#give rates as COMSOL equations
#all variables used here have to be defined as COMSOL params as seen before
#C1_rate='rho_act*coverage*A*\
#        exp(-\
#            max(\
#                (Ga_CHO+(alpha1+'+str(nCHO)+')*(phiM-phi)*F_const)/RT+alpha1*(7+log10(max(cOH_at_0,0.0)*Lmol))*'+str(np.log(10.))+', \
#                (Ga_CHOH+(alpha2+'+str(nCHOH)+')*(phiM-phi)*F_const)/RT+alpha2*(7+log10(max(cOH_at_0,0.0)*Lmol))*'+str(np.log(10.))+\
#            ')\
#        )'#/F_const/'+str(species['C1']['zeff'])
#C2_rate='rho_act*coverage^2*A*\
#        exp(-\
#            max(\
#                (Ga_OCCOH+(alpha1+'+str(nOCCOH)+')*(phiM-phi)*F_const)/RT+alpha1*(7+log10(max(cOH_at_0,0.0)*Lmol))*'+str(np.log(10.))+', \
#                (Ga_OCCO+'+str(nOCCO)+'*(phiM-phi)*F_const)/RT'\
#            ')\
#        )'sigma_max=5.0 #smoothness factor for maximum, the larger it is the smoother the function is approximated


comsol_args['parameter']['OH_min']=['1e-30 [mol/m^3]','Minimal OH- concentration allowed in the evaluation of the rates'] #of the rate coverage with which the Langmuir isotherm will be scaled'] #0.44


comsol_args['variables']={}

comsol_args['variables']['coverage_only_CO']=['Kads_CO*cp[[CO]]*Lmol/(1.+cp[[CO]]*Lmol*Kads_CO)*max_coverage','CO Coverage according to Langmuir isotherm if there is no other species']
#comsol_args['variables']['coverage_only_H']=['sqrt(Kads_H*Lmol)/(1.+sqrt(Lmol*Kads_H))*max_coverage_H','H Coverage according to Langmuir isotherm assuming H2O as proton donor']
comsol_args['parameter']['phiM_ref_she']=['-7*log(10.)*RT/F_const','Reference potential vs. SHE']

comsol_args['variables']['cOH_at_0']=['comp1.at0(0,cp[[OH-]])','OH- concentration at the electrode']

if educt=='CO2':
    comsol_args['variables']['coverage']=['Kads_CO2*cp[[CO2]]*Lmol/(1.+cp[[CO2]]*Lmol*Kads_CO2)*max_coverage','CO2 Coverage according to Langmuir isotherm = CO coverage (assuming no barrier between the states).']
elif educt=='CO':
    comsol_args['variables']['coverage']=['Kads_CO*cp[[CO]]*Lmol/(1.+cp[[CO]]*Lmol*Kads_CO)*max_coverage','CO Coverage according to Langmuir isotherm']
    #competetive adsorption with H2
    #comsol_args['variables']['coverage']=['(coverage_only_CO-coverage_only_CO*coverage_only_H)/(1-coverage_only_H*coverage_only_CO)','CO coverage according to Langmuir isotherm including competitive ads.']

#comsol_args['variables']['coverage_H']=['Kads_H*cp[[HCO3-]]*Lmol/(1.+cp[[HCO3-]]*Lmol*Kads_H)*max_coverage_H','H Coverage according to Langmuir isotherm']
#comsol_args['variables']['coverage_H']=['(1-coverage)*coverage_only_H','H coverage according to Langmuir isotherm including competitive ads.']


#comsol_args['variables']['jCHO']=['rho_act*coverage*A*(max(cOH_at_0,OH_min)*Lmol)^(alpha_CHO)*'+\
#                         'exp(-'+\
#                            '(Ga_CHO+(alpha_CHO+n_CHO)*(phiM-phi)*F_const)/RT+alpha_CHO*(7)*log(10)'+\
#                         ')','rate of CHO']
#comsol_args['variables']['jCHOH']=['rho_act*coverage*A*(max(cOH_at_0,OH_min)*Lmol)^(alpha_CHO)*'+\
#                         'exp(-'+\
#                            '(Ga_CHOH+(alpha_CHOH+n_CHOH)*(phiM-phi)*F_const)/RT+alpha_CHOH*(7)*log(10)'+\
#                         ')','rate of CHOH']
#comsol_args['variables']['jOCCOH']=['rho_act*coverage^2*A*(max(cOH_at_0,OH_min)*Lmol)^(alpha_CHO)*'+\
#                        'exp(-'+\
#                            '(Ga_OCCOH+(alpha_OCCOH+n_OCCOH)*(phiM-phi)*F_const)/RT+alpha_OCCOH*(7)*log(10)'+\
#                        ')','rate of OCCOH']
#comsol_args['variables']['jOCCO']=['rho_act*coverage^2*A*(max(cOH_at_0,OH_min)*Lmol)^(alpha_CHO)*'+\
#                        'exp(-'+\
#                            '(Ga_OCCO+(alpha_OCCO+n_OCCO)*(phiM-phi)*F_const)/RT+alpha_OCCO*(7)*log(10)'+\
#                        ')','rate of OCCO']
boundary_thickness=7.93E-05 #in m
system['boundary thickness']=boundary_thickness
#'delta_phi_inf_min_iR
comsol_args['variables']['delta_phi']=['0.0','Electrolyte potential drop']
comsol_args['variables']['jCHO']=['SA*rho_act*coverage*A*'+\
                         'exp('+\
                            '-(Ga_CHO+(alpha_CHO+n_CHO)*(phiM-phiM_ref_she+delta_phi)*F_const)/RT+alpha_CHO*(14+log10(max(cOH_at_0,OH_min)*Lmol))*log(10)'+\
                         ')','rate of CHO']
comsol_args['variables']['jCHOH']=['SA*rho_act*coverage*A*'+\
                         'exp('+\
                            '-(Ga_CHOH+DG_CHOH+(alpha_CHOH+n_CHOH)*(phiM-phiM_ref_she+delta_phi)*F_const)/RT+alpha_CHOH*(14+log10(max(cOH_at_0,OH_min)*Lmol))*log(10)'+\
                         ')','rate of CHOH']
comsol_args['variables']['jOCCOH']=['SA*rho_act*coverage^2*A*'+\
                        'exp('+\
                            '-(Ga_OCCOH+DG_OCCOH+(alpha_OCCOH+n_OCCOH)*(phiM-phiM_ref_she+delta_phi)*F_const)/RT+alpha_OCCOH*(14+log10(max(cOH_at_0,OH_min)*Lmol))*log(10)'+\
                        ')','rate of OCCOH']
comsol_args['variables']['jOCCO']=['SA*rho_act*coverage^2*A*'+\
                        'exp('+\
                            '-(Ga_OCCO+(alpha_OCCO+n_OCCO)*(phiM-phiM_ref_she+delta_phi)*F_const)/RT+alpha_OCCO*(14+log10(max(cOH_at_0,OH_min)*Lmol))*log(10)'+\
                        ')','rate of OCCO']
#comsol_args['variables']['jH']=['SA*rho_act*coverage_H^2*A*'+\
#                        'exp('+\
#                            '-(Ga_H+(alpha_H+n_H)*(phiM)*F_const)/RT+alpha_H*(7+log10(max(cOH_at_0,OH_min)*Lmol))*log(10)'+\
#                        ')','rate of H']
#additional comsol outputs
#name, equation, unit
comsol_args['outputs']=[\
        ['jCHO','jCHO','mol/m^2/s'],\
        ['jCHOH','jCHOH','mol/m^2/s'],\
        ['jOCCOH','jOCCOH','mol/m^2/s'],\
#        ['jH','jH','mol/m^2/s'],\
        ['jOCCO','jOCCO','mol/m^2/s'],\
        ['coverage','coverage','']]
        #last one in list is the name of the file, first one is variable name
#        ['coverage_H','coverage_H','']] #last one in list is the name of the file, first one is variable name

method=0 #method2 with stationary solver only working one, so far...

#H2_rate='jH'
if method==0: 
    C1_rate='jCHOH'
#    C1_rate='min('+\
#                'jCHO, '+\
#                'jCHOH'+\
#            ')'#/F_const/'+str(species['C1']['zeff'])
    C2_rate='jOCCOH'
#    C2_rate='min('+\
#                'jOCCOH, '+\
#                'jOCCO'+\
#            ')'
elif method==1:
    comsol_args['parameter']['sigma_min']=['10.0','Smooth the calculation of the maximum by exponential summation. The larger this is, the smoother the maximum function is around the maximum']
    C1_rate='-1./sigma_min*'+\
            'log('+\
                'exp(sigma_min*jCHO)+'+\
                'exp(sigma_min*jCHOH)'+\
            ')'
    C2_rate='-1./sigma_min*'+\
            'log('+\
                'exp(sigma_min*jOCCO)+'+\
                'exp(sigma_min*jOCCOH)'+\
            ')'                
elif method==2:
    C1_rate='rho_act*coverage*A*'+\
            'exp(-'+\
                'max('+\
                    '(Ga_CHO+(alpha_CHO+n_CHO)*(phiM-phi)*F_const)/RT+alpha_CHO*(7+log10(max(cOH_at_0,OH_min)*Lmol))*log(10), '+\
                    '(Ga_CHOH+(alpha_CHOH+n_CHOH)*(phiM-phi)*F_const)/RT+alpha_CHOH*(7+log10(max(cOH_at_0,OH_min)*Lmol))*log(10)'+\
                ')'+\
            ')'#/F_const/'+str(species['C1']['zeff'])
    C2_rate='rho_act*coverage^2*A*'+\
            'exp(-'+\
                'max('+\
                    '(Ga_OCCOH+(alpha_OCCOH+n_OCCOH)*(phiM-phi)*F_const)/RT+alpha_OCCOH*(7+log10(max(cOH_at_0,OH_min)*Lmol))*log(10), '+\
                    '(Ga_OCCO+(alpha_OCCO+n_OCCO)*(phiM-phi)*F_const)/RT+alpha_OCCO*(7+log10(max(cOH_at_0,OH_min)*Lmol))*log(10)'+\
                ')'+\
            ')'
elif method==3:
    comsol_args['parameter']['avoid_div_zero']=['0.0','Avoid division by zero by adding small value on top of numerator and denominator']
    comsol_args['parameter']['sigma_max']=['5.0','Smooth the calculation of the maximum by exponential summation. The larger this is, the smoother the maximum function is around the maximum']
    C1_rate='rho_act*coverage*A*'+\
            '('+\
                '(1.+avoid_div_zero)/'+\
                '('+\
                    '((1+avoid_div_zero)/(max(cOH_at_0,OH_min)*Lmol+avoid_div_zero))^(sigma_max*alpha_CHO)*'+\
                        'exp(-'+\
                            'sigma_max*('+\
                                '(Ga_CHO+(alpha_CHO+n_CHO)*(phiM-phi)*F_const)/RT+alpha_CHO*7*log(10)'+\
                            ')'+\
                        ')+'+\
                    '((1+avoid_div_zero)/(max(cOH_at_0,OH_min)*Lmol+avoid_div_zero))^(sigma_max*alpha_CHOH)*'+\
                        'exp(-'+\
                            'sigma_max*('+\
                                '(Ga_CHOH+(alpha_CHOH+n_CHOH)*(phiM-phi)*F_const)/RT+alpha_CHOH*7*log(10)'+\
                            ')'+\
                        ')+'+\
                    'avoid_div_zero'+\
                ')'+\
            ')^(1./sigma_max)'
    
    C2_rate='rho_act*coverage^2*A*'+\
            '('+\
                '(1.+avoid_div_zero)/'+\
                '('+\
                    '((1+avoid_div_zero)/(max(cOH_at_0,OH_min)*Lmol+avoid_div_zero))^(sigma_max*alpha_OCCOH)*'+\
                        'exp(-'+\
                            'sigma_max*('+\
                                '(Ga_OCCOH+(alpha_OCCOH+n_OCCOH)*(phiM-phi)*F_const)/RT+alpha_OCCOH*7*log(10)'+\
                            ')'+\
                        ')+'+\
                    ''+\
                        'exp(-'+\
                            'sigma_max*('+\
                                '(Ga_OCCO+n_OCCO)*(phiM-phi)*F_const)/RT'+\
                            ')'+\
                        ')+'+\
                    'avoid_div_zero'+\
                ')'+\
            ')^(1./sigma_max)'

#electrode_reactions['C1']['rates']=[C1_rate,'0.0']
#electrode_reactions['C2']['rates']=[C2_rate,'0.0']

species['C1']['flux-equation']=C1_rate
species['C2']['flux-equation']=C2_rate
#species['H2']['flux-equation']=H2_rate


if not nobuffer:
    visc=viscosity(species['HCO3-']['bulk concentration']/10**3), #Pa*s at 25C of KHCO3 solution
#system['electrolyte viscosity']=visc[0]

comsol_args['bin_path']='/Applications/COMSOL53a/Multiphysics/bin/comsol'

###########################################################################
#BOUNDARY CONDITIONS FOR PBE
###########################################################################

potentials=[-1.0] #,-0.75,-0.5,-0.25,0.0]
results=[]
for potential in potentials:
    descriptors={'phiM':list(np.linspace(0.0,-1.2,nphi))}
    system['phiM']=potential

    #'potential','gradient','robin'
    pb_bound={
    #        'potential': {'wall':'zeta'},
    #        'gradient': {'bulk':0.0}}
        'potential':{'bulk':0.0},
        'wall':system['phiM']}


    ###########################################################################
    #SETUP AND RUN
    ###########################################################################
    if nobuffer:
        tp=Transport(
            species=species,
            electrode_reactions=electrode_reactions,
            system=system,
            pb_bound=pb_bound,
            comsol_args=comsol_args,
            descriptors=descriptors,
            nx=nx)
    else:
        tp=Transport(
            species=species,
            electrode_reactions=electrode_reactions,
            electrolyte_reactions=electrolyte_reactions,
            system=system,
            pb_bound=pb_bound,
            comsol_args=comsol_args,
            descriptors=descriptors,
            nx=nx)
    
    
    tp.set_calculator('comsol') #odespy') #--bdf')
    #tp.set_calculator('Crank-Nicolson--LF')
    #tp.set_initial_concentrations('Gouy-Chapman')
    
    c=Calculator(transport=tp,tau_jacobi=1e-5,ntout=1,dt=1e-1,tmax=10,mode='stationary',desc_method='internal-cont') #time-dependent')
    #scale_pb_grid
    c.run() #1.0)
    tp.save() #saves all data to pickle files to enable restart or plotting later
    
#    p=Plot(transport=tp)
#    p.plot(large_plots=['concentrations_reaction','desc_current_density'],\
#            small_plots=['potential','concentrations_electrolyte','current_density','pH'])
    #p.plot(large_plots=['electrode flux'])
    #results.append(potential,tp.species['CO2']['flux'])
#plt.plot(results,'-o')
#plt.show()
    #p.add_plot('polarization')
#    p.plot(cout)
 #   p.plot_polarization()
    ###########################################################################

