from catint.transport import Transport
from catint.calculator import Calculator 
from catint.plot import Plot
import numpy as np
import sys
from units import *
from read_data import read_data

pH_i=13.0
nobuffer=True 
#False #True #False #True #False #True #False #True #False #True #False #True #False #True #False #True 

nx=200
nphi=20 #520 #1040 #260 #130

SA=1.
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
    'H2':   {   'reaction':            '2 H2O + 2 e- -> H2 + 2 OH-'}}
    #'H2':           {   'reaction':             '2 HCO3- + 2 e- -> H2 + 2 CO32-'},
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
    'phiPZC': -0.75, #+unit_R*298.14/unit_F*pH_i*np.log(10.), #value at SHE: https://www.sciencedirect.com/science/article/pii/S002207280300799X
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
H_i = 0.00078*system['pressure']*1000.
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
    'OH-':              {   'symbol':               'OH^-',
                            'name':                 'hydroxyl',
                            'diffusion':            5.273e-9,
                            'bulk concentration':   OHm_i},
    'H+':               {   'symbol':               'H^+',
                            'name':                 'hydronium',
                            'diffusion':            9.311e-9,   #CRC handbook, IONIC CONDUCTIVITY AND DIFFUSION AT INFINITE DILUTION
                            'bulk concentration':   Hm_i},
    'H2':               {   'symbol':               'H_2',
                            'name':                 'hydrogen',
                            'diffusion':            4.50e-009,
                            'bulk concentration':   H_i},
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
comsol_args['parameter']['Ga_H']=[str(1.13*unit_F)+'[J/mol]','H Activation Energy']
comsol_args['parameter']['Ga_H_back']=[str(0.2*unit_F)+'[J/mol]','H Activation Energy']
#comsol_args['parameter']['eVToJmol']=[str(eVTokcal*1000*calToJ)+'[J/eV/mol]','eV to J/mol Conversion factor']
comsol_args['parameter']['alpha_H']=['0.45','Butler-Volmer Parameter'] #std 0.5
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





comsol_args['parameter']['OH_min']=['1e-30 [mol/m^3]','Minimal OH- concentration allowed in the evaluation of the rates'] #of the rate coverage with which the Langmuir isotherm will be scaled'] #0.44


comsol_args['boundary_variables']={}
comsol_args['global_variables']={}

comsol_args['parameter']['phiM_ref_she']=['-7*log(10.)*RT/F_const','Reference potential vs. SHE']

comsol_args['boundary_variables']['cOH_at_0']=['cp[[OH-]]','OH- concentration at the electrode']
#comp1.at0(0,cp[[OH-]])','OH- concentration at the electrode']
comsol_args['boundary_variables']['pH_at_0']=['14+log10(max(cp[[OH-]],OH_min)*Lmol)','pH at the electrode']
#14+log10(max(cOH_at_0,OH_min)*Lmol)','pH at the electrode']
#comsol_args['boundary_variables']['pH_at_0']=['0.0','pH at the electrode']


boundary_thickness=7.93E-05 #in m
system['boundary thickness']=boundary_thickness
#'delta_phi_inf_min_iR
comsol_args['boundary_variables']['delta_phi']=['0.0','Electrolyte potential drop']
comsol_args['boundary_variables']['jH']=['rho_act*A*'+\
                         'exp('+\
                            '-(Ga_H+alpha_H*(phiM-phiM_ref_she)*F_const)/RT'+\
                         ')','rate of H']
comsol_args['boundary_variables']['jH_back']=['rho_act*A*('+\
                         'exp('+\
                            '-(Ga_H+alpha_H*(phiM-phiM_ref_she)*F_const)/RT'+\
                         ')+'+\
                         'exp('+\
                            '(Ga_H_back+(1-alpha_H)*(phiM-phiM_ref_she)*F_const)/RT+log(10)*(14-pH_at_0)'+\
                         '))','rate of H back reaction']

method=0 #method2 with stationary solver only working one, so far...

H2_rate='jH'
species['H2']['flux-equation']=H2_rate


if not nobuffer:
    visc=viscosity(species['HCO3-']['bulk concentration']/10**3), #Pa*s at 25C of KHCO3 solution
#system['electrolyte viscosity']=visc[0]

comsol_args['bin_path']='/Applications/COMSOL53a/Multiphysics/bin/comsol'

#descriptor method
comsol_args['desc_method']='internal-cont'
comsol_args['par_name']='phiM'
comsol_args['model_type']='tp_dilute_species'
comsol_args['solver']='parametric'

###########################################################################
#BOUNDARY CONDITIONS FOR PBE
###########################################################################


potentials=[-1.0] #,-0.75,-0.5,-0.25,0.0]
results=[]
for potential in potentials:
    descriptors={'phiM':list(np.linspace(0.0,-1.5,nphi))}
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
    
    c=Calculator(transport=tp,tau_jacobi=1e-5,ntout=1,dt=1e-1,tmax=10,mode='stationary') #time-dependent')

    c.run() #1.0)
    tp.save() #saves all data to pickle files to enable restart or plotting later
    
    ###########################################################################

