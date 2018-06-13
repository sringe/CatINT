import sys
sys.path.insert(0,'/scratch/users/sringe/transport/catint2')
sys.path.insert(0,'/scratch/users/sringe/transport/catmap')
from shutil import copyfile as copy
from catint.transport import Transport
from catint.calculator import Calculator
from catint.plot import Plot
from catint.catmap_wrapper import CatMAP
import numpy as np
import sys
from units import *
from read_data import read_data

only_catmap=True

pH_i=7.0 #12.2 #12.2 #6.8
nobuffer=True #False #True #False #True #False #True 

educt='CO2' #CO2 or CO

nx=200
nflux_comsol=100
grid_factor=200
mix_scf=0.5
nphi=15

tau_scf=0.01

RF=1

min_desc_delta=0.2
max_desc_delta=0.2

grid_factor_domain=grid_factor
grid_factor_bound=grid_factor

include_protons=False

use_elreac=True
if nobuffer:
    use_elreac=False

###########################################################################
#REACTIONS
###########################################################################
#reactants:             (first line are the educts, second the products)
#constant:              (dimensionless (mol/m^3))
#rates:                 (forward and backward rates)

if use_elreac:
    electrolyte_reactions=['phosphate-base']
    if include_protons:
        electrolyte_reactions+=['phosphate-acid']

electrode_reactions={
    #'H2':           {   'reaction':            '2 H2O + 2 e- -> H2 + 2 OH-'},
    #'H2':           {   'reaction':             '2 HCO3- + 2 e- -> H2 + 2 CO32-'},
    'H2':           {   'reaction':            '2 H2O + 2 e- -> H2 + 2 OH-'},
    'CH4':          {   'reaction':            'CO + 5 H2O + 6 e- -> CH4 + 6 OH-'},
    'CH3CH2OH':     {   'reaction':            '2 CO + 7 H2O + 8 e- -> CH3CH2OH + 8 OH-'}
#    'HCOOH':        {    'reaction':            'CO2 + 2 H2O + 2 e- ->  HCOOH + 2 OH-'}
#    'CH2O':         {   'reaction':            'CO2 + 3 H2O + 4 e- -> CH2O + 4 OH-'}}
   # 'C2H4':         {   'reaction':            '2 CO2 + 8 H2O + 12 e- -> C2H4 + 12 OH-'}}
    }

#if educt=='CO2':
#    electrode_reactions['CO']={   'reaction':            'CO2 + H2O + 2 e- ->  CO + 2 OH-'}

#reactions=\
#    {
#    'buffe':           {   'reactants':            [['CO2','H2O'],['H2CO3']],
#                            'constant':             2.63e-3},                               #KH
#    'buffer-acid':      {   'reactants':            [['CO2','H2O'],['HCO3-','H+']],
#                            'constant':             (4.44e-7)*1000.0},                      #K1a
#    'buffer-base':      {   'reactants':            [['CO2','OH-'],['HCO3-']],
#                            'constant':             (4.44e7)/1000.0,                        #K1b
#                            'rates':                [(5.93e3)/1000.0,(5.93e3)/(4.44e7)]},   #"k1f, k1r"
#    'buffer-base2':     {   'reactants':            [['HCO3-','OH-'],['CO32-','H2O']],
#                            'constant':             (4.66e3)/1000.0,
#                            'rates':                [(1.0e8)/1000.0,(1.0e8)/(4.66e3)]},     #"k2f,k2r"
#    'buffer2':          {   'reactants':            [['CO2','CO32-','H2O'],['HCO3-','HCO3-']],
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
    'Stern capacitance': 20, #std: 20
    'bulk_pH':pH_i,
    'potential drop':'Stern', #Stern', #either Stern or full
    #'ion radius': 6.0
#    'Henry constants':{
#        'CO2': 0.03419,
#        'CO':9.5e-4
#        }
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
#CO2_i = 0.03419*system['pressure']*1000. #initial CO2(aq) bulk concentrations at t=0 and Pressure P in [mol/m3] units
#                        #from Henry constant (29.41 atm/M
CO_i = 9.7e-4*system['pressure']*1000.
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
#HCO3m_i=electrolyte_reactions['buffer-base']['constant']*CO2_i*OHm_i
#CO32m_i=electrolyte_reactions['buffer-base2']['constant']*HCO3m_i*OHm_i
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
#pH_i=14+np.log10(OHm_i/1000.)
#print 'pH',pH_i
#CO2_i = 1000.
#CO_i = 0.0

Clm_i=0.0 #0.2*1000.-OHm_i

Hm_i=10**(-pH_i)*1000.0
if not include_protons:
    Hm_i=0.0

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
    K_i = Clm_i+OHm_i-Hm_i #+0.1/1000.
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
                            'bulk_concentration':   K_i},
    'OH-':              {   'symbol':               'OH^-',
                            'name':                 'hydroxyl',
                            'diffusion':            5.273e-9,
                            'bulk_concentration':   OHm_i},
#    'Cl-':              {   'symbol':               'Cl^-',
#                            'name':                 'chloride',
#                            'diffusion':            2.03e-9, #http://www.aqion.de/site/194
#                            'bulk concentration':   Clm_i},
    'H2':               {   'symbol':               'H_2',
                            'name':                 'hydrogen',
                            'diffusion':            4.50e-009},
    'CO':               {   'symbol':               'CO',
                            'name':                 'carbon monoxide',
                            'diffusion':            2.03e-9,
                            'bulk_concentration':   CO_i},
    'CH4':              {   'symbol':               'CH_4',
                            'name':                 'methane',
                            'diffusion':            1.49e-009},
#    'CH2O':             {   'symbol':               'CH_2O',
#                            'name':                 'formaldehyde',
#                            'diffusion':            1e-9},
#   # 'C2H4':             {   'symbol':               'C_2H_4',
#   #                         'name':                 'ethylene',
#   #                         'diffusion':            1.87e-009},
#    'HCOOH':            {   'symbol':               'HCOO^-',
#                            'name':                 'formate',
#                            'diffusion':            1.454e-009},
    'CH3CH2OH':             {     'name':           'ethanol', #assuming EtOH
                            'symbol':               'CH_3CH_2OH',
                            'diffusion':            0.84e-009},
    }

if include_protons:
    species.update(\
    {
    'H+':               {   'symbol':               'H^+',
                            'name':                 'hydronium',
                            'diffusion':            9.311e-9,   #CRC handbook, IONIC CONDUCTIVITY AND DIFFUSION AT INFINITE DILUTION
                            'bulk concentration':   Hm_i}
    })

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
comsol_args['parameter']['e0']=['1[C]','electronic charge']

system['active site density']=4.1612542339231805e-07

comsol_args['parameter']['RF']=[RF,'Roughness Factor']
comsol_args['parameter']['grid_factor_domain']=[str(grid_factor_domain),'Grid factor']
comsol_args['parameter']['grid_factor_bound']=[str(grid_factor_bound),'Grid factor']
comsol_args['parameter']['grid_factor']=[str(grid_factor),'Grid factor']
comsol_args['nflux']=nflux_comsol

###########################################################################
#RATE EQUATIONS/FLUXES
###########################################################################

species['H2']['flux']='catmap' #H2_rate
species['CO']['flux']='catmap' #CO_rate
species['CH3CH2OH']['flux']='catmap' #CO_rate
species['CH4']['flux']='catmap' #CO_rate

boundary_thickness=7.93E-05 #in m

if not nobuffer:
    visc=viscosity(species['HCO3-']['bulk concentration']/10**3), #Pa*s at 25C of KHCO3 solution
system['boundary thickness']=boundary_thickness
#system['electrolyte viscosity']=visc[0]

#descriptor method
#comsol_args['desc_method']='external' #internal-cont'
#comsol_args['model_type']='tp_dilute_species'
#comsol_args['solver']='parametric'

comsol_args['par_method']='internal'
#comsol_args['desc_method']='internal-cont'

###########################################################################
#BOUNDARY CONDITIONS FOR PBE
###########################################################################

catmap_args={}
catmap_args['n_inter']='automatic'
catmap_args['min_desc_delta']=min_desc_delta
catmap_args['max_desc_delta']=max_desc_delta
catmap_args['desc_method']='automatic'

potentials=[-1.0] #,-0.75,-0.5,-0.25,0.0]
results=[]

for potential in potentials:
    descriptors={'phiM':list(np.linspace(-0.4,-1.0,nphi))}
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
            catmap_args=catmap_args,
            model_name='CO2R',
            descriptors=descriptors,
            nx=nx)
    else:
        tp=Transport(
            species=species,
            electrode_reactions=electrode_reactions,
            electrolyte_reactions=electrolyte_reactions,
            system=system,
            pb_bound=pb_bound,
            catmap_args=catmap_args,
            comsol_args=comsol_args,
            model_name='CO2R',
            descriptors=descriptors,
            nx=nx)
    
    
    tp.set_calculator('comsol') #odespy') #--bdf')
    
    if only_catmap:
        cm=CatMAP(transport=tp,model_name='CO2R')
        for p in descriptors['phiM']:
            print '!!! now running p = '+str(p)
            tp.system['phiM']=p
            cm.run()
    else:
        c=Calculator(transport=tp,tau_scf=tau_scf,ntout=1,dt=1e-1,tmax=10,mix_scf=mix_scf)
        c.run()
#        tp.save() #saves all data to pickle files to enable restart or plotting later
    
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

