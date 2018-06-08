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

pH_i=7.0 #6.8
nobuffer=False #True #False #True #False #True 

educt='CO2' #CO2 or CO

nx=400 #200
nflux_comsol=10
grid_factor=200
nphi=1

RF=1

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
    #'H2':           {   'reaction':            '2 H2O + 2 e- -> H2 + 2 OH-'},
    #'H2':           {   'reaction':             '2 HCO3- + 2 e- -> H2 + 2 CO32-'},
    'H2':           {   'reaction':            '2 H2O + 2 e- -> H2 + 2 OH-'},
    'CO':           {   'reaction':             'CO2 + H2O + 2 e- -> CO + 2 OH-'},
    'CH4':          {   'reaction':            'CO2 + 6 H2O + 8 e- -> CH4 + 8 OH-'},
    'CH3CH2OH':     {   'reaction':            '2 CO2 + 9 H2O + 12 e- -> CH3CH2OH + 12 OH-'}}
#    'CH2O':         {   'reaction':            'CO2 + 3 H2O + 4 e- -> CH2O + 4 OH-'}}
   # 'C2H4':         {   'reaction':            '2 CO2 + 8 H2O + 12 e- -> C2H4 + 12 OH-'}}


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
    'bulk_pH':pH_i
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
CO2_i = 0.03419*system['pressure']*1000. #initial CO2(aq) bulk_concentrations at t=0 and Pressure P in [mol/m3] units
#                        #from Henry constant (29.41 atm/M
#CO_i = 9.5e-4*system['pressure']*1000.
#CO32m_i = ((2*bic_i+electrolyte_reactions['buffer2']['constant']*CO2_i)-\
#            (np.sqrt((2*bic_i+electrolyte_reactions['buffer2']['constant']*CO2_i)**2\
#            -4.0*(bic_i)**2)))/2  #initial (CO3)2- bulk_concentrations at t=0 [mol/m3]

## Initial composition of the bulk electrolyte at t=0
#HCO3m_i = bic_i-CO32m_i #initial HCO3- bulk_concentrations at t=0 [mol/m3]
#K_i = bic_i #initial K+ bulk_concentrations at t=0 [mol/m3]
#OHm_i = HCO3m_i/electrolyte_reactions['buffer-base']['constant']/CO2_i #initial OH- bulk_concentrations at t=0 [mol/m3]
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
#pH_i=14+np.log10(OHm_i/1000.)
#print 'pH',pH_i
#CO2_i = 1000.
CO_i = 0.0

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
    K_i = OHm_i-Hm_i #+0.1/1000.
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
#bulk_concentrations [mol/m^3]
#zeff number of electrons need to derive product from CO2
#req number of educt molecules needed to derive product
#Henry              [atm/M] (from http://butane.chem.uiuc.edu/pshapley/GenChem1/L23/web-L23.pdf)
species=\
    {
    'K':                {   'symbol':               'K^+',
                            'name':                 'potassium',
                            'diffusion':            1.957e-9,
                            'bulk_concentration':   K_i},
    'CO2':              {   'symbol':               'CO_2',
                            'name':                 'carbon dioxide',
                            'diffusion':            1.91e-9,
                            'bulk_concentration':   CO2_i},
    'OH-':              {   'symbol':               'OH^-',
                            'name':                 'hydroxyl',
                            'diffusion':            5.273e-9,
                            'bulk_concentration':   OHm_i},
    'H+':               {   'symbol':               'H^+',
                            'name':                 'hydronium',
                            'diffusion':            9.311e-9,   #CRC handbook, IONIC CONDUCTIVITY AND DIFFUSION AT INFINITE DILUTION
                            'bulk_concentration':   Hm_i},
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
#    'HCOO-':            {   'symbol':               'HCOO^-',
#                            'name':                 'formate',
#                            'diffusion':            1.454e-009},
    'CH3CH2OH':             {     'name':           'ethanol', #assuming EtOH
                            'symbol':               'CH_3CH_2OH',
                            'diffusion':            0.84e-009},
    }

if not nobuffer:
    species['CO32-']={   'symbol':               'CO_3^{2-}',
                            'name':                 'carboxylate',
                            'diffusion':            9.23e-010,
                            'bulk_concentration':   CO32m_i}
    species['HCO3-']={   'symbol':               'HCO_3^-',
                            'name':                 'bicarbonate',
                            'diffusion':            1.185e-009,
                            'bulk_concentration':   HCO3m_i}

#else:
#    species['Cl-']={    'symbol':               'Cl^-',
#                        'name':                 'chlorine',
#                        'diffusion':            1.185e-009,
#                        'bulk_concentration':   0.1/1000.}
    
###########################################################################

###########################################################################
#Additional COMSOL variable definitions
###########################################################################

comsol_args={}
comsol_args['parameter']={}
comsol_args['parameter']['e0']=['1[C]','electronic charge']

system['active site density']=4.1612542339231805e-07

comsol_args['parameter']['RF']=[RF,'Roughness Factor']
comsol_args['parameter']['grid_factor']=[str(grid_factor),'Grid factor']
comsol_args['nflux']=nflux_comsol

###########################################################################
#RATE EQUATIONS/FLUXES
###########################################################################

species['H2']['flux']='catmap' #H2_rate
species['CO']['flux']='catmap' #CO_rate
species['CH3CH2OH']['flux']='catmap' #CO_rate
species['CH4']['flux']='catmap' #CO_rate
species['CO2']['flux']='catmap' #CO_rate

boundary_thickness=7.93E-05 #in m

if not nobuffer:
    visc=viscosity(species['HCO3-']['bulk_concentration']/10**3), #Pa*s at 25C of KHCO3 solution
system['boundary thickness']=boundary_thickness
#system['electrolyte viscosity']=visc[0]

#descriptor method
#comsol_args['desc_method']='external' #internal-cont'
#comsol_args['model_type']='tp_dilute_species'
#comsol_args['solver']='parametric'
###########################################################################
#BOUNDARY CONDITIONS FOR PBE
###########################################################################

potentials=[-1.0] #,-0.75,-0.5,-0.25,0.0]
results=[]

for potential in potentials:
    descriptors={'phiM':list(np.linspace(-0.4144,-0.4144,nphi))}
    system['phiM']=potential

    #'potential','gradient','robin'
    pb_bound={
    #        'potential': {'wall':'zeta'},
    #        'gradient': {'bulk':0.0}}
        'potential':{'bulk':0.0},
        'wall':system['phiM']}

    catmap_args={}
    catmap_args['n_inter']='automatic'

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
        for p in descriptors['phiM']: #np.linspace(-0.9344,-0.9344,1):
            print '!!! now running p = '+str(p)
            tp.system['phiM']=p
            cm.run()
    else:
        c=Calculator(transport=tp,tau_scf=1e-6,ntout=1,dt=1e-1,tmax=10)
        c.run()
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

