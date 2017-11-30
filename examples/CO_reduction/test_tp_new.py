from catint.transport import Transport
from catint.calculator import Calculator 
from catint.plot import Plot
import numpy as np
import sys
from units import *

###########################################################################
#REACTIONS
###########################################################################
#reactants:             (first line are the educts, second the products)
#constant:              (dimensionless (mol/m^3))
#rates:                 (forward and backward rates)



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
    'educts': ['CO'],#,'OH'],  #the educt which is converted to products
    'products': [['C1','C2']],#,['H2']],
    'temperature':  298,     #K
    'pressure':     1.,      #atm
    'water viscosity':  8.90e-004, #Pa*s at 25C
    #calculate the electrolyte viscosity. This will be used to rescale the diffusion coefficients
    #according to Einstein-Stokes relation: D_in_electrolyte = D_in_water * mu0/mu
    'epsilon': 78.36,
   # 'exclude species': ['H+'], #exclude this species from PNP equations
    'migration': True
    }
###########################################################################

###########################################################################
#INITIAL CONCENTRATIONS
###########################################################################
#set up the initial concentrationss from this constants:
#CO2_i = 0.03419*system['pressure']*1000. #initial CO2(aq) bulk concentrations at t=0 and Pressure P in [mol/m3] units
CO_i = 9.5e-4*system['pressure']*1000.
                        #from Henry constant (29.41 atm/M
#CO32m_i = ((2*bic_i+reactions['buffer2']['constant']*CO2_i)-\
#            (np.sqrt((2*bic_i+reactions['buffer2']['constant']*CO2_i)**2\
#            -4.0*(bic_i)**2)))/2  #initial (CO3)2- bulk concentrations at t=0 [mol/m3]
# Initial composition of the bulk electrolyte at t=0
#HCO3m_i = bic_i-CO32m_i #initial HCO3- bulk concentrations at t=0 [mol/m3]
OHm_i = 1e-7/1000 #HCO3m_i/reactions['buffer-base']['constant']/CO2_i #initial OH- bulk concentrations at t=0 [mol/m3]
K_i = OHm_i #bic_i #initial K+ bulk concentrations at t=0 [mol/m3]
pH_i = 14+np.log10(OHm_i/1000.0) #initial pH (in log. arg must be conc in M)
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
                            'diffusion':            1.957e-009,
                            'bulk concentration':   K_i},
#    'CO2':              {   'symbol':               'CO_2',
#                            'name':                 'carbon dioxide',
#                            'diffusion':            1.91e-009,
#                            'bulk concentration':   CO2_i},
#    'CO32-':            {   'symbol':               'CO_3^{2-}',
#                            'name':                 'carboxylate',
#                            'diffusion':            9.23e-010,
#                            'bulk concentration':   CO32m_i},
#    'HCO3-':            {   'symbol':               'HCO_3^-',
#                            'name':                 'bicarbonate',
#                            'diffusion':            1.185e-009,
#                            'bulk concentration':   HCO3m_i},
    'OH-':              {   'symbol':               'OH^-',
                            'name':                 'hydroxyl',
                            'diffusion':            5.273e-009,
                            'bulk concentration':   OHm_i},
#    'H+':               {   'symbol':               'H^+',
#                            'name':                 'hydronium',
#                            'bulk concentration':   10**(-pH_i)},
#    'H2':               {   'symbol':               'H_2',
#                            'name':                 'hydrogen',
#                            'diffusion':            4.50e-009,
#                            'zeff':                 2.0,
#                            'req':                  [1]},
    'CO':               {   'symbol':               'CO',
                            'name':                 'carbon monoxide',
                            'diffusion':            2.03e-009,
                            'zeff':                 2.0,
                            'bulk concentration':   CO_i,
                            'req':                  [1]},
    'C1':               {   'symbol':               'CH4',
                            'name':                 'C1-methane',
                            'diffusion':            1.49e-009,
                            'zeff':                 8.0,
                            'req':                  [1]},
    'C2':               {   'symbol':               'C_2H_5OH',
                            'name':                 'C2-ethanol',
                            'diffusion':            0.84e-009,
                            'zeff':                 12.0,
                            'req':                  [2]}
    }
###########################################################################

###########################################################################
#Additional COMSOL variable definitions
###########################################################################

comsol_params={}
comsol_params['A1']=['1.e13[1/s]','Exponential prefactor']
comsol_params['A2']=['1.e13[m^3/s/mol]','Exponential prefactor']
comsol_params['Ga_CHO']=['1.11746219[eV]','CHO Activation Energy']
comsol_params['Ga_CHOH']=['2.37467774[eV]','CHOH Activation Energy']
comsol_params['Ga_OCCO']=['0.578959276[eV]','OCCO Activation Energy']
comsol_params['Ga_OCCOH']=['1.10495851[eV]','OCCOH Activation Energy']
comsol_params['eVToJmol']=[str(eVTokcal*1000*calToJ)+'[J/eV/mol]','eV to J/mol Conversion factor']
comsol_params['alpha1']=['0.5','Butler-Volmer Parameter']
comsol_params['alpha2']=['2.0','Butler-Volmer Parameter']
comsol_params['e0']=['1[C]','electronic charge']


###########################################################################
#RATE EQUATIONS/FLUXES
###########################################################################

#give rates as COMSOL strings
CO_rate=''
OHm_rate=''
C1_rate='[[CO]]*A1*exp(-max(Ga_CHO+alpha1*V*e0, Ga_CHOH+alpha2*V*e0)*eVToJmol/RT)'#/F_const/'+str(species['C1']['zeff'])
C2_rate='[[CO]]^2*A2*exp(-max(Ga_OCCOH+alpha1*V*e0, Ga_OCCO)*eVToJmol/RT)' #/F_const/'+str(species['C2']['zeff'])
formulas=[C1_rate,C2_rate]
for product,formula in zip(system['products'][0],formulas):
    CO_rate+='-'+str(species[product]['req'][0])+'*'+formula #+'/'+str(unit_F*species[product]['zeff'])+'*'+formula
    OHm_rate+='+'+formula

#OHm_rate+='-'+

#CO_rate='-'+str(species['C1']['req'][0])+'/'+str(species['C1']['zeff'])+'[[CO]]*A*exp(-max([Ga_CHO+0.5*V, Ga_CHOH+2*V])/RT*'+str(eVToJmol)+')-[[CO]]**2*A*exp(-max([Ga_OCCOH+0.5*V, Ga_OCCO])/RT*'+str(eVToJmol)+')'

data_fluxes={}
data_fluxes['C1']=C1_rate
data_fluxes['C2']=C2_rate
data_fluxes['CO']=CO_rate
data_fluxes['OH-']=OHm_rate

for key in species:
    if key in data_fluxes: #'flux' in species[key]:
        species[key]['flux']=data_fluxes[key]

boundary_thickness=7.93E-05 #in m

#visc=viscosity(species['HCO3-']['bulk concentration']/10**3), #Pa*s at 25C of KHCO3 solution
system['boundary thickness']=boundary_thickness
#system['electrolyte viscosity']=visc[0]

###########################################################################
#BOUNDARY CONDITIONS FOR PBE
###########################################################################

#'potential','gradient','robin'
pb_bound={
#        'potential': {'wall':'zeta'},
#        'gradient': {'bulk':0.0}}
    'potential':{'bulk':0.0},
    'gradient': {'wall':0.0}} 

###########################################################################
#SETUP AND RUN
###########################################################################
tp=Transport(
    species=species,
#    reactions=reactions,
    system=system,
    pb_bound=pb_bound,
    comsol_params=comsol_params,
    nx=40)


tp.set_calculator('comsol') #odespy') #--bdf')
#tp.set_calculator('Crank-Nicolson--LF')
#tp.set_initial_concentrations('Gouy-Chapman')

c=Calculator(transport=tp,tau_jacobi=1e-5,ntout=1)
#scale_pb_grid
cout=c.run(dt=1e-3,tmax=20.0) #1.0)

p=Plot(transport=tp)
p.plot(cout)

###########################################################################

