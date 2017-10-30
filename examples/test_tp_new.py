from transport.transport import Transport
from transport.calculator import Calculator 
from transport.plot import Plot
from read_data import read_data
import numpy as np
import sys

###########################################################################
#REACTIONS
###########################################################################
#reactants:             (first line are the educts, second the products)
#constant:              (dimensionless (mol/m^3))
#rates:                 (forward and backward rates)

reactions=\
    {
    'buffe':           {   'reactants':            [['CO2','H2O'],['H2CO3']],
                            'constant':             2.63e-3},                               #KH
    'buffer-acid':      {   'reactants':            [['CO2','H2O'],['HCO3-','H+']],
                            'constant':             (4.44e-7)*1000.0},                      #K1a
    'buffer-base':      {   'reactants':            [['CO2','OH-'],['HCO3-']],
                            'constant':             (4.44e7)/1000.0,                        #K1b
                            'rates':                [(5.93e3)/1000.0,(5.93e3)/(4.44e7)]},   #"k1f, k1r"
    'buffer-base2':     {   'reactants':            [['HCO3-','OH-'],['CO32-','H2O']],
                            'constant':             (4.66e3)/1000.0,
                            'rates':                [(1.0e8)/1000.0,(1.0e8)/(4.66e3)]},     #"k2f,k2r"
    'buffer2':          {   'reactants':            [['CO2','CO32-','H2O'],['HCO3-','HCO3-']],
                            'constant':             9.52e3}                                 #K3
    }
###########################################################################

###########################################################################
#READ DATA FILE
###########################################################################

data_fluxes,boundary_thickness,viscosity,bic_i=read_data()

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
    'exclude species': ['H+'], #exclude this species from PNP equations
    'migration': True
    }
###########################################################################

###########################################################################
#INITIAL CONCENTRATIONS
###########################################################################
#set up the initial concentrationss from this constants:
CO2_i = 0.03419*system['pressure']*1000. #initial CO2(aq) bulk concentrations at t=0 and Pressure P in [mol/m3] units
                        #from Henry constant (29.41 atm/M
CO32m_i = ((2*bic_i+reactions['buffer2']['constant']*CO2_i)-\
            (np.sqrt((2*bic_i+reactions['buffer2']['constant']*CO2_i)**2\
            -4.0*(bic_i)**2)))/2  #initial (CO3)2- bulk concentrations at t=0 [mol/m3]
# Initial composition of the bulk electrolyte at t=0
HCO3m_i = bic_i-CO32m_i #initial HCO3- bulk concentrations at t=0 [mol/m3]
K_i = bic_i #initial K+ bulk concentrations at t=0 [mol/m3]
OHm_i = HCO3m_i/reactions['buffer-base']['constant']/CO2_i #initial OH- bulk concentrations at t=0 [mol/m3]
pH_i = 14+np.log10(OHm_i/1000.0) #initial pH (in log. arg must be conc in M)
###########################################################################


###########################################################################
#SPECIES DATA
###########################################################################
#(the charges are determined automatically from the species name)
#diffusion          [m/s]       (infinite dilution in water at 25C)
#bulk concentrations [mol/m^3]
#zeff number of electrons need to derive product from CO2
#Henry              [atm/M] (from http://butane.chem.uiuc.edu/pshapley/GenChem1/L23/web-L23.pdf)
species=\
    {
    'K':                {   'symbol':               'K^+',
                            'name':                 'potassium',
                            'diffusion':            1.957e-009,
                            'bulk concentration':   K_i},
    'CO2':              {   'symbol':               'CO_2',
                            'name':                 'carbon dioxide',
                            'diffusion':            1.91e-009,
                            'bulk concentration':   CO2_i},
    'CO32-':            {   'symbol':               'CO_3^{2-}',
                            'name':                 'carboxylate',
                            'diffusion':            9.23e-010,
                            'bulk concentration':   CO32m_i},
    'HCO3-':            {   'symbol':               'HCO_3^-',
                            'name':                 'bicarbonate',
                            'diffusion':            1.185e-009,
                            'bulk concentration':   HCO3m_i},
    'OH-':              {   'symbol':               'OH^-',
                            'name':                 'hydroxyl',
                            'diffusion':            5.273e-009,
                            'bulk concentration':   OHm_i},
    'H+':               {   'symbol':               'H^+',
                            'name':                 'hydronium',
                            'bulk concentration':   10**(-pH_i)},
    'H2':               {   'symbol':               'H_2',
                            'name':                 'hydrogen',
                            'diffusion':            4.50e-009,
                            'zeff':                 2.0},
    'CO':               {   'symbol':               'CO',
                            'name':                 'carbon monoxide',
                            'diffusion':            2.03e-009,
                            'zeff':                 2.0},
    'CH4':              {   'symbol':               'CH_4',
                            'name':                 'methane',
                            'zeff':                 8.0,
                            'diffusion':            1.49e-009},
    'C2H4':             {   'symbol':               'C_2H_4',
                            'name':                 'ethylene',
                            'zeff':                 12.0,
                            'diffusion':            1.87e-009},
    'HCOO-':            {   'symbol':               'HCOO^-',
                            'name':                 'formate',
                            'zeff':                 2.0,
                            'diffusion':            1.454e-009},
    'etol':             {   'name':                 'ethanol',
                            'symbol':               'C_2H_5OH',
                            'zeff':                 12.0,
                            'diffusion':            0.84e-009},
    'propol':           {   'name':                 'n-propanol',
                            'symbol':               'C_3H_7OH',
                            'zeff':                 18.0,
                            'diffusion':            1.3e-009},
    'allyl':            {   'name':                 'allylalcohol',
                            'symbol':               'C_3H_5OH',
                            'zeff':                 16.0,
                            'diffusion':            1.1e-009},
    'metol':            {   'name':                 'Methanol',
                            'symbol':               'CH_3OH',
                            'zeff':                 6.0,
                            'diffusion':            0.84e-009},
    'acet':             {   'name':                 'acetate',
                            'symbol':               'CH_3COO^-',
                            'zeff':                 8.0,
                            'diffusion':            1.089e-009},
    'etgly':            {   'name':                 'ethyleneglycol',
                            'symbol':               'C_2H_4OHOH',
                            'zeff':                 10.0,
                            'diffusion':            1.102e-009},      #at 0.02 fraction of dlycol
    'unknown':          {   'name':                 'unknown',
                            'symbol':               'C_2H_2O_2',
                            'zeff':                 6.0, #assumed: glyoxal
                            'diffusion':            0.0}
    }
###########################################################################

###########################################################################
#MODIFY DATA FROM FILES
###########################################################################


potential=str(-1.16953) #-0.95526)
total_flux=0.0
for key in species:
    if key!='unknown' and key in data_fluxes: #'flux' in species[key]:
        species[key]['flux']=data_fluxes[key][potential]
        total_flux+=data_fluxes[key][potential]

if total_flux<1.0:
    species['unknown']['flux']=1.-total_flux
else:
    species['unknown']['flux']=0.0

print species['HCO3-']['bulk concentration']
visc=viscosity(species['HCO3-']['bulk concentration']/10**3), #Pa*s at 25C of KHCO3 solution
system['boundary thickness']=boundary_thickness
system['current density']=data_fluxes['current_density']['-0.95526']
system['electrolyte viscosity']=visc[0]
###########################################################################
#BOUNDARY CONDITIONS FOR PBE
###########################################################################

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
    reactions=reactions,
    system=system,
    pb_bound=pb_bound,
    nx=500)


tp.set_calculator('odespy') #--bdf')
#tp.set_calculator('Crank-Nicolson--LF')
#tp.set_initial_concentrations('Gouy-Chapman')

c=Calculator(transport=tp,tau_jacobi=1e-5,ntout=1)
#scale_pb_grid
cout=c.run(dt=1e-3,tmax=20.0) #1.0)

p=Plot(transport=tp)
p.plot(cout)

###########################################################################

