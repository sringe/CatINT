from catint.transport import Transport
from catint.calculator import Calculator 
from catint.plot import Plot
import numpy as np
import sys
from units import *

pH_i=13.0
nobuffer=True #False #True #False #True #False #True #False #True #False #True #False #True #False #True #False #True #False #True

use_elreac=True
if nobuffer:
    use_elreac=False
###########################################################################
#REACTIONS
###########################################################################
#reactants:             (first line are the educts, second the products)
#constant:              (dimensionless (mol/m^3))
#rates:                 (forward and backward rates)

electrolyte_reactions=\
    {
#    'buffe':            {   'reaction':             'CO2 + H2O <-> H2CO3', 
#                            'constant':             2.63e-3},                               #KH
    'buffer-acid':      {   'reaction':            'CO2 + H2O <-> HCO3- + H+', 
#                            'constant':             (4.44e-7)*1000.0,
    #                        'rates':                [3.7e-2,3.7e-2/(4.44e-7*1000.0)]},                      #K1a (=Kc=K0)
                            'constant':             0.00138951310, #Schulz2006
                             'rates':               [3.71e-2,2.67e4/1000.]}, #from Schulz2006
    'buffer-acid2':     {   'reaction':             'HCO3- <-> CO32- + H+',
                            'constant':             1.1888e-06, #Schulz2006
                            'rates':                [59.44,5e10/1000.]}, #from Schulz2006
    'buffer-base':      {   'reaction':            'CO2 + OH- <-> HCO3-', 
#                            'constant':             (4.44e7)/1000.0,                        #K1b (=Kc!=K0, since unit conversion factor is missing)
#                            'rates':                [(5.93e3)/1000.0,(5.93e3)/(4.44e7)]},   #"k1f, k1r"
                            'constant':             22966.014418, #Schulz2006
                            'rates':                [2.23e3/1000.,9.71e-5]}, #Schulz2006
    'buffer-base2':     {   'reaction':            'HCO3- + OH- <-> CO32- + H2O', 
#                            'constant':             (4.66e3)/1000.0,
#                            'rates':                [(1.0e8)/1000.0,(1.0e8)/(4.66e3)]},     #"k2f,k2r"
                            'constant':             19.60784, #Schulz2006
                            'rates':                [6e9/1000.,3.06e5]}, #Schulz2006
#    'buffer2':          {   'reaction':            'CO2 + CO32- + H2O <->  2 HCO3-', 
#                            'constant':             9.52e3}                                 #K3
    'self-dissociation of water':            {   'reaction':             'H2O <-> OH- + H+',
                            'constant':             1e-14,
                            'rates':                [1.3e8*1e-14,1.3e8]} # from https://en.wikipedia.org/wiki/Self-ionization_of_water
    }
electrode_reactions={
    'H2-1':   {'reaction': '2 H2O + 2 e- -> H2 + 2 OH-'},
    'H2-2':   {'reaction': '2 H+ + 2 e- -> H2'}}

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
    'phiPZC': 0.0,
    'Stern capacitance': 20 #std: 20
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
#zeff number of electrons need to derive product from educt
#req number of educt molecules needed to derive product
#Henry              [atm/M] (from http://butane.chem.uiuc.edu/pshapley/GenChem1/L23/web-L23.pdf)
species=\
    {
    'K':                {   'symbol':               'K^+',
                            'name':                 'potassium',
                            'diffusion':            1.957e-009,
                            'kind':                 'electrolyte',
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
                            'kind':                 'electrolyte',
                            'flux':                 1e-8,
                            'zeff':
                            'bulk concentration':   OHm_i},
#    'H+':               {   'symbol':               'H^+',
#                            'name':                 'hydronium',
#                            'bulk concentration':   10**(-pH_i)},
    'H2':               {   'symbol':               'H_2',
                            'name':                 'hydrogen',
                            'diffusion':            4.50e-009,
                            'flux':                 {'kind':'educt','values':['OH-'],'reqs':[]},\
#                            'req':                  [1]},
    }}
###########################################################################


boundary_thickness=7.93E-05 #in m

#visc=viscosity(species['HCO3-']['bulk concentration']/10**3), #Pa*s at 25C of KHCO3 solution
system['boundary thickness']=boundary_thickness
#system['electrolyte viscosity']=visc[0]

###########################################################################
#BOUNDARY CONDITIONS FOR PBE
###########################################################################

descriptors={'pH':range(1,15,1)}
system['phiM']=-0.5

pb_bound={
    'potential':{'bulk':0.0},
    'wall':system['phiM']}


###########################################################################
#SETUP AND RUN
###########################################################################
tp=Transport(
    species=species,
    electrode_reactions=electrode_reactions,
    electrolyte_reactions=electrolyte_reactions,
    system=system,
    pb_bound=pb_bound,
    comsol_params=comsol_params,
    descriptors=descriptors,
    nx=40)


tp.set_calculator('comsol') #odespy') #--bdf')
#tp.set_calculator('Crank-Nicolson--LF')
#tp.set_initial_concentrations('Gouy-Chapman')

c=Calculator(transport=tp,tau_jacobi=1e-5,ntout=1,dt=1e-1,tmax=10)
#scale_pb_grid
cout=c.run() #1.0)

p=Plot(transport=tp)
p.plot(cout)
    #results.append(potential,tp.species['CO2']['flux'])
#plt.plot(results,'-o')
#plt.show()
    #p.add_plot('polarization')
#    p.plot(cout)
 #   p.plot_polarization()
    ###########################################################################

