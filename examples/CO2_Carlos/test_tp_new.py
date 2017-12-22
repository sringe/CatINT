from catint.transport import Transport
from catint.calculator import Calculator 
from catint.plot import Plot
from read_data import read_data
import numpy as np
import sys

###########################################################################
#REACTIONS
###########################################################################
#important: the species used here must be also listed in the species dictionary later
#reaction:              define chemical reaction as string
#constant:              (dimensionless (mol/m^3))
#rates:                 (forward and backward rates)

electrolyte_reactions=\
    {
#    'buffe':            {   'reaction':             'CO2 + H2O <-> H2CO3', 
#                            'constant':             2.63e-3},                               #KH
#    'buffer-acid':      {   'reaction':            'CO2 + H2O <-> HCO3- + H+', 
#                            'constant':             (4.44e-7)*1000.0},                      #K1a (=Kc=K0)
    'buffer-base':      {   'reaction':            'CO2 + OH- <-> HCO3-', 
                            'constant':             (4.44e7)/1000.0,                        #K1b (=Kc!=K0, since unit conversion factor is missing)
                            'rates':                [(5.93e3)/1000.0,(5.93e3)/(4.44e7)]},   #"k1f, k1r"
    'buffer-base2':     {   'reaction':            'HCO3- + OH- <-> CO32- + H2O', 
                            'constant':             (4.66e3)/1000.0,
                            'rates':                [(1.0e8)/1000.0,(1.0e8)/(4.66e3)]},     #"k2f,k2r"
    'buffer2':          {   'reaction':            'CO2 + CO32- + H2O <->  2 HCO3-', 
                            'constant':             9.52e3}                                 #K3
    }
#reactions=electrolyte_reactions
#rates can be given as reaction rate or current densities
#if the rates are now given here, they are calculated from the species fluxes 
electrode_reactions=\
    {
    'H2':           {   'reaction':            '2 H2O + 2 e- -> H2 + 2 OH-'},
    'CO':           {   'reaction':            'CO2 + H2O + 2 e- ->  CO + 2 OH-'},
    'CH4':          {   'reaction':            'CO2 + 6 H2O + 8 e- -> CH4 + 8 OH-'},
    'C2H4':         {   'reaction':            '2 CO2 + 8 H2O + 12 e- -> C2H4 + 12 OH-'},
    'metol':        {   'reaction':            'CO2 + 5 H2O + 6 e- ->  metol + 6 OH-'},
    'acet':         {   'reaction':            '2 CO2 + 6 H2O + 8 e- -> acet + 8 OH-'}, #!!! actually we use Hacet here
    'etol':         {   'reaction':            '2 CO2 + 9 H2O + 12 e- -> etol + 12 OH-'},
    'etgly':        {   'reaction':            '2 CO2 + 8 H2O + 10 e- -> etgly + 10 OH-'},
    'allyl':        {   'reaction':            '3 CO2 + 11 H2O + 16 e- -> allyl + 16 OH-'},
    'propol':       {   'reaction':            '3 CO2 + 13 H2O + 18 e- -> propol + 18 OH-'},
    'HCOO-':        {   'reaction':            'CO2 + 2 H2O + 2 e- -> HCOO- + 2 OH-'},  #consider HCOOH here
    'unknown':      {   'reaction':            '2 CO2 + 4 H2O + 6 e- -> unknown + 6 OH-'} 
    }

#K = a(CO32-)*a(H2O)/(a(HCO3-)*a(OH-)) = c(CO32-)[mol/l]/c(HCO3-)[mol/l]/c(OH-)[mol/l]*unit_conversion = Kc*unit_conversion

###########################################################################

###########################################################################
#READ DATA FILE
###########################################################################

data_fluxes,boundary_thickness,viscosity,bic_i=read_data()

###########################################################################
#SYSTEM VARIABLES
###########################################################################
system=\
    {
    'temperature':  298,     #K
    'pressure':     1.,      #atm
    'water viscosity':  8.90e-004, #Pa*s at 25C
    #calculate the electrolyte viscosity. This will be used to rescale the diffusion coefficients
    #according to Einstein-Stokes relation: D_in_electrolyte = D_in_water * mu0/mu
    'exclude species':['H+','H2O'], #here all the species from the reactions defined above which should be not 
    #considered for the transport or rate calculations (activity of 1 e.g. for water), should be listed here
    #for equation definitions
    'epsilon': 78.36,
    'Stern capacitance': 18., #in muF/cm^2
    'migration': False,
    'electrode reactions': True,
    'electrolyte reactions': True
    }
###########################################################################

###########################################################################
#INITIAL CONCENTRATIONS
###########################################################################
#set up the initial concentrationss from this constants:
CO2_i = 0.03419*system['pressure']*1000. #initial CO2(aq) bulk concentrations at t=0 and Pressure P in [mol/m3] units
                        #from Henry constant (29.41 atm/M
CO32m_i = ((2*bic_i+electrolyte_reactions['buffer2']['constant']*CO2_i)-\
            (np.sqrt((2*bic_i+electrolyte_reactions['buffer2']['constant']*CO2_i)**2\
            -4.0*(bic_i)**2)))/2  #initial (CO3)2- bulk concentrations at t=0 [mol/m3]
# Initial composition of the bulk electrolyte at t=0
HCO3m_i = bic_i-CO32m_i #initial HCO3- bulk concentrations at t=0 [mol/m3]
K_i = bic_i #initial K+ bulk concentrations at t=0 [mol/m3]
OHm_i = HCO3m_i/electrolyte_reactions['buffer-base']['constant']/CO2_i #initial OH- bulk concentrations at t=0 [mol/m3]
pH_i = 14+np.log10(OHm_i/1000.0) #initial pH (in log. arg must be conc in M)
###########################################################################


###########################################################################
#SPECIES DATA
###########################################################################
#(the charges are determined automatically from the species name)
#diffusion          [m/s]       (infinite dilution in water at 25C)
#bulk concentrations [mol/m^3]
#zeff number of electrons need to derive product from CO2
#if the fluxes of other reactant should be summed up to get the flux of another reactant
#put a dictionary as shown below. values are the reactants which contribute to the flux,
#reqs are the reaction equivalents
#flux can be equation as string or dictionary (see above) or the value as float
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
#                            'flux':                 {'kind':'educt','values':['C2H4','CH4','CO','HCOO-','etol','propol','allyl','metol','acet','etgly','unknown'],\
#                                                    'reqs':'number_of_catoms'}},
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
#                            'flux':                 {'kind':'product','values':['C2H4','CH4','CO','HCOO-','etol','propol','allyl','metol','acet','etgly','unknown','H2'],\
#                                                    'reqs':'zeff'}},
   # 'H+':               {   'symbol':               'H^+',
   #                         'name':                 'hydronium',
   #                         'bulk concentration':   10**(-pH_i)},
    'H2':               {   'symbol':               'H_2',
                            'name':                 'hydrogen',
                            'diffusion':            4.50e-009},
    'CO':               {   'symbol':               'CO',
                            'name':                 'carbon monoxide',
                            'diffusion':            2.03e-009},
    'CH4':              {   'symbol':               'CH_4',
                            'name':                 'methane',
                            'diffusion':            1.49e-009},
    'C2H4':             {   'symbol':               'C_2H_4',
                            'name':                 'ethylene',
                            'diffusion':            1.87e-009},
    'HCOO-':            {   'symbol':               'HCOO^-',
                            'name':                 'formate',
                            'diffusion':            1.454e-009},
    'etol':             {   'name':                 'ethanol',
                            'symbol':               'C_2H_5OH',
                            'diffusion':            0.84e-009},
    'propol':           {   'name':                 'n-propanol',
                            'symbol':               'C_3H_7OH',
                            'diffusion':            1.3e-009},
    'allyl':            {   'name':                 'allylalcohol',
                            'symbol':               'C_3H_5OH',
                            'diffusion':            1.1e-009},
    'metol':            {   'name':                 'Methanol',
                            'symbol':               'CH_3OH',
                            'diffusion':            0.84e-009},
    'acet':             {   'name':                 'acetate',
                            'symbol':               'CH_3COO^-',
                            'diffusion':            1.089e-009},
    'etgly':            {   'name':                 'ethyleneglycol',
                            'symbol':               'C_2H_4OHOH',
                            'diffusion':            1.102e-009},      #at 0.02 fraction of dlycol
    'unknown':          {   'name':                 'unknown', #glyoxal
                            'symbol':               'C_2H_2O_2',
                            'diffusion':            0.0}
    }
###########################################################################

###########################################################################
#MODIFY DATA FROM FILES
###########################################################################

#set up descriptors:
#this works but it does not adapt the fluxes! there is currently no 
#descriptor based data input implemented!
#descriptors={'phiM':[[key for key in data_fluxes['CO']][0]]}
descriptors={'phiM':[-1.16953]} #[[key for key in data_fluxes['CO']][0]]}
print 'descriptors',descriptors
#set potential to the value you want
potential=str(-1.16953) #-0.95526)
total_flux=0.0
for key in species:
    if key!='unknown' and key in data_fluxes: #'flux' in species[key]:
        electrode_reactions[key]['current density']=[data_fluxes[key][potential],0.0]
        total_flux+=data_fluxes[key][potential]

#unknown species only if flux of others is smaller than 1
if total_flux<1.0:
    electrode_reactions['unknown']['current density']=[1.-total_flux,0.0]
else:
    electrode_reactions['unknown']['current density']=[0.0,0.0]



visc=viscosity(species['HCO3-']['bulk concentration']/10**3), #Pa*s at 25C of KHCO3 solution
system['boundary thickness']=boundary_thickness
#system['current density']=data_fluxes['current_density']['-0.95526']
system['electrolyte viscosity']=visc[0]

###########################################################################
#BOUNDARY CONDITIONS FOR PBE
###########################################################################

system['phiM']=potential
#by default Robin BCs will be taken at the electrode and Dirichlet in the bulk

###########################################################################
#SETUP AND RUN
###########################################################################


tp=Transport(
    species=species,
    electrolyte_reactions=electrolyte_reactions,
    electrode_reactions=electrode_reactions,
    system=system,
    descriptors=descriptors,
    nx=40)


tp.set_calculator('comsol') #odespy') #--bdf')
#tp.set_calculator('Crank-Nicolson--LF')
#tp.set_initial_concentrations('Gouy-Chapman')

c=Calculator(transport=tp,tau_jacobi=1e-5,ntout=1,dt=1e-1,tmax=10.0,mode='stationary')
#scale_pb_grid
c.run()

#plots:
#plots: 'concentrations','potential','efield','current_density' (are plotted in a separate figure for all descriptors)
#descriptor_plots: same as above (the results requested will be plotted into a single figure for all descriptors)
p=Plot(transport=tp,logscale=False) #transport=tp,plots=['concentrations','potential','efield','current_density'],descriptor_plots=['_potential'])
p.plot()

tp.save() #save all results for later use
###########################################################################

