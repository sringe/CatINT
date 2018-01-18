
from catint.transport import Transport
from catint.calculator import Calculator 
from catint.plot import Plot
import numpy as np
import sys
from units import *
from read_data import read_data


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
   # 'exclude species': ['H+'], #exclude this species from PNP equations
    'migration': True,
    'electrode reactions': False,
    'electrolyte reactions': True, #False,
    'phiPZC': 0.0,
    'Stern capacitance': 200
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
                        #from Henry constant (29.41 atm/M
CO_i = 9.5e-4*system['pressure']*1000.
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
#req number of educt molecules needed to derive product
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
    }
###########################################################################

boundary_thickness=7.93E-05 #in m

visc=viscosity(species['HCO3-']['bulk concentration']/10**3), #Pa*s at 25C of KHCO3 solution
system['boundary thickness']=boundary_thickness
system['electrolyte viscosity']=visc[0]

###########################################################################
#BOUNDARY CONDITIONS FOR PBE
###########################################################################

potentials=[-1.0] #,-0.75,-0.5,-0.25,0.0]
results=[]
for potential in potentials:
    descriptors={'phiM':[-1.0]}
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
    tp=Transport(
        species=species,
        electrolyte_reactions=electrolyte_reactions,
        system=system,
        pb_bound=pb_bound,
        descriptors=descriptors,
        nx=40)
    
    
    tp.set_calculator('comsol') #odespy') #--bdf')
    #tp.set_calculator('Crank-Nicolson--LF')
    #tp.set_initial_concentrations('Gouy-Chapman')
    
    c=Calculator(transport=tp,tau_jacobi=1e-5,ntout=1,dt=1e-1,tmax=10)
    #scale_pb_grid
    c.run() #1.0)
    
    p=Plot(transport=tp)
    p.plot()
    #results.append(potential,tp.species['CO2']['flux'])
#plt.plot(results,'-o')
#plt.show()
    #p.add_plot('polarization')
#    p.plot(cout)
 #   p.plot_polarization()
    ###########################################################################
