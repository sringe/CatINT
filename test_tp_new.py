from transport import Transport

#t=Transport(integrator='Crank-Nicolson')

###########################################################################
#REACTIONS
###########################################################################
#reactants:             (first line are the educts, second the products)
#constant:              (dimensionless (mol/m^3))
#rates:                 (forward and backward rates)

reactions=\
    {
    'buffer':           {   'reactants':            [['CO2','H2O'],['H2CO3']],
                            'constant':             2.63e-003},
    'buffer-acid':      {   'reactants':            [['CO2','H2O'],['HCO3-','H+']],
                            'constant':             (4.44e-007)*1000.0},
    'buffer-base':      {   'reactants':            [['CO2','OH-'],['HCO3-']],
                            'constant':             (4.44e007)/1000.0,
                            'rates':                [(5.93e003)/1000.0,(5.93e003)/(4.44e007)]},
    'buffer-base2':     {   'reactants':            [['HCO3-','OH-'],['CO32-','H2O']],
                            'constant':             (4.66e003)/1000.0,
                            'rates':                [(1.0e008)/1000.0,(1.0e008)/(4.66e003)]},
    'buffer2':          {   'reactants':            [['CO2','CO32-','H2O'],['HCO3-','HCO3-']],
                            'constant':             9.52e003}
    }
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
    'electrolyte viscosity':    xx/1000. #Pa*s at 25C of KHCO3 solution
    }
###########################################################################

###########################################################################
#INITIAL CONCENTRATIONS
###########################################################################
#set up the initial concentrationss from this constants:
bic_i=1.0 #mol/l
bic_i*=1000. #convert to mol/m^3
CO2_i = 0.03419*P*1000. #initial CO2(aq) bulk concentrations at t=0 and Pressure P in [mol/m3] units
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
#Henry              [atm/M] (from http://butane.chem.uiuc.edu/pshapley/GenChem1/L23/web-L23.pdf)
species=\
    {
    'K':                {   'name':                 r'K^+',
                            'diffusion':            1.957e-009,
                            'bulk concentration':   K_i},
    'CO2':              {   'name':                 r'CO_2',
                            'diffusion':            1.91e-009,
                            'bulk concentration':   CO2_i},
    'CO32-':            {   'name':                 r'CO_3^{2-}',
                            'diffusion':            9.23e-010,
                            'bulk concentration':   CO32m_i},
    'HCO3-':            {   'name':                 r'HCO_3^-',
                            'diffusion':            1.185e-009,
                            'bulk concentration':   HCO3m_i},
    'OH-':              {   'name':                 r'OH^-',
                            'diffusion':            5.273e-009,
                            'bulk concentraiton':   OHm_i},
    'H+':               {   'name':                 r'H^+',
                            'bulk concentration':   10**(-pH_i)},
    'H2':               {   'name':                 r'H2',
                            'diffusion':            4.50e-009,
                            'Henry':                1282.05},
    'CO':               {   'name':                 r'CO',
                            'diffusion':            2.03e-009,
                            'Henry':                1052.63},
    'CH4':              {   'name':                 r'CH_4',
                            'diffusion':            1.49e-009},
    'C2H4':             {   'name':                 r'C_2H_4',
                            'diffusion':            1.87e-009},
    'HCOO-':            {   'name':                 r'HCOO^-',
                            'diffusion':            1.454e-009},
    'etol':             {   'name':                 'Ethanol',
                            'diffusion':            0.84e-009},
    'propol':           {   'name':                 'n-Propanol',
                            'diffusion':            1.3e-009},
    'allyl':            {   'name':                 'Allylalcohol',
                            'diffusion':            1.1e-009},
    'metol':            {   'name':                 'Methanol',
                            'diffusion':            0.84e-009},
    'acet':             {   'name':                 'Acetate',
                            'diffusion':            1.089e-009},
    'etgly':            {   'name':                 'Ethyleneglycol',
                            'diffusion':            1.102e-009}      #at 0.02 fraction of dlycol
    }
###########################################################################

###########################################################################
#SETUP AND RUN
###########################################################################
t=Transport(
    species=species,
    reactions=reactions,
    system=system,
    integrator='odeint', #,
    integrator_correction=False, #True,
    dt=1e-9,
    tmax=1e-6)

#t=Transport(integrator='FTCS')
t.run()
t.plot()
###########################################################################

