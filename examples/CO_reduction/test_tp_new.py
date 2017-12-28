from catint.transport import Transport
from catint.calculator import Calculator 
from catint.plot import Plot
import numpy as np
import sys
from units import *
from read_data import read_data

pH_i=6.8
nobuffer=False #True #False #True #False #True #False #True

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
                            'constant':             (4.44e-7)*1000.0,
                            'rates':                [3.7e-2,3.7e-2/(4.44e-7*1000.0)]},                      #K1a (=Kc=K0)
    'buffer-base':      {   'reaction':            'CO2 + OH- <-> HCO3-', 
                            'constant':             (4.44e7)/1000.0,                        #K1b (=Kc!=K0, since unit conversion factor is missing)
                            'rates':                [(5.93e3)/1000.0,(5.93e3)/(4.44e7)]},   #"k1f, k1r"
    'buffer-base2':     {   'reaction':            'HCO3- + OH- <-> CO32- + H2O', 
                            'constant':             (4.66e3)/1000.0,
                            'rates':                [(1.0e8)/1000.0,(1.0e8)/(4.66e3)]},     #"k2f,k2r"
#    'buffer2':          {   'reaction':            'CO2 + CO32- + H2O <->  2 HCO3-', 
#                            'constant':             9.52e3}                                 #K3
    'self-dissociation of water':            {   'reaction':             'H2O <-> OH- + H+',
                            'constant':             1e-14,
                            'rates':                [1.3e8*1e-14,1.3e8]} # from https://en.wikipedia.org/wiki/Self-ionization_of_water
    }
electrode_reactions={
    'C1':   {'reaction': 'CO + 5 H2O + 6 e- -> C1 + 6 OH-'}, #methane
    'C2':   {'reaction': '2 CO + 7 H2O + 8 e- -> C2 + 8 OH-'}} #ethanol

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
    'phiPZC': 0.0,
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

##1) option: initialize with CO2_i and OHm_i
## Initial composition of the bulk electrolyte at t=0
#HCO3m_i = bic_i-CO32m_i #initial HCO3- bulk concentrations at t=0 [mol/m3]
#K_i = bic_i #initial K+ bulk concentrations at t=0 [mol/m3]
#OHm_i = HCO3m_i/electrolyte_reactions['buffer-base']['constant']/CO2_i #initial OH- bulk concentrations at t=0 [mol/m3]
#pH_i = 14+np.log10(OHm_i/1000.0) #initial pH (in log. arg must be conc in M)

##2) option: initialize with HCO3m_i and OHm_i
print 'CO2 before',CO2_i
OHm_i=10**(pH_i-14.)*1000.0
HCO3m_i=0.1*1000.
CO32m_i=electrolyte_reactions['buffer-base2']['constant']*HCO3m_i*OHm_i
CO2_i=HCO3m_i/OHm_i/electrolyte_reactions['buffer-base']['constant']
print 'CO2 after',CO2_i

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
    K_i = OHm_i
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
#                            'diffusion':            4.50e-009,
#                            'zeff':                 2.0,
#                            'req':                  [1]},
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
###########################################################################

###########################################################################
#Additional COMSOL variable definitions
###########################################################################


comsol_params={}
comsol_params['A']=['1.e13[1/s]','Exponential prefactor']
#comsol_params['A2']=['1.e13[1/s]','Exponential prefactor']
comsol_params['Ga_CHO']=[str(1.11746219*unit_F)+'[J/mol]','CHO Activation Energy']
comsol_params['Ga_CHOH']=[str(2.37467774*unit_F)+'[J/mol]','CHOH Activation Energy']
comsol_params['Ga_OCCO']=[str(0.578959276*unit_F)+'[J/mol]','OCCO Activation Energy']
comsol_params['Ga_OCCOH']=[str(1.10495851*unit_F)+'[J/mol]','OCCOH Activation Energy']
#comsol_params['eVToJmol']=[str(eVTokcal*1000*calToJ)+'[J/eV/mol]','eV to J/mol Conversion factor']
comsol_params['alpha_CHO']=['0.5','Butler-Volmer Parameter'] #std 0.5
comsol_params['alpha_CHOH']=['0.5','Butler-Volmer Parameter'] #std 2.0
comsol_params['alpha_OCCOH']=['0.5','Butler-Volmer Parameter'] #std 0.5
comsol_params['alpha_OCCO']=['0.5','Butler-Volmer Parameter'] #std 0.5
comsol_params['n_CHO']=['1','Butler-Volmer Parameter'] #std 0.5
comsol_params['n_CHOH']=['1','Butler-Volmer Parameter'] #std 2.0
comsol_params['n_OCCOH']=['0','Butler-Volmer Parameter'] #std 0.5
comsol_params['n_OCCO']=['0','Butler-Volmer Parameter'] #std 0.5
comsol_params['e0']=['1[C]','electronic charge']
comsol_params['erho']=['80.3e-6[C/cm^2]','surface density of active sites x elementary charge']
comsol_params['Lmol']=['1[l/mol]','conversion factor']
#active site density. singh paper: 7.04e-6[mol/m^2]
#here: 3x3 Cu211 cell as example. area=6.363x7.794*1e-20, active sites=3 (step top sites, 1.004495558139274e-05), 9 (all top sites, 3.013486674417822e-05)
comsol_params['rho_act']=['1.004495558139274e-05[mol/m^2]','Density of Active Sites'] #from Singh paper: 7.04e-6
comsol_params['Ga_CO_ads']=[str(-0.3*unit_F)+'[J/mol]','Adsorption barrier for CO on Cu211']
comsol_params['Kads']=['exp(-Ga_CO_ads/RT)','Equilibrium constant for CO adsorption']
comsol_params['max_coverage']=['0.44','Maximal coverage with which the Langmuir isotherm will be scaled'] #0.44

###########################################################################
#RATE EQUATIONS/FLUXES
###########################################################################



#give rates as COMSOL equations
#all variables used here have to be defined as COMSOL params as seen before
#C1_rate='rho_act*coverage*A*\
#        exp(-\
#            max(\
#                (Ga_CHO+(alpha1+'+str(nCHO)+')*(phiM-phi)*F_const)/RT+alpha1*(7+log10(max([[OH-]],0.0)*Lmol))*'+str(np.log(10.))+', \
#                (Ga_CHOH+(alpha2+'+str(nCHOH)+')*(phiM-phi)*F_const)/RT+alpha2*(7+log10(max([[OH-]],0.0)*Lmol))*'+str(np.log(10.))+\
#            ')\
#        )'#/F_const/'+str(species['C1']['zeff'])
#C2_rate='rho_act*coverage^2*A*\
#        exp(-\
#            max(\
#                (Ga_OCCOH+(alpha1+'+str(nOCCOH)+')*(phiM-phi)*F_const)/RT+alpha1*(7+log10(max([[OH-]],0.0)*Lmol))*'+str(np.log(10.))+', \
#                (Ga_OCCO+'+str(nOCCO)+'*(phiM-phi)*F_const)/RT'\
#            ')\
#        )'sigma_max=5.0 #smoothness factor for maximum, the larger it is the smoother the function is approximated

comsol_params['OH_min']=['1e-30 [mol/m^3]','Minimal OH- concentration allowed in the evaluation of the rates'] #of the rate coverage with which the Langmuir isotherm will be scaled'] #0.44


comsol_variables={}
comsol_variables['coverage']=['Kads*[[CO]]*Lmol/(1.+[[CO]]*Lmol*Kads)*max_coverage','CO Coverage according to Langmuir isotherm']
#comsol_variables['jCHO']=['rho_act*coverage*A*(max([[OH-]],OH_min)*Lmol)^(alpha_CHO)*'+\
#                         'exp(-'+\
#                            '(Ga_CHO+(alpha_CHO+n_CHO)*(phiM-phi)*F_const)/RT+alpha_CHO*(7)*log(10)'+\
#                         ')','rate of CHO']
#comsol_variables['jCHOH']=['rho_act*coverage*A*(max([[OH-]],OH_min)*Lmol)^(alpha_CHO)*'+\
#                         'exp(-'+\
#                            '(Ga_CHOH+(alpha_CHOH+n_CHOH)*(phiM-phi)*F_const)/RT+alpha_CHOH*(7)*log(10)'+\
#                         ')','rate of CHOH']
#comsol_variables['jOCCOH']=['rho_act*coverage^2*A*(max([[OH-]],OH_min)*Lmol)^(alpha_CHO)*'+\
#                        'exp(-'+\
#                            '(Ga_OCCOH+(alpha_OCCOH+n_OCCOH)*(phiM-phi)*F_const)/RT+alpha_OCCOH*(7)*log(10)'+\
#                        ')','rate of OCCOH']
#comsol_variables['jOCCO']=['rho_act*coverage^2*A*(max([[OH-]],OH_min)*Lmol)^(alpha_CHO)*'+\
#                        'exp(-'+\
#                            '(Ga_OCCO+(alpha_OCCO+n_OCCO)*(phiM-phi)*F_const)/RT+alpha_OCCO*(7)*log(10)'+\
#                        ')','rate of OCCO']
comsol_variables['jCHO']=['rho_act*coverage*A*'+\
                         'exp('+\
                            '-(Ga_CHO+(alpha_CHO+n_CHO)*(phiM-phi)*F_const)/RT+alpha_CHO*(7+log10(max([[OH-]],OH_min)*Lmol))*log(10)'+\
                         ')','rate of CHO']
comsol_variables['jCHOH']=['rho_act*coverage*A*'+\
                         'exp('+\
                            '-(Ga_CHOH+(alpha_CHOH+n_CHOH)*(phiM-phi)*F_const)/RT+alpha_CHOH*(7+log10(max([[OH-]],OH_min)*Lmol))*log(10)'+\
                         ')','rate of CHOH']
comsol_variables['jOCCOH']=['rho_act*coverage^2*A*'+\
                        'exp('+\
                            '-(Ga_OCCOH+(alpha_OCCOH+n_OCCOH)*(phiM-phi)*F_const)/RT+alpha_OCCOH*(7+log10(max([[OH-]],OH_min)*Lmol))*log(10)'+\
                        ')','rate of OCCOH']
comsol_variables['jOCCO']=['rho_act*coverage^2*A*'+\
                        'exp('+\
                            '-(Ga_OCCO+(alpha_OCCO+n_OCCO)*(phiM-phi)*F_const)/RT+alpha_OCCO*(7+log10(max([[OH-]],OH_min)*Lmol))*log(10)'+\
                        ')','rate of OCCO']
#additional comsol outputs
#name, equation, unit
comsol_outputs=[\
        ['jCHO','jCHO','mol/m^2/s'],\
        ['jCHOH','jCHOH','mol/m^2/s'],\
        ['jOCCOH','jOCCOH','mol/m^2/s'],\
        ['jOCCO','jOCCO','mol/m^2/s']] #last one in list is the name of the file, first one is variable name

method=0 #method2 with stationary solver only working one, so far...

if method==0: 
    C1_rate='min('+\
                'jCHO, '+\
                'jCHOH'+\
            ')'#/F_const/'+str(species['C1']['zeff'])
    C2_rate='min('+\
                'jOCCOH, '+\
                'jOCCO'+\
            ')'
elif method==1:
    comsol_params['sigma_min']=['10.0','Smooth the calculation of the maximum by exponential summation. The larger this is, the smoother the maximum function is around the maximum']
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
                    '(Ga_CHO+(alpha_CHO+n_CHO)*(phiM-phi)*F_const)/RT+alpha_CHO*(7+log10(max([[OH-]],OH_min)*Lmol))*log(10), '+\
                    '(Ga_CHOH+(alpha_CHOH+n_CHOH)*(phiM-phi)*F_const)/RT+alpha_CHOH*(7+log10(max([[OH-]],OH_min)*Lmol))*log(10)'+\
                ')'+\
            ')'#/F_const/'+str(species['C1']['zeff'])
    C2_rate='rho_act*coverage^2*A*'+\
            'exp(-'+\
                'max('+\
                    '(Ga_OCCOH+(alpha_OCCOH+n_OCCOH)*(phiM-phi)*F_const)/RT+alpha_OCCOH*(7+log10(max([[OH-]],OH_min)*Lmol))*log(10), '+\
                    '(Ga_OCCO+(alpha_OCCO+n_OCCO)*(phiM-phi)*F_const)/RT+alpha_OCCO*(7+log10(max([[OH-]],OH_min)*Lmol))*log(10)'+\
                ')'+\
            ')'
elif method==3:
    comsol_params['avoid_div_zero']=['0.0','Avoid division by zero by adding small value on top of numerator and denominator']
    comsol_params['sigma_max']=['5.0','Smooth the calculation of the maximum by exponential summation. The larger this is, the smoother the maximum function is around the maximum']
    C1_rate='rho_act*coverage*A*'+\
            '('+\
                '(1.+avoid_div_zero)/'+\
                '('+\
                    '((1+avoid_div_zero)/(max([[OH-]],OH_min)*Lmol+avoid_div_zero))^(sigma_max*alpha_CHO)*'+\
                        'exp(-'+\
                            'sigma_max*('+\
                                '(Ga_CHO+(alpha_CHO+n_CHO)*(phiM-phi)*F_const)/RT+alpha_CHO*7*log(10)'+\
                            ')'+\
                        ')+'+\
                    '((1+avoid_div_zero)/(max([[OH-]],OH_min)*Lmol+avoid_div_zero))^(sigma_max*alpha_CHOH)*'+\
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
                    '((1+avoid_div_zero)/(max([[OH-]],OH_min)*Lmol+avoid_div_zero))^(sigma_max*alpha_OCCOH)*'+\
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

electrode_reactions['C1']['rates']=[C1_rate,'0.0']
electrode_reactions['C2']['rates']=[C2_rate,'0.0']

boundary_thickness=7.93E-05 #in m

if not nobuffer:
    visc=viscosity(species['HCO3-']['bulk concentration']/10**3), #Pa*s at 25C of KHCO3 solution
system['boundary thickness']=boundary_thickness
#system['electrolyte viscosity']=visc[0]

###########################################################################
#BOUNDARY CONDITIONS FOR PBE
###########################################################################

potentials=[-1.0] #,-0.75,-0.5,-0.25,0.0]
results=[]
for potential in potentials:
    descriptors={'phiM':list(np.linspace(0.0,-1.2,130))}
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
            comsol_params=comsol_params,
            comsol_variables=comsol_variables,
            comsol_outputs=comsol_outputs,
            descriptors=descriptors,
            nx=200)
    else:
        tp=Transport(
            species=species,
            electrode_reactions=electrode_reactions,
            electrolyte_reactions=electrolyte_reactions,
            system=system,
            pb_bound=pb_bound,
            comsol_params=comsol_params,
            comsol_variables=comsol_variables,
            comsol_outputs=comsol_outputs,
            descriptors=descriptors,
            nx=200)
    
    
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

