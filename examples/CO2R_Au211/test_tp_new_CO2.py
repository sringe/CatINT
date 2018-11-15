import sys
sys.path.insert(0,'/scratch/users/sringe/transport/catint2')
sys.path.insert(0,'/scratch/users/sringe/transport/catmap')
from shutil import copyfile as copy
from catint.transport import Transport
from catint.calculator import Calculator
from catint.plot import Plot
from catint.catmap_wrapper import CatMAP
from catint.comsol_reader import Reader
import numpy as np
from units import *
from read_data import read_data
from extrapolate_surface_conc import extrapolate

transport_mode='comsol' #None #extrapolate'
#can be one of the following:
#   None            only catmap
#   'comsol'        iterative catmap-comsol
#   'extrapolate'  use extrapolated log(c_surface) -> potential curves and run pure catmap with these

pH_i=6.8
nobuffer=False #True #False #True #False #True #False #True 

educt='CO2' #CO2 or CO

nx=200
dflux_comsol=0.01
grid_factor=100
mix_scf=0.1
nphi=None #40
dphi=0.05

include_ramp_comsol=['PZC','CS'] #,'reactions']

tau_scf=0.03 #required accuracy of current density

RF=1

min_desc_delta=0.2
max_desc_delta=0.2

grid_factor_domain=100 #grid_factor
grid_factor_bound=200 #grid_factor

include_protons=False

#put here a results folder with which the surface concentrations should be initialized
init_folder=None #'try8_w_tp_cdl_au111_eps6_beta_0.5_Ga_0.0/comsol_results_id000_0006_0001' #None #'try8_w_tp_cdl_au111_eps6_beta_0.5_Ga_0.0/comsol_results_id000_0011_0001'
#put here a catmap-comsol transport calculation which is used to extrapolate transport to other potentials
extrapol_folder=['try8_w_tp_cdl_au111_eps6_beta_0.5_Ga_0.0'] #None #['try8_w_tp_cdl_au111_eps6_beta_0.5_Ga_0.0'] #try8_w_tp_cdl_au111_eps6_beta_0.5_Ga_0.0']
#None #['try8_w_tp_cdl_au111_eps6_beta_0.5_Ga_0.0','try8_w_tp_cdl_au111_eps6_beta_0.5_Ga_0.0_2'] #try7_w_tp_cdl_au111_eps6_beta_0.5_Ga_0.0','try7_w_tp_cdl_au111_eps6_beta_0.5_Ga_0.0_2']


use_elreac=True
if nobuffer:
    use_elreac=False

###########################################################################
#REACTIONS
###########################################################################
if use_elreac:
    electrolyte_reactions=['bicarbonate-base']
    if include_protons:
        electrolyte_reactions+=['bicarbonate-acid','water-diss']

electrode_reactions={
    #'H2':           {   'reaction':            '2 H2O + 2 e- -> H2 + 2 OH-'},
    #'H2':           {   'reaction':             '2 HCO3- + 2 e- -> H2 + 2 CO32-'},
#    'H2':           {   'reaction':            '2 H2O + 2 e- -> H2 + 2 OH-'},
    'CO':           {   'reaction':             'CO2 + H2O + 2 e- -> CO + 2 OH-'},
#    'CH4':          {   'reaction':            'CO2 + 6 H2O + 8 e- -> CH4 + 8 OH-'},
#    'CH3CH2OH':     {   'reaction':            '2 CO2 + 9 H2O + 12 e- -> CH3CH2OH + 12 OH-'},
#    'HCOOH':        {    'reaction':            'CO2 + 2 H2O + 2 e- ->  HCOOH + 2 OH-'}
#    'CH2O':         {   'reaction':            'CO2 + 3 H2O + 4 e- -> CH2O + 4 OH-'}}
   # 'C2H4':         {   'reaction':            '2 CO2 + 8 H2O + 12 e- -> C2H4 + 12 OH-'}}
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
    'epsilon': 78.36,
    #'epsilon_func': 'Booth',
#    'exclude species': ['CO32-','HCO3-'], #exclude this species from PNP equations
    'migration': True,
    'electrode reactions': True,
    'electrolyte reactions': use_elreac, #False,
    'phiPZC': 0.2, #ModernAspects of Electrochemistry Books/, value in water
    'bulk_pH':pH_i,
    'init_folder':init_folder,
    'potential drop':'Stern', #either Stern or full
    'Stern capacitance': 30, #std: 20, Journal of Electroanalytical Chemistry 414 (1996) 209-220
    'Stern epsilon':2 #value or Booth
    }
###########################################################################

###########################################################################
#READ DATA FILE
###########################################################################

data_fluxes,boundary_thickness,viscosity,bic_i=read_data()

OHm_i=10**(pH_i-14.)*1000.0
Hm_i=10**(-pH_i)*1000.0
if not include_protons:
    Hm_i=0.0

###########################################################################


###########################################################################
#SPECIES DATA
###########################################################################

species=\
    {
    'K+':             {'bulk_concentration':   'charge_neutrality',\
                        'MPB_radius':           2*4.1e-10},
    #'Cl-':            {'bulk_concentration':    (1.45-0.09)*1000.},
    #'HCO3-':            {'bulk_concentration':  91.0944666093},
    #'CO32-':            {'bulk_concentration':  0.0267841528009},
    #'K+':               {'bulk_concentration':  91.1480980107,\
    #                      'MPB_radius':         6.62e-10},
    'CO2':              {'bulk_concentration':   'Henry'},
    'OH-':              {'bulk_concentration':   OHm_i},
    #'H+':               {'bulk_concentration':  Hm_i},
#    'HCO3-':            {'bulk_concentration':  0.1*1000.},
    'CO':               {'bulk_concentration':0.0},
#    'CO2':              {}
    }

if include_protons:
    species.update(\
    {
    'H+':               {'bulk_concentration':   Hm_i}
    })

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
comsol_args['solver_settings']={}
comsol_args['solver_settings']['direct']={}
comsol_args['solver_settings']['direct']['nliniterrefine']=True
comsol_args['solver_settings']['ramp']={}
comsol_args['solver_settings']['ramp']['names']=include_ramp_comsol
comsol_args['solver_settings']['ramp']['dramp']=dflux_comsol
#comsol_args['solver_settings']['solver_sequence']='tds_elstat'

###########################################################################
#RATE EQUATIONS/FLUXES
###########################################################################

species['CO']['flux']='catmap' #CO_rate
species['CO2']['flux']='catmap' #CO2_rate

boundary_thickness=7.93E-05 #in m

#if not nobuffer:
#    visc=viscosity(species['HCO3-']['bulk_concentration']/10**3), #Pa*s at 25C of KHCO3 solution
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

phimin=-0.5 #-0.5 #-0.7
phimax=-2.0

for potential in potentials:
    if nphi is not None:
        descriptors={'phiM':list(np.linspace(phimin,phimax,nphi))} #-0.0592*6.8,-2.0,nphi))}
    elif dphi is not None:
        descriptors={'phiM':list(np.linspace(phimin,phimax,-(phimax-phimin)/dphi+1))}
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

    
    if transport_mode != 'comsol':
        cm=CatMAP(transport=tp,model_name='CO2R')
        if transport_mode=='extrapolate':
            #for extrapolation of transport to high overpotential regime
            if transport_mode == 'extrapolate':
                extra=extrapolate(tp=tp,extrapol_folder=extrapol_folder)
                #extra.plot()
        for pot in descriptors['phiM']:
            print '!!! now running pot = '+str(pot)
            tp.system['phiM']=pot
            if transport_mode=='extrapolate':
                #set the surface concentrations according to extrapolated functions
                for sp in tp.species:
                    tp.species[sp]['surface_concentration']=10**extra.extrapol_func[sp](pot)
                #set voltage drop (= phi-phi0) according to extrapolated function
                tp.system['potential']=[extra.extrapol_func['voltage_diff_drop'](pot)]
                tp.system['surface_pH']=extra.extrapol_func['surface_pH'](pot) #lambda x: extra.extrapol_func['OH-'](pot)-3.+14.
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

