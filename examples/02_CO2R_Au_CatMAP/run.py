import numpy as np
from catint.transport import Transport
from catint.calculator import Calculator
from catint.units import unit_NA

pH = 6.8
system=\
    {
    #ENVIRONMENTAL CONDITIONS
    'temperature':  298,         #K
    'pressure':     1.013,       #bar
    'bulk_pH': pH,
    #MASS TRANSPORT
    'boundary thickness': 8.E-05, #m
    #ELECTROSTATICS
    'epsilon': 78.36,            #eps_0
    'migration': True,
    #REACTIONS
    'electrode reactions': True,
    'electrolyte reactions': True,
    #CHARGING
    'charging_scheme':'comsol',
    'phiM':-0.5,
    'phiPZC': 0.16,              #V vs. SHE
    'Stern capacitance': 20.,    #micro F/cm2
    #KINETICS
    'potential drop':'Stern',
    'active site density': 9.61e-05/unit_NA*(1e10)**2, #mol sites/m^2
    #INITIALIZATION
   # 'init_folder':init_folder,
    }
electrolyte_reactions=['bicarbonate-base','water-diss',{'additional_cell_reactions':'bicarbonate-acid'}]
electrode_reactions={
    'CO': {'reaction': 'CO2 + H2O + 2 e- -> CO + 2 OH-'},
}
species=\
    {
    'K+':               {'bulk_concentration':   'charge_neutrality',
                         #'MPB_radius':           2*3.5e-10},
                         'MPB_radius':           2*4.1e-10},
    'CO2':              {'bulk_concentration':   'Henry'},
    'OH-':              {'bulk_concentration':   10**(pH-14.)*1000.0}, #mol/m^3
    'CO':               {'bulk_concentration':   0.0}
    }
phimin=-0.5
phimax=-2.0
dphi=0.01
descriptors={'phiM':list(np.linspace(phimin,phimax,int(-(phimax-phimin)/dphi+1.)))}
comsol_args={}
#parameter
comsol_args['parameter']={}   
comsol_args['parameter']['grid_factor']=[str(100),'Grid factor']
comsol_args['parameter']['grid_factor_domain']=[str(100),'Grid factor']
comsol_args['parameter']['grid_factor_bound']=[str(200),'Grid factor']
#solver_settings
comsol_args['solver_settings']={}
comsol_args['solver_settings']['direct']={}
comsol_args['solver_settings']['direct']['nliniterrefine']=True
comsol_args['solver_settings']['ramp']={}
comsol_args['solver_settings']['ramp']['names']=['PZC','CS']
comsol_args['solver_settings']['ramp']['dramp']=0.01
#par_method
comsol_args['par_method']='internal'

#SOLVER SEQUENCE
#comsol_args['solver_settings']['solver_sequence']='tds_elstat'
#OTHER PARAMETER
#comsol_args['parameter']['RF']=[1,'Roughness Factor']
catmap_args={}
#CATMAP DESCRIPTOR RAMPING
catmap_args['desc_method']='automatic'
#catmap_args['min_desc']=0.0
catmap_args['min_desc_delta']=0.2
catmap_args['max_desc_delta']=0.2
#INTERACTIONS
catmap_args['n_inter']='automatic'
species['CO']['flux']='catmap' #CO production rate
species['CO2']['flux']='catmap' #CO2 consumption rate
nx=200
tp=Transport(
    species=species,
    electrode_reactions=electrode_reactions,
    electrolyte_reactions=electrolyte_reactions,
    system=system,
    catmap_args=catmap_args,
    comsol_args=comsol_args,
    model_name='CO2R',
    descriptors=descriptors,
    nx=nx)
tp.set_calculator('comsol')
c=Calculator(transport=tp,tau_scf=0.008,ntout=1,dt=1e-1,tmax=10,mix_scf=0.02)
c.run()
