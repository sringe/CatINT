from catmap.model import ReactionModel
import sys
from catmap import analyze
import os

os.system('rm *.log *.pkl')

mkm_file=sys.argv[1]

model = ReactionModel(setup_file = mkm_file, max_log_line_length=0)
model.output_variables+=['consumption_rate','production_rate', 'free_energy', 'selectivity', 'interacting_energy','turnover_frequency']
model.run()
possible_prdt = ('H2_g',) 
#possible_prdt = ('H2_g','CO_g','CH4_g') #,'CH3CH2OH_g',)
model.output_labels['my_selectivity'] = possible_prdt
prdt_idx = {}
print 'searching',[prdt for prdt in possible_prdt]
print 'output labels',model.output_labels['production_rate']
for prdt in possible_prdt:
    prdt_idx[prdt] = model.output_labels['production_rate'].index(prdt)
my_selectivity = []
for descri, rate in model.production_rate_map:
    total = 1.e-99
    for prdt in possible_prdt:
        total += rate[prdt_idx[prdt]]
    selectivity = [rate[prdt_idx[prdt]]/total for prdt in possible_prdt]
    my_selectivity.append([descri, selectivity])
model.my_selectivity_map = my_selectivity

from catmap import analyze
vm = analyze.VectorMap(model)

vm.plot_variable = 'rate'
vm.log_scale = True
vm.min = 1e-10
vm.max = 1e6
fig = vm.plot(save=False)
fig.savefig('rate.pdf')

vm.plot_variable = 'production_rate'
vm.log_scale = True
vm.min = 1e-10
vm.max = 1e6
fig = vm.plot(save=False)
fig.savefig('production_rate.pdf')

idx = [i for i in range(len(model.interacting_energy_map)) if model.interacting_energy_map[i][0][0] == 0.0][0]
descrip, coverages = model.coverage_map[idx]
rxn_parameters = model.scaler.get_rxn_parameters(descrip)
#rate_constants = model.solver.get_rate_constants(rxn_parameters,coverages)
#kfs, krs, dkfs, dkrs = model.rate_constants(rxn_parameters,coverages,
#    model._gas_energies,model._site_energies,
#    model.temperature,model.interaction_response_function,
#    model._mpfloat,model._matrix,model._math.exp)
model.solver.get_interacting_energies(rxn_parameters)
all_ads = model.adsorbate_names + model.transition_state_names
N_ads = len(all_ads)
energies = rxn_parameters[:N_ads]
eps_vector = rxn_parameters[N_ads:]
cvg = coverages + [0]*len(model.transition_state_names)
print 'the energies'
print model.interaction_function(cvg,energies,eps_vector,model.thermodynamics.adsorbate_interactions.interaction_response_function,False,False)

def plot_fed(corr):
    ma = analyze.MechanismAnalysis(model)
    ma.surface_colors = ['k','b','r','yellow','green','orange','cyan']
    ma.label_args['size'] = 14
    ma.energy_type = 'free_energy' #can also be free_energy/potential_energy
    ma.include_labels = True #way too messy with labels
    ma.pressure_correction = corr #assume all pressures are 1 bar (so that energies are the same as from DFT)
    ma.coverage_correction = False
    ma.energy_type = 'interacting_energy'
    ma.include_labels = True
    if not corr:
        fig = ma.plot(save='FED.pdf',plot_variants=[0.0],method=2) #desc_val[0]])
    else:
        fig = ma.plot(save='FED_pressure_corrected.pdf',plot_variants=[0.0],method=2) #desc_val[0]])
#plot_fed(True)
plot_fed(False)
