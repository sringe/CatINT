from catmap.model import ReactionModel
import sys
from catmap import analyze

mkm_file=sys.argv[1]

model = ReactionModel(setup_file = mkm_file, max_log_line_length=0)
model.output_variables+=['consumption_rate','production_rate', 'free_energy', 'selectivity', 'interacting_energy','turnover_frequency']
model.run()
def plot_fed(corr):
    ma = analyze.MechanismAnalysis(model)
    ma.surface_colors = ['k','b','r','yellow','green','orange','cyan']
    ma.label_args['size'] = 14
    ma.energy_type = 'free_energy' #can also be free_energy/potential_energy
    ma.include_labels = True #way too messy with labels
    ma.pressure_correction = corr #assume all pressures are 1 bar (so that energies are the same as from DFT)
    ma.coverage_correction = False
    ma.include_labels = True
    if not corr:
        fig = ma.plot(save='FED.pdf',plot_variants=[0.0]) #desc_val[0]])
    else:
        fig = ma.plot(save='FED_pressure_corrected.pdf',plot_variants=[0.0]) #desc_val[0]])
#plot_fed(True)
plot_fed(False)
