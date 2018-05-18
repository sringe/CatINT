import sys
sys.path.insert(0,'/scratch/users/sringe/transport/catmap')
from catmap.model import ReactionModel
#import sys
from string import Template
import os
import numpy as np
from glob import glob

if len(sys.argv)>2:
    if sys.argv[-1]=='--cleanup':
        files=glob('*.png')+glob('*.pdf')+glob('*.pkl')+glob('*.log')
        for f in files:
            os.remove(f)

include_rate_control = False

pH = ['0','4','7','10','14']
i=0
j = pH[i]
print 'pH = '+j    
mkm_template = Template(open('catmap_'+sys.argv[1]+'_template.mkm').read())
mkm_text = mkm_template.substitute(pH_new = j)
mkm_file = 'catmap_'+sys.argv[1]+'.mkm'
with open(mkm_file,'w') as f:
    f.write(mkm_text)
print 'Reading',mkm_file
model = ReactionModel(setup_file = mkm_file)
model.output_variables+=['production_rate', 'free_energy', 'selectivity', 'interacting_energy']
model.run()

from catmap import analyze
vm = analyze.VectorMap(model)
#vm.plot_variable = 'rate'
#vm.descriptor_labels = ['U vs. SHE (V)']
#vm.log_scale = False
#vm.min = -100
#vm.max = 100
#fig = vm.plot(save=False)
#fig.savefig('rate'+j+'.png')
vm.descriptor_labels = ['U vs. SHE (V)']
vm.plot_variable = 'rate'
vm.log_scale = True
vm.min = 1e-10
vm.max = 1e6
fig = vm.plot(save=False)
fig.savefig('rate.pdf')

#vm.plot_variable = 'production_rate'
#vm.descriptor_labels = ['U vs. SHE (V)']
#vm.log_scale = True
#vm.min = -10 #1e-10
#vm.max = 6 #1e6
#fig = vm.plot(save=False)
#fig.savefig('production_rate'+j+'.png')
vm.descriptor_labels = ['U vs. SHE (V)']
vm.plot_variable = 'production_rate'
vm.log_scale = True
vm.min = 1e-10
vm.max = 1e6
fig = vm.plot(save=False)
fig.savefig('production_rate.pdf')

vm = analyze.VectorMap(model)
vm.log_scale = False
vm.plot_variable = 'coverage'
vm.descriptor_labels = ['coverage (ML)']
vm.min = 0
vm.max = 1
fig = vm.plot(save=False)
fig.savefig('coverage'+j+'.png')

ma = analyze.MechanismAnalysis(model)
ma.energy_type = 'interacting_energy'
label_size = 10
ma.kwarg_dict = {
        'VTA':{'label_positions':None,'label_args':{'color':'black','rotation':45,'ha':'center','size':label_size}},
        'VHA':{'label_positions':None,'label_args':{'color':'blue','rotation':45,'ha':'center','size':label_size}},
        'VTB':{'label_positions':None,'label_args':{'color':'darkviolet','rotation':45,'ha':'center','size':label_size}},
        'VHB':{'label_positions':None,'label_args':{'color':'firebrick','rotation':45,'ha':'center','size':label_size}},
        }

ma.subplots_adjust_kwargs = {'top': 0.87, 'bottom':0.22}
ma.include_labels = True
ma.label_args['size'] = 12
ma.pressure_correction = False
ma.coverage_correction = False
fig = ma.plot(save=False, plot_variants = [0.0])
ax = fig.add_subplot(111)
ax.set_ylim([-0.5,2.5])
fig.savefig('FED'+j+'.png')


