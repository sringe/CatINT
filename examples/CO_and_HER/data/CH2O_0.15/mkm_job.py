from catmap import ReactionModel
import sys
from string import Template
import os


include_rate_control = False
FED = True
#for s in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
for s in [1.0]:
    mkm_template = Template(open('COR_template.mkm').read())
    mkm_text = mkm_template.substitute(strength=s)
    with open('COR.mkm','w') as f:
        f.write(mkm_text)
    model = ReactionModel(setup_file = 'COR.mkm')
    model.output_variables+=['production_rate', 'free_energy', 'selectivity', 'interacting_energy', 'interaction_matrix']
    if include_rate_control:
        model.output_variables += ['rate_control']
    model.run()

    possible_prdt = ('CH4_g', 'H2_g', 'CH3CH2OH_g', 'CH2O_g', 'CH4O2_g')
    model.output_labels['my_selectivity'] = possible_prdt
    prdt_idx = {}
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
    fig.savefig('rate'+str(s)+'.pdf')

    vm.plot_variable = 'production_rate'
    vm.log_scale = True
    vm.min = 1e-10
    vm.max = 1e6
    fig = vm.plot(save=False)
    fig.savefig('production_rate'+str(s)+'.pdf')

    vm = analyze.VectorMap(model)
    vm.log_scale = False
    vm.unique_only = False
    vm.plot_variable = 'coverage'
    vm.min = 0
    vm.max = 1
    fig = vm.plot(save=False)
    fig.savefig('coverage'+str(s)+'.pdf')

    vm = analyze.VectorMap(model)
    vm.log_scale = True
    vm.unique_only = False
    vm.plot_variable = 'coverage'
    vm.min = 1e-20
    vm.max = 1
    fig = vm.plot(save=False)
    fig.savefig('coverageLog'+str(s)+'.pdf')

    vm = analyze.VectorMap(model)
    vm.log_scale = False
    vm.plot_variable = 'my_selectivity'
    vm.include_labels = ['CH4_g', 'H2_g', 'CH3CH2OH_g', 'CH2O_g', 'CH4O2_g']
    vm.min = 0
    vm.max = 1
    fig = vm.plot(save=False)
    fig.savefig('my_selectivity'+str(s)+'.pdf')

    vm = analyze.VectorMap(model)
    vm.plot_variable = 'my_selectivity'
    vm.include_labels = ['CH4_g', 'H2_g', 'CH3CH2OH_g', 'CH2O_g', 'CH4O2_g']
    vm.min = 1e-10
    vm.max = 1
    vm.log_scale = True
    fig = vm.plot(save=False)
    fig.savefig('my_selectivityLog'+str(s)+'.pdf')

    if include_rate_control:
        mm = analyze.MatrixMap(model)
        mm.plot_variable = 'rate_control'
        mm.log_scale = False
        mm.min = -2
        mm.max = 2
        mm.plot(save='rate_control.pdf')
    
    os.system('cp C2_pH0.pkl C2_pH0_'+str(s)+'.pkl')

    ma = analyze.MechanismAnalysis(model)
    ma.energy_type = 'interacting_energy'
    label_size = 10
    ma.kwarg_dict = {
            'C2_via_OCCHO-p':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
            'C2_via_OCCOH-p':{'label_positions':None,'label_args':{'color':'b','rotation':45,'ha':'center','size':label_size}},
            'C2_via_OCCO-SH':{'label_positions':None,'label_args':{'color':'g','rotation':45,'ha':'center','size':label_size}},
            'C2_via_OCCOH-p':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
            'CH4_via_CHO-ele':{'label_positions':None,'label_args':{'color':'b','rotation':45,'ha':'center','size':label_size}},
            'CH4_via_COH-ele':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
            'HER_Heyrovsky':{'label_positions':None,'label_args':{'color':'g','rotation':45,'ha':'center','size':label_size}},
            }

    if FED:
        ma.subplots_adjust_kwargs = {'top': 0.87, 'bottom':0.22}
        ma.surface_colors = ['k','b','r','yellow','green','orange','cyan']
        ma.include_labels = True
        ma.label_args['size'] = 14
        ma.pressure_correction = False
        ma.coverage_correction = False
        fig = ma.plot(save=False, plot_variants = [0.0, -0.5, -0.7])

        ax = fig.add_subplot(111)
        delta = 0.05
        y0 = 0.92
        x0 = 0.85
        x1 = 0.9
        y1 = 0.9
        for surf,col in zip(ma.surface_names, ma.surface_colors):
            ax.annotate('____',[x0,y0],xycoords='axes fraction',color=col)
            ax.annotate(surf,[x1,y1],xycoords='axes fraction',color=col)
            y0-= delta
            y1-= delta
        fig.savefig('FED_int_'+str(s)+'.pdf')
