from catint.plot import Plot
from catint.transport import Transport


tp=Transport()
p=Plot(transport=tp,init_from_file='calc_std_settings')
p.plot(large_plots=['concentrations_reaction','desc_current_density'],\
        small_plots=['potential','concentrations_electrolyte','current_density','pH'])
