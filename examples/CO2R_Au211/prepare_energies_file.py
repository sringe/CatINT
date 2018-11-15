import sys
import numpy as np
from my_io import replace_at_string

surface='Au-fcc111'

zero_field=False

sigma_params={\
    'Au-fcc111':{
#fitting parameters, in brackets the number of converged points
#    'CO2_t':[-0.00089829, -0.01080444,  0.27591133],               #wo-extrapol-eps6
        'CO2_t':[-1.89651545e-04,  2.11273194e-02,  6.19034594e-01+0.4],    #w-extrapol-eps6 (5)
#    'CO2_t':[-6.40273008e-05  2.66094980e-02  7.45267835e-01]      #w-extrapol-eps10 (4)
#    'CO2_t':[-6.04210199e-06,  2.47639012e-02,  7.78606901e-01],   #w-extrapol-eps80 (2)
        'CO_t':[-0.00102957, -0.03303381,  0.4962946 ],                  #eps6-bridge (3) - R
#    'CO_t':[ 3.64651593e-05 -1.11837285e-02  7.15122431e-01],      #eps6 (7)
        'COOH_t':[-3.61959351e-04, -2.96221404e-03,  6.89265032e-01-0.2]     #eps6 (3) - R
    },
    'Au-fcc211':{
        'COOH_t':[-2.05088705e-04, -1.86659496e-03,  5.28198285e-01],  #eps6 (3) - R
#        'CO_t':[-3.13816399e-04, -1.26773916e-02,  4.31694279e-01],    #eps6 (10)
        'CO_t':[-0.00115802, -0.01070715,  0.44634201],                  #eps6-bridge (3) - R
    }}

sigma_params=sigma_params[surface]

#assume same sigma dependence of transition state as COOH
sigma_params['COOH-H2O-ele_t']=[]
for val in sigma_params['COOH_t']:
    sigma_params['COOH-H2O-ele_t'].append(val)
#COOH to CO barrier
sigma_params['COOH-H2O-ele_t'][-1]=0.0 #0.4
if zero_field:
    for s in sigma_params:
        sigma_params[s][0]=0.0
        sigma_params[s][1]=0.0
vibrations={
        'gas':{
            'CO2_g':[24.1,70.7,635.8,640.5,1312.2,2361.2],
            'CO_g':[89.8,127.2,2145.5],
            'H2O_g':[103.0,180.6,245.1,1625.9,3722.8,3830.3],
            'H2_g':[3.8,13.5,4444.5]},
        'Au-fcc111':{
            #CO2 at -0.9 e charge
            'CO2_t':[(114.3+159.4)/2.,(180.2+187.0)/2.,(208.6+217.3)/2.,(239.4+262.0)/2.,\
                (301.0+311.0)/2.,(503.5+517.6)/2.,(558.2+566.3)/2.,(1173.8+1178.3)/2.,\
                (1872.9+1906.8)/2.],
            'COOH_t':[(89.1+115.4)/2.,(145.3+198.8)/2.,(237.0+248.3)/2.,(258.5+272.6)/2.,\
                (295.4+311.5)/2.,(517.1+524.5)/2.,(638.6+655.0)/2.,(792.6+795.7)/2.,\
                (1028.1+1037.3)/2.,(1340.9+1347.8)/2.,(1654.9+1661.4)/2.,(3545.0+3557.7)/2.],
            'CO_t':[(122.4+136.9)/2.,(144.8+166.3)/2.,(183.5+196.1)/2.,(213.5+241.5)/2.,\
                    (2071.9+2074.6)/2.]}
            }


vibrations_new=vibrations[surface]
vibrations_new.update(vibrations['gas'])
vibrations=vibrations_new

free_en_corrs={
        'CO2_g':0.33}
#raw energies relative to CO2
raw_energies={
        'H2_g':0.0,
        'CO2_g':0.0,
        'H2O_g':0.0,
        'CO_g':0.057936319999996044-free_en_corrs['CO2_g'],
        'COOH_t':sigma_params['COOH_t'][-1]-free_en_corrs['CO2_g'],
        'CO_t':sigma_params['CO_t'][-1]-free_en_corrs['CO2_g'],
        'CO2_t':sigma_params['CO2_t'][-1]-free_en_corrs['CO2_g'],
        'COOH-H2O-ele_t':sigma_params['COOH-H2O-ele_t'][-1],
        }
for key in list(set(raw_energies.keys()+vibrations.keys())):
    rkey=key.split('_')[0]

    if '_t' in key:
        #first work on mkm file
        search_str=\
            'species_definitions[\''+key+'\'][\'sigma_params\']'
        print 'search1',search_str
        replace_str=\
            ['species_definitions[\''+key+'\'][\'sigma_params\']=[{},{},{}]'.format(sigma_params[key][0],sigma_params[key][1],sigma_params[key][2])]
        replace_at_string(search_str,'catmap_CO2R_template.mkm',replace_str)

    #2nd work on energies file

    if '_t' in key:
        site='Au\t211'
    elif '_g' in key:
        site='None\tgas'

    if key in vibrations:
        vibs=sorted(vibrations[key]) #better safe than sorry, sort the vibrations here, the largest should be in the end by default in ase
    else:
        vibs=[]

    search_str=\
        site+'\t'+rkey+'\t'
    print 'search',search_str
    print 'test',key
    print 'test2',raw_energies.keys()
    replace_str=\
        [site+'\t'+rkey+'\t'+str(raw_energies[key])+'\tfcc\t'+str(vibs)+'\t[]\tmy ads calcs']
    replace_at_string(search_str,'catmap_CO2R_energies.txt',replace_str)
