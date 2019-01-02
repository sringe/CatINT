import sys
import numpy as np
from my_io import replace_at_string
from units import Rydberg

#free_en_corrs={'CO2_g':0.33,'CO_g':0.0,'H2_g':0.0,'H2O_g':0.0}
free_en_corrs={'CO2_g':0.33,'CO_g':0.0,'H2_g':0.0,'H2O_g':0.0} #09}

abinitio_energies = {
         'CO_gas': -46.04587566*Rydberg+free_en_corrs['CO_g'],
         'H2_gas': -2.42121318*Rydberg+free_en_corrs['H2_g'],
         'CO2_gas': -80.15782693*Rydberg+free_en_corrs['CO2_g'],
         'H2O_gas': -36.47522813*Rydberg+free_en_corrs['H2O_g'],
         }

abinitio_energies['H']=abinitio_energies['H2_gas']/2.
abinitio_energies['O']=abinitio_energies['H2O_gas']-abinitio_energies['H2_gas']
abinitio_energies['C']=abinitio_energies['CO2_gas']-abinitio_energies['O']*2.

abinitio_energies['CO_gas']-=(abinitio_energies['C']+abinitio_energies['O'])
abinitio_energies['H2_gas']-=(2.*abinitio_energies['H'])
abinitio_energies['CO2_gas']-=(abinitio_energies['C']+abinitio_energies['O']*2.)
abinitio_energies['H2O_gas']-=(abinitio_energies['H']*2.+abinitio_energies['O'])

for a in abinitio_energies:
    abinitio_energies[a]=round(abinitio_energies[a],10)

surface='Au-fcc211'
ts='COOH-H2O-ele'

zero_field=False

alpha_factor={'CO':0.4,'CO2':0.7,'COOH':1.}
dipole_factor={'COOH':-20.,'CO':0.3,'CO2':0.3}
shift={'COOH':0.1,'CO':0.07,'CO2':0.2}

sigma_params={\
    'Au-fcc111':{
#fitting parameters, in brackets the number of converged points
#    'CO2_t':[-0.00089829, -0.01080444,  0.27591133],               #wo-extrapol-eps6
#    'CO2_t':[-6.40273008e-05  2.66094980e-02  7.45267835e-01]      #w-extrapol-eps10 (4)
#    'CO2_t':[-6.04210199e-06,  2.47639012e-02,  7.78606901e-01],   #w-extrapol-eps80 (2)
#    'CO_t':[ 3.64651593e-05 -1.11837285e-02  7.15122431e-01],      #eps6 (7)

#take these:
        #'COOH_t':[-3.61959351e-04, -2.96221404e-03,  6.89265032e-01],     #eps6 (3) - R, 111
        'COOH_t':[-4.61907388e-04,-2.84367353e-04,5.61124442e-01],  #eps6-highpw-211
#        'COOH_t':[-4.76480730e-04 -6.32919405e-03  7.35459652e-01], #eps6-highpw
        'CO2_t':[-1.89651545e-04,  2.11273194e-02,  6.19034594e-01],    #w-extrapol-eps6 (5)
 #       'CO2_t':[-2.16059516e-04,  1.60712400e-02,  6.07344189e-01],    #w-extrapol-eps6-2x2
#        'CO2_t': 
        'CO_t':[-0.00102957, -0.03303381,  0.4962946 ],                  #eps6-bridge (3) - R
    },
    'Au-fcc211':{
#        'COOH_t':[-4.53877629e-04, -8.81491164e-04,  5.75756161e-01],  #eps6-highpw
        'COOH_t':[-4.53877629e-04*alpha_factor['COOH'], -8.81491164e-04*dipole_factor['COOH'],  5.75756161e-01+shift['COOH']],
#        'CO_t':[-3.13816399e-04, -1.26773916e-02,  4.31694279e-01],    #eps6 (10)
        'CO2_t':[-6.58975689e-04*alpha_factor['CO2'],  2.85121677e-02*dipole_factor['CO2'],  9.02204236e-01+shift['CO2']],
       # 'CO_t':[-0.00115802, -0.01070715,  0.44634201],                  #eps6-bridge (3) - R
       'CO_t':[-0.00115802*alpha_factor['CO'], -0.01070715*dipole_factor['CO'],  0.44634201+shift['CO']]
    }}

sigma_params=sigma_params[surface]

#assume same sigma dependence of transition state as COOH
sigma_params[ts+'_t']=[]
for val in sigma_params['COOH_t']:
    sigma_params[ts+'_t'].append(val)
#COOH to CO barrier
sigma_params[ts+'_t'][-1]=0.0 #0.4
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

vibrations['Au-fcc211']=vibrations['Au-fcc111']

vibrations_new=vibrations[surface]
vibrations_new.update(vibrations['gas'])
vibrations=vibrations_new

#raw energies relative to CO2
raw_energies={
        'H2_g':abinitio_energies['H2_gas'],
        'CO2_g':abinitio_energies['CO2_gas'],
        'H2O_g':abinitio_energies['H2O_gas'],
        'CO_g':abinitio_energies['CO_gas'],
        'COOH_t':sigma_params['COOH_t'][-1]-free_en_corrs['CO2_g'],
        'CO_t':sigma_params['CO_t'][-1]-free_en_corrs['CO2_g'],
        'CO2_t':sigma_params['CO2_t'][-1]-free_en_corrs['CO2_g'],
        ts+'_t':sigma_params[ts+'_t'][-1],
        }
print raw_energies
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
