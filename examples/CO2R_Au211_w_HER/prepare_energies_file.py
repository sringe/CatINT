import sys
import numpy as np
from my_io import replace_at_string
from units import *


DeltaG=-0.27319*eVToJ #J (solvation energy)
WaterDensity=55.5 #mol/L
rho_l=WaterDensity*10**3*unit_NA #liquid water particle density in [1/m^3]
PaToatm=9.86923e-6
kT=unit_kB*unit_T
#E = P*V
#E*P
Gen=np.log(kT*np.exp((DeltaG-kT*np.log(1/rho_l))/kT)*9.86923e-6)*0.0257
Pvap=kT*np.exp((DeltaG-kT*np.log(1/rho_l))/(kT))*PaToatm
#Pvap=kT [1/m^3]*np.exp((DeltaG-kT*np.log([1/m^3]/rho_l))/(kT))
#np.log(Pvap/kT)/kT/[1/m^3])*kT
Gen=0.0257*np.log(3.1690*1000./1e5)
free_en_corrs={'CO2_g':0.33,'CO_g':0.0,'H2_g':0.09,'H2O_g':Gen,\
        'CO2_t':0.25*1.5,'COOH_t':0.25,'CO_t':-0.1}
#free_en_corrs['CO_g']=free_en_corrs['CO2_g']-0.33
#free_en_corrs={'CO2_g':0.45,'CO_g':0.0,'H2_g':0.0,'H2O_g':0.0} #09}
#free_en_corrs={'CO2_g':0.41,'CO_g':-0.18,'H2_g':0.09,'H2O_g':-0.21} #09}

#beef, highpw
abinitio_energies = {
         'CO_gas': -46.05684028*Rydberg+free_en_corrs['CO_g'],
         'H2_gas': -2.42319084*Rydberg+free_en_corrs['H2_g'],
         'CO2_gas': -80.17847584*Rydberg+free_en_corrs['CO2_g'],
         'H2O_gas': -36.48757865*Rydberg+free_en_corrs['H2O_g'],
         }

#rpbe-d3, highpw
#abinitio_energies={
#        'CO_gas':-45.61372421*Rydberg+free_en_corrs['CO_g'],
#        'H2_gas':-2.35720879*Rydberg+free_en_corrs['H2_g'],
#        'CO2_gas':-79.42951519*Rydberg+free_en_corrs['CO2_g'],
#        'H2O_gas':-36.11827578*Rydberg+free_en_corrs['H2O_g']}

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

#new implicit data from (211) facet
#sigma_params['CO2_t']=[-3.69386461e-04,  3.34776775e-02,  8.88053574e-01] #[-4.40195996e-04,  2.55686722e-02,  8.21548063e-01] #(111)-BEEF
#sigma_params['COOH_t']=[-3.92928033e-04,  1.51464537e-03,  5.70060077e-01] #[-4.61907388e-04, -2.84367353e-04,  5.62678756e-01] #211
#sigma_params['CO_t']=[-3.67141924e-04, -1.22885293e-02,  4.37721603e-01] #[-9.11793749e-05,  4.30663723e-03,  4.91072027e-01] #[0.0018451327711406771, 0.03420296578328402, -0.45936375172574206] #0.00218647, 0.04507817, -0.3058627089720931]

#new explicit/implicit data from (100) facet. using implicit slopes and explicit data points
#sigma_params['CO2_t']=[-0.00188943, -0.03172613, 1.2113668269074331] #implicit-explicit
#sigma_params['CO2_t']=[-0.00188943, -0.03172613,  0.41581459] #implicit
#sigma_params['CO2_t']=[0.05118888,1.9795144707714756] #explicit (linear fit)
#sigma_params['CO2_t']=[-0.0004612829823,0.00855365031,1.102204236] #old (wrong) params that gave the correct slope

#sigma_params['CO2_t']=[-2.86600929e-04,  2.97720125e-02,  7.37600203e-01] #BEEF-vdW, surfpar, 211
sigma_params['CO2_t']=[-2.17559621e-04,  2.30636099e-02,  5.65804645e-01] #BEEF-vdW, surfpar, 100

#sigma_params['COOH_t']=-0.000115225255, 0.00267094987, 0.1635177832156014 #implicit-explicit
#sigma_params['COOH_t']=[0.000115225255, 0.00267094987, 0.7347896980000002] #implicit

#sigma_params['COOH_t']=[-9.02956820e-05,  2.26896383e-03,  2.53214079e-01] #BEEF-vdW, surfpar, 211
sigma_params['COOH_t']=[-2.30004216e-04,-1.03411773e-03,  3.44172799e-01] #BEEF-vdW, surfpar, 100

#sigma_params['CO_t']=[-9.11793749e-05, 0.00430663723, -0.19944434712717363] #implicit-explicit
#sigma_params['CO_t']=[0.002525766808954512, 0.07345291335810258, -0.08036375172574228] #explicit (polynomial fit)
#sigma_params['CO_t']=[-9.11793749e-05, 0.00430663723, 0.49107202699999997] #implicit

#sigma_params['CO_t']=[-1.89106972e-04, -9.42574086e-03,  3.87255672e-01] #BEEF-vdW, surfpar, 211
sigma_params['CO_t']=[-1.02922269e-04,  1.97610024e-03,  5.04130607e-01] #BEEF-vdW, surfpar, 100

#assume same sigma dependence of transition state as COOH
sigma_params[ts+'_t']=[]
for val in sigma_params['COOH_t']:
    sigma_params[ts+'_t'].append(val)
#COOH to CO barrier at 0 V vs. SHE
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

print sigma_params['COOH_t'][-1],free_en_corrs['CO2_g'],free_en_corrs['H2_g']/2.
#raw energies relative to CO2
raw_energies={
        'H2_g':abinitio_energies['H2_gas'],
        'CO2_g':abinitio_energies['CO2_gas'],
        'H2O_g':abinitio_energies['H2O_gas'],
        'CO_g':abinitio_energies['CO_gas'],
        'COOH_t':sigma_params['COOH_t'][-1]+free_en_corrs['COOH_t']-free_en_corrs['CO2_g']-free_en_corrs['H2_g']/2.,
        'CO_t':sigma_params['CO_t'][-1]+free_en_corrs['CO_t']+free_en_corrs['H2O_g']-free_en_corrs['CO2_g']-free_en_corrs['H2_g'],
        'CO2_t':sigma_params['CO2_t'][-1]+free_en_corrs['CO2_t']-free_en_corrs['CO2_g'],
        ts+'_t':sigma_params[ts+'_t'][-1],
        }
print raw_energies
print sigma_params['CO2_t'][-1],free_en_corrs['CO2_g']
for key in list(set(raw_energies.keys()+vibrations.keys())):
    rkey=key.split('_')[0]

    if '_t' in key:
        #first work on mkm file
        search_str=\
            'species_definitions[\''+key+'\'][\'sigma_params\']'
        if len(sigma_params[key])==3:
            replace_str=\
                ['species_definitions[\''+key+'\'][\'sigma_params\']=[{},{},{}]'.format(sigma_params[key][0],sigma_params[key][1],sigma_params[key][2])]
        elif len(sigma_params[key])==2:
            replace_str=\
                ['species_definitions[\''+key+'\'][\'sigma_params\']=[{},{}]'.format(sigma_params[key][0],sigma_params[key][1])]
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
