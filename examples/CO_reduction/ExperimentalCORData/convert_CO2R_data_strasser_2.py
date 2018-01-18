import numpy as np
from glob import glob
import sys
from units import *

#converts the product rates

#fnames=glob('hori_FE_CO2R_*.csv')
fnames=glob('Strasser_prate_KHCO3_01_*NHE.csv')

#molecular weight [g/mol] / gas density at 25 C [kg/m^3] -> L/mol
Vm={'H2':2.016/8.23E-2,
    'CH4':16.043/6.567E-1,
    'C2H4':28.05/1.153,   #https://encyclopedia.airliquide.com/ethylene
    'CO':28.01/1.145}

nel={'H2':2,
    'CH4':8,
    'C2H4':12,
    'CO':2}

nprod={'H2':1,
       'CH4':1,
       'C2H4':1,
       'CO':1}

alldata=[]
header=''
for fname in fnames:
    header+='#U vs SHE (V),'
    header+=fname.split('.')[0].split('_')[-2]+','
    print 'Converting ',fname
    data=np.loadtxt(fname,delimiter=',')
    x=data[:,0]
    y=data[:,1]
    sp=fname.split('.')[0].split('_')[-2]
    print 'sp:',sp, nel[sp],Vm[sp],nprod[sp]
    y=y/60./60./Vm[sp]*unit_F/10.*nel[sp]/nprod[sp]
    x = list(x) + [0]*(15 - len(x))
    y = list(y) + [0]*(15 - len(y))
    x=np.array(x)
    y=np.array(y)
    alldata.append(x)
    alldata.append(y)
alldata=map(list, zip(*alldata))
print 'alldata'
print alldata
#    if len(pdata)!=len(data):
#        print 'Something wrong here , ', len(pdata), len(data)
#    np.savetxt(fname.split('.')[0].split('_')[-1],alldata,delimiter=',')
header=header[:-1]
print header
print len(alldata[0])
np.savetxt('strasser_jpart_from_prate_KHCO3_test_CO2R_SHE.csv',alldata,delimiter=',',fmt="%s",header=header)
