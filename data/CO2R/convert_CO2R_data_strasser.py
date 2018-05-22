import numpy as np
from glob import glob
import sys

#fnames=glob('hori_FE_CO2R_*.csv')
fnames=glob('Strasser_KHCO3_02_*NHE.csv')
jdata=np.loadtxt('Strasser_KHCO3_02_j_NHE.csv',delimiter=',')

alldata=[]
header=''
header+='#U vs SHE (V),'
alldata.append([str(x) for x,y in jdata])
for fname in fnames:
    if '_j_' in fname:
        continue
    header+=fname.split('.')[0].split('_')[-2]+','
    print 'Converting ',fname
    data=np.loadtxt(fname,delimiter=',')
    pdata=[]
    for xj,dj in jdata:
        #print '>> ',xj
        closest=np.infty
        i=-1
        index=None
        for xd,dd in data:
            i+=1
            dist=abs(xd-xj)
            if dist<closest and dist<2e-2:
                closest=dist
                index=i
        print 'index',xj, index
        if index is None:
            pdata.append('')
        else:
            print xj, dj, data[index,1],-dj*data[index,1]/100.
            pdata.append(str(-dj*data[index,1]/100.))
            #    print '>> >> appending ', xd
    print 'adding pdata=',pdata
    alldata.append(pdata)
alldata=map(list, zip(*alldata))
#    if len(pdata)!=len(data):
#        print 'Something wrong here , ', len(pdata), len(data)
#    np.savetxt(fname.split('.')[0].split('_')[-1],alldata,delimiter=',')
header=header[:-1]
print header
print len(alldata[0])
np.savetxt('strasser_jpart_KHCO3_test_CO2R_SHE.csv',alldata,delimiter=',',fmt="%s",header=header)
