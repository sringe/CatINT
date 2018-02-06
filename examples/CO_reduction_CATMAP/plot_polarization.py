import numpy as np
import matplotlib.pyplot as plt
import sys
from glob import glob

colorlist={}
colorlist['CO']='r'
colorlist['H2']='b'
colorlist['CO2']='k'
colorlist['H']='olive'
symbols=['x','o']

def read_data(files):
    pdata={}
    iarg=-1
    for arg2 in files:
        iarg+=1
        data=np.loadtxt(arg2)
        x=data[:,0]
        y=data[:,1]
        species=arg2.split('/')[-1].split('_')[-1].split('.')[0]
        if k==0:
            symbol='o'
        else:
            symbol='-'
        if species not in pdata:
            pdata[species]=[]
        pdata[species].append([x[0],y[0]])
    for sp in pdata:
        pdata[sp].sort(key=lambda x: x[0])
        pdata[sp]=np.array(pdata[sp])
    return pdata

k=-1
for arg in sys.argv[1:]:
    k+=1
    print 'Reading ',arg
    results_folder=arg+'/catmap_output'
    
    ax1=plt.subplot('211')
    ax2=plt.subplot('212')
    
    ax1.set_title('Polarization')
    pdata=read_data(glob(results_folder+'/*/j*'))
    cdata=read_data(glob(results_folder+'/*/cov*'))
    for isp,sp in enumerate(pdata):
        x=pdata[sp][:,0]
        y=pdata[sp][:,1]
        ax1.semilogy(x,y,'-o',color=colorlist[sp],label=sp) #,label=arg.split('/')[-1])
    for isp,sp in enumerate(cdata):
        x=cdata[sp][:,0]
        y=cdata[sp][:,1]
        ax2.plot(x,y,'-o',color=colorlist[sp],label=sp)
    ax1.legend()
#    data=np.loadtxt('data.txt')
#    ax1.plot(data[:,0],data[:,1],'-')
    ax2.set_title('Coverages')
#    iarg=-1
#    for arg2 in glob(results_folder+'/*/cov*'):
#        data=np.loadtxt(arg2)
#        x=data[:,0]
#        y=data[:,1]
#        species=arg2.split('/')[-1].split('_')[-1].split('.')[0]
#        if k==0:
#            ax2.semilogy(x,y,'o',color=colorlist[species],label=arg2.split('/')[-1])
#        else:
#            ax2.semilogy(x,y,'-',color=colorlist[species])
    ax2.legend()
plt.show()
