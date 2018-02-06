import numpy as np
import matplotlib.pyplot as plt
import sys
from glob import glob
from catint.experimental import EXPDATA

exp=EXPDATA()

j_log_plot=True

fig=plt.figure(figsize=(7, 10))

colorlist={}
colorlist['CO']='r'
colorlist['H2']='b'
colorlist['CO2']='k'
colorlist['CH4']='olive'
colorlist['CH3CH2OH']='orange'
colorlist['H']='lightblue'
colorlist['CHOH']='y'
colorlist['OCCOH']='darkred'
colorlist['CHO']='0.75'
colorlist['OCCO']='g'
colorlist['COOH']='darkblue'
symbols=['x','o']

products=['CO','CH4','CH3CH2OH','H2']

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
    results_folder=arg+'/catmap_output'
    
    ax1=plt.subplot('211')
    ax2=plt.subplot('212')
    
    ax1.set_title('Polarization')
    pdata=read_data(glob(results_folder+'/*/j*'))
    cdata=read_data(glob(results_folder+'/*/cov*'))
    for isp,sp in enumerate(pdata):
        if sp not in products:
            continue
        x=pdata[sp][:,0]
        y=pdata[sp][:,1]
        color=exp.get_color(sp)
        if color is None:
            color=colorlist[sp]
        if j_log_plot:
            func=ax1.semilogy
        else:
            func=ax1.plot
        func(x,y,'-',color=color,label=sp) #,label=arg.split('/')[-1])
    for isp,sp in enumerate(cdata):
        x=cdata[sp][:,0]
        y=cdata[sp][:,1]
        color=exp.get_color(sp)
        if color is None:
            color=colorlist[sp]
        ax2.plot(x,y,'-',color=color,label=sp)
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
exp.plot_data(reference=['hori','jaramillo'],ax=ax1,species=['H$_2$','CO','CH$_4$','C2-sum'],pH=['6.8','7.2'],system=['pc-Cu'],scale='RHE',only_points=True,take_log=j_log_plot)
plt.show()
