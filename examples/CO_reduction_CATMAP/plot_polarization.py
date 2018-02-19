from itertools import cycle
import numpy as np
import matplotlib.pyplot as plt
import sys
from glob import glob
from catint.experimental import EXPDATA

exp=EXPDATA()

j_log_plot=True

fig=plt.figure(figsize=(6, 7))

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
symbols={
        'CO':'x',
        'H2':'o',
        'COOH':'d',
        'HCOOH':'3',
        'CO2':'D',
        'CH4':'p',
        'CH3CH2OH':'1',
        'CHO':'2'
        }

products=['CO','CH4','CH3CH2OH','H2','HCOOH']

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

linestyles=cycle(['-','--',':'])

k=-1
for arg in sys.argv[1:]:
    linestyle=next(linestyles)
    k+=1
    results_folder=arg+'/catmap_output'
    
    ax1=plt.subplot('211')
    ax2=plt.subplot('212')
    
    ax1.set_title('Polarization')
    pdata=read_data(glob(results_folder+'/*/j*'))
    cdata=read_data(glob(results_folder+'/*/cov*'))
    for isp,sp in enumerate(pdata):
        print 'checking',sp
        if sp not in products:
            continue
        x=pdata[sp][:,0]
        y=pdata[sp][:,1]
        color=exp.get_color(sp)
        if color is None:
            color=colorlist[sp]
        if sp in symbols:
            symbol=symbols[sp]
        else:
            symbol=''
#        symbol=''
        if j_log_plot:
            func=ax1.semilogy
        else:
            func=ax1.plot
        if k==0:
            func(x,y,linestyle+symbol,color=color,label=sp) #,label=arg.split('/')[-1])
        else:
            func(x,y,linestyle,color=color)
    for isp,sp in enumerate(cdata):
        x=cdata[sp][:,0]
        y=cdata[sp][:,1]
        color=exp.get_color(sp)
        if color is None:
            color=colorlist[sp]
        if sp in symbols:
            symbol=symbols[sp]
        else:
            symbol=''
#        symbol=''
        print 'the color',sp,color
        if k==0:
            ax2.plot(x,y,linestyle+symbol,color=color,label=sp)
        else:
            ax2.plot(x,y,linestyle+symbol,color=color)
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
print 'plotting exp'
exp.plot_data(reference=['hori','jaramillo'],ax=ax1,species=['H$_2$','CO','CH$_4$','C2-sum','HCOOH'],pH=['6.8','7.2'],system=['pc-Cu'],scale='RHE',only_points=True,take_log=j_log_plot)
ax1.set_ylim([1e-12,1e4])
ax1.set_xlim([-1.2,0])
ax2.set_xlim([-1.2,0])
plt.show()
