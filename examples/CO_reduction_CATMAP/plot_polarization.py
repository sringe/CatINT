import numpy as np
import matplotlib.pyplot as plt
import sys
from glob import glob

results_folder=sys.argv[1]+'/catmap_output'

ax1=plt.subplot('211')
ax2=plt.subplot('212')

ax1.set_title('Polarization')
iarg=-1
for arg in glob(results_folder+'/*/j*'):
    iarg+=1
    data=np.loadtxt(arg)
    x=data[:,0]
    y=data[:,1]
    #if any([yy<0 for yy in y]):
    #    y=-y
    #y=abs(y)
    if iarg==0:
        ax1.plot(x,y,'.',label=arg.split('/')[-1])
    else:
        ax1.plot(x,y,'.')
ax1.legend()
data=np.loadtxt('data.txt')
ax1.plot(data[:,0],data[:,1],'-')
ax2.set_title('Coverages')
iarg=-1
for arg in glob(results_folder+'/cov*'):
    data=np.loadtxt(arg)
    x=data[:,0]
    y=data[:,1]
    if iarg==0:
        ax2.plot(x,y,'.',label=arg.split('/')[-1])
    else:
        ax2.plot(x,y,'.')
ax2.legend()
plt.show()
