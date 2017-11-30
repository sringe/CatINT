from scipy.optimize import fmin
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
from glob import glob


files=glob("results_x*")

toplot=[]
for file in files:
    data=np.loadtxt(file)
    fi=interp1d(data[:,0],data[:,1])
    xval=float(file.split('_')[1].split('.')[0][1:])
    xval2=float(file.split('_')[2].split('.')[0][2:])
    
    def f(xx):
        data=[]
        print xx
        if len(xx)>1:
            for x in xx:
                if x<min(fi.x):
                    data.append(1e10*(x-min(fi.x))**2+fi(min(fi.x)))
                elif x>max(fi.x):
                    data.append(1e10*(x-max(fi.x))**2+fi(max(fi.x)))
                else:
                    data.append(fi(x))
        else:
            if type(xx)==float:
                x=xx
            else:
                x=xx[0]
            if x<min(fi.x):
                return 1e10*(x-min(fi.x))**2+fi(min(fi.x))
            elif x>max(fi.x):
                return 1e10*(x-max(fi.x))**2+fi(max(fi.x))
            else:
                return fi(x)
        return data

    xfine=np.linspace(min(fi.x)-1e-6,max(fi.x)+1e-6,1000)
#    plt.plot(xfine,f(xfine),'-')
#    plt.plot(data[:,0],data[:,1],'o')
#    plt.show()
    init_guess=10**np.average(np.log10(data[:,0]))
    xopt, fopt, iter, funcalls, warnflag=fmin(f,np.average(data[:,0]),full_output=True, disp=False)
    toplot.append([xval,xopt[0],xval2])
    print xval, xopt[0]
toplot=np.array(toplot)
plt.scatter(toplot[:,0],toplot[:,1],c=toplot[:,2])
plt.autoscale()
plt.show()
