import numpy as np
from scipy.interpolate import interp1d

#read some co2 input data:
def read_data():
    data={} #collections.defaultdict(list)
    fname='data/Poly_Cu_full.txt'
#    n_pot = sum(1 for line in open(fname))-6
    with open(fname,'r') as of:
        for i,line in enumerate(of):
            if i==1:
                keys=line.split()
                for key in keys:
                    data[key]={}
            if i==5:
                boundary_thickness=float(line)
            elif i==4:
                bic_conc=float(line)
            if i>5:
                data_tmp=map(float,line.split())
                j=-1
                for key,d in zip(keys,data_tmp):
                    j+=1
                    if j==1:
                        data[key][str(data_tmp[0])]=-d/1000*100**2
                    elif j>1:
                        data[key][str(data_tmp[0])]=-d*data_tmp[1]/1000*100**2
    fname='data/viscosity.txt'
    data_vis=np.loadtxt(fname)
    viscosity=interp1d(data_vis[:,0],data_vis[:,1]/1000.)
    return data, boundary_thickness,viscosity,bic_conc*1000.
