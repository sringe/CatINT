import re
import matplotlib.pyplot as plt
from glob import glob
import numpy as np
import sys

products=['CO'] #,'CO2','OH']
states=['COOH-H2O-ele','COOH','CO2']

for product in products:
    for state in states:
        pair=product+'_'+state
        files=glob('CO2R_results/catmap_output/desc_*/rc_'+pair+'.tsv')
        data=[]
        for f in files:
            print 'f',f
            pot=float(re.findall('desc\_(-?[0-9.]+)',f)[0])
            x,y=np.loadtxt(f)[0]
            data.append([x,y])
            print pot,y
        data=np.array(sorted(data))
        plt.plot(data[:,0],data[:,1],'-',label=pair)
plt.legend()
plt.show()
