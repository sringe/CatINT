import numpy as np
import sys
import matplotlib.pyplot as plt

for arg in sys.argv[1:]:
    data=np.loadtxt(arg,delimiter=',')
    plt.semilogy(data[:,0],data[:,1],'-o')
    plt.semilogy(data[:,0],data[:,2],'-o')
    plt.semilogy(data[:,0],data[:,3],'-o')
    plt.semilogy(data[:,0],data[:,4],'-o')
    plt.semilogy(data[:,0],data[:,5],'-o')
plt.show()
