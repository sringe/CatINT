import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt('data/buffer_dependence.txt')
x=data[:,0]
y1=data[:,1]/6.kk
y2=data[:,2]

plt.bar(x,y1)
plt.bar(x,y2)
plt.show()
