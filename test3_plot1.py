import matplotlib.pyplot as plt
import numpy as np
import math 
L = np.arange(0, 2*np.pi, 0.01)
plt.figure(figsize=(16, 9),dpi=300)
plt.xlim((0,2*math.pi))
F=5*np.cos(L)
plt.plot(L, F,linewidth=0.5,color='b')
x1=[]
y1=[]
with open('Test3_1_Finite_difference.txt', 'r') as f:
    data = f.readlines()
    for line in data:
        odom = line.split(' ')
        x1.append(float(odom[0]))
        y1.append(float(odom[1]))
plt.scatter(x1,y1,s=3, c='#000000',alpha=0.7)
x2=[]
y2=[]
with open('Test3_1_Fourier_spectra.txt', 'r') as f:
    data = f.readlines()
    for line in data:
        odom = line.split(' ')
        x2.append(float(odom[0]))
        y2.append(float(odom[1]))
plt.scatter(x2,y2,s=3,c='r',alpha=0.7)
plt.show()