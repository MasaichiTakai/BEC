from tokenize import PseudoExtras
import numpy.fft as fft
import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib import cm
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.colors import Normalize
    
x=np.linspace(-10,10,1000)
y=np.linspace(-10,10,1000)
xx, yy = np.meshgrid(x, y, sparse=False, indexing='ij')
V=0.1*(xx**2+yy**2)
my_g = 0.5*np.sqrt(4*1000)
psi_array = 0.3*np.ones((1000,1000))#np.sqrt(np.maximum(0, my_g - V)/1000)
for i in range(1000):
    for j in range(1000):
        psi_array[i][j]/=(1+3*np.exp((-(i-500)**2-(j-500)**2)/4000))

fig = plt.figure(figsize = (8, 8))
ax = fig.add_subplot(111, projection="3d")
ax.plot_surface(xx,yy,psi_array, cmap = 'jet',rstride=2, cstride=2,
norm=Normalize(vmin=0, vmax=0.34))


ax.view_init(elev=90, azim=0)
ax.set_zlim(0,0.4)
ax.grid(False)

plt.show()