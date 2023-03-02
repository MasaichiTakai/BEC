import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import stats
 
#凝縮体の位相と位相欠陥の場所を出力

def findVortex(g):
        #-3から３までの値を取る
        vortex=np.zeros((200,200))
        minlist=[]
        step=30.0/200
        for i in range(200):
            for j in range(200):
                if np.linalg.norm(np.array([i,j])-np.array([100,100]))<50:
                    around=[g[i-1][j], g[i-1][j+1], g[i-1][j-1],
                            g[i][j],   g[i][j+1],   g[i][j-1],
                            g[i+1][j], g[i+1][j+1], g[i+1][j-1]]
                    #それぞれ第１象限から第４まで
                    flag1,flag2,flag3,flag4=0,0,0,0

                    for k in range(9):
                        if       -np.pi < around[k] and around[k] <= -np.pi*0.5:flag1=1
                        elif -np.pi*0.5 < around[k] and around[k] <= 0:flag2=1
                        elif          0 < around[k] and around[k] <= np.pi*0.5:flag3=1
                        elif  np.pi*0.5 < around[k] and around[k] <= np.pi:flag4=1
                    if flag1==1 and flag2==1 and flag3==1 and flag4==1:
                        if (vortex[i-1][j-1]!=1 and vortex[i-1][j]!=1 and vortex[i-1][j+1]!=1
                            and  vortex[i][j-1]!=1 and vortex[i][j]!=1 and vortex[i][j+1]!=1
                            and vortex[i+1][j-1]!=1 and vortex[i+1][j]!=1 and vortex[i+1][j+1]!=1):
                            vortex[i][j]=1
                            minlist.append([-15+step*j,15-step*i])
  
        return vortex

fol="text.npy"
ang=np.angle(np.load(fol))
vortex=findVortex(ang)
fig, ax = plt.subplots(2,1)
ax[0].imshow(vortex,vmin=0,vmax=1,cmap='jet')
ax[1].imshow(ang,vmin=-np.pi,vmax=np.pi,cmap='jet')
plt.show()