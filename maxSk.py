import matplotlib.pyplot as plt
import numpy as np
import random

def calcmax(arr):
        arr=np.asarray(arr)
        lenx=arr.shape[0]
        leny=arr.shape[1]

        randmin=0.1
        randmax=1-randmin
        rand=[[random.randint(int(randmin*lenx),int(randmax*lenx)),random.randint(int(randmin*leny),int(randmax*leny))]for i in range(1000)]
        minlist=[]
        while len(rand)>0:
            [i,j]=rand.pop(0)
            flag=0
            while flag==0:
                if   arr[i][j]>arr[i][j-1]:j=j-1
                elif arr[i][j]>arr[i][j+1]:j=j+1
                elif arr[i][j]>arr[i-1][j]:i=i-1
                elif arr[i][j]>arr[i+1][j]:i=i+1
                else :
                    flag=1
                    minlist.append([i,j])
        #インデックス半径40以下のもののみ考える
        centerx=int(lenx/2)
        centery=int(leny/2)
        arrtmp=minlist
        minlist=[]
        distr=80
        for i in range(len(arrtmp)):
            dist=np.sqrt( (centerx-arrtmp[i][0])**2 + (centery-arrtmp[i][1])**2 )
            if dist<distr:minlist.append([arrtmp[i][0],arrtmp[i][1]])
        sk=np.zeros((lenx,leny))
        for i in range(len(minlist)):
            sk[minlist[i][0]][minlist[i][1]]=1
        minlist=[]
        step=30/lenx
        for i in range(lenx):
            for j in range(leny):
                if sk[i][j]==1:minlist.append([-15+step*j,15-step*i])
        
        #skには最小値のとこに１が格納されて、minlistにはその最小値の座標が格納されてる
        
        #フーリエ変換用
        Sk=np.zeros((500,500),dtype=complex)
        maxk=4.2
        dk=2*maxk/500
        for i in range(500):
            ky=maxk-i*dk
            for j in range(500):
                kx=-maxk+j*dk
                for k in range(len(minlist)):
                    tmp=minlist[k][0]*kx+minlist[k][1]*ky
                    Sk[i][j]+=np.exp(1j*tmp)
        Sk=np.abs(Sk)/len(minlist)
        
        Sk2=Sk[250:251,300:500]
        return Sk2

#omega=0の時の各aでのSkの最大値を求める
#a=1.2~2.0までの間0.01ずつでの最大値を取るインデックスを出力

V0=2.0
kappa=np.array([1.3,1.4,1.5,1.6,1.7,1.8,1.9])
a=[]
for k in kappa:
    gammay=1.03
    arg=np.pi/6
    k1=k*np.array([1.0,0])
    k2=k*np.array([-np.cos(arg),np.sin(arg)])
    k3=k*np.array([-np.cos(arg),-np.sin(arg)])
    x=np.linspace(-15,15,200)
    y=np.linspace(-15,15,200)
    xx, yy = np.meshgrid(x, y, sparse=False, indexing='ij')
    V = -V0*(np.sin(k1[0]*xx+k1[1]*yy)**2 + np.sin(k2[0]*xx+k2[1]*yy)**2 + np.sin(k3[0]*xx+k3[1]*yy)**2 )
    b=calcmax(V)
    tmp=0
    for i in range(1,200):
        if b[0][tmp]<b[0][i]: tmp=i
    a.append(300+tmp)

np.save("skdata",a)

        
