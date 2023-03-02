
from .wave_function import WaveFunction2D
from .parameter_object import ParameterObject
import matplotlib .pyplot as plt
import random

import numpy as np
import numpy.fft as fft


class ImaginaryTimeStepper:
    '''This class is used to solve the Gross-Pitaevski equation (GPE).

    Attributes:
        psi_0: A WaveFunction2D instance that has the initial wavefunction in it.
        parameterObject: A ParameterObject instance that contains all setting and parameters for the simulation.
        dt: A float that is an alias for ParameterObject.dt
        epsilon_iteration_step_limit: A float that is an alias for ParameterObject.epsilon_limit
        maxIterations: A float that is an alias for ParameterObject.maxIterations
        psi_n: A WavveFunction2D instance that contains the current wave function at time step t_n
        n: An integer that refers to the current time step
        dataM: A DataManager instance that is handling the saving of the time steps to the disk.
        globalAttributes: A dictionary that has the contents of ParameterObject but the right format for the DataManager.
        '''
    def __init__(self, psi_0, parameterObject):
        '''Initializes an instance of this class.

        Arguments:
            psi_0: A WaveFunction2D instance that has the initial wavefunction in it.
            parameterObject: A ParameterObject instance that contains all setting and parameters for the simulation.
        '''

        self.psi_0 = psi_0
        self.paramObj = parameterObject
        self.resultPsi=psi_0

        # set up aliases for parameterObject
        self.dt = self.paramObj.dt
        self.epsilon_iteration_step_limit = self.paramObj.epsilon_limit
        self.maxIterations = self.paramObj.maxIterations

        # 各時間ステップでの波動関数は全てここに保存される
        self.psi_n = self.psi_0
        self.n = 0
        self.count=1

    #その時点での波動関数の二乗を返す
    def returnFrame(self):
        '''Returns |psi|**2 of the current time step.
        '''
        return np.abs(self.psi_n.psi_array)**2
    #その時点での波動関数を返す
    def returnPsi(self):
        '''Returns psi of the current time step.
        '''
        return self.psi_n.psi_array
    #PDEでの計算に使われるアルファの計算。BEEPで呼び出される
    def calcAlpha(self):
        # This follows [1] and [2]. A stabilization parameter is used for faster convergence.
        b_ = self.paramObj.V + self.paramObj.beta2*np.abs(self.psi_n.psi_array)**2
        bmin = np.min(b_)
        bmax = np.max(b_)
        alpha = 0.5 * (bmax + bmin)
        if self.count==1:
            #print(1/alpha)
            self.count=2
        return alpha
    #次の時刻の波動関数をアルファとGを用いて計算。BEEPで呼び出される
    #@jit
    def calcNextPsi_m(self, G_m, alpha):
        # set up aliases for grid parameters
        a, b, c, d = (-15,15,-15,15)
        M, N = (200,200)

        #  計算後の結果はここに保存
        psi_m_next = WaveFunction2D(self.paramObj)
        
        # set up a grid in fourier space
        p = np.arange(-M//2, M//2, 1)
        q = np.arange(-N//2, N//2, 1)
        pp, qq = np.meshgrid(p, q, indexing='ij')

        # set up 'derivation' constants in fourier space
        lambda_q = 2*qq*np.pi/(d-c)
        my_p = 2*pp*np.pi/(b-a)

        # shift 'derivation' constants, will give wrong results otherwise
        my_p = np.fft.fftshift(my_p)
        lambda_q = np.fft.fftshift(lambda_q)

        # 本丸
        GPdt=self.dt#*-1j
        psi_m_hat = self.psi_n.psi_hat_array + GPdt * G_m.psi_hat_array
        psi_m_hat *= 2 / (2 + GPdt*(2*alpha + my_p**2 + lambda_q**2))

        # return the result
        psi_m_next.setPsiHat(psi_m_hat)
        return psi_m_hat
    #波動関数の規定状態を計算  
    #空間微分はFFTで、時間微分はセミインプリシットで。dtは小さく
    def BFFP(self,fol,vMax,aperiod,energy0):

        # set up epsilon
        epsilon_iteration_step = 1
        epsilon_sum = 0
        psi_max_old = np.zeros((200,200))
        #計算前後での変化量を見るために計算前の波動関数を保存する変数
        psi_max = self.psi_0.psi_array
        self.n = 0
        f=open(fol+"/file.txt","w")

        G_m = WaveFunction2D(self.paramObj)

    
        # 虚時間発展
        # conditions for loop are exit conditions
        while epsilon_iteration_step > self.epsilon_iteration_step_limit and self.n < self.maxIterations:
            
            # Aを計算
            alpha = self.calcAlpha()

            # 線形の部分は反復が一回で良い
            # 波動関数のFTT、角運動量、Gを計算
            self.psi_n.calcFFT()
            self.psi_n.calcL()
            G_m.setPsi(self.psi_n.calcG_m(self.psi_n, alpha))
            G_m.calcFFT()

            # 次の波動関数を計算
            self.n += 1
            self.psi_n.setPsiHat(self.calcNextPsi_m(G_m, alpha))
            self.psi_n.calcIFFT()
            self.psi_n.norm()

            # calculate epsilon
            psi_max_old = psi_max
            psi_max = self.psi_n.psi_array
            epsilon_iteration_step = np.max(np.abs(psi_max_old - psi_max)) / self.dt #最も変化したところを出力
            epsilon_sum += epsilon_iteration_step

            # 変化量がある値を超えるとその時の波動関数を保存する（ここでは省略）
            if epsilon_sum > self.paramObj.epsilon_threshold:
                f.write("frame is {}, and epsilon is {}.\n".format(self.n,format(epsilon_iteration_step,"1E")))
                #print("frame is {}, and epsilon is {}.".format(self.n,format(epsilon_iteration_step,"1E")))
                epsilon_sum = 0
            
            if self.n%100000==0:
                self.plotD(self.psi_n,f,fol)
                #np.save(fol+"/text",self.psi_n.psi_array)

        #ループ終了後の操作
        #f.write("last frame is {}, and epsilon is {}.\n".format(self.n,format(epsilon_iteration_step,"1E")))
        #print("frame is {}, and epsilon is {}.".format(self.n,format(epsilon_iteration_step,"1E")))
        np.save(fol+"/text",self.psi_n.psi_array)
        
        skdata=np.load("skdata.npy")
        aperiod=float(aperiod)
        v1=int(skdata[int((aperiod-1.2)*100)])-250
        v2=[250-int(v1/2*np.sqrt(3)), 250+int(v1/2)]
        v3=[v2[0], 250-int(v1/2)]
        v1=[250,v1+250]
        
        if epsilon_iteration_step <= self.epsilon_iteration_step_limit:
            self.plotD(self.psi_n,f,fol)
        
        #[ri,ru,lu,median]=self.calcSk(np.abs(self.psi_n.psi_array)**2,f,fol,vMax,v1,v2,v3)
        E=1#self.psi_n.calcEnergy()/energy0
        #print("last energy  is ",E)
        f.write("last Energy={}\n".format(E))
        f.close()
        #print("{},{},{},{},{},{},{},{}".format(vMax,aperiod,round(E,5),round(ri,5),round(ru,5),round(lu,5),round(median,5),round(np.log10(epsilon_iteration_step),5)))
        
    def plotD(self,psi_arr,f,fol):
        #波動関数をプロットする
        # E=psi_arr.calcEnergy()
        # f.write("Energy={}\n".format(E))
        plt.figure(figsize=(10,8))
        min = np.min(np.abs(psi_arr.psi_array)**2)
        max = np.max(np.abs(psi_arr.psi_array)**2)
        plt.title(self.n)
        plt.imshow(np.abs(psi_arr.psi_array)**2,vmin=min,vmax=max,cmap='jet')
        plt.colorbar()
        filename=fol+"/"+str(self.n/1000000)+".png"
        plt.savefig(filename)

    def calcSk(self,arr,f,fol,vMax,v1,v2,v3):
        arr=np.asarray(arr)
        lenx=200
        leny=200
        vMax=float(vMax)
        # randmin=0.2
        # randmax=1-randmin
        # rand=[[random.randint(int(randmin*lenx),int(randmax*lenx)),random.randint(int(randmin*leny),int(randmax*leny))]for i in range(10000)]
        # minlist=[]
        # while len(rand)>0:
        #     [i,j]=rand.pop(0)
        #     flag=0
        #     while flag==0:
        #         if i-1==-1 or j-1==-1 or i+1==200 or j+1==200:break

        #         if   arr[i][j]>arr[i][j-1]:j=j-1
        #         elif arr[i][j]>arr[i][j+1]:j=j+1
        #         elif arr[i][j]>arr[i-1][j]:i=i-1
        #         elif arr[i][j]>arr[i+1][j]:i=i+1
        #         else :
        #             flag=1
        #             minlist.append([i,j])
        # #インデックス半径40以下のもののみ考える
        # centerx=int(lenx/2)
        # centery=int(leny/2)
        # arrtmp=minlist
        # minlist=[]
        # distr=40
        # for i in range(len(arrtmp)):
        #     dist=np.sqrt( (centerx-arrtmp[i][0])**2 + (centery-arrtmp[i][1])**2 )
        #     if dist<distr:minlist.append([arrtmp[i][0],arrtmp[i][1]])
        # sk=np.zeros((lenx,leny))
        # for i in range(len(minlist)):
        #     sk[minlist[i][0]][minlist[i][1]]=1
        # minlist=[]
        # step=30/lenx
        # for i in range(lenx):
        #     for j in range(leny):
        #         if sk[i][j]==1:minlist.append([-15+step*j,15-step*i])
        
        
        # fig, ax = plt.subplots()
        # ax.imshow(sk,vmin=0,vmax=1,cmap='jet')
        # plt.savefig(fol+"/sk1.png")

        minlist=self.findVortex(fol,arr,vMax)
        

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
        plt.figure(figsize=(10,8))
        max=np.max(Sk)
        min=np.min(Sk)
        
        plt.imshow(Sk,vmin=min,vmax=max,cmap='jet')
        plt.savefig(fol+"/sk2.png")
        med=(Sk[v1[0]][v1[1]] + Sk[v2[0]][v2[1]] + Sk[v3[0]][v3[1]] )/3.0
        f.write("right Sk={}\n".format(Sk[v1[0]][v1[1]]))
        f.write("upright Sk={}\n".format(Sk[v2[0]][v2[1]]))
        f.write("upleft Sk={}\n".format(Sk[v3[0]][v3[1]]))
        f.write("median of Sk={}\n".format(med))
        return [Sk[v1[0]][v1[1]],Sk[v2[0]][v2[1]],Sk[v3[0]][v3[1]],med]
    

    def BFFPE(self,fol):

        # set up epsilon
        epsilon_iteration_step = 1
        epsilon_sum = 0
        psi_max_old = np.zeros((200,200))
        #計算前後での変化量を見るために計算前の波動関数を保存する変数
        psi_max = self.psi_0.psi_array
        self.n = 0
        f=open(fol+"/file.txt","w")

        G_m = WaveFunction2D(self.paramObj)

    
        # 虚時間発展
        # conditions for loop are exit conditions
        while epsilon_iteration_step > self.epsilon_iteration_step_limit and self.n < self.maxIterations:
            
            # Aを計算
            alpha = self.calcAlpha()

            # 線形の部分は反復が一回で良い
            # 波動関数のFTT、角運動量、Gを計算
            self.psi_n.calcFFT()
            self.psi_n.calcL()
            G_m.setPsi(self.psi_n.calcG_m(self.psi_n, alpha))
            G_m.calcFFT()

            # 次の波動関数を計算
            self.n += 1
            self.psi_n.setPsiHat(self.calcNextPsi_m(G_m, alpha))
            self.psi_n.calcIFFT()
            self.psi_n.norm()

            # calculate epsilon
            psi_max_old = psi_max
            psi_max = self.psi_n.psi_array
            epsilon_iteration_step = np.max(np.abs(psi_max_old - psi_max)) / self.dt #最も変化したところを出力
            epsilon_sum += epsilon_iteration_step

            # 変化量がある値を超えるとその時の波動関数を保存する（ここでは省略）
            if epsilon_sum > self.paramObj.epsilon_threshold:
                epsilon_sum = 0
            

        #ループ終了後の操作
        
        E=self.psi_n.calcEnergy()
        #print("last energy at omega=0  is ",E)
        f.write("last Energy={}\n".format(E))
        f.close()
        return E

    def findVortex(self,fol,arr,vMax):
        #-3から３までの値を取る
        g=np.angle(np.load(fol+"/text.npy"))
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
        fig, ax = plt.subplots()
        ax.imshow(vortex,vmin=0,vmax=1,cmap='jet')
        plt.savefig(fol+"/sk1.png")
        np.save("aiueo",minlist)
        return minlist
