
import numpy as np
import numpy.fft as fft

from .parameter_object import ParameterObject

# implementation of the simpson integration rule for the calculation of the observables
def simpson(y, x):
    '''Simpson rule for any amount of intervals.
    Only works for equal intervals.
    '''
    N = len(x)
    if N%2==1:
        s = simpson_even(y, x)
    else:
        s = simpson_even(y[:-1], x[:-1])
        # trapeziod rule for last iterval
        s += (x[-1]-x[-2]) * (y[-1]+y[-2])/2
    return s

#@jit
def simpson_even(y, x):
    '''Simpson rule for even amount of intervals.
    Only works for equal intervals.
    '''
    N = len(x)
    a = x[0]
    b = x[-1]
    h = (b-a)/(N-1)
    
    s = 2*np.sum(y[2:-1:2]) + 4*np.sum(y[1:-1:2])
    s += y[0] + y[-1]
    return s * h/3


class WaveFunction2D:
    '''This class represents a complex wave function in 2 dimensions and contains all neccessary functions to manipulate it.
    Especially this class also contains functionality for the fourier transform psi_hat of the wave function and several other "components".

    Attributes:
        paramObj: A ParameterObject instance.
        psi_array: A 2D numpy array that contains complex values of the wavefunction.
        フーリエ変換
        psi_hat_array: A 2D numpy array that contains complex values of the fourier transform of the wavefunction.
        角運動量演算子
        L_psi_array: A 2D numpy array that contains complex values of the angular momentum operator applied to the wavefunction.
        ナブラ^2の演算子
        nabla_psi_array: A 2D numpy array that contains complex values of the Nabla operator applied to the wavefunction.
        E: The energy expectation value of the wavefunction.
        L_expectation: The angular momentum expectation value of the wavefunction.
        Nabla_expectation: The kinetic energy expectation value of the wavefunction.
    '''
    def __init__(self, parameterObject):
        
        self.paramObj = parameterObject

        self.psi_array = np.zeros((200,200)) + (0+0j)
        self.psi_hat_array = np.zeros((200,200)) + (0+0j)
        self.L_psi_array = np.zeros((200,200)) + (0+0j)
        self.nabla_psi_array = np.zeros((200,200)) + (0+0j)

        self.E = None
        self.L_expectation = None
        self.Nabla_expectation = None
    
    #波動関数を設定するだけ。time steeperで呼び出される
    def setPsi(self, array):
        self.psi_array = array
    #FFTされた波動関数
    def setPsiHat(self, array):
        self.psi_hat_array = array
    #初期値の設定。トーマスフェルミ
    def initPsi_0(self):
        gamma_y=1.0
        my_g = 0.5*np.sqrt(4*self.paramObj.beta2*gamma_y*1.0)
        self.psi_array = np.sqrt(np.maximum(0, my_g - self.paramObj.V)/self.paramObj.beta2)
        self.norm()
    #規格化。毎ステップ行う
    def norm(self):
        self.psi_array /= np.sqrt( np.sum(np.abs(self.psi_array[1:-1, 1:-1])**2) * self.paramObj.dx * self.paramObj.dy )
        return self.psi_array

    def getNorm(self):
        return np.sqrt( np.sum(np.abs(self.psi_array[1:-1, 1:-1])**2) * self.paramObj.dx * self.paramObj.dy )
    #波動関数のFFTを計算。計算結果はself.psi_hat_arrayに入れる
    def calcFFT(self):
        self.psi_hat_array = np.fft.fft2(self.psi_array) /(np.prod((256,256)))
    #波数空間上の波動関数をユークリッド空間に戻す
    def calcIFFT(self):
        self.psi_array = np.fft.ifft2(self.psi_hat_array) * (np.prod((256,256)))
    #エネルギーの期待値を返す
    def calcEnergy(self):
        # set up aliases
        x = self.paramObj.x
        y = self.paramObj.y

        # add each term of the GPE energy [1] to avoid mistakes
        #dE = 0.5*np.abs(self.nabla_psi_array)**2 
        dE = self.paramObj.Vopt * np.abs(self.psi_array)**2
        #dE += self.paramObj.beta2/2 * np.abs(self.psi_array)**4
        #dE -= self.paramObj.omega * (np.conjugate(self.psi_array)*self.L_psi_array).real

        # 2d integrate by using the simpson rule twice on different axes
        self.E = simpson(simpson(dE, y), x)
        return self.E
    #角運動量の期待値を返す
    def calcL_expectation(self):
        # set up aliases
        x = self.paramObj.x
        y = self.paramObj.y

        # set up integrant
        dL = np.conjugate(self.psi_array)*self.L_psi_array

        # 2d integrate by using the simpson rule twice on different axes
        self.L_expectation = simpson(simpson(dL, y), x)
        return self.L_expectation
    #運動エネルギー
    def calcNabla_expectation(self):
        '''Calculates the expectation value for the kinetic energy.
        '''
        # set up aliases
        x = self.paramObj.x
        y = self.paramObj.y

        # set up integrant
        dN = 0.5*np.abs(self.nabla_psi_array)**2

        # 2d integrate by using the simpson rule twice on different axes
        self.Nabla_expectation = simpson(simpson(dN, y), x)
        return self.Nabla_expectation
    #波動関数の勾配を返す
    def calcNabla(self):
        """Calculates the spatial derivative operator nabla applied to the wavefunction.
        Retruns:
            The Wavefunction after the nabla operator was applied.
        """
        # set up aliases
        a, b, c, d = (-15,15,-15,15)
        M, N = (200,200)
        
        # set up grid in fourier space
        p = np.arange(-M//2, M//2, 1)
        q = np.arange(-N//2, N//2, 1)
        pp, qq = np.meshgrid(p, q, indexing='ij')

        # set up 'derivation' constants in fourier space
        lambda_q = 2*qq*np.pi/(d-c)
        my_p = 2*pp*np.pi/(b-a)

        # shift 'derivation' constants, will give wrong results otherwise
        my_p = np.fft.fftshift(my_p)
        lambda_q = np.fft.fftshift(lambda_q)

        #do spatial derivation by fourier transforming twice
        n_psi = WaveFunction2D(self.paramObj)
        n_psi.setPsiHat((my_p + lambda_q) * self.psi_hat_array)
        n_psi.calcIFFT()
        self.nabla_psi_array = n_psi.psi_array

        return self.nabla_psi_array

    #@jit
    def calcL(self):
        # set some aliases
        a, b, c, d = (-15,15,-15,15)
        M, N = (200,200)

        # calculating D_x(Psi) and D_y(Psi) first to later Calculate L
        # set up grid in fourier space
        p = np.arange(-M//2, M//2, 1)
        q = np.arange(-N//2, N//2, 1)
        pp, qq = np.meshgrid(p, q, indexing='ij')

        # set up 'derivation' constants in fourier space
        lambda_q = 2*qq*np.pi/(d-c)
        my_p = 2*pp*np.pi/(b-a)

        # shift 'derivation' constants, will give wrong results otherwise
        my_p = np.fft.fftshift(my_p)
        lambda_q = np.fft.fftshift(lambda_q)

        # x偏微分での波動関数のFFT
        Dx_psi = WaveFunction2D(self.paramObj)
        Dx_psi.setPsiHat(my_p * self.psi_hat_array)
        Dx_psi.calcIFFT()

        # y偏微分での波動関数のFFT
        Dy_psi = WaveFunction2D(self.paramObj)
        Dy_psi.setPsiHat(lambda_q * self.psi_hat_array)
        Dy_psi.calcIFFT()

        # adding Dx and Dy up to L
        xx, yy = np.meshgrid(self.paramObj.x, self.paramObj.y, sparse=False, indexing='ij')

        self.L_psi_array = xx * Dy_psi.psi_array - yy * Dx_psi.psi_array

    #PDEの計算で使う（差分法のとこ）
    def calcG_m(self, psi_m, alpha):
        
        # adding each term of G_m seperately to avoid mistakes.
        # using formula from [1]
        g = alpha * psi_m.psi_array + 0j
        g -= self.paramObj.V * psi_m.psi_array 
        g -= self.paramObj.beta2 * np.abs(self.psi_array)**2 * psi_m.psi_array 
        g += self.paramObj.omega * psi_m.L_psi_array
        return g
