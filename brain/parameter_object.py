import numpy as np
from enum import IntEnum

class PotentialChoice(IntEnum):
    '''A simple enumerator to aviod confusion.
    '''
    HARMONIC = 0
    HARMONIC_QUARTIC = 1
    HARMONIC_OPTIC = 2
    NOT_IMPLEMENTED = -1

class Psi0Choice(IntEnum):
    '''A simple enumerator to aviod confusion.
    '''
    THOMAS_FERMI = 0
    GAUSS = 1
    NOT_IMPLEMENTED = -1


class ParameterObject:

    #生成するときに引数を省略したらここに記載している値を用いる
    def __init__(self, resolutionX = 200, resolutionY = 200,
    x_low = -15, x_high = 15, y_low = -15, y_high = 15,
    beta2 = 1000, omega = 0.9,
    epsilon_limit=1e-10, epsilon_threshold=1, dt=0.005, maxIterations=30_000,
    ):
        '''Initializes the instance with all the given parameters.
        Contains all standard parameters.
        '''

        self.resolutionX = resolutionX
        self.resolutionY = resolutionY

        # set the bounddaries of the image box
        self.x_low = x_low
        self.x_high = x_high
        self.y_low = y_low
        self.y_high = y_high

        # calculate the spatial step and make a 2D coordinate array
        self.dx = (self.x_high - self.x_low)/self.resolutionX
        self.x = np.linspace(self.x_low, self.x_high, self.resolutionX)

        self.dy = (self.y_high - self.y_low)/self.resolutionY
        self.y = np.linspace(self.y_low, self.y_high, self.resolutionY)

        # constants for the BEC itself
        self.beta2 = beta2
        self.omega = omega

        # numerical parameters
        self.epsilon_limit = epsilon_limit
        self.epsilon_threshold = epsilon_threshold
        self.dt = dt
        self.maxIterations = maxIterations
        self.Vopt=0
    
    #周期ポテンシャル
    def initVperiodic(self, V0 = 1, kappa = np.pi):
        '''This function initializes a harmonic + periodic (optic) potential V(r) ~ r^2 + sin^2(k*r).
        '''
        xx, yy = np.meshgrid(self.x, self.y, sparse=False, indexing='ij')
        gammay=1.0
        arg=np.pi/6
        k1=kappa*np.array([1.0,0])
        k2=kappa*np.array([-np.cos(arg),np.sin(arg)])
        k3=kappa*np.array([-np.cos(arg),-np.sin(arg)])
        xx, yy = np.meshgrid(self.x, self.y, sparse=False, indexing='ij')
        self.Vopt = V0*(np.sin(k1[0]*xx+k1[1]*yy)**2 + np.sin(k2[0]*xx+k2[1]*yy)**2 + np.sin(k3[0]*xx+k3[1]*yy)**2 )
        self.V = 0.5 * (xx**2 + gammay*yy**2) + self.Vopt


    #追試
    def initV(self):
        '''Initializes the right potential based on potential_choice.
        '''
        gammay=1
        gammax=1
        xx, yy = np.meshgrid(self.x, self.y, sparse=False, indexing='ij')
        a1=[0,2.2]
        a2=[1.1*np.sqrt(3),-1.1]
        self.V=0
        for i in np.arange (-5,6,1):
            for j in range(-5,6,1):
                rx=i*a1[0]+j*a2[0]
                ry=i*a1[1]+j*a2[1]
                self.V+=6*np.exp(-((xx-rx)**2+(yy-ry)**2)/0.325**2)

        self.V +=0.5 * (gammax*xx**2 + gammay*yy**2)
