#親ディレクトリをインポートするため
import sys
sys.path.append("..")
from brain import ParameterObject, WaveFunction2D, ImaginaryTimeStepper
import numpy as np
import os


#### initialize parameters
#resolution  解像度
res_x = 200
res_y = 200
x_low = -15
x_high = 15
y_low = -15
y_high = 15
beta2 = 1000
epsilon_limit=1e-8
#threshold 閾値
epsilon_threshold=1
dt=0.01
maxIterations=100000
omega = 0.77

#p.initV()
#psi0.setPsi(np.load("data/text.npy"))

# a=1.6
# paramarr=[[0.02,"a"]]
# paramarr=np.array(paramarr)
# for j in range(paramarr.shape[0]):
#     p.initVperiodic(float(paramarr[j][0]), a)
#     psi0 = WaveFunction2D(p)
#     psi0.initPsi_0()
#     i = ImaginaryTimeStepper(psi0, p)
#     i.BFFP(paramarr[j][1])


if __name__ == '__main__':
    args = sys.argv
    foldername="v="+args[1]+"a="+args[2]
    os.makedirs(foldername,exist_ok=True)
    # p = ParameterObject(res_x, res_y, x_low, x_high, y_low, y_high,
    #                 beta2, 0    , epsilon_limit, epsilon_threshold, dt, maxIterations)
    # p.initVperiodic(float(args[1]),float(args[2]) )
    # psi0 = WaveFunction2D(p)
    # psi0.initPsi_0()
    # i = ImaginaryTimeStepper(psi0, p)
    energy0=1#i.BFFPE(foldername)

    p = ParameterObject(res_x, res_y, x_low, x_high, y_low, y_high,
                    beta2, omega, epsilon_limit, epsilon_threshold, dt, maxIterations)

    # p.initVperiodic(float(args[1]),float(args[2]) )
    p.initV( )
    psi0 = WaveFunction2D(p)
    psi0.initPsi_0()
    #psi0.setPsi(np.load(foldername+"/text.npy"))
    i = ImaginaryTimeStepper(psi0, p)
    i.BFFP(foldername,args[1],args[2],energy0)
