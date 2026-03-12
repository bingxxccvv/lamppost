import numpy as np

data = np.loadtxt("/data/home/bingbing/mywork/raytransfer-my/taskcode/photons13/photons4trf_a9.66446e-01.i3.49448e-02.Mdl_0.00000e+00.Q_2.50000e-01.dat")

# 每个半径84行
block = 84

last = data[-block:]
prev = data[-2*block:-block]

new = np.zeros_like(last)

# r 用指数外推
r_last = last[:,0]
r_prev = prev[:,0]

r_new = r_last**2 / r_prev
new[:,0] = r_new

# 其他列线性外推
for i in range(1,11):
    new[:,i] = 2*last[:,i] - prev[:,i]

data_new = np.vstack([data,new])

np.savetxt("/data/home/bingbing/mywork/raytransfer-my/taskcode/photons13/photons4trf_a9.66446e-01.i3.49448e-02.Mdl_0.00000e+00.Q_2.50000e-01fix.dat",data_new,fmt="%.15Le")
