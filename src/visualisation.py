"""!
Opening file docstring contains LATEX parsing - \f$\beta\f$

"""

# def foo():
#     '''!
#     random docs
#     \f$\alpha\f$

import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import glob
import os
import numpy as np
Nchkpnts = len(os.listdir('out/')) - 1
grid_res = h5py.File('out/1.chkpnt', 'r')['c'][...].shape[0]

c = np.zeros((Nchkpnts, grid_res, grid_res))
c_prev = np.zeros((Nchkpnts, grid_res, grid_res))
dt = np.zeros((Nchkpnts))
for i in range(1, Nchkpnts):
    data = h5py.File('out/'+str(i)+'.chkpnt', 'r')
    test = data['c'][...]
    c[i, :, :] = test
    dt[:] = data['dt'][...]

ims = []
fig, ax = plt.subplots()

for k in range(0, len(c)):
    im = []
    shw = ax.imshow(c[k, :, :])
    im.append(shw)
    im.append(plt.text((len(c[0, :, :])-1)/2, -3*len(c[0, :, :]) /
                       15, s=("T = " + str(k+1)), ha="center", weight="bold"))
    ims.append(im)
fig.colorbar(shw)

final_an = anim.ArtistAnimation(fig, ims, interval=500, repeat=False)
plt.gca().set_aspect('equal', adjustable='box')

plt.axis('off')
plt.show()


# # f_c = np.zeros((grid_size,grid_size,len(t))) #Bulk free energy density
# # a = [2,6,9]#polynomial coefficients
# # for t in range(len(time)):
# #   for i in range(grid_size):
# #     for j in range(grid_size):
# #       f_c(i,j,t) = 0
# #       for n in range(len(a)):
# #         f_c(i,j,t) += a[n]*c[i,j,t]**(n-1)

# # F = np.zeros((len(t))) #Free energy functional


# # for t in range(len(time)):
# #   for i in range(grid_size):
# #     for j in range(grid_size):
# #       F(t) += f_c(i,j,t) + 0.5*K*np.linalg.norm(np.gradient(c(i,j,t)))*dA


# # # plt.text(count/2-0.5, -(count/15+count/25), s = dat.initial, ha = "center")
# # # plt.text(count/2-0.5, -(count/15+2*count/25), s = "N = " + str(dat.n) + ", T = " + str(dat.t) + ", J = " + str(dat.j) + ", B = " + str(dat.b), ha = "center")
# # # plt.text(count/2-0.5, -(count/15+3*count/25), s = "Compiled from: " + dat.file, ha = "center")
# # # plt.title(dat.title, weight = "bold")
# # # plt.xlabel("Time (" + mag.tunits +")")
# # # plt.ylabel("Magnetisation (" + mag.eunits +")")
# # # plt.title("Average spin versus time", weight = "bold")
# # grid_size = 3
