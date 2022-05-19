"""!
Opening file docstring contains LATEX parsing - \f$\beta\f$

"""

# def foo():
#     '''!
#     random docs
#     \f$\alpha\f$

from tkinter import N
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import glob
import os
import numpy as np

def read_metadata(filename):
  with open(filename) as f:

    grid_params = f.readline().split()[1:]
    sys_params = f.readline().split()[1:]
    sys_params = sys_params + f.readline().split()[1:]

  chkpnt_times = np.genfromtxt(filename, skip_header=4)

  return np.array(grid_params, dtype = int), np.array(sys_params, dtype = float), np.array(chkpnt_times[:,1], dtype = float)


Nchkpnts = len(os.listdir('out/')) -1
grid_res = h5py.File('out/1.chkpnt', 'r')['c'][...].shape[0]
grid_params,sys_params,t_array = read_metadata('out/metadata.dat')
dimension = grid_params[0]
gridres = 2**grid_params[1]
L,A,M,K,p0,p1 = sys_params
c = np.zeros((Nchkpnts,grid_res,grid_res))
c_prev = np.zeros((Nchkpnts,grid_res,grid_res))
dt = np.zeros((Nchkpnts))

for i in range(1,Nchkpnts):
    data = h5py.File('out/'+str(i)+'.chkpnt', 'r')
    test= data['c'][...]
    c[i,:,:] = test
    dt[i] = data['dt'][...]

ims =[]
fig, ax = plt.subplots()

for k in range(0, Nchkpnts):
    im = []
    shw = ax.imshow(c[k,:,:])
    im.append(shw)
    im.append(plt.title('Evolution of the concentration', weight = "bold"))
    im.append(plt.text(grid_res*0.425,grid_res*1.1 , s = ("T = " + str(t_array[k])), weight = "bold"))
    ims.append(im)

fig.colorbar(shw)
final_an = anim.ArtistAnimation(fig, ims, interval = 5, repeat = True)
final_an.save('Conc_Evolution.mp4', fps = 10)
plt.gca().set_aspect('equal', adjustable='box')
plt.axis('off')

plt.clf()

f_c = np.zeros((Nchkpnts,grid_res,grid_res)) #Bulk free energy density 
for t in range(Nchkpnts):
  for i in range(grid_res):
    for j in range(grid_res):
        f_c[t,i,j] = A*((c[t,i,j]-p0)**2)*(c[t,i,j]-p1)**2 
F = np.zeros((Nchkpnts)) #Free energy functional


c_halo = np.zeros((Nchkpnts, grid_res+2, grid_res+2))
c_halo[:, 1:-1, 1:-1] = c
c_halo[:,0,:] = c_halo[:,-2,:]
c_halo[:,-1,:] = c_halo[:,1,:]
c_halo[:,:,0] = c_halo[:,:,-2]
c_halo[:,:,-1] = c_halo[:,:,1]
for t in range(Nchkpnts):
  grad_c_x = np.gradient(c_halo[t,:,:], axis=0)[1:-1,1:-1]
  grad_c_y = np.gradient(c_halo[t,:,:], axis=1)[1:-1,1:-1]
  grad_c_sqd = (grad_c_x ** 2 + grad_c_y ** 2)
  for i in range(grid_res):
    for j in range(grid_res):
       F[t] += f_c[t,i,j] + 0.5*K*grad_c_sqd[i,j]
F = F/(grid_res**2)

plt.plot(t_array,F)

plt.title('Free Energy Functional', weight = "bold")
plt.ylabel('Energy', weight = "bold")
plt.xlabel('Time', weight = "bold")
plt.savefig('Free_Energy.png')

