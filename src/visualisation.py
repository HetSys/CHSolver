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
import numpy as np

def read_metadata(filename):
  '''! Reads the grid parameters, system parameters and checkpoint times from the metadata.dat file in the out directory.
  '''
  with open(filename) as f:
    grid_params = f.readline().split()[1:]
    sys_params = f.readline().split()[1:]
    sys_params = sys_params + f.readline().split()[1:]
  chkpnt_times = np.genfromtxt(filename, skip_header=4)
  return np.array(grid_params, dtype = int), np.array(sys_params, dtype = float), np.array(chkpnt_times[:,1], dtype = float)

def read_hdf5():
  '''! Reads concentration at current timestep (c), concentration at previous timestep (c_prev) and the corresponding timestep (dt) from the collection of HDF5 checkpoint files in the out/ directory.
  '''
  c = np.zeros((Nchkpnts,grid_res,grid_res))
  c_prev = np.zeros((Nchkpnts,grid_res,grid_res))
  dt = np.zeros((Nchkpnts))
  for i in range(1,Nchkpnts):
      data = h5py.File('out/'+str(i)+'.chkpnt', 'r')
      test= data['c'][...]
      c[i,:,:] = test
      dt[i] = data['dt'][...]
  return c, c_prev, dt

grid_params,sys_params,t_array = read_metadata('out/metadata.dat')
dimension = grid_params[0]
grid_res = 2**grid_params[1]
L,A,M,K,p0,p1 = sys_params
Nchkpnts = len(t_array)
c, c_prev, dt = read_hdf5()

def find_nearest_t_index(t_array, t):
    '''! Produces an animation for the evolution of species concentration, outputted as an mp4 to the main directory.
    @param animation_fps Number of frames per second for the animation
    '''
    t_array = np.asarray(t_array)
    index = (np.abs(t_array - t)).argmin()
    return index 

def plot_conc_evol(animation_fps = 10 , ti = 0, tf = t_array[-1]):
  '''! Produces an animation for the evolution of species concentration, outputted as an mp4 to the main directory.
  @param animation_fps Number of frames per second for the animation
  @param ti Time to intialise the visualiation (will use time closest to that specified)
  @param tf Time to finish the visualiation (will use time closest to that specified)
  '''
  ti_index = find_nearest_t_index(t_array, ti)
  tf_index = find_nearest_t_index(t_array, tf)
  ims = []
  fig, ax = plt.subplots()
  for k in range(ti_index, tf_index+1):
      im = []
      ax.set_xticks([])
      ax.set_yticks([])
      shw = ax.imshow(c[k,:,:])
      im.append(shw)
      im.append(plt.title('Concentration Evolution', weight = "bold", fontsize = 14))
      im.append(plt.text(grid_res/2,grid_res*1.06 , s = ("t = " + str(f"{t_array[k]:.3f}")), weight = "bold", fontsize = 12, ha = 'center'))
      ims.append(im)
  fig.colorbar(shw)
  ax.set_axis_off()
  fig.add_axes(ax)
  final_an = anim.ArtistAnimation(fig, ims, interval = 5, repeat = True)
  final_an.save('Conc_Evolution.gif', writer = anim.PillowWriter(fps = animation_fps))
  return None

def plot_free_energy(ti = 0, tf = t_array[-1]):
  '''! Produces a plot of Free energy versus time, outputted as a png to the main directory.
  @param ti Time to intialise the visualiation (will use time closest to that specified)
  @param tf Time to finish the visualiation (will use time closest to that specified)
  '''
  ti_index = find_nearest_t_index(t_array, ti)
  tf_index = find_nearest_t_index(t_array, tf)
  f_c = np.zeros((Nchkpnts,grid_res,grid_res)) #Bulk free energy density 
  for t in range(Nchkpnts):
    for i in range(grid_res):
      for j in range(grid_res):
          f_c[t,i,j] = A*((c[t,i,j]-p0)**2)*(c[t,i,j]-p1)**2 

  F = np.zeros((Nchkpnts)) #Free energy functional
  ## Setting up an alternate c array to account for PBCs
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
        F[t] += f_c[t,i,j] + 0.5*K*grad_c_sqd[i,j]## dV = 1? so no term in here
  F = F/(grid_res**2) ## Divided by total volume, since each dV is 1, just equal to number of cells

  plt.plot(t_array[ti_index:tf_index+1],F[ti_index:tf_index+1])
  plt.title('Free Energy Functional', weight = "bold")
  plt.ylabel('Energy', weight = "bold")
  plt.xlabel('Time', weight = "bold")
  plt.savefig('Free_Energy.png')
  return None


plot_conc_evol(animation_fps = 10, ti = 0, tf = 27)
plt.clf() ## Clear previous figure
plot_free_energy(ti = 0, tf = 24)