"""!
Opening file docstring contains LATEX parsing - \f$\beta\f$

"""

# def foo():
#     '''!
#     random docs
#     \f$\alpha\f$

import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np
from .dataclass import CHData

def find_nearest_t_index(t_array:np.array, t):
    '''! Produces an animation for the evolution of species concentration, outputted as an mp4 to the main directory.
    @param animation_fps Number of frames per second for the animation
    '''
    t_array = np.asarray(t_array)
    index = (np.abs(t_array - t)).argmin()
    return index 

def plot_conc_evol(data_obj:CHData, animation_fps = 10 , ti = 0, tf=-1):
  '''! Produces an animation for the evolution of species concentration, outputted as an mp4 to the main directory.
  @param animation_fps Number of frames per second for the animation
  @param ti Time to intialise the visualiation (will use time closest to that specified)
  @param tf Time to finish the visualiation (will use time closest to that specified)
  '''

  t_array = data_obj.T
  grid_res = 2 ** data_obj.grid_level
  c = data_obj.C

  if tf==-1:
    tf = t_array[-1]

  ti_index = find_nearest_t_index(t_array, ti)
  tf_index = find_nearest_t_index(t_array, tf)
  ims = []
  fig, ax = plt.subplots()
  #plt.axis("off")
  for k in range(ti_index, tf_index+1):
      im = []
      ax.set_xticks([])
      ax.set_yticks([])
      shw = ax.imshow(c[k,:,:], interpolation="bicubic")
      im.append(shw)
      im.append(plt.title('Concentration Evolution', weight = "bold", fontsize = 14))
      im.append(plt.text(0.5, -0.05 , s = ("t = " + str(f"{t_array[k]:.3f}")), weight = "bold", fontsize = 12, ha = 'center', transform=ax.transAxes))
      ims.append(im)
  fig.colorbar(shw)
  ax.set_axis_off()
  fig.add_axes(ax)
  final_an = anim.ArtistAnimation(fig, ims, interval = 5, repeat = True, blit=False)
  final_an.save('Conc_Evolution.gif', writer = anim.PillowWriter(fps = animation_fps))

def plot_free_energy(data_obj:CHData, ti=0, tf=-1):
  '''! Produces a plot of Free energy versus time, outputted as a png to the main directory.
  @param ti Time to intialise the visualiation (will use time closest to that specified)
  @param tf Time to finish the visualiation (will use time closest to that specified)
  '''
  t_array = data_obj.T
  grid_res = 2 ** data_obj.grid_level
  c = data_obj.C

  L = data_obj.L
  A = data_obj.A
  M = data_obj.M
  K = data_obj.K
  p0 = data_obj.p0
  p1 = data_obj.p1

  num_checkpoints = t_array.shape[0]

  
  if tf==-1:
    tf = t_array[-1]

  ti_index = find_nearest_t_index(t_array, ti)
  tf_index = find_nearest_t_index(t_array, tf)
  f_c = np.zeros_like(c) #Bulk free energy density 
  for t in range(num_checkpoints):
    for i in range(grid_res):
      for j in range(grid_res):
          f_c[t,i,j] = A*((c[t,i,j]-p0)**2)*(c[t,i,j]-p1)**2 

  F = np.zeros((num_checkpoints)) #Free energy functional
  ## Setting up an alternate c array to account for PBCs
  c_halo = np.zeros((num_checkpoints, grid_res+2, grid_res+2))
  print(c_halo.shape)
  c_halo[:, 1:-1, 1:-1] = c
  print(c_halo.shape)
  c_halo[:,0,:] = c_halo[:,-2,:]
  print(c_halo.shape)
  c_halo[:,-1,:] = c_halo[:,1,:]
  print(c_halo.shape)
  c_halo[:,:,0] = c_halo[:,:,-2]
  print(c_halo.shape)
  c_halo[:,:,-1] = c_halo[:,:,1]
  print(c_halo.shape)
  for t in range(num_checkpoints):
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


#plot_conc_evol(animation_fps = 10, ti = 0, tf = 27)
#plt.clf() ## Clear previous figure
#plot_free_energy(ti = 0, tf = 24)
