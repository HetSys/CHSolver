import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np

def find_previous_t_index(t_array:np.array, t):
    '''! Finds the index of the closest prior time to t within the t_array.
    @param t_array array of all outputted times.
    @param t desired start time
    '''
    t_array = np.asarray(t_array)
    index = (np.abs(t_array - t)).argmin()
    if t_array[index] > t:
      index += -1
    return index


def plot_conc_evol(data_obj, animation_fps = 10 , ti = 0, tf=-1, metadata = True, filetype = 'mp4'):
  '''! Produces an animation for the evolution of species concentration, outputted as an mp4 to the main directory.
  @param animation_fps Number of frames per second for the animation. Default is 10.
  @param ti Time to intialise the visualiation (will use closest time prior to that specified). Default is 0
  @param tf Time to finish the visualiation (will use closest time prior to that specified). Default is the final time recorded.
  @param metadata Show metadata at top of the visual. Default is True

  '''

  t_array = data_obj.T
  c = data_obj.C
  L = data_obj.L
  A = data_obj.A
  M = data_obj.M
  K = data_obj.K
  p0 = data_obj.p0
  p1 = data_obj.p1

  if tf==-1:
    tf = t_array[-1]

  ti_index = find_previous_t_index(t_array, ti)
  tf_index = find_previous_t_index(t_array, tf)
  ims = []
  fig, ax = plt.subplots()
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
  if metadata == True:
    plt.suptitle('L = '+str(L)+', A = '+str(A)+', M = '+str(M)+ ', K = '+str(K)+', p0 = '+str(p0)+', p1 = '+str(p1), fontsize = 8)
  final_an = anim.ArtistAnimation(fig, ims, interval = 5, repeat = True, blit=False)
  if filetype == 'gif':
    final_an.save('Conc_Evolution.gif', writer = anim.PillowWriter(fps = animation_fps))
  elif filetype == 'mp4':
    final_an.save('Conc_Evolution.mp4', writer = "ffmpeg", fps = animation_fps)
  else:
    print('Unable to use specified filetype, defaulting to mp4.')
    final_an.save('Conc_Evolution.mp4', writer = "ffmpeg", fps = animation_fps)

  return final_an

def plot_conc_snapshot(data_obj, t, metadata = True):
  '''! Produces a snapshort of species concentration at a specified time, outputted as an png to the main directory.
@param t Snapshot time to visualise (will use closest time prior to that specified). Default is 0
@param metadata Show metadata at top of the visual. Default is True
'''

  t_array = data_obj.T
  c = data_obj.C
  L = data_obj.L
  A = data_obj.A
  M = data_obj.M
  K = data_obj.K
  p0 = data_obj.p0
  p1 = data_obj.p1
  t_index = find_previous_t_index(t_array, t)

  fig, ax = plt.subplots()
  shw = ax.imshow(c[t_index,:,:], interpolation="bicubic")
  ax.set_title('Concentration at time '+ str(f"{t_array[t_index]:.3f}"), weight = "bold", fontsize = 14)
  fig.colorbar(shw)
  ax.axis('off')
  t_= t_array[t_index]
  t_dec = f"{(t_-int(t_)):.3f}"[2:]
  t_str = (str(int(t_))+'_'+t_dec)
  if metadata == True:
    plt.suptitle('L = '+str(L)+', A = '+str(A)+', M = '+str(M)+ ', K = '+str(K)+', p0 = '+str(p0)+', p1 = '+str(p1), fontsize = 8)
  plt.savefig('Conc_snapshot_t_'+t_str+'.png')
  return fig, ax

def plot_free_energy(data_obj, ti=0, tf=-1, metadata = True):
  '''! Produces a plot of Free energy (F) versus time, outputted as a png to the main directory.
  \f$F = \int_V[f(c) + \frac{1}{2}\kappa(\nabla c)^2]\ dV$, where $f(c) = A(c-p_0)^2(c-p_1)^2\f$
  @param ti Time to intialise the visualiation (will use closest time prior to that specified). Default is 0
  @param tf Time to finish the visualiation (will use closest time prior to that specified). Default is the final time recorded.
  @param metadata Show metadata at top of the visual. Default is True
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

  ti_index = find_previous_t_index(t_array, ti)
  tf_index = find_previous_t_index(t_array, tf)
  f_c = np.zeros_like(c) #Bulk free energy density
  for t in range(num_checkpoints):
    for i in range(grid_res):
      for j in range(grid_res):
          f_c[t,i,j] = A*((c[t,i,j]-p0)**2)*(c[t,i,j]-p1)**2

  F = np.zeros((num_checkpoints)) #Free energy functional
  ## Setting up an alternate c array to account for PBCs
  c_halo = np.zeros((num_checkpoints, grid_res+2, grid_res+2))
  c_halo[:, 1:-1, 1:-1] = c
  c_halo[:,0,:] = c_halo[:,-2,:]
  c_halo[:,-1,:] = c_halo[:,1,:]
  c_halo[:,:,0] = c_halo[:,:,-2]
  c_halo[:,:,-1] = c_halo[:,:,1]
  for t in range(num_checkpoints):
    grad_c_x = np.gradient(c_halo[t,:,:], axis=0)[1:-1,1:-1]
    grad_c_y = np.gradient(c_halo[t,:,:], axis=1)[1:-1,1:-1]
    grad_c_sqd = (grad_c_x ** 2 + grad_c_y ** 2)
    for i in range(grid_res):
      for j in range(grid_res):
        F[t] += f_c[t,i,j] + 0.5*K*grad_c_sqd[i,j]## dV = 1? so no term in here
  F = F/(grid_res**2) ## Divided by total volume, since each dV is 1, just equal to number of cells

  fig, ax = plt.subplots()
  ax.plot(t_array[ti_index:tf_index+1],F[ti_index:tf_index+1])
  ax.set_title('Free Energy Functional', weight = "bold")
  ax.set_ylabel('Energy', weight = "bold")
  ax.set_xlabel('Time', weight = "bold")
  if metadata == True:
    plt.suptitle('L = '+str(L)+', A = '+str(A)+', M = '+str(M)+ ', K = '+str(K)+', p0 = '+str(p0)+', p1 = '+str(p1), fontsize = 8)
  plt.savefig('Free_Energy.png')
  return fig, ax

if __name__ == "__main__":
  from dataclass import CHData
  dat = CHData()

  # Read outputs from default loc and plot visualisation
  dat.read_outputs("out")
  plot_conc_evol(dat, filetype = 'mp4')
  plt.clf()
  plot_free_energy(dat)
  plt.clf()
  # plot_conc_snapshot(dat,0.5)
