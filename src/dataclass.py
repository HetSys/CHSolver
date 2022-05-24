import json
import os
import numpy as np
import h5py

class CHData():
    '''!
    Class to hold data for file I/O, binding to F90 solvers, and visualisation
    '''
    fname = None
    input_data = {}
    _data = {}


    def __init__(self, fname="input-data.json"):
        '''!
        Loads data from JSON file given by fname
        Will not save old files before opening a new one
        '''
        
        self.L = 1.0
        self.A = 1.0
        self.M = 0.25
        self.K = 0.0004
        self.p0 = -1.0
        self.p1 = 1.0
        self.T = np.linspace(0, 0.1, 3)


        self.grid_type = "r"
        self.grid_level = 7

        self.fname = fname

        self.cwd = os.getcwd()
        self.filepath = self.cwd + os.sep + fname

        if os.path.isfile(self.filepath):
            self._read_jsonfile()
        else:
            self.save_rundata("default")
            



    def _read_jsonfile(self):
        try:
            self._data = json.loads(open(self.filepath).read())
        except json.decoder.JSONDecodeError:
            raise RuntimeError(f"Error in reading JSON file {self.fname}, likely due to reading a blank file")


    @property
    def run_names(self) -> list:
        '''!
        Get names of all available runs
        '''
        keys = list(self._data.keys())
        return keys

    def read_rundata(self, run_name: str) -> None:
        '''!
        Search for a run with name matching run_name
        Throws KeyError if name not found

        '''
        try:
            self.input_data = self._data[run_name]
        except KeyError:
            print(f"Run {run_name} not found.")
            print("Currently recognised runs:")
            print(self.run_names)
            raise  # Raise the caught KeyError again

        self.L = self.input_data["L"]
        self.A = self.input_data["A"]
        self.M = self.input_data["M"]
        self.K = self.input_data["K"]
        self.p0 = self.input_data["p0"]
        self.p1 = self.input_data["p1"]
        self.grid_level = self.input_data["grid_level"]
        self.grid_type = self.input_data["grid_type"]

        if "T" in self.input_data.keys():
            self.T = self.input_data["T"]

    def print_rundata(self):
        print(f"L = {self.L}")
        print(f"A = {self.A}")
        print(f"M = {self.M}")
        print(f"K = {self.K}")
        print(f"p0 = {self.p0}")
        print(f"p1 = {self.p1}")
        print(f"Grid type = {self.grid_type}")
        print(f"Grid level = {self.grid_level}")
        print(f"Output Times = {self.T}")

    def save_rundata(self, run_name: str) -> None:
        '''!
        Save current input data as a new run in JSON file.
        Will overwrite if run with the same name exists
        '''
        self.input_data["L"] = self.L
        self.input_data["A"] = self.A
        self.input_data["M"] = self.M
        self.input_data["K"] = self.K
        self.input_data["p0"] = self.p0
        self.input_data["p1"] = self.p1
        self.input_data["grid_level"] = self.grid_level
        self.input_data["grid_type"] = self.grid_type
        self.input_data["T"] = list(self.T)

        if os.path.isfile(self.filepath):
            # Update internal _data with any changes to the file
            self._read_jsonfile()

        self._data[run_name] = self.input_data
        self._save_jsonfile()
    
    def _save_jsonfile(self):
        with open(self.filepath, 'w') as f:
                json.dump(self._data, f, indent=2)

    def solve(self, cmd_args=""):
        exe = "./chsolver"
        CH_cmds = f" -l {self.L} -a {self.A} -m {self.M} -k {self.K} -0 {self.p0} -1 {self.p1}"
        T_cmds = " -t {" + ":".join([str(t) for t in self.T]) + "}"
        other_cmds = f" -L {self.grid_level} -i {self.grid_type}"

        full_cmd = exe + CH_cmds + T_cmds + other_cmds + " " + cmd_args
        error_code = os.system(full_cmd)

        if error_code != 0:
            # Error has occurred
            print(full_cmd)
            raise FortranSourceError("Error occurred in solving backend")

        self.read_outputs("out")

    def read_outputs(self, outdir):
        '''!
        Reads output metadata and checkpoint files from outdir
        '''

        grid_params, sys_params, checkpoint_times = _read_metadata(os.getcwd() + os.sep + outdir + os.sep + "metadata.dat")
        
        self.grid_dimensionality = grid_params[0]
        self.grid_level = grid_params[1]

        self.L, self.A, self.M, self.K, self.p0, self.p1 = sys_params

        self.T = checkpoint_times

        grid_res = 2**self.grid_level

        c, c_prev, dt = _read_hdf5_files(self.T.shape[0], grid_res, outdir)

        self.C = c


def _read_metadata(filename):
  '''! Reads the grid parameters, system parameters and checkpoint times from the metadata.dat file in the out directory.
  '''
  with open(filename) as f:
    grid_params = np.array(f.readline().split()[1:], dtype=int)
    sys_params = f.readline().split()[1:]
    sys_params = np.array(sys_params + f.readline().split()[1:], dtype=float)
  chkpnt_times = np.genfromtxt(filename, skip_header=4, dtype=float)[:, 1]



  return grid_params, sys_params, chkpnt_times

def _read_hdf5_files(num_checkpoints, grid_res, outdir):
  '''! Reads concentration at current timestep (c), concentration at previous timestep (c_prev) and the corresponding timestep (dt) from the collection of HDF5 checkpoint files in the out/ directory.
  '''
  c = np.zeros((num_checkpoints,grid_res,grid_res))
  c_prev = np.zeros((num_checkpoints,grid_res,grid_res))
  dt = np.zeros((num_checkpoints))
  for i in range(1,num_checkpoints):
      data = h5py.File(outdir + os.sep + str(i) + '.chkpnt', 'r')
      test= data['c'][...]
      c[i,:,:] = test
      dt[i] = data['dt'][...]
  return c, c_prev, dt

class FortranSourceError(BaseException):
    pass

if __name__ == "__main__":
    CHData() # Generate input-data.json if it does not exist