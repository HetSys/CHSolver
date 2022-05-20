import json
import glob
import os
import numpy as np
import h5py

class Json_handler():
    '''!
    Singleton handler for JSON I/O

    '''
    _instance = None  # Start with no instances of this class
    _input_fname = "input-data.json"  # fname of JSON file

    def __new__(cls):
        if (cls._instance is None):
            # Instance does not exist already, create one
            cls._instance = super().__new__(cls)
            cls._instance._load_data()
        return cls._instance

    def _load_data(self) -> None:
        '''!
        Loads data from self._input_fname JSON file
        '''
        self.cwd = os.getcwd()
        # Recursive search to find file
        self.filepath = glob.glob((self.cwd + os.sep + "**" + os.sep +
                                   self._input_fname), recursive=True)[0]
        self._data = json.loads(open(self.filepath).read())

    def save_data(self) -> None:
        '''!
        Save local run data back to json file
        '''
        with open(self.filepath, 'w') as f:
            json.dump(self._data, f)

    @property
    def run_names(self) -> list:
        '''!
        Get names of all available runs
        '''
        return list(self._data.keys())

    def get_rundata(self, run_name: str) -> dict:
        '''!
        Search for a run with name matching run_name
        Throws KeyError if name not found

        '''
        try:
            return self._data[run_name]
        except KeyError:
            print(f"Run {run_name} not found.")
            print("Currently recognised runs:")
            print(self.get_run_names())
            raise  # Raise the caught KeyError again

    def set_rundata(self, run_name: str, rundata: dict):
        '''!
        Write new/overwrite data associated with run_name

        '''
        self._data[run_name] = rundata
        self.save_data()


def generate_json_file():
    '''!
    Generates a fresh json file with sample input configs

    '''
    cwd = os.getcwd()
    input_fname = "input-data.json"
    # Recursive search to find file
    filepath = glob.glob((cwd + os.sep + "**" + os.sep +
                          input_fname), recursive=True)
    if filepath:
        data = json.loads(open(filepath[0]).read())
    else:
        data = {}

    data["default"] = {
        "grid_type": "r",
        "grid_level": 7,
        "A": 1.0,
        "L": 1.0,
        "M": 0.25,
        "K": 0.0004,
        "p0": -1.0,
        "p1": 1.0,
        "T" : tuple(np.linspace(0, 0.1, 3))
    }

    with open(Json_handler._input_fname, 'w') as f:
        json.dump(data, f, indent=2)

def _read_metadata(filename):
  '''! Reads the grid parameters, system parameters and checkpoint times from the metadata.dat file in the out directory.
  '''
  with open(filename) as f:
    grid_params = np.array(f.readline().split()[1:], dtype=int)
    sys_params = f.readline().split()[1:]
    sys_params = np.array(sys_params + f.readline().split()[1:], dtype=float)
  chkpnt_times = np.genfromtxt(filename, skip_header=4, dtype=float)[:, 1]



  return grid_params, sys_params, chkpnt_times

def _read_hdf5_files(num_checkpoints, grid_res):
  '''! Reads concentration at current timestep (c), concentration at previous timestep (c_prev) and the corresponding timestep (dt) from the collection of HDF5 checkpoint files in the out/ directory.
  '''
  c = np.zeros((num_checkpoints,grid_res,grid_res))
  c_prev = np.zeros((num_checkpoints,grid_res,grid_res))
  dt = np.zeros((num_checkpoints))
  for i in range(1,num_checkpoints):
      data = h5py.File('out/'+str(i)+'.chkpnt', 'r')
      test= data['c'][...]
      c[i,:,:] = test
      dt[i] = data['dt'][...]
  return c, c_prev, dt
