from .DataIO import Json_handler, _read_hdf5_files, _read_metadata
import os


class Data_class():
    '''!
    Class to hold data for file I/O, binding to F90 solvers, and visualisation
    '''
    _jh = Json_handler()

    def __init__(self):
        '''!
        Define attribute namespaces

        '''

    def open_jsonfile(self, fname):
        '''!
        Call the json handler to open a file given by fname
        '''
        self._jh._load_data(fname)
    
    def get_runnames(self):
        return self._jh.run_names

    def get_rundata(self, run_name: str) -> None:
        '''!
        Retrieve input run data from JSON file
        '''
        self.input_data = self._jh.get_rundata(run_name)

        self.L = self.input_data["L"]
        self.A = self.input_data["A"]
        self.M = self.input_data["M"]
        self.K = self.input_data["K"]
        self.p0 = self.input_data["p0"]
        self.p1 = self.input_data["p1"]
        self.grid_level = self.input_data["grid_level"]
        self.grid_init = self.input_data["grid_init"]

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
        self.input_data["grid_init"] = self.grid_init

        self._jh.set_rundata(run_name, self.input_data)
    
    def save_jsonfile(self):
        self._jh.save_data()

    def read_outputs(self, outdir):
        '''!
        Reads output metadata and checkpoint files from outdir
        '''

        grid_params, sys_params, checkpoint_times = _read_metadata(os.get_cwd() + os.sep + outdir + os.sep + "metadata.dat")
        
        self.grid_dimensionality = grid_params[0]
        self.grid_level = grid_params[1]

        self.L, self.A, self.M, self.K, self.p0, self.p1 = sys_params

        self.T = checkpoint_times

        grid_res = 2**self.grid_level

        c, c_prev, dt = _read_hdf5_files(self.T.shape[0], grid_res, outdir)

        self.C = c

