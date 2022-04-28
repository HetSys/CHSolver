import json
import glob
import os


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
        "grid_level": 2,
        "A": 1.0,
        "L": 1.0,
        "M": 1.0,
        "K":1.0,
        "p0": 0.0,
        "p1": 1.0
    }

    with open(Json_handler._input_fname, 'w') as f:
        json.dump(data, f, indent=2)
