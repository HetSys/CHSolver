import json
import glob
import os

class Json_handler():
  '''!
  Singleton handler for JSON I/O

  '''
  _instance = None # Start with no instances of this class
  _input_fname = "input-data.json" # fname of JSON file

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
    self.filepath = glob.glob("".join([self.cwd, os.sep, "**", os.sep, self._input_fname]), recursive = True)[0] # Recursive search to find file
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

  def get_rundata(self, run_name:str) -> dict:
    try:
      return self._data[run_name]
    except KeyError:
      print(f"Run {run_name} not found.")
      print(f"Currently recognised runs:")
      [print(key) for key in self.get_run_names()]
      return {}

  def set_rundata(self, run_name:str, rundata:dict):
    self._data[run_name] = rundata


