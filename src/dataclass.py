from json_handler import Json_handler




      


class Data_class():
  '''!
  Class to hold data for file I/O, binding to F90 solvers, and visualisation
  '''

  def __init__(self):
    '''!
    Define attribute namespaces
    
    '''

    self.input_data = {
      "grid_type"  : "r",
      "grid_level" : 2,
      "A"          : 1,
      "L"          : 1,
      "M"          : 1,
      "p0"         : 0,
      "p1"         : 0 
    }