from DataIO import Json_handler


class Data_class():
    '''!
    Class to hold data for file I/O, binding to F90 solvers, and visualisation
    '''
    _jh = Json_handler()

    def __init__(self):
        '''!
        Define attribute namespaces

        '''

        self.input_data = self._jh.get_rundata("default")

    def get_rundata(self, run_name: str) -> None:
        '''!
        Retrieve input run data from JSON file
        '''
        self.input_data = self._jh.get_rundata(run_name)

    def save_rundata(self, run_name: str) -> None:
        '''!
        Save current input data as a new run in JSON file.
        Will overwrite if run with the same name exists
        '''
        self._jh.set_rundata(run_name, self.input_data)
