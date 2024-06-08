from abc import ABC, abstractmethod

class Transducers(ABC):

    @abstractmethod
    def get_sensitivity(self):
        pass

    @abstractmethod
    def get_pressure_resp(self):
        pass

class ElectroDynamic(Transducers):
    def __init__(
            self,
            Sd = None,
            Mms = None,
            Rec = None,
            Lec = None,
            Qts = None,
            Qes = None,
            Qms = None,
            Cms = None,
            Rms = None,
            fs = None,
            Bl = None,
    ):
        self.Sd = Sd
        self.Mms = Mms 
        self.Rec = Rec 
        self.Lec = Lec 
        self.Qts = Qts 
        self.Qes = Qes 
        self.Qms = Qms 
        self.Cms = Cms 
        self.Rms = Rms 
        self.fs = fs 
        self.Bl = Bl 

    def imp_to_ts(self, f_z, z, print_ts, ):
        pass

    def ts_to_imp(self, f_min, f_max, ):
        pass

    def print_ts(self):
        print(79*'-')
        print(f'Sd = {self.Sd}')
        print(f'Mms = {self.Mms}')
        print(f'Rec = {self.Rec}')
        print(f'Lec = {self.Lec}')
        print(f'Qts = {self.Qts}')
        print(f'Qes = {self.Qes}')
        print(f'Qms = {self.Qms}')
        print(f'Cms = {self.Cms}')
        print(f'Rms = {self.Rms}')
        print(f'fs = {self.fs}')
        print(f'Bl = {self.Bl}')
        print(79*'-')