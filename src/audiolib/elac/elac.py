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
            Sd,
            Mms,
            Rec,
            Lec,
            Qts,
            Qes,
            Qms,
            Cms,
            Rms,
            fs,
            Bl,
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

    def imp_to_ts(self, f_z, z, ):
        pass

    def ts_to_imp(self, f_min, f_max, ):
        pass
