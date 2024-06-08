import numpy as np

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

    def imp_to_ts(self, f_z, z, print_ts=True, ):
        z_prime = np.gradient(z, f_z, )
        # z has to points of slope 0: fs and when it gets into Lec slope
        # Pick first to obtain fs
        self.Rec = z[0]
        idx_fs = np.argmin(z_prime)
        z_max = z[idx_fs]
        self.fs = f_z[idx_fs] # Hz
        print(self.fs)
        r0 = z_max / self.Rec
        Z_at_f1_f2 = np.sqrt(r0)*self.Rec
        idx_f1 = np.argmin(np.abs(z[:idx_fs] - Z_at_f1_f2))
        f1 = f_z[idx_f1]
        # Limit f2 search frequency range to [fs:(2*fs-f1)] to avoid Zmax@Lec:
        idx_limit_high_freq_f2 = int(2*idx_fs - idx_f1)
        idx_f2 = np.argmin(
            np.abs(
                z[idx_fs:idx_limit_high_freq_f2] - Z_at_f1_f2
            )
        ) + idx_fs
        f2 = f_z[idx_f2]
        if print_ts:
            print(f'fs : {np.round(self.fs)}, f1: {np.round(f1)}, f2: {np.round(f2)}')


    def get_pressure_resp(self):
        pass

    def get_sensitivity(self):
        pass

    def plot_z_params():
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
