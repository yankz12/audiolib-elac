import audiolib.plotting as al_plt
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.backend_bases import MouseButton
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

    def imp_to_ts(self, f_z, z, print_ts=True, plot_params=True):
        # Error on fs is present with derivation due to noise
        z_prime = np.gradient(z, f_z, )
        # z has to points of slope 0: fs and when it gets into Lec slope
        # Pick first to obtain fs
        self.Rec = z[0]
        self.fs = self._manual_pick_fs(f_z, z, )
        print(self.fs)
        idx_fs = np.argwhere(f_z==self.fs)
        z_max = z[idx_fs]
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
        if plot_params:
            self.plot_z_params(f_z, z, f1, f2, r0)

    def _manual_pick_fs(self, f_z, z):
        fig, ax = al_plt.plot_rfft_freq(f_z, z, xscale='log', )
        ax.set_title('Manually hover over fs and fix with "Space"-Button.')
        fs_selection = plt.ginput(
            n=1,
            timeout=0,
            show_clicks='true',
            mouse_add=None,
            mouse_pop=None,
            mouse_stop=MouseButton.RIGHT,
        )
        plt.close(fig)
        return fs_selection[0][0]

    def get_pressure_resp(self):
        pass

    def get_sensitivity(self):
        pass

    def plot_z_params(self, f_z, z, f1, f2, r0):
        v_Rec = self.Rec*np.ones(len(f_z))
        v_z_f1_f2 = np.sqrt(r0)*self.Rec*np.ones(len(f_z))
        fig, ax = al_plt.plot_rfft_freq(
            f_z,
            z,
            xscale = 'log',
            yscale='lin',
        )
        ax.axvline(x=f1, ymin=0, ymax=100, linestyle='--')
        ax.axvline(x=f2, ymin=0, ymax=100, linestyle='--')
        ax.axvline(x=self.fs, ymin=0, ymax=100, linestyle='--')
        ax.plot(f_z, v_z_f1_f2, linestyle='--')
        ax.plot(f_z, v_Rec, linestyle='--')
        plt.show(block=False)
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
