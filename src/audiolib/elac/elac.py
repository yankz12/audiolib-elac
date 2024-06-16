import audiolib.plotting as al_plt
import audiolib.tools as al_tls
import matplotlib.pyplot as plt
import numpy as np
import warnings

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
            f_z = None,
            z = None,
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
        param_list = [
            Mms,
            Rec,
            Lec,
            Qts,
            Qes,
            Qms,
            Cms,
            Rms,
            fs ,
            Bl ,
        ]
        non_none_param_idcs = [
            i for i in range(len(param_list)) if param_list[i] != None
        ] # Get idcs of input params which are not None
        imp_input_given = (f_z is not None) and (z is not None)
        if imp_input_given and (len(non_none_param_idcs) == 0):
            self.f_z = f_z
            self.z = z
        if imp_input_given and (len(non_none_param_idcs) > 0):
            self.f_z = f_z
            self.z = z
            warnings.warn(
                'Impedance curve and TS-Params given: Will overwrite given ' 
                + 'TS-params by inherent TS-parameter calculation via |Z|.'
            )
        if not imp_input_given and (len(non_none_param_idcs) > 0):
            self.Mms = Mms
            self.Rec = Rec
            self.Lec = Lec
            self.Qes = Qes
            self.Qms = Qms
            self.fs = fs
            self.Bl = Bl
            self._update_dependent_params()

    def imp_to_ts(self, plot_params=True):
        # TODO: Add Mms calculation from added mass method
        # TODO: Add Lec calculation from imaginary part average divided by omega
        self.Rec = self.z[0]
        self.fs, self._z_max = self._manual_pick_fs(self.f_z, self.z, )
        idx_fs = al_tls.closest_idx_to_val(arr=self.f_z, val=self.fs)
        self._r0 = self._z_max / self.Rec
        Z_at_f1_f2 = np.sqrt(self._r0)*self.Rec
        idx_f1 = al_tls.closest_idx_to_val(arr=self.z[:idx_fs], val=Z_at_f1_f2)
        self._f1  = self.f_z[idx_f1]
        # Limit f2 search frequency range to [fs:(2*fs)] to avoid Zmax@Lec:
        idx_limit_high_freq_f2 = int(2*idx_fs)
        print(f'high limit: {self.f_z[idx_limit_high_freq_f2]} Hz')
        idx_f2 = idx_fs + al_tls.closest_idx_to_val(
            arr = self.z[idx_fs:idx_limit_high_freq_f2],
            val = Z_at_f1_f2,
        )
        self._f2 = self.f_z[idx_f2]
        self.Qms = self.fs*np.sqrt(self._r0) / (self._f2 - self._f1)
        self.Qes = self.Qms / (self._r0 - 1)
        self.Qts = self.Qms*self.Qes / (self.Qms + self.Qes)

        if plot_params:
            self.plot_z_params()

    def ts_to_imp(self, ):
        # TODO: Define proper frequency range if f_z not given
        omega = self.f_z*2*np.pi
        z_ms = self.Rms + 1j*omega*self.Mms + (1 / (1j*omega*self.Cms))
        z_ls = self.Rec + 1j*omega*self.Lec + (self.Bl**2 / z_ms)

    def _manual_pick_fs(self, f_z, z):
        fig, ax = al_plt.plot_rfft_freq(f_z, z, xscale='log', )
        ax.set_title(r'Manually hover over f$_s$ and Z$_{max}$ and select ' + 
                     r'with "Space"-Button.')
        ax.set_ylabel(r'|Z| [$\Omega$]')
        fs_selection = plt.ginput(
            n=1,
            timeout=0,
            show_clicks='true',
            mouse_add=None,
            mouse_pop=None,
            mouse_stop=MouseButton.RIGHT,
        )
        plt.close(fig)
        fs = fs_selection[0][0]
        zmax = fs_selection[0][1]
        return fs, zmax

    def get_pressure_resp(self):
        pass

    def get_sensitivity(self):
        pass

    def plot_z_params(self, ):
        v_Rec = self.Rec*np.ones(len(self.f_z))
        v_z_f1_f2 = np.sqrt(self._r0)*self.Rec*np.ones(len(self.f_z))
        _, ax = al_plt.plot_rfft_freq(
            self.f_z,
            self.z,
            xscale = 'log',
            yscale='lin',
        )
        ax.axvline(x=self._f1, ymin=0, ymax=100, linestyle='--', color='r', label=r'f$_1$')
        ax.axvline(x=self._f2, ymin=0, ymax=100, linestyle='--', color='cyan', label=r'f$_2$')
        ax.axvline(x=self.fs, ymin=0, ymax=100, linestyle='--', color='k', label=r'f$_s$')
        ax.plot(self.f_z, v_z_f1_f2, linestyle='--', label=r'$\sqrt{r_0} R_{ec}$')
        ax.plot(self.f_z, v_Rec, linestyle='--', label=r'$R_{ec}$')
        ax.set_ylabel(r'|Z| [$\Omega$]')
        ax.legend()
        plt.show(block=False)

    def _update_dependent_params(self):
        # Update dependent variables when their parameters are changed/set.
        # TODO: Define dependent and independent variables, like:
        # Dependent: Qts = Qms*Qes/(Qms+Qes), Cms = ..., Rms = 
        # Independent: Lec, Rec, Sd, etc.
        self.Qts = self.Qms * self.Qes / (self.Qms + self.Qes)
        self.Cms = 1 / ((2*np.pi*self.fs)**2 * self.Mms )
        self.Rms = 1 / (2*np.pi*self.fs*self.Qms*self.Cms)

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
