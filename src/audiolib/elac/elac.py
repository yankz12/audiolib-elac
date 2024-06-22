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
    """
    Class to characterize electrodynamic transducers. Complete Thiele-Small 
    parameters are automatically calculated when added-mass-measurement is
    given. Alternatively, Thiele-Small parameters can be passed during 
    initilization in order to calculate modeled impedance curve from parameters.
    
    TODO: Update Parameters with added mass inputs
    Parameters
    ----------
    f_z : list or np.array
        frequency vector of z, measured on infinite Baffle
        If given with z, will re-calculate TS-params and overwrite any input
        TS-params of this object.
    z : list or np.array
        Electrical input impedance, measured on infinite Baffle
        If given with f_z, will re-calculate TS-params and overwrite any input
        TS-params of this object.
    Sd : float
        Radiating surface of transducer
    Mms : float
        Moving mass of transducer
    Rec : float
        Electrical Resistance of coil
    Lec : float
        Electrical Inductance of coil
    Qes : float
        Electrical quality factor
    Qms : float
        Mechanical quality factor
    Qts : float
        Total quality factor, resulting from Qes and Qms
    Cms : float
        Compliance of mechanical suspension
    Rms : float
        Resistance of mechanical suspension
    fs : float
        Resonance frequency
    Bl : float
        Force factor

    TODO: Examples
    --------
    """
    def __init__(
            self,
            f_z = None,
            z = None,
            f_z_added_mass = None,
            z_added_mass = None,
            added_mass = None,
            Sd = None,
            fs = None,
            Rec = None,
            Lec = None,
            Qes = None,
            Qts = None,
            Mms = None,
            Cms = None,
            Rms = None,
            Qms = None,
            Bl = None,
    ):
        self.Sd = Sd
        param_list = [
            fs ,
            Rec,
            Lec,
            Qes,
            Mms,
            Cms,
            Rms,
            Qms,
            Qts,
            Bl ,
        ]
        not_none_param_idcs = [
            i for i in range(len(param_list)) if param_list[i] != None
        ] # Get idcs of input params which are not None
        imp_input_given = (f_z is not None) and (z is not None)
        added_mass_imp_input_given = (
            (f_z_added_mass is not None) and (z_added_mass is not None)
        )
        if imp_input_given:
            self.f_z = f_z
            self.z = z
            if (len(not_none_param_idcs) > 0) and added_mass_imp_input_given:
                warnings.warn(
                    'Impedance curve and TS-Params given: Will overwrite given ' 
                    + 'TS-params by inherent TS-parameter calculation via |Z|.'
                )
                self.imp_to_ts()
                # self.Mms = self._calc_mms_via_added_mass(
                #     f_z = 
                # )
        if not imp_input_given and (len(not_none_param_idcs) > 0):
            self.Mms = Mms
            self.Rec = Rec
            self.Lec = Lec
            self.Qes = Qes
            self.Qms = Qms
            self.fs = fs
            self.Bl = Bl
            self._update_dependent_params(plot_params=False)

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
        self._update_dependent_params()


        if plot_params:
            self.plot_z_params()

    def ts_to_imp(self, freq_range=None, freq_resolution=None, ):
        """
        Return complex-valued impedance curve, calculated from TS-parameters

        Parameters
        ----------
        freq_range : list or array [low_end, high_end]
            Frequency range of returned impedance curve, optional.
            Not necessary if f_z and z is given to object beforehand
        freq_resolution :
            Frequency resolution of returned impedance curve, optional.
            Not necessary if f_z and z is given to object beforehand
        
        Returns
        ------- 
        f_es : np.array
            Frequency vector of z
        z_es : np.array, complex
            Modeled electrical input impedance, derived from TS-parameters
        """
        if self.f_z is not None:
            f_es = self.f_z
            omega_es = f_es*2*np.pi
        else:
            f_es = np.arange(freq_range[0], freq_range[1], freq_resolution)
            omega_es = 2*np.pi*f_es
        z_ms = self.Rms + 1j*omega_es*self.Mms + (1 / (1j*omega_es*self.Cms))
        z_es = self.Rec + 1j*omega_es*self.Lec + (self.Bl**2 / z_ms)
        return f_es, z_es,

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

    def _calc_mms_via_added_mass(self, f_z, z added_mass, ):
        """
        Parameters
        ----------
        added_mass : float
            Weight of added mass in grams

        Returns
        -------
        Mms : float
            Moving mass, derived by Mms = added_mass / ((fs / fs_new)**2 - 1)
        """
        # TODO: Write this, modify imp_to_ts in order to re-use its fs calc
        fs_new, _ = self._manual_pick_fs(self.f_z, self.z, )
        Mms = added_mass / ((self.fs / fs_new)**2 - 1)
        pass

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
        # TODO: Implement user feedback: input Qts (e.g. datasheet) vs.
        #   calculated Qts fromt his method
        self.Qts = self.Qms * self.Qes / (self.Qms + self.Qes)
        self.Cms = 1 / ((2*np.pi*self.fs)**2 * self.Mms )
        self.Rms = 1 / (2*np.pi*self.fs*self.Qms*self.Cms)
        self.Qes = self.Qms / (self._r0 - 1)
        self.Bl = np.sqrt(self.Rms*(self._z_max - self.Rec))

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
