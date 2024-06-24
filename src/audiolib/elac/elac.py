import audiolib.plotting as al_plt
import audiolib.tools as al_tls
import matplotlib.pyplot as plt
import numpy as np
import warnings

from matplotlib.backend_bases import MouseButton
from abc import ABC, abstractmethod

# Frequency region in which to derive Lec from
LEC_CALC_FREQ_RANGE = [500, 2000]

class Transducers(ABC):
    def __init__(self, c, rho, ):
        self.c = c
        self.rho = rho

    @abstractmethod
    def get_sensitivity(self):
        pass

    @abstractmethod
    def get_pressure_resp(self):
        pass

    @property
    def c(self):
        return self._c
    
    @c.setter
    def c(self, c):
        self._c = c

    @property
    def rho(self):
        return self._rho
    
    @rho.setter
    def rho(self, rho):
        self._rho = rho

class ElectroDynamic(Transducers):
    """
    Class to characterize electrodynamic transducers. Complete Thiele-Small 
    parameters are automatically calculated when added-mass-measurement is
    given. Alternatively, Thiele-Small parameters can be passed during 
    initilization in order to calculate modeled impedance curve from parameters.
    
    TODO: Update Parameters with added mass inputs
    Parameters
    ----------
    name : str
        Name of the transducer/model
    c : float
        Speed of sound during impedance measurement [m/s]
    rho : float
        Air density during impedance measurement [kg/m³]
    f_z : list or np.array, obsolete if all TS-parameters are given
        frequency vector of z, on infinite Baffle.
        If given with z, will re-calculate TS-params and overwrite any input
        TS-params of this object.
    z : list or np.array, obsolete if all TS-parameters are given
        Electrical input impedance, on infinite Baffle.
        If given with f_z, will re-calculate TS-params and overwrite any input
        TS-params of this object.
    f_z_added_mass : list or np.array, optional
        frequency vector of z of added mass driver, on infinite Baffle.
        Will not calculate Mms if not given.
    z_added_mass : list or np.array, optional
        Electrical input impedance of added mass driver, on infinite Baffle.
        Will not calculate Mms if not given.
    added_mass : float, optional
        Mass in kg
    Lec_estimation_range : array of len 2, obsolete if all TS-parameters given
        Frequency range in which to estimate Lec
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
            name = None,
            c = None,
            rho = None,
            f_z = None,
            z_abs = None,
            z_rad = None,
            f_z_added_mass = None,
            z_added_mass = None,
            added_mass = None,
            Lec_estimation_range = None,
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
        self.name = name
        self.Lec_estimation_range = Lec_estimation_range
        self.Sd = Sd
        self._c = c
        self._rho = rho
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
        imp_input_given = (
            (f_z is not None) and (z_abs is not None) and (z_rad is not None)
        )
        added_mass_imp_input_given = not None in [
            f_z_added_mass,
            z_added_mass,
            added_mass,
        ]
        if imp_input_given:
            # TODO: Perform format parsing on z: imag/real, abs/rad?
            self.f_z = f_z
            self.z_abs = z_abs
            self.z_rad = z_rad
            if (len(not_none_param_idcs) > 0):
                warnings.warn(
                    'Impedance curve and TS-Params given: Will overwrite given ' 
                    + 'TS-params by inherent TS-parameter calculation via |Z|.'
                )
            if added_mass_imp_input_given:
                self.f_z_added_mass = f_z_added_mass
                self.z_abs_added_mass = z_added_mass
                self.added_mass = added_mass
            else:
                warnings.warn(
                    'Mms not calculated: missing added mass frequency vector ' 
                    + 'and/or added mass impedance vector and/or added weight.'
                )
            self.imp_to_ts()
        if not imp_input_given and (len(not_none_param_idcs) > 0):
            self.Mms = Mms
            self.Rec = Rec
            self.Lec = Lec
            self.Qes = Qes
            self.Qms = Qms
            self.fs = fs
            self.Bl = Bl
            self._update_dependent_ts_params(plot_params=False)

    def imp_to_ts(self, added_mass_available=True, plot_params=False):
        self.Rec = self.z_abs[0]
        self.fs, self._z_max = self._manual_pick_fs(self.f_z, self.z_abs, 'Free-Air')
        if added_mass_available:
            self.Mms = self._calc_mms_via_added_mass()
        idx_fs = al_tls.closest_idx_to_val(arr=self.f_z, val=self.fs)
        self._r0 = self._z_max / self.Rec
        Z_at_f1_f2 = np.sqrt(self._r0)*self.Rec
        idx_f1 = al_tls.closest_idx_to_val(arr=self.z_abs[:idx_fs], val=Z_at_f1_f2)
        self._f1  = self.f_z[idx_f1]
        # Limit f2 search frequency range to [fs:(2*fs)] to avoid Zmax@ high f:
        idx_limit_high_freq_f2 = int(2*idx_fs)
        idx_f2 = idx_fs + al_tls.closest_idx_to_val(
            arr = self.z_abs[idx_fs:idx_limit_high_freq_f2],
            val = Z_at_f1_f2,
        )
        self._f2 = self.f_z[idx_f2]
        self.Lec = self._estimate_Lec()
        self._update_dependent_ts_params()
        if plot_params:
            self.plot_z_params()

    def ts_to_imp(self, freq_range=None, freq_resolution=None, ):
        """
        Return complex-valued impedance curve, calculated from TS-parameters

        Parameters
        ----------
        freq_range : list or array [low_end, high_end]
            Frequency range of returned impedance curve, optional.
            Not necessary if f_z, z_abs and z_rad is given to object beforehand
        freq_resolution :
            Frequency resolution of returned impedance curve, optional.
            Not necessary if f_z , z_abs and z_rad is given to object beforehand
        
        Returns
        ------- 
        f_es : np.array
            Frequency vector of z
        z_es : np.array, complex
            Modeled electrical input impedance, derived from TS-parameters
        """
        if self.f_z is not None:
            f_es = np.array(self.f_z)
            omega_es = f_es*2*np.pi
        else:
            f_es = np.arange(freq_range[0], freq_range[1], freq_resolution)
            omega_es = 2*np.pi*f_es
        z_ms = self.Rms + 1j*omega_es*self.Mms + (1 / (1j*omega_es*self.Cms))
        z_es = self.Rec + 1j*omega_es*self.Lec + (self.Bl**2 / z_ms)
        return f_es, z_es,

    def _manual_pick_fs(self, f_z, z, z_explanation, ):
        fig, ax = al_plt.plot_rfft_freq(f_z, z, xscale='log', )
        ax.set_title(f'{z_explanation}:\nManually hover over f$_s$ and Z$_{{max}}$ '+ 
                     fr'and select with "Space"-Button.')
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

    def _estimate_Lec(self):
        z_complex = self.z_abs*np.exp(1j*np.array(self.z_rad))
        omega = 2*np.pi*np.array(self.f_z)
        low_idx = al_tls.closest_idx_to_val(
            arr=self.f_z, val=self.Lec_estimation_range[0]
        )
        high_idx = al_tls.closest_idx_to_val(
            arr=self.f_z, val=self.Lec_estimation_range[1]
        )
        Lec = np.mean(
            z_complex.imag[low_idx:high_idx] / omega[low_idx:high_idx]
        )
        return Lec

    def _calc_mms_via_added_mass(self, ):
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
        fs_new, _ = self._manual_pick_fs(
            self.f_z_added_mass,
            self.z_abs_added_mass,
            'Added Mass',
        )
        Mms = self.added_mass / ((self.fs / fs_new)**2 - 1)
        return Mms

    def get_pressure_resp(self):
        pass

    def get_sensitivity(self):
        pass

    def plot_z_params(self, ):
        """
        Plot the parameters, from which the TS-parameters are calculated 
        (r0, f1, f2, fs, Rec) on top of the input impedance plot. This is to
        verify the integrity of the calculation.

        No input parameters, no returns, just opens a plot.
        """
        v_Rec = self.Rec*np.ones(len(self.f_z))
        v_z_f1_f2 = np.sqrt(self._r0)*self.Rec*np.ones(len(self.f_z))
        fig, ax = al_plt.plot_rfft_freq(
            self.f_z,
            self.z_abs,
            xscale = 'log',
            yscale='lin',
        )
        ax.axvline(x=self._f1, ymin=0, ymax=100, linestyle='--', color='r', label=r'f$_1$')
        ax.axvline(x=self._f2, ymin=0, ymax=100, linestyle='--', color='cyan', label=r'f$_2$')
        ax.axvline(x=self.fs, ymin=0, ymax=100, linestyle='--', color='k', label=r'f$_s$')
        ax.plot(self.f_z, v_z_f1_f2, linestyle='--', label=r'$\sqrt{r_0} R_{ec}$')
        ax.plot(self.f_z, v_Rec, linestyle='--', label=r'$R_{ec}$')
        ax.axvspan(
            xmin=self.Lec_estimation_range[0],
            xmax=self.Lec_estimation_range[1],
            alpha=.5,
            color='b',
            label='$L_{{ec}}$ est. range'
        )
        ax.set_ylabel(r'|Z| [$\Omega$]')
        ax.legend()
        plt.show(block=False)
        return fig, ax

    def _update_dependent_ts_params(self):
        """
        Update dependent variables when their parameters are changed/set.
        Dependent variables are e.g. Qts = Qms*Qes / (Qms+Qes), since it depends
        on other TS-parameters, which have to be calculated beforehand.
        TODO: Implement user feedback: input Qts (e.g. datasheet) vs.
          calculated Qts fromt his method
        """
        self.Qms = self.fs*np.sqrt(self._r0) / (self._f2 - self._f1)
        self.Qes = self.Qms / (self._r0 - 1)
        self.Qts = self.Qms * self.Qes / (self.Qms + self.Qes)
        self.Cms = 1 / ((2*np.pi*self.fs)**2 * self.Mms)
        self.Rms = 1 / (2*np.pi*self.fs*self.Qms*self.Cms)
        self.Bl = np.sqrt(self.Rms*(self._z_max - self.Rec))
        self.Vas = (self.Sd*self.c)**2*self.rho/((2*np.pi*self.fs)**2*self.Mms)

    def print_ts(self):
        # TODO: Implement standardized unit-SI-prefix-conversion in tools.
        print('\n')
        print(79*'-')
        print(f'{self.name}: Thiele-Small-Parameters')
        print(f'  fs = {np.round(self.fs, 2)} Hz')
        print(f'  Sd = {self.Sd} m²')
        print(f'  Vas = {np.round(self.Vas*1000, 3)} L')
        print(f'  Rec = {np.round(self.Rec, 2)} Ohm')
        print(f'  Lec = {np.round(self.Lec*1000, 2)} mH') 
        print(f'  Qts = {np.round(self.Qts, 3)}')
        print(f'  Qes = {np.round(self.Qes, 3)}')
        print(f'  Qms = {np.round(self.Qms, 3)}')
        print(f'  Mms = {np.round(self.Mms*1000, 2)} g')
        print(f'  Cms = {np.round(self.Cms*1e6, 2)} µm/N')
        print(f'  Rms = {np.round(self.Rms, 2)} kg/s')
        print(f'  Bl = {np.round(self.Bl, 2)} Tm')
        print(79*'-')
        print('\n')

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, c):
        self._c = c
        self._update_dependent_ts_params()

    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, rho):
        self._rho = rho
        self._update_dependent_ts_params()