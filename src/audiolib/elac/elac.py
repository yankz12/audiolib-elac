import audiolib.plotting as al_plt
import audiolib.tools as al_tls
import matplotlib.pyplot as plt
import numpy as np
import warnings

from matplotlib.backend_bases import MouseButton
from abc import ABC, abstractmethod

# TODO: Implement more beautiful Unit- and prefix-conversion in audiolib.tools
STANDARD_UNITS = {
    'fs' : 'Hz',
    'Sd' : 'cm²',
    'Vas' : 'L',
    'Rec' : 'Ohm',
    'Lec' : 'mH',
    'Qts' : '',
    'Qes' : '',
    'Qms' : '',
    'Mms' : 'g',
    'Cms' : 'µm/N',
    'Rms' : 'kg/s',
    'Bl' : 'Tm',
}

STANDARD_PREFIX_CONV = {
    'fs' : 1,
    'Sd' : 1e4, # mm² to cm²
    'Vas' : 1000, # m³ to L
    'Rec' : 1,
    'Lec' : 1000, # H to mH
    'Qts' : 1,
    'Qes' : 1,
    'Qms' : 1,
    'Mms' : 1000, # kg to g
    'Cms' : 1e6, # m/N to µm/N
    'Rms' : 1,
    'Bl' : 1,
}

class Transducers(ABC):
    def __init__(self, c, rho, ):
        self.c = c
        self.rho = rho

    @abstractmethod
    def get_baffled_sens(self):
        pass

    @abstractmethod
    def get_baffled_p_resp(self):
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
    Class to characterize electrodynamic transducers. Able to derive TS-Params
    from measured electrical input impedance and able to calculate modeled input
    impedance from given TS-Params. Check README for more practical application
    explanation.
    
    Parameters
    ----------
    name : str, optional, defaults to empty string
        Name of the transducer/model
    c : float, optional
        Speed of sound during impedance measurement [m/s]
    rho : float,  optional
        Air density during impedance measurement [kg/m³]
    f_z : list or np.array, obsolete if all TS-parameters are given
        frequency vector of z, on infinite Baffle.
        If given with z_abs and z_rad, will re-calculate TS-params and overwrite
        any input TS-params of this instance.
    z_abs : list or np.array, obsolete if all TS-parameters are given
        Absolute values of Electrical input impedance, on infinite Baffle.
        If given with f_z and z_rad, will re-calculate TS-params and overwrite
        any input TS-params of this instance.
    z_rad : list or np.array, obsolete if all TS-parameters are given
        Phase of Electrical input impedance in radians, on infinite Baffle.
        If given with f_z and z_abs, will re-calculate TS-params and overwrite
        any input TS-params of this instance.
    added_element_method : str, optional
        Either 'mass' or 'volume', depending on the procedure used for TS
        parameter estimation
    f_z_added : list or np.array, optional
        frequency vector of z of added mass or volume driver, on infinite Baffle.
        Will not calculate Mms, Cms and Vas if not given
    z_added : list or np.array, optional
        Absolute values of Electrical input impedance of added mass or volume
        measurement, on infinite Baffle. Will not calculate Mms if not given.
    added_element_quantity : float, optional
        Quantity of Element that has been added for measurement:
        Mass [kg] or Volume [m³]
    Lec_estimation_range : array of len 2, obsolete if all TS-parameters given
        Frequency range in which to estimate Lec as a simple inductance.
        TODO: Implement semi-inductance model.
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

    Examples
    --------
    See test script in this repository under ../audiolib-elac/test/test_z.py
    Dataset is a measured impedance of a Dayton Audio LS10-44 woofer.
    """

    def __init__(
            self,
            name = '',
            c = None,
            rho = None,
            f_z = None,
            z_abs = None,
            z_rad = None,
            added_element_method = None,
            added_element_quantity = None,
            f_z_added = None,
            z_added = None,
            Lec_estimation_range = None,
            Sd = None,
            Vas = None,
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
        self._Lec_estimation_range = Lec_estimation_range
        self.Sd = Sd
        self.Vas = Vas
        self.fs = fs
        self.Rec = Rec
        self.Lec = Lec
        self.Qes = Qes
        self.Qts = Qts
        self.Mms = Mms
        self.Cms = Cms
        self.Rms = Rms
        self.Qms = Qms
        self.Bl = Bl
        self._c = c
        self._rho = rho

        imp_input_given = (
            (f_z is not None) and (z_abs is not None) and (z_rad is not None)
        )
        added_imp_input_given = not None in [
            added_element_method,
            f_z_added,
            z_added,
            added_element_quantity,
        ]
        if added_imp_input_given:
            self.added_element_method = added_element_method
            self.f_z_added = f_z_added
            self.z_abs_added = z_added
            self.added_element_quantity = added_element_quantity

        not_none_ts_params = [
            key for key,val in self.ts_parameter_dict.items() if val is not None and key != 'Sd'
        ] # List alls TS-params, that are not None
        # TODO: make Sd exeption more beautiful or independent from dict
        if imp_input_given:
            # TODO: Perform format parsing on z: imag/real, abs/rad?
            self.f_z = f_z
            self.z_abs = z_abs
            self.z_rad = z_rad
            if (len(not_none_ts_params) > 0):
                warnings.warn(
                    'Impedance curve and TS-Params given: Will overwrite given ' +
                    'TS-params by inherent TS-parameter calculation via ' +
                    'given impedance curve.'
                )
            if added_element_method == 'mass':
                self._calc_cms = self._calc_cms_via_added_mass
                self._calc_mms = self._calc_mms_via_added_mass
                self._calc_vas = self._calc_vas_via_added_mass
                self.imp_to_ts(added_mass = True)
            elif added_element_method == 'volume':
                self._calc_cms = self._calc_cms_via_added_vol
                self._calc_mms = self._calc_mms_via_added_vol
                self._calc_vas = self._calc_vas_via_added_vol
                self.imp_to_ts(added_vol = True)
            else:
                raise ValueError(
                    'Faulty added_element_method ' + 
                    f'{added_element_method}. Has to be either "mass" or' +
                    '"volume"!'
                )
        
        if not imp_input_given and (len(not_none_ts_params) > 0):
            self.Mms = Mms
            self.Rec = Rec
            self.Lec = Lec
            self.Qes = Qes
            self.Qms = Qms
            self.fs = fs
            self.Bl = Bl
            self._update_dependent_ts_params(plot_params=False)

    def _get_z_curve_params(
        self,
        f_z = None,
        z_abs = None,
        manual_pick_title_text = '',
    ):
        Rec = z_abs[0]
        fs, z_max = self._manual_pick_fs(f_z, z_abs, manual_pick_title_text, )
        idx_fs = al_tls.closest_idx_to_val(arr=f_z, val=fs)
        idx_limit_high_freq_f2 = int(2*idx_fs)
        r0 = z_max / Rec
        Z_at_f1_f2 = np.sqrt(r0)*Rec
        idx_f1 = al_tls.closest_idx_to_val(arr=z_abs[:idx_fs], val=Z_at_f1_f2)
        f1  = f_z[idx_f1]
        # Limit f2 search frequency range to [fs:(2*fs)] to avoid Zmax@ high f:
        idx_f2 = idx_fs + al_tls.closest_idx_to_val(
            arr = z_abs[idx_fs:idx_limit_high_freq_f2],
            val = Z_at_f1_f2,
        )
        f2 = f_z[idx_f2]
        return fs, Rec, r0, f1, f2, z_max, 

    def imp_to_ts(
        self,
        added_mass = None,
        added_vol = None,
        plot_params = False,
    ):
        self.fs, self.Rec, self._r0, self._f1, self._f2, self._z_max = self._get_z_curve_params(
            f_z = self.f_z,
            z_abs = self.z_abs,
            manual_pick_title_text = 'Free Air',
        )
        self.Qms = self.fs*np.sqrt(self._r0) / (self._f2 - self._f1)
        self.Qes = self.Qms / (self._r0 - 1)
        self.Qts = self.Qms * self.Qes / (self.Qms + self.Qes)

        # Order of calculations is different for each "added" method:
        if added_mass:
            self.Mms = self._calc_mms()
            self.Cms = self._calc_cms()
            self.Vas = self._calc_vas()
        if added_vol:
            self.Vas = self._calc_vas()
            self.Cms = self._calc_cms()
            self.Mms = self._calc_mms()

        self.Lec = self._estimate_Lec()
        self.Rms = 1 / (2*np.pi*self.fs*self.Qms*self.Cms)
        self.Bl = np.sqrt(self.Rms*(self._z_max - self.Rec))
        # self._update_dependent_ts_params()
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
        freq_resolution : float
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
        """
        Interactive pyplot Figure to manually choose the point of (fs, Z_max).
        Use the mouse to zoom into the plot. Hover over the point of choice and
        confirm by hitting the "Space"-Button.
        """
        fig, ax = al_plt.plot_rfft_freq(f_z, z, xscale='log', )
        ax.set_title(
            f'{z_explanation}:' +
            f'\nManually hover over f$_s$ and Z$_{{max}}$ ' + 
            f'and confirm with "Space"-Button.' + 
            f'\nYou can use left mouse button to zoom etc.'
        )
        plt.tight_layout()
        ax.set_ylabel(r'|Z| [$\Omega$]')
        cursor = al_plt.BlittedCursor(ax) # Crosshair cursor on plot
        fig.canvas.mpl_connect('motion_notify_event', cursor.on_mouse_move)
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
            arr = self.f_z, val = self._Lec_estimation_range[0]
        )
        high_idx = al_tls.closest_idx_to_val(
            arr = self.f_z, val = self._Lec_estimation_range[1]
        )
        Lec = np.mean(
            z_complex.imag[low_idx:high_idx] / omega[low_idx:high_idx]
        )
        return Lec

    def _calc_mms_via_added_mass(self, ):
        fs_new, _ = self._manual_pick_fs(
            self.f_z_added,
            self.z_abs_added,
            'Added Mass',
        )
        print(fs_new)
        Mms = self.added_element_quantity / ((self.fs / fs_new)**2 - 1)
        return Mms

    def _calc_cms_via_added_mass(self, ):
        return 1 / ((2*np.pi*self.fs)**2 * self.Mms)

    def _calc_vas_via_added_mass(self, ):
        return self.rho*self.c**2 * self.Cms * self.Sd**2

    def _calc_vas_via_added_vol(self, ):
        fs_added, Rec, r0, f1, f2, z_max = self._get_z_curve_params(
            f_z = self.f_z_added,
            z_abs = self.z_abs_added,
            manual_pick_title_text = 'Added volume',
        )
        Qmsc = fs_added * np.sqrt(r0) / (f2 - f1) # sc for speaker closed
        Qesc = Qmsc / (r0 - 1)
        return ( (fs_added * Qesc / (self.fs * self.Qes)) - 1) * self.added_element_quantity

    def _calc_cms_via_added_vol(self):
        return self.Vas / (self.rho * self.c**2 * self.Sd**2)

    def _calc_mms_via_added_vol(self, ):
        # TODO: Implement via 
        # Vas = Vb * ( (fs_added_vol * Q_added_vol / (fs * Qes) ) - 1 )
        # Cms = Vas / (rho * c**2 * Sd**2);
        # Mms = 1 / (4*pi**2 * f_added_vol**2 * Cms)
        return 1 / ((2*np.pi*self.fs)**2 * self.Cms)

    def get_baffled_p_resp(self):
        pass

    def get_baffled_sens(self):
        pass

    def plot_z_params(self, ):
        """
        Plot the parameters, from which the TS-parameters are calculated 
        (r0, f1, f2, fs, Rec) on top of the input impedance plot. This is to
        verify the integrity of the calculation.

        Returns:
        fig : 
        """
        v_Rec = self.Rec*np.ones(len(self.f_z))
        v_z_f1_f2 = np.sqrt(self._r0)*self.Rec*np.ones(len(self.f_z))
        fig, ax_mag, ax_arg = al_plt.plot_mag_phase(
            freq_h=self.f_z,
            magnitude=self.z_abs,
            phase_deg=np.array(self.z_rad)/np.pi*180,
            xscale = 'log',
            yscale_mag='lin',
        )
        ax_mag.axvline(x=self._f1, ymin=0, ymax=100, linestyle='--', color='r', label=r'f$_1$')
        ax_mag.axvline(x=self._f2, ymin=0, ymax=100, linestyle='--', color='cyan', label=r'f$_2$')
        ax_mag.axvline(x=self.fs, ymin=0, ymax=100, linestyle='--', color='k', label=r'f$_s$')
        ax_mag.plot(self.f_z, v_z_f1_f2, linestyle='--', label=r'$\sqrt{r_0} R_{ec}$')
        ax_mag.plot(self.f_z, v_Rec, linestyle='--', label=r'$R_{ec}$')
        ax_mag.axvspan(
            xmin=self._Lec_estimation_range[0],
            xmax=self._Lec_estimation_range[1],
            alpha=.5,
            color='b',
            label='$L_{{ec}}$ est. range'
        )
        ax_arg.axvspan(
            xmin=self._Lec_estimation_range[0],
            xmax=self._Lec_estimation_range[1],
            alpha=.5,
            color='b',
            label='$L_{{ec}}$ est. range'
        )
        ax_mag.set(
            ylabel=r'|Z| [$\Omega$]',
            xlabel='',
        )
        ax_mag.legend()
        plt.show(block=False)
        return fig, ax_mag, ax_arg, 

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
        self.Rms = 1 / (2*np.pi*self.fs*self.Qms*self.Cms)
        self.Bl = np.sqrt(self.Rms*(self._z_max - self.Rec))


        # if self.c is None or self.rho is None:
        #     warnings.warn('c and/or rho not defined: Unable to calculate Vas!')
        # else:
        #     self.Vas = (self.Sd*self.c)**2*self.rho/(
        #         (2*np.pi*self.fs)**2*self.Mms
        #     )

    def print_ts(self):
        print('\n')
        print(79*'-')
        print(f'{self.name}: Thiele-Small-Parameters')
        for key, val in self.ts_parameter_dict.items():
            val = (
                'N/A' if val is None else str(
                    np.round(val*STANDARD_PREFIX_CONV[key], 2)
                )
            )
            print('  ' + key + ' = ' + val + f' {STANDARD_UNITS[key]}')
        print(79*'-')
        print('\n')

    @property
    def ts_parameter_dict(self,):
        param_dict = {
            'fs' : self.fs,
            'Sd' : self.Sd,
            'Vas' : self.Vas,
            'Rec' : self.Rec,
            'Lec' : self.Lec,
            'Qts' : self.Qts,
            'Qes' : self.Qes,
            'Qms' : self.Qms,
            'Mms' : self.Mms,
            'Cms' : self.Cms,
            'Rms' : self.Rms,
            'Bl' : self.Bl,
        }
        return param_dict

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self, c):
        raise ValueError(
            'Cannot set speed of sound in existing ' +
            'instance of ElectroDynamic: Creation of new instance necessary!'
        )

    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, rho):
        raise ValueError(
            'Cannot set fluid density in existing ' +
            'instance of ElectroDynamic: Creation of new instance necessary!'
        )

    @property
    def Lec_estimation_range(self):
        return self._Lec_estimation_range

    @Lec_estimation_range.setter
    def Lec_estimation_range(self, Lec_estimation_range):
        self._Lec_estimation_range = Lec_estimation_range
        print('\nRecalculating Lec from new estimation range ... ')
        self._Lec = self._estimate_Lec()
        print(f'New Lec: {np.round(self._Lec*1000, 2)} mH')
