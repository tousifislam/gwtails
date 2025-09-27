#############################################################################
##
##      Filename: tails.py
##
##      Author: Tousif Islam
##
##      Created: 12-07-2023
##
##      Modified: 09-27-2025
##
#############################################################################

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
matplotlib.rcParams['mathtext.fontset'] ='stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral' 
matplotlib.rcParams['axes.linewidth'] = 1 #set the value globally
plt.rcParams["figure.figsize"] = (6,6)
plt.rcParams['font.size'] = '18'

from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import gwtools

class PostMergerAmplitudeFit:
    """
    Class to fit post-merger late-tile tail data
    """
    def __init__(self, t=None, f=None, filename=None, qinput=1000, throw_junk_until_indx=None, 
                 qnm_fit_window=None, tail_fit_window=None, crossterm_fit_window=None,
                 fit_tail_envelop=False, m1_to_M_scale=True):
        """
        Parameters:
        -----------
        t : array-like, optional
            Time array. Default is None - in that case you should provide a filename 
            that contains the strain/psi4.
        f : array-like, optional  
            Function array - it could be strain h_lm or psi4_lm. Default is None - 
            in that case you should provide a filename that contains the strain/psi4.
        
        filename : str, optional
            Name of the waveform data file (Required); should be in .txt or .dat format;
            columns should have the following format: time, h_real, h_imag, ..
            Default is None - in that case you should provide (time 't' and strain/psi4 'f').
        qinput : float, default 1000
            Mass ratio values with qinput>=1.
        
        throw_junk_until_indx : int, optional
            Last index until which data should be discarded before applying qnm or tail fits; 
            Default is None.
        tail_fit_window : list, optional
            Provide a window for fitting tail coefficients only; e.g. [200,800] (in M);
            Default is None (in that case no tail fit is performed).
        qnm_fit_window : list, optional
            Provide a window for fitting qnm coefficients only; e.g. [200,800] (in M);
            Default is None (in that case no qnm fit is performed).
        crossterm_fit_window : list, optional
            Provide a window for fitting cross-term coefficients only; e.g. [100,200] (in M);
            Default is None (in that case no cross-term fit is performed).
        fit_tail_envelop : bool, default False
            Whether the data has oscillatory tail behaviour - in which case, we fit the tail envelop.
        
        m1_to_M_scale : bool, default True
            Whether a mass-scale transformation from m1 (mass of the primary) to total mass M is needed.
        """
        # Handle data input: either from file or direct provision
        if filename is not None:
            self.filename = filename
            self.raw_time, self.raw_h = self._read_data_from_file()
        elif t is not None and f is not None:
            self.raw_time, self.raw_h = t, f
        else:
            raise ValueError("Either provide 'filename' or both 't' and 'f' arrays")
        print("..... time and strain/psi4 data is loaded")
        
        # Mass-scale conversion and data cleaning
        if m1_to_M_scale:
            m1toM = 1 / (1 + 1/qinput)
            self.raw_time, self.raw_h = self.raw_time * m1toM, self.raw_h * m1toM
            print("..... waveform mass-scale changed from m1 to M")
        
        # Remove junk data from start
        self.throw_junk_until_indx = throw_junk_until_indx
        if throw_junk_until_indx is not None:
            self.raw_time = self.raw_time[throw_junk_until_indx:]
            self.raw_h = self.raw_h[throw_junk_until_indx:]
            print("..... junk data is removed from the start")
        
        # Set merger time coordinate and interpolate
        self.peak_time = self._get_peak_time(self.raw_time, self.raw_h)
        self.t_interp, self.h_interp = self._get_interp_data()
        print("..... strain/psi4 is cast onto a common time grid")
                
        # Store and perform tail fits
        self.tail_fit_window = tail_fit_window
        self.fit_tail_envelop = fit_tail_envelop
        if tail_fit_window is not None:
            print("..... fitting the tails with a decaying power law")
            if fit_tail_envelop:
                print('.. tail data has oscillations - fitting tail envelop instead')
                self.popt_tail, self.pcov_tail, self.t_envelop, self.A_envelop = self._fit_tail_envelop_data()
            else:
                self.popt_tail, self.pcov_tail = self._fit_tail_data()
        
        # Store and perform QNM fits
        self.qnm_fit_window = qnm_fit_window
        if qnm_fit_window is not None:
            print("..... fitting QNM with the 220 mode")
            self.popt_qnm, self.pcov_qnm = self._fit_qnm_data()
        
        # Store and perform oscillatory intermediate fits
        self.oscillatory_fit_window = crossterm_fit_window
        if tail_fit_window is not None and qnm_fit_window is not None:
            print("..... fitting oscillatory intermediate part")
            if crossterm_fit_window is None:
                self.oscillatory_fit_window = [qnm_fit_window[1]+15, tail_fit_window[0]-15]
            self.popt_cross_terms, self.pcov_cross_terms = self._fit_cross_terms_data()
        
    def _read_data_from_file(self):
        """Reads data directly from file (filename should be absolute path)."""
        data = np.loadtxt(self.filename, unpack=True)
        print(f'.......... Data read. Shape: {data.shape}')
        return data[0], data[1] - 1j*data[2]
    
    def plot_raw_amplitude(self):
        """Plots raw amplitudes without time shift or junk removal."""
        plt.figure(figsize=(6,6))
        plt.semilogy(self.raw_time, abs(self.raw_h), color='royalblue')
        plt.xlabel('raw time')
        plt.ylabel('raw amplitude')
    
    def _get_peak_via_quadratic_fit(self, t, func):
        """Finds peak using quadratic fit around maximum value."""
        index = np.clip(np.argmax(func), 2, len(t) - 3)
        
        # Setup quadratic fit around peak
        t_local = t[index-2:index+3] - t[index]
        f_local = func[index-2:index+3]
        
        # Solve normal equations for quadratic fit
        X = np.vstack([np.ones(5), t_local, t_local**2])
        XTX_inv = np.linalg.inv(X @ X.T)
        coefs = XTX_inv @ X @ f_local
        
        # Return peak time and amplitude
        peak_offset = -coefs[1] / (2 * coefs[2])
        peak_amp = coefs[0] - coefs[1]**2 / (4 * coefs[2])
        return t[index] + peak_offset, peak_amp
    
    def _get_peak_time(self, t, modes):
        """Finds peak time using quadratic fit around maximum amplitude."""
        return self._get_peak_via_quadratic_fit(t, abs(modes)**2)[0]
    
    def _get_interp_data(self):
        """Interpolate amplitudes after time shift or junk removal."""
        t_transform = self.raw_time - self.peak_time
        t_interp = np.arange(-50,max(t_transform)-100,0.1)
        h_interp = gwtools.interpolate_h(t_transform, self.raw_h, t_interp)
        return t_interp, h_interp
    
    def plot_interpolated_amplitude(self):
        """Plots interpolated amplitudes after time shift or junk removal."""
        plt.figure(figsize=(6,6))
        plt.semilogy(self.t_interp, abs(self.h_interp), color='royalblue')
        plt.xlabel('time')
        plt.ylabel('amplitude')
        
    def freq(time, h, s=0.1):
        """Calculate frequency from phase derivative."""
        phase_interp = gwtools.interpolant_h(time, gwtools.phase(h), s=s)
        return abs(splev(time, phase_interp, der=1))
    
    def _cut_data(self, tini, tout):
        """Extract data within time window."""
        mask = (self.t_interp >= tini) & (self.t_interp <= tout)
        return self.t_interp[mask], abs(self.h_interp[mask])
    
    def qnm_fit_func(self, t, Aqnm, tau):
        """QNM exponential decay function."""
        return Aqnm * np.exp(-t/tau)
    
    def tail_fit_func(self, t, Atail, c, n):
        """Tail power-law decay function."""
        return Atail / (t + c)**n
    
    def qnm_and_tail_fit_func(self, t, phi_tail, omega):
        """Combined QNM + tail + cross-terms fitting function."""
        qnm_sq = self.qnm_fit_func(t, *self.popt_qnm)**2
        tail_sq = self.tail_fit_func(t, *self.popt_tail)**2
        
        # Cross terms
        qnm_amp = self.popt_qnm[0] * np.exp(-t/self.popt_qnm[1])
        tail_amp = self.popt_tail[0] * (t + self.popt_tail[1])**(-self.popt_tail[2])
        cross = 2 * qnm_amp * tail_amp * np.cos(phi_tail + omega*t)
        
        return np.sqrt(qnm_sq + tail_sq + cross)
    
    def _fit_tail_data(self):
        """Fit tail data with power-law decay."""
        xdata, ydata = self._cut_data(*self.tail_fit_window)
        return curve_fit(self.tail_fit_func, xdata, ydata, maxfev=10000)
    
    def _fit_tail_envelop_data(self):
        """Fit envelope of oscillatory tail data."""
        xdata, ydata = self._cut_data(*self.tail_fit_window)
        peaks = find_peaks(np.log10(ydata))[0]
        x_env, y_env = xdata[peaks], ydata[peaks]
        popt, pcov = curve_fit(self.tail_fit_func, x_env, y_env, maxfev=10000)
        return popt, pcov, x_env, y_env
    
    def _fit_qnm_data(self):
        """Fit QNM data with exponential decay."""
        xdata, ydata = self._cut_data(*self.qnm_fit_window)
        return curve_fit(self.qnm_fit_func, xdata, ydata, maxfev=10000)
    
    def _fit_cross_terms_data(self):
        """Fit oscillatory intermediate region."""
        xdata, ydata = self._cut_data(*self.oscillatory_fit_window)
        return curve_fit(self.qnm_and_tail_fit_func, xdata, ydata, maxfev=10000)
        
    def _plot_qnm_fit(self, t_plot=None, xmax=None):
        """Plot QNM fit comparison."""
        if t_plot is None:
            t_plot = np.arange(10, self.tail_fit_window[1], 0.1)
        
        y_min = 10**np.floor(np.log10(min(abs(self.h_interp))))
        
        plt.semilogy(self.t_interp, abs(self.h_interp), color='grey', lw=3, alpha=0.4, label='Data')
        plt.semilogy(t_plot, self.qnm_fit_func(t_plot, *self.popt_qnm), '--', label='QNM Fit')
        if xmax is None:
            plt.xlim(0, 1020)
        else:
            plt.xlim(0, xmax)
        plt.ylim(y_min, 1e1)
        plt.xlabel('t')
        plt.ylabel('Amplitude')
        plt.legend(fontsize=14)
        plt.show()
    
    def _plot_tail_fits(self, t_plot=None, xmax=None):
        """Plot tail fit comparison."""
        if t_plot is None:
            t_plot = np.arange(10, self.tail_fit_window[1], 0.1)
        
        y_min = 10**np.floor(np.log10(min(abs(self.h_interp))))
        
        plt.semilogy(self.t_interp, abs(self.h_interp), color='grey', lw=4, alpha=0.4, label='Data')
        plt.semilogy(t_plot, self.tail_fit_func(t_plot, *self.popt_tail), '--', label='Tail Fit')
        if xmax is None:
            plt.xlim(-100, 1020)
        else:
            plt.xlim(-100, xmax)
        plt.ylim(y_min, 1e1)
        plt.xlabel('time')
        plt.ylabel('Amplitude')
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.show()
    
    def _plot_all_fits(self, t_plot=None, xmax=None):
        """Plot all fits in a two-panel comparison."""
        if t_plot is None:
            t_plot = np.arange(10, self.tail_fit_window[1], 0.1)
        
        y_min = 10**np.floor(np.log10(min(abs(self.h_interp))))
        
        plt.figure(figsize=(14, 6))
        
        # Left panel: Individual fits
        plt.subplot(121)
        plt.semilogy(self.t_interp, abs(self.h_interp), color='grey', lw=4, alpha=0.4, label='Data')
        plt.semilogy(t_plot, self.tail_fit_func(t_plot, *self.popt_tail), '--', label='Tail Fit')
        plt.semilogy(t_plot, self.qnm_fit_func(t_plot, *self.popt_qnm), '--', label='QNM Fit')
        if xmax is None:
            plt.xlim(-100, 1020)
        else:
            plt.xlim(-100, xmax)
        plt.ylim(y_min, 1e1)
        plt.xlabel('time')
        plt.ylabel('Amplitude')
        plt.legend(fontsize=14)
        
        # Right panel: Combined fit
        plt.subplot(122)
        plt.semilogy(self.t_interp, abs(self.h_interp), color='grey', lw=4, alpha=0.4, label='Data')
        plt.semilogy(t_plot, self.qnm_and_tail_fit_func(t_plot, *self.popt_cross_terms), 
                    '--', color='g', label='QNM+Tail Fit')
        if xmax is None:
            plt.xlim(-100, 1020)
        else:
            plt.xlim(-100, xmax)
        plt.ylim(y_min, 1e1)
        plt.xlabel('time')
        plt.ylabel('Amplitude')
        plt.legend(fontsize=14)
        
        plt.tight_layout()
        plt.show()
        
    def _print_all_fits(self):
        """Print all fitted parameters."""
        if self.tail_fit_window is not None:
            print(f'Atail : {self.popt_tail[0]:.9e}')
            print(f'c     : {self.popt_tail[1]:.9e}')
            print(f'n     : {self.popt_tail[2]:.9e}')
        
        if self.qnm_fit_window is not None:
            print(f'Aqnm  : {self.popt_qnm[0]:.9f}')
            print(f'tau   : {self.popt_qnm[1]:.9f}')
        
        if self.tail_fit_window is not None and self.qnm_fit_window is not None:
            if hasattr(self, 'popt_cross_terms'):
                print(f'phi_tail : {self.popt_cross_terms[0]:.9f}')
                print(f'omega    : {self.popt_cross_terms[1]:.9f}')