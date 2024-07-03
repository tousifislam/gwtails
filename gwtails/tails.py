#############################################################################
##
##      Filename: tails.py
##
##      Author: Tousif Islam
##
##      Created: 12-07-2023
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
    Class to fit post-merger data
    """
    def __init__(self, filename, qinput=1000, throw_junk_until_indx=None, 
                 qnm_fit_window=None, tail_fit_window=None, crossterm_fit_window=None,
                 fit_tail_envelop=False):
        """
        filename: name of the waveform data file (Required); should be in .txt or .dat format;
                  columns should have the following format: time, h_real, h_imag, ..
                  
        qinput: mass ratio values with qinput>=1
        throw_junk_until_indx: last index until which data should be discarded before applying qnm or tail fits; 
                               Default is None
                               
        tail_fit_window: provide a window for fitting tail coefficients only;
                         e.g. [200,800] (in M);
                         Default is None (in that case no tail fit is performed)
        qnm_fit_window:  provide a window for fitting qnm coefficients only;
                         e.g. [200,800] (in M);
                         Default is None (in that case no qnm fit is performed)
        crossterm_fit_window: provide a window for fitting cross-term coefficients only;
                              e.g. [100,200] (in M);
                              Default is None (in that case no cross-term fit is performed)

        fit_tail_envelop: Whether the data has oscillatory tail behaviour - in which case, we fit the tail envelop;
                          Default is False.
        """
        self.filename = filename
        self.raw_time, self.raw_h = self._read_data_from_file()
        m1toM = 1 / (1+1/qinput)
        self.raw_time = self.raw_time * m1toM
        self.raw_h = self.raw_h * m1toM
        
        self.throw_junk_until_indx = throw_junk_until_indx
        if self.throw_junk_until_indx is not None:
            self.raw_time = self.raw_time[self.throw_junk_until_indx:]
            self.raw_h = self.raw_h[self.throw_junk_until_indx:]
            
        self.peak_time = self._get_peak_time(self.raw_time, self.raw_h)
        self.t_interp, self.h_interp = self._get_interp_data()
        
        self.tail_fit_window = tail_fit_window
        self.fit_tail_envelop = fit_tail_envelop
        if self.tail_fit_window is not None:
            if self.fit_tail_envelop is True:
                print('.. tail data has oscillations - fitting tail envelop instead')
                self.popt_tail, self.pcov_tail, self.t_envelop, self.A_envelop = self._fit_tail_envelop_data()
            else:
                self.popt_tail, self.pcov_tail = self._fit_tail_data()
            
        self.qnm_fit_window = qnm_fit_window
        if self.qnm_fit_window is not None:
            self.popt_qnm, self.pcov_qnm = self._fit_qnm_data()
        
        self.oscillatory_fit_window = crossterm_fit_window
        if self.tail_fit_window is not None and self.qnm_fit_window is not None:
            if self.oscillatory_fit_window is None:
                self.oscillatory_fit_window = [self.qnm_fit_window[1]+15, self.tail_fit_window[0]-15]
            self.popt_cross_terms, self.pcov_cross_terms = self._fit_cross_terms_data()
        
    def _read_data_from_file(self):
        data = np.loadtxt(self.filename, unpack=True)
        print('Data read. Shape of data : (%d,%d)'%(data.shape[0],data.shape[1]))
        raw_time = data[0]
        raw_h = data[1] - 1j*data[2]
        return raw_time, raw_h
    
    def plot_raw_amplitude(self):
        plt.figure(figsize=(6,6))
        plt.semilogy(self.raw_time, abs(self.raw_h), color='royalblue')
        plt.xlabel('raw time')
        plt.ylabel('raw amplitude')
        
    def _get_peak_via_quadratic_fit(self, t, func):
        index = np.argmax(func)
        index = max(2, min(len(t) - 3, index))
        testTimes = t[index-2:index+3] - t[index]
        testFuncs = func[index-2:index+3]
        xVecs = np.array([np.ones(5),testTimes,testTimes**2.])
        invMat = np.linalg.inv(np.array([[v1.dot(v2) for v1 in xVecs] \
            for v2 in xVecs]))
        yVec = np.array([testFuncs.dot(v1) for v1 in xVecs])
        coefs = np.array([yVec.dot(v1) for v1 in invMat])
        return t[index] - coefs[1]/(2.*coefs[2]), coefs[0] - coefs[1]**2./4/coefs[2]
    
    def _get_peak_time(self, t, modes):
        normSqrVsT = abs(modes)**2
        return self._get_peak_via_quadratic_fit(t, normSqrVsT)[0]
    
    def _get_interp_data(self):
        t_transform = self.raw_time - self.peak_time
        t_interp = np.arange(-50,max(t_transform)-100,0.1)
        h_interp = gwtools.gwtools.interpolate_h(t_transform, self.raw_h, t_interp)
        return t_interp, h_interp
    
    def plot_interpolated_amplitude(self):
        plt.figure(figsize=(6,6))
        plt.semilogy(self.t_interp, abs(self.h_interp), color='royalblue')
        plt.xlabel('time')
        plt.ylabel('amplitude')
        
    def freq(time, h, s=0.1):
        tmp = gwtools.interpolant_h(time, gwtools.phase(h), s=s)
        f = splev(time, tmp, der=1)
        return abs(f)

    def _cut_data(self, tini, tout):
        input_time = self.t_interp 
        input_h = self.h_interp
        indx = np.logical_and((input_time>=tini),(input_time<=tout))
        twindow = input_time[indx]
        Awindow = abs(input_h[indx])
        return twindow, Awindow
    
    def qnm_fit_func(self, t, Aqnm, tau):
        return Aqnm * np.exp(-t/tau)

    def tail_fit_func(self, t, Atail, c, n):
        return Atail / (t+c)**n
    
    def qnm_and_tail_fit_func(self, t, phi_tail, omega):
        qnm_terms = self.qnm_fit_func(t, *self.popt_qnm)**2
        tail_terms =  self.tail_fit_func(t, *self.popt_tail)**2
        cross_terms =  2 * self.popt_qnm[0] * np.exp(-t/self.popt_qnm[1]) * self.popt_tail[0] * ((t+self.popt_tail[1])**(-self.popt_tail[2])) * np.cos(phi_tail + omega*t)
        return np.sqrt(qnm_terms + tail_terms + cross_terms)
    
    def _fit_tail_data(self):
        xdata, ydata = self._cut_data(tini=self.tail_fit_window[0], tout=self.tail_fit_window[1])
        popt, pcov = curve_fit(self.tail_fit_func, xdata, ydata, maxfev=10000)
        return popt, pcov
    
    def _fit_tail_envelop_data(self):
        xdata, ydata = self._cut_data(tini=self.tail_fit_window[0], tout=self.tail_fit_window[1])
        peaks = find_peaks(np.log10(ydata))[0]
        xdata_envelop = xdata[peaks]
        ydata_envelop = ydata[peaks]
        popt, pcov = curve_fit(self.tail_fit_func, xdata_envelop, ydata_envelop, maxfev=10000)
        return popt, pcov, xdata_envelop, ydata_envelop
    
    def _fit_qnm_data(self):
        xdata, ydata = self._cut_data(tini=self.qnm_fit_window[0], tout=self.qnm_fit_window[1])
        popt, pcov = curve_fit(self.qnm_fit_func, xdata, ydata, maxfev=10000)
        return popt, pcov
    
    def _fit_cross_terms_data(self):
        xdata, ydata = self._cut_data(tini=self.oscillatory_fit_window[0], tout=self.oscillatory_fit_window[1])
        popt, pcov = curve_fit(self.qnm_and_tail_fit_func, xdata, ydata, maxfev=10000)
        return popt, pcov
    
    def _plot_qnm_fit(self):
        tcommon = np.arange(10,1000,0.1)
        plt.semilogy(self.t_interp, abs(self.h_interp), color='grey', lw=3, alpha=0.4, label='Data')
        plt.semilogy(tcommon, self.qnm_fit_func(tcommon, *self.popt_qnm), '--', label='QNM Fit')
        plt.xlim(0,1020)
        plt.ylim(10**np.floor(np.log10(min(abs(self.h_interp)))),1e1)
        plt.xlabel('t')
        plt.ylabel('Amplitude')
        plt.legend(fontsize=14)
        plt.show()
        
    def _plot_tail_fits(self):
        tcommon = np.arange(10,1000,0.1)
        plt.semilogy(self.t_interp, abs(self.h_interp), color='grey', lw=4, alpha=0.4, label='Data')
        plt.semilogy(tcommon, self.tail_fit_func(tcommon, *self.popt_tail), '--', label='Tail Fit')
        plt.xlim(-100,1000)
        plt.ylim(10**np.floor(np.log10(min(abs(self.h_interp)))),1e1)
        plt.xlabel('time')
        plt.ylabel('Amplitude')
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.show()
        
    def _plot_all_fits(self):
        tcommon = np.arange(10,1000,0.1)

        plt.figure(figsize=(14,6))
        plt.subplot(121)
        plt.semilogy(self.t_interp, abs(self.h_interp), color='grey', lw=4, alpha=0.4, label='Data')
        plt.semilogy(tcommon, self.tail_fit_func(tcommon, *self.popt_tail), '--', label='Tail Fit')
        plt.semilogy(tcommon, self.qnm_fit_func(tcommon, *self.popt_qnm), '--', label='QNM Fit')
        plt.xlim(-100,1000)
        plt.ylim(10**np.floor(np.log10(min(abs(self.h_interp)))),1e1)
        plt.xlabel('time')
        plt.ylabel('Amplitude')
        plt.legend(fontsize=14)
        plt.subplot(122)
        plt.semilogy(self.t_interp, abs(self.h_interp), color='grey', lw=4, alpha=0.4, label='Data')
        plt.semilogy(tcommon, self.qnm_and_tail_fit_func(tcommon, *self.popt_cross_terms), '--', color='g', label='QNM+Tail Fit')
        plt.xlim(-100,1000)
        plt.ylim(10**np.floor(np.log10(min(abs(self.h_interp)))),1e1)
        plt.xlabel('time')
        plt.ylabel('Amplitude')
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.show()
        
    def _print_all_fits(self):
        if self.tail_fit_window is not None:
            print('Atail : %.9e'%self.popt_tail[0])
            print('c : %.9e'%self.popt_tail[1])
            print('n : %.9e'%self.popt_tail[2])
        if self.qnm_fit_window is not None:
            print('Aqnm : %.9f'%self.popt_qnm[0])
            print('tau : %.9f'%self.popt_qnm[1])
        if self.tail_fit_window is not None and self.qnm_fit_window is not None:
            if self.oscillatory_fit_window is None:
                print('phi_tail : %.9f'%self.popt_cross_terms[0])
                print('omega : %.9f'%self.popt_cross_terms[1])