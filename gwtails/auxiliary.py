#############################################################################
##
##      Filename: auxiliary.py
##
##      Author: Tousif Islam
##
##
##      Modified: 09-27-2025
##
#############################################################################

import numpy as np

from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import gwtools

def compute_tail_exponent(t, f, use_filter=True, window_length=800, polyorder=4):
    """
    Compute tail exponent p from time series data.
    
    Parameters:
    -----------
    t : array-like
        Time array
    f : array-like  
        Strain/psi4 amplitude data
    use_filter : bool, default True
        Whether to apply Savgol filter to the computed exponent
    window_length : int, default 800
        Length of filter window (must be odd)
    polyorder : int, default 4
        Order of polynomial used in filtering
        
    Returns:
    --------
    p : array-like
        Tail exponent (filtered if use_filter=True)
    """
    from scipy.signal import savgol_filter
    
    # Compute tail exponent
    p = np.gradient(np.log(abs(f))) / np.gradient(np.log(t))
    
    # Apply filter if requested
    if use_filter:
        # Ensure window_length is odd and valid
        if window_length % 2 == 0:
            window_length += 1
        window_length = min(window_length, len(p))
        
        p = savgol_filter(p, window_length, polyorder)
    
    return p

def compute_frequency(t, h, method='spline', s=0.01, window_length=800, polyorder=4):
    """
    Computes orbital frequency of a given gravitational wave time-series.
    
    Parameters:
    -----------
    t : array-like
        Time array
    h : array-like
        Waveform mode array or psi4 array
    method : str, default 'spline'
        Method for computing frequency: 'spline' or 'savgol'
    s : float, default 0.01
        Smoothing parameter for spline interpolation
    window_length : int, default 800
        Window length for Savgol filter (used if method='savgol')
    polyorder : int, default 4
        Polynomial order for Savgol filter
        
    Returns:
    --------
    freq : array-like
        Orbital frequency
    """
    from scipy.interpolate import splev
    from scipy.signal import savgol_filter
    
    if method == 'spline':
        # Use spline interpolation with smoothing
        phase_interp = gwtools.interpolant_h(t, gwtools.phase(h), s=s)
        return abs(splev(t, phase_interp, der=1))
    
    elif method == 'savgol':
        # Use Savgol filter on phase derivative
        phase = gwtools.phase(h)
        freq_raw = np.gradient(phase, t)
        
        # Ensure window_length is odd and valid
        if window_length % 2 == 0:
            window_length += 1
        window_length = min(window_length, len(freq_raw))
        
        return abs(savgol_filter(freq_raw, window_length, polyorder))
    
    else:
        raise ValueError("method must be 'spline' or 'savgol'")

def find_tail_start(t, h, freq_threshold=1e-3, consecutive_zeros=5, 
                    method='spline', merger_time=10, **freq_kwargs):
    """
    Find the start of the tail by detecting consecutive frequency zeros.
    
    Parameters:
    -----------
    t : array-like
        Time array
    h : array-like
        Strain/psi4 data
    freq_threshold : float, default 1e-4
        Threshold below which absolute frequency is considered zero
    consecutive_zeros : int, default 5
        Number of consecutive zero frequencies required
    method : str, default 'spline'
        Method for frequency computation ('spline' or 'savgol')
    merger_time : float, default 10
        Time before which data is ignored (ensures tail is after merger)
    **freq_kwargs : dict
        Additional arguments for compute_frequency function
        
    Returns:
    --------
    tail_start_time : float
        Time at which tail starts
    tail_start_index : int
        Index at which tail starts (in original array)
    frequency : array
        Computed normalized frequency array
    zero_mask : array
        Boolean mask showing where frequency is considered zero
    """
    
    # Filter data to only consider times after merger
    post_merger_mask = t > merger_time
    t_post = t[post_merger_mask]
    h_post = h[post_merger_mask]
    
    if len(t_post) == 0:
        print(f"No data found after merger time t={merger_time}")
        return None, None, None, None
    
    # Compute frequency on post-merger data
    frequency = compute_frequency(t_post, h_post, method=method, **freq_kwargs)
    
    # Normalize frequency by maximum value
    frequency = frequency / max(frequency)
    
    # Create mask for "zero" frequencies using absolute value
    zero_mask = np.abs(frequency) < freq_threshold
    
    # Find consecutive zero sequences
    def find_consecutive_groups(mask, min_length):
        """Find groups of consecutive True values with minimum length."""
        # Find transitions
        diff = np.diff(np.concatenate(([False], mask, [False])).astype(int))
        starts = np.where(diff == 1)[0]
        ends = np.where(diff == -1)[0]
        
        # Filter by minimum length
        valid_groups = []
        for start, end in zip(starts, ends):
            if end - start >= min_length:
                valid_groups.append((start, end))
        
        return valid_groups
    
    # Find groups of consecutive zeros
    zero_groups = find_consecutive_groups(zero_mask, consecutive_zeros)
    
    if len(zero_groups) == 0:
        print(f"No sequences of {consecutive_zeros} consecutive zero frequencies found after t={merger_time}")
        return None, None, frequency, zero_mask
    
    # Take the first group as the tail start
    tail_start_index_post = zero_groups[0][0]
    tail_start_time = t_post[tail_start_index_post]
    
    return tail_start_time

def calculate_chi_square(data, model, sigma):
    """Calculate chi-square statistic."""
    chi2 = np.sum((data - model)**2 / sigma**2)
    return chi2

def reduced_chi2(data, model, sigma, n_params):
    """Calculate reduced chi-square (chi2 / degrees of freedom)."""
    chi2 = np.sum((data - model)**2 / sigma**2)
    n_data = len(data)
    dof = n_data - n_params  # degrees of freedom
    return chi2 / dof
