#############################################################################
##
##      Filename: bayesian_tails.py
##
##      Author: Tousif Islam
##
##      Created: 12-07-2023
##
##      Modified: 09-27-2025
##
#############################################################################

import emcee
import numpy as np
import corner

import emcee
import numpy as np

class TailFitMCMC:
    """
    Unified Bayesian fitting class for tail decay function: Atail * (t + ctail)**ptail
    Supports fixing any combination of parameters.

    Example:
    
        # Create fit instance
        fitter = gwtails.TailFitMCMC(t=t, # post-merger time
                                     A=abs(h[(2,2)]), # post-merger amplitue
                                     t_tail_window=[1200, 8000], 
                                     lsq_params=None, # no initial guess
                                     log_Atail_range=(0, 20), 
                                     ctail_range=(0, 2500), 
                                     ptail_range=(-15, -1),
                                     fixed_params=None,
                                     percentage_err_in_data=2.5)
        
        # Perform fitting
        samples, percentiles, log_likelihood = fitter.fit(n_walkers=32, n_steps=2000, burn_in=500)
        
        # quick corner plot
        fitter.plot_corner();

    """
    
    def __init__(self, t, A, t_tail_window, lsq_params=None, 
                 log_Atail_range=(-5, 15), ctail_range=(0, 2500), ptail_range=(-15, -1),
                 fixed_params=None, percentage_err_in_data=2.5):
        """
        Initialize the class.
        
        Parameters:
            t : array-like
                Time array
            A : array-like  
                Amplitude data
            t_tail_window : list
                [t_start, t_end] for fitting window
            lsq_params : list or array, optional
                [Atail, ctail, ptail] from least squares fit to use as initial guess
            log_Atail_range : tuple, optional
                (min, max) range for log10(Atail) prior, default (-5, 15)
            ctail_range : tuple, optional
                (min, max) range for ctail parameter prior, default (0, 2500)
            ptail_range : tuple, optional
                (min, max) range for ptail parameter prior, default (-15, -1)
            fixed_params : dict, optional
                Dictionary specifying fixed parameters, e.g., {'ptail': -6}, {'ctail': 0}, 
                {'ptail': -6, 'ctail': 0}, or None for all free parameters
            percentage_err_in_data : float, optional
                Percentage error in data for likelihood calculation, default 2.5%
        """
        self.t = np.array(t)
        self.A = np.array(A)
        self.t_tail_window = t_tail_window
        self.lsq_params = lsq_params
        self.percentage_err_in_data = percentage_err_in_data
        
        # Prior ranges
        self.log_Atail_range = log_Atail_range
        self.ctail_range = ctail_range
        self.ptail_range = ptail_range
        
        # Handle fixed parameters
        self.fixed_params = fixed_params if fixed_params is not None else {}
        
        # Determine free parameters
        self.param_names = ['log10_Atail', 'ctail', 'ptail']
        self.free_params = [p for p in self.param_names if p not in self.fixed_params]
        self.ndim = len(self.free_params)
        
        # Extract data within window
        mask = (self.t >= t_tail_window[0]) & (self.t <= t_tail_window[1])
        self.t_fit = self.t[mask]
        self.A_fit = self.A[mask]
        
        # Initialize results as None
        self.samples = None
        self.percentiles = None
        self.sampler = None
        self.log_likelihood_samples = None
    
    def _unpack_params(self, free_param_values):
        """Convert free parameter values to full parameter set."""
        params = {}
        free_idx = 0
        
        for param in self.param_names:
            if param in self.fixed_params:
                params[param] = self.fixed_params[param]
            else:
                params[param] = free_param_values[free_idx]
                free_idx += 1
        
        return params['log10_Atail'], params['ctail'], params['ptail']
    
    def tail_model(self, t, free_param_values):
        """Tail decay model."""
        log10_Atail, ctail, ptail = self._unpack_params(free_param_values)
        Atail = 10**log10_Atail
        return Atail * (t + ctail)**ptail
    
    def log_likelihood(self, free_param_values):
        """Log likelihood function."""
        try:
            model = self.tail_model(self.t_fit, free_param_values)
            if np.any(~np.isfinite(model)) or np.any(model <= 0):
                return -np.inf
            
            sigma = (self.percentage_err_in_data / 100.0) * self.A_fit  # Relative error
            chi2 = np.sum((self.A_fit - model)**2 / sigma**2)
            return -0.5 * chi2
        except:
            return -np.inf
    
    def log_prior(self, free_param_values):
        """Log prior function."""
        params_dict = {}
        free_idx = 0
        
        for param in self.param_names:
            if param not in self.fixed_params:
                params_dict[param] = free_param_values[free_idx]
                free_idx += 1
        
        # Check bounds for free parameters only
        if 'log10_Atail' in params_dict:
            if not (self.log_Atail_range[0] < params_dict['log10_Atail'] < self.log_Atail_range[1]):
                return -np.inf
        
        if 'ctail' in params_dict:
            if not (self.ctail_range[0] < params_dict['ctail'] < self.ctail_range[1]):
                return -np.inf
        
        if 'ptail' in params_dict:
            if not (self.ptail_range[0] < params_dict['ptail'] < self.ptail_range[1]):
                return -np.inf
        
        return 0.0
    
    def log_probability(self, free_param_values):
        """Log posterior probability."""
        lp = self.log_prior(free_param_values)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(free_param_values)
    
    def _get_initial_guess(self):
        """Get initial guess for free parameters."""
        if self.lsq_params is not None:
            Atail_lsq, ctail_lsq, ptail_lsq = self.lsq_params
            log10_Atail_lsq = np.log10(abs(Atail_lsq))
            all_params = {'log10_Atail': log10_Atail_lsq, 'ctail': ctail_lsq, 'ptail': ptail_lsq}
            print("..... Using least-squares parameters as initial guess")
        else:
            log_A_median = np.log10(1/np.median(self.A_fit))
            all_params = {'log10_Atail': log_A_median, 'ctail': 50, 'ptail': -5}
            print("..... Using data-driven initial guess")
        
        # Extract only free parameters
        initial_guess = [all_params[p] for p in self.free_params]
        
        # Print initial guess
        guess_str = ", ".join([f"{p}={all_params[p]:.3f}" for p in self.param_names 
                               if p not in self.fixed_params])
        fixed_str = ", ".join([f"{p}={self.fixed_params[p]:.3f} (fixed)" 
                              for p in self.fixed_params])
        
        if fixed_str:
            print(f".......... Initial guess: {guess_str}, {fixed_str}")
        else:
            print(f".......... Initial guess: {guess_str}")
        
        return initial_guess
    
    def _initialize_walkers(self, initial_guess, n_walkers):
        """Initialize MCMC walkers."""
        # Get bounds and widths for free parameters
        ranges = {
            'log10_Atail': self.log_Atail_range,
            'ctail': self.ctail_range,
            'ptail': self.ptail_range
        }
        
        bounded_guess = initial_guess.copy()
        scatter = []
        
        for i, param in enumerate(self.free_params):
            param_range = ranges[param]
            # Ensure within bounds
            if not (param_range[0] < bounded_guess[i] < param_range[1]):
                bounded_guess[i] = (param_range[0] + param_range[1]) / 2
            # Calculate scatter as 10% of range width
            scatter.append(0.1 * (param_range[1] - param_range[0]))
        
        scatter = np.array(scatter)
        pos = []
        max_attempts = 1000
        
        for i in range(n_walkers):
            attempts = 0
            while attempts < max_attempts:
                candidate = bounded_guess + scatter * np.random.randn(self.ndim)
                
                # Check if within bounds
                valid = True
                for j, param in enumerate(self.free_params):
                    param_range = ranges[param]
                    if not (param_range[0] < candidate[j] < param_range[1]):
                        valid = False
                        break
                
                if valid:
                    pos.append(candidate)
                    break
                attempts += 1
            
            # Fallback to uniform sampling
            if attempts >= max_attempts:
                candidate = np.array([
                    np.random.uniform(ranges[param][0], ranges[param][1]) 
                    for param in self.free_params
                ])
                pos.append(candidate)
        
        return np.array(pos)
    
    def fit(self, n_walkers=32, n_steps=2000, burn_in=500):
        """
        Perform Bayesian fitting using MCMC.
        
        Parameters:
        -----------
        n_walkers : int, default 32
            Number of MCMC walkers
        n_steps : int, default 2000
            Number of MCMC steps
        burn_in : int, default 500
            Number of burn-in steps to discard
        """
        # Print fitting message
        fixed_str = ", ".join([f"{p}={v}" for p, v in self.fixed_params.items()])
        if fixed_str:
            print(f"..... Fitting tail data with the function : Atail * (t + ctail)**ptail ({fixed_str} fixed)")
        else:
            print("..... Fitting tail data with the function : Atail * (t + ctail)**ptail")
        
        # Get initial guess
        initial_guess = self._get_initial_guess()
        
        # Initialize walkers
        pos = self._initialize_walkers(initial_guess, n_walkers)
        
        # Run MCMC
        self.sampler = emcee.EnsembleSampler(
            n_walkers, self.ndim, 
            lambda params, *args: self.log_probability(params)
        )
        
        print("..... Running MCMC")
        self.sampler.run_mcmc(pos, n_steps, progress=True)
        
        # Get samples after burn-in
        self.samples = self.sampler.get_chain(discard=burn_in, flat=True)
        
        # Get log likelihood values after burn-in
        self.log_likelihood_samples = self.sampler.get_log_prob(discard=burn_in, flat=True)
        
        # Calculate percentiles
        self._calculate_percentiles()
        
        print("..... Done!")
        
        return self.samples, self.percentiles, self.log_likelihood_samples
    
    def _calculate_percentiles(self):
        """Calculate parameter percentiles from MCMC samples."""
        self.percentiles = {}
        
        print("..... Bayesian Results:")
        
        # Process free parameters
        for i, param in enumerate(self.free_params):
            mcmc = np.percentile(self.samples[:, i], [16, 50, 84])
            q = np.diff(mcmc)
            self.percentiles[param] = {
                'median': mcmc[1],
                'upper_err': q[1], 
                'lower_err': q[0]
            }
            
            if param == 'log10_Atail':
                # Report both log and linear scale
                Atail_mcmc = 10**mcmc
                Atail_q = np.diff(Atail_mcmc)
                self.percentiles['Atail'] = {
                    'median': Atail_mcmc[1],
                    'upper_err': Atail_q[1],
                    'lower_err': Atail_q[0]
                }
                avg_err = (q[0] + q[1]) / 2
                avg_Atail_err = (Atail_q[0] + Atail_q[1]) / 2
                print(f".......... log10_Atail = {mcmc[1]:.3f} +/- {avg_err:.3f}")
                print(f".......... Atail = {Atail_mcmc[1]:.3e} +/- {avg_Atail_err:.3e}")
            else:
                avg_err = (q[0] + q[1]) / 2
                print(f".......... {param} = {mcmc[1]:.3f} +/- {avg_err:.3f}")
        
        # Store fixed parameters
        for param, value in self.fixed_params.items():
            self.percentiles[param] = {'median': value, 'upper_err': 0, 'lower_err': 0}
            print(f".......... {param} = {value:.3f} (fixed)")
        
        # Print log likelihood statistics
        if self.log_likelihood_samples is not None:
            log_like_mcmc = np.percentile(self.log_likelihood_samples, [16, 50, 84])
            log_like_q = np.diff(log_like_mcmc)
            avg_log_like_err = (log_like_q[0] + log_like_q[1]) / 2
            print(f".......... log_likelihood = {log_like_mcmc[1]:.2f} +/- {avg_log_like_err:.2f}")
    
    def get_best_fit_params(self):
        """Get best-fit parameters from percentiles."""
        if self.percentiles is None:
            raise ValueError("Must run fit() first")
        
        return {
            'Atail': self.percentiles['Atail']['median'],
            'ctail': self.percentiles['ctail']['median'],
            'ptail': self.percentiles['ptail']['median']
        }
    
    def predict(self, t_pred=None):
        """
        Predict model values using best-fit parameters.
        
        Parameters:
        -----------
        t_pred : array-like, optional
            Time points for prediction. If None, uses fitting time range.
        """
        if self.percentiles is None:
            raise ValueError("Must run fit() first")
        
        if t_pred is None:
            t_pred = self.t_fit
        
        best_params = self.get_best_fit_params()
        return best_params['Atail'] * (t_pred + best_params['ctail'])**best_params['ptail']
    
    def plot_corner(self, **kwargs):
        """
        Create a corner plot of the MCMC samples.
        
        Parameters:
        -----------
        **kwargs : optional
            Additional keyword arguments to pass to corner.corner()
        
        Returns:
        --------
        fig : matplotlib.figure.Figure
            The corner plot figure
        """
        if self.samples is None:
            raise ValueError("Must run fit() first")
        
        try:
            import corner
        except ImportError:
            raise ImportError("corner package is required for plotting. Install with: pip install corner")
        
        # Create labels for free parameters only
        label_map = {
            'log10_Atail': r'$\rm log_{10}(A_{\rm tail})$',
            'ctail': r'$c_{\rm tail}$',
            'ptail': r'$p_{\rm tail}$'
        }
        labels = [label_map[param] for param in self.free_params]
        
        # Create corner plot
        fig = corner.corner(self.samples, labels=labels, **kwargs)
        
        return fig