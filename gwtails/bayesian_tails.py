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

class TailFitMCMC:
    """
    Bayesian fitting class for tail decay function: Atail / (t + c)**n
    """
    
    def __init__(self, t, A, t_tail_window, lsq_params=None, 
                 log_Atail_range=(-5, 15), c_range=(0, 2500), n_range=(1, 15)):
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
                [Atail, c, n] from least squares fit to use as initial guess
            log_Atail_range : tuple, optional
                (min, max) range for log10(Atail) prior, default (-5, 15)
            c_range : tuple, optional
                (min, max) range for c parameter prior, default (0, 2500)
            n_range : tuple, optional
                (min, max) range for n parameter prior, default (1, 15)
        """

        # Simulation time and post-merger amplitude
        self.t = np.array(t)
        self.A = np.array(A)

        # fit window
        self.t_tail_window = t_tail_window
        self.lsq_params = lsq_params
        
        # Prior ranges
        self.log_Atail_range = log_Atail_range
        self.c_range = c_range
        self.n_range = n_range
        
        # Extract data within window
        mask = (self.t >= t_tail_window[0]) & (self.t <= t_tail_window[1])
        self.t_fit = self.t[mask]
        self.A_fit = self.A[mask]
        
        # Initialize results as None
        self.samples = None
        self.percentiles = None
        self.sampler = None
    
    def tail_model(self, t, params):
        """Tail decay model with log-scale Atail."""
        log_Atail, c, n = params
        Atail = 10**log_Atail
        return Atail / (t + c)**n
    
    def log_likelihood(self, params):
        """Log likelihood function."""
        try:
            model = self.tail_model(self.t_fit, params)
            if np.any(~np.isfinite(model)) or np.any(model <= 0):
                return -np.inf
            
            # Use relative errors based on data magnitude
            sigma = 0.1 * self.A_fit  # 10% relative error
            chi2 = np.sum((self.A_fit - model)**2 / sigma**2)
            return -0.5 * chi2
        except:
            return -np.inf
    
    def log_prior(self, params):
        """Log prior function with user-defined bounds."""
        log_Atail, c, n = params
        
        # Use user-defined prior ranges
        if (self.log_Atail_range[0] < log_Atail < self.log_Atail_range[1] and      
            self.c_range[0] < c < self.c_range[1] and             
            self.n_range[0] < n < self.n_range[1]):                 
            return 0.0
        return -np.inf
    
    def log_probability(self, params):
        """Log posterior probability."""
        lp = self.log_prior(params)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(params)
    
    def _get_initial_guess(self):
        """Get initial guess for parameters."""
        if self.lsq_params is not None:
            Atail_lsq, c_lsq, n_lsq = self.lsq_params
            log_Atail_lsq = np.log10(abs(Atail_lsq))
            initial_guess = [log_Atail_lsq, c_lsq, n_lsq]
            print("..... Using least-squares parameters as initial guess")
        else:
            # Data-driven initial guess
            log_A_median = np.log10(1/np.median(self.A_fit))
            initial_guess = [log_A_median, 50, 5]
            print("..... Using data-driven initial guess")
        
        print(f".......... Initial guess: log_Atail={initial_guess[0]:.3f}, c={initial_guess[1]:.3f}, n={initial_guess[2]:.3f}")
        return initial_guess
    
    def _initialize_walkers(self, initial_guess, n_walkers):
        """Initialize MCMC walkers."""
        # Dimension of the fit function
        ndim = 3
        
        # Ensure initial guess is within bounds, if not, center it
        bounded_guess = initial_guess.copy()
        if not (self.log_Atail_range[0] < bounded_guess[0] < self.log_Atail_range[1]):
            bounded_guess[0] = (self.log_Atail_range[0] + self.log_Atail_range[1]) / 2
        if not (self.c_range[0] < bounded_guess[1] < self.c_range[1]):
            bounded_guess[1] = (self.c_range[0] + self.c_range[1]) / 2
        if not (self.n_range[0] < bounded_guess[2] < self.n_range[1]):
            bounded_guess[2] = (self.n_range[0] + self.n_range[1]) / 2
        
        # Adaptive scatter based on prior ranges (use 10% of range width)
        log_Atail_width = self.log_Atail_range[1] - self.log_Atail_range[0]
        c_width = self.c_range[1] - self.c_range[0]
        n_width = self.n_range[1] - self.n_range[0]
        
        scatter = np.array([
            0.1 * log_Atail_width,  # 10% of log_Atail range
            0.1 * c_width,          # 10% of c range  
            0.1 * n_width           # 10% of n range
        ])
        
        pos = []
        max_attempts = 1000  # Prevent infinite loops
        
        for i in range(n_walkers):
            attempts = 0
            while attempts < max_attempts:
                candidate = bounded_guess + scatter * np.random.randn(ndim)
                # Ensure within user-defined prior bounds
                if (self.log_Atail_range[0] < candidate[0] < self.log_Atail_range[1] and 
                    self.c_range[0] < candidate[1] < self.c_range[1] and 
                    self.n_range[0] < candidate[2] < self.n_range[1]):
                    pos.append(candidate)
                    break
                attempts += 1
            
            # If we couldn't find a valid candidate, use uniform sampling within bounds
            if attempts >= max_attempts:
                candidate = np.array([
                    np.random.uniform(self.log_Atail_range[0], self.log_Atail_range[1]),
                    np.random.uniform(self.c_range[0], self.c_range[1]), 
                    np.random.uniform(self.n_range[0], self.n_range[1])
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
        # Print initial message
        print("..... Fitting tail data with the function : Atail / (t + c)**n")
        
        # Get initial guess
        initial_guess = self._get_initial_guess()
        
        # Initialize walkers
        pos = self._initialize_walkers(initial_guess, n_walkers)
        
        # Run MCMC
        ndim = 3
        self.sampler = emcee.EnsembleSampler(
            n_walkers, ndim, 
            lambda params, *args: self.log_probability(params)
        )
        
        print("..... Running MCMC")
        self.sampler.run_mcmc(pos, n_steps, progress=True)
        
        # Get samples after burn-in
        self.samples = self.sampler.get_chain(discard=burn_in, flat=True)
        
        # Calculate percentiles
        self._calculate_percentiles()
        
        # Print completion message
        print("..... Done!")
        
        return self.samples, self.percentiles
    
    def _calculate_percentiles(self):
        """Calculate parameter percentiles from MCMC samples."""
        self.percentiles = {}
        labels = ['log_Atail', 'c', 'n']
        
        print("..... Bayesian Results:")
        for i, label in enumerate(labels):
            mcmc = np.percentile(self.samples[:, i], [16, 50, 84])
            q = np.diff(mcmc)
            self.percentiles[label] = {
                'median': mcmc[1],
                'upper_err': q[1], 
                'lower_err': q[0]
            }
            
            if label == 'log_Atail':
                # Report both log and linear scale
                Atail_mcmc = 10**mcmc
                Atail_q = np.diff(Atail_mcmc)
                self.percentiles['Atail'] = {
                    'median': Atail_mcmc[1],
                    'upper_err': Atail_q[1],
                    'lower_err': Atail_q[0]
                }
                # Use symmetric +/- format
                avg_err = (q[0] + q[1]) / 2
                avg_Atail_err = (Atail_q[0] + Atail_q[1]) / 2
                print(f".......... log_Atail = {mcmc[1]:.3f} +/- {avg_err:.3f}")
                print(f".......... Atail = {Atail_mcmc[1]:.3e} +/- {avg_Atail_err:.3e}")
            else:
                # Use symmetric +/- format
                avg_err = (q[0] + q[1]) / 2
                print(f".......... {label} = {mcmc[1]:.3f} +/- {avg_err:.3f}")
    
    def get_best_fit_params(self):
        """Get best-fit parameters from percentiles."""
        if self.percentiles is None:
            raise ValueError("Must run fit() first")
        
        return {
            'Atail': self.percentiles['Atail']['median'],
            'c': self.percentiles['c']['median'],
            'n': self.percentiles['n']['median']
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
        return best_params['Atail'] / (t_pred + best_params['c'])**best_params['n']
        