Quick Start
===========

This guide shows the basic usage of ``gwtails`` for analyzing post-merger
late-time tails in eccentric BBH mergers.

Least-squares fitting with PostMergerAmplitudeFit
--------------------------------------------------

The :class:`~gwtails.PostMergerAmplitudeFit` class performs least-squares fitting
of QNM ringdown, power-law tails, and cross-term components.

**From a data file:**

.. code-block:: python

   import gwtails

   fit = gwtails.PostMergerAmplitudeFit(
       filename="path/to/waveform_data.dat",
       qinput=1000,
       qnm_fit_window=[20, 70],
       tail_fit_window=[200, 1000],
   )

   # Print fitted parameters
   fit._print_all_fits()

   # Plot all fits
   fit._plot_all_fits()

**From arrays:**

.. code-block:: python

   import numpy as np
   import gwtails

   t = np.load("time.npy")
   h = np.load("strain.npy")

   fit = gwtails.PostMergerAmplitudeFit(
       t=t,
       f=h,
       qinput=1000,
       qnm_fit_window=[20, 70],
       tail_fit_window=[200, 1000],
       crossterm_fit_window=[30, 400],
   )

   # Plot individual fits
   fit._plot_qnm_fit()
   fit._plot_tail_fits()

Bayesian MCMC fitting with TailFitMCMC
----------------------------------------

The :class:`~gwtails.TailFitMCMC` class provides Bayesian parameter estimation
for the tail decay function :math:`A_{\mathrm{tail}} (t + c_{\mathrm{tail}})^{p_{\mathrm{tail}}}`.

.. code-block:: python

   import gwtails

   fitter = gwtails.TailFitMCMC(
       t=t,                          # post-merger time array
       A=abs(h),                     # post-merger amplitude
       t_tail_window=[1200, 8000],
       log_Atail_range=(0, 20),
       ctail_range=(0, 2500),
       ptail_range=(-15, -1),
       percentage_err_in_data=2.5,
   )

   # Run MCMC
   samples, percentiles, log_like = fitter.fit(
       n_walkers=32,
       n_steps=2000,
       burn_in=500,
   )

   # Get best-fit parameters
   params = fitter.get_best_fit_params()
   print(params)

   # Corner plot
   fitter.plot_corner()

**Fixing parameters:**

.. code-block:: python

   # Fix the tail exponent
   fitter = gwtails.TailFitMCMC(
       t=t,
       A=abs(h),
       t_tail_window=[1200, 8000],
       fixed_params={"ptail": -6},
   )

   samples, percentiles, log_like = fitter.fit()

Utility functions
-----------------

.. code-block:: python

   import gwtails

   # Compute tail exponent from time series
   p = gwtails.compute_tail_exponent(t, amplitude)

   # Compute orbital frequency
   freq = gwtails.compute_frequency(t, h, method="spline")

   # Find where the tail starts
   tail_start = gwtails.find_tail_start(t, h, freq_threshold=1e-3)
