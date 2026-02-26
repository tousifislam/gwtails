gwtails
=======

.. image:: ../../gwtails_logo.png
   :align: center
   :width: 400px

**gwtails** is a Python package to analyze post-merger late-time tails in
eccentric binary black hole (BBH) mergers.

This code is based on the papers:

* `Phenomenology and origin of late-time tails in eccentric binary black hole mergers <https://arxiv.org/abs/2407.04682>`_ — Islam, Faggioli, Khanna, Field, van de Meent, Buonanno (2024)
* `Bayesian analysis of late-time tails in spin-aligned eccentric binary black hole mergers <https://arxiv.org/abs/2511.21898>`_ — Islam, Faggioli, Khanna (2025)

Features
--------

* **Least-squares fitting** of QNM, tail, and cross-term components via :class:`~gwtails.PostMergerAmplitudeFit`
* **Bayesian MCMC fitting** with flexible parameter fixing via :class:`~gwtails.TailFitMCMC`
* **Auxiliary utilities** for computing tail exponents, frequencies, and finding tail start times
* Built-in visualization for fit diagnostics and corner plots

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   quickstart
   api/index
   notebooks/index
   citation
   contributing
