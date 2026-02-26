"""
gwtails: Post-merger late-time tail analysis for eccentric binary black hole mergers.
"""

__version__ = "0.1.0"
__author__ = "Tousif Islam"

from .tails import PostMergerAmplitudeFit
from .auxiliary import (
    compute_tail_exponent,
    compute_frequency,
    find_tail_start,
    calculate_chi_square,
    reduced_chi2,
)
from .bayesian_tails import TailFitMCMC

__all__ = [
    "PostMergerAmplitudeFit",
    "TailFitMCMC",
    "compute_tail_exponent",
    "compute_frequency",
    "find_tail_start",
    "calculate_chi_square",
    "reduced_chi2",
]
