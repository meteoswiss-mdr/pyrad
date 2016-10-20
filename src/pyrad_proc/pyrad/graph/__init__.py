"""
==================================
Plotting (:mod:`pyrad.graph`)
==================================

.. currentmodule:: pyrad.graph

Functions to plot graphics.

Plots
=====

.. autosummary::
    :toctree: generated/

    plot_ppi
    plot_rhi
    plot_cappi
    plot_bscope
    plot_quantiles
    plot_timeseries
    plot_timeseries_comp
    plot_sun_retrieval_ts
    get_colobar_label

"""

from .plots import plot_ppi, plot_rhi, plot_cappi, plot_bscope, plot_quantiles
from .plots import plot_timeseries, plot_timeseries_comp
from .plots import plot_sun_retrieval_ts

from .plots import get_colobar_label

__all__ = [s for s in dir() if not s.startswith('_')]
