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

    plot_surface
    plot_latitude_slice
    plot_longitude_slice
    plot_latlon_slice
    plot_ppi
    plot_ppi_map
    plot_rhi
    plot_bscope
    plot_time_range
    plot_rhi_profile
    plot_along_coord
    plot_field_coverage
    plot_density
    plot_cappi
    plot_quantiles
    plot_histogram
    plot_histogram2
    plot_antenna_pattern
    plot_timeseries
    plot_timeseries_comp
    plot_monitoring_ts
    plot_scatter_comp
    plot_intercomp_scores_ts
    plot_sun_hits
    plot_sun_retrieval_ts
    get_colobar_label

"""

from .plots import plot_ppi, plot_ppi_map, plot_rhi, plot_cappi, plot_bscope
from .plots import plot_histogram, plot_histogram2, plot_density, plot_scatter
from .plots import plot_timeseries, plot_timeseries_comp, plot_monitoring_ts
from .plots import plot_sun_hits, plot_sun_retrieval_ts, plot_antenna_pattern
from .plots import plot_intercomp_scores_ts, plot_scatter_comp, plot_quantiles
from .plots import plot_along_coord, plot_field_coverage
from .plots import plot_surface, plot_latitude_slice, plot_longitude_slice
from .plots import plot_latlon_slice, plot_time_range

from .plots import get_colobar_label

__all__ = [s for s in dir() if not s.startswith('_')]
