"""
======================================================
Dataset processing (:mod:`pyrad.proc`)
======================================================

.. currentmodule:: pyrad.proc

Initiate the dataset processing.

Auxiliary functions
====================

.. autosummary::
    :toctree: generated/

    get_process_func
    process_raw
    process_save_radar
    process_point_measurement

Echo classification and filtering
=================================

.. autosummary::
    :toctree: generated/

    process_echo_id
    process_echo_filter
    process_cdf
    process_filter_snr
    process_filter_visibility
    process_outlier_filter
    process_hydroclass

Phase processing and attenuation correction
===========================================

.. autosummary::
    :toctree: generated/

    process_estimate_phidp0
    process_correct_phidp0
    process_smooth_phidp_single_window
    process_smooth_phidp_double_window
    process_kdp_leastsquare_single_window
    process_kdp_leastsquare_double_window
    process_phidp_kdp_Maesaka
    process_phidp_kdp_lp
    process_attenuation

Monitoring, calibration and noise correction
============================================

.. autosummary::
    :toctree: generated/

    process_correct_bias
    process_correct_noise_rhohv
    process_rhohv_rain
    process_zdr_rain
    process_sun_hits
    process_selfconsistency_kdp_phidp
    process_selfconsistency_bias
    process_monitoring
    process_estimate_phidp0
    process_time_avg
    process_weighted_time_avg
    process_time_avg_flag
    process_colocated_gates
    process_intercomp
    process_intercomp_time_avg

Retrievals
==========

.. autosummary::
    :toctree: generated/

    process_signal_power
    process_snr
    process_l
    process_cdr
    process_rainrate
    process_wind_vel
    process_windshear

Trajectory functions
====================

.. autosummary::
    :toctree: generated/

    process_trajectory
    process_traj_atplane
    process_traj_antenna_pattern

COSMO data
==========

.. autosummary::
    :toctree: generated/

    process_cosmo
    process_cosmo_lookup_table
    process_cosmo_coord


"""

from .process_aux import get_process_func, process_raw, process_save_radar
from .process_aux import process_point_measurement
from .process_traj import process_trajectory, process_traj_atplane, \
    process_traj_antenna_pattern

from .process_echoclass import process_echo_id, process_echo_filter
from .process_echoclass import process_filter_snr, process_filter_visibility
from .process_echoclass import process_outlier_filter, process_hydroclass
from .process_echoclass import process_cdf

from .process_phase import process_correct_phidp0
from .process_phase import process_smooth_phidp_single_window
from .process_phase import process_smooth_phidp_double_window
from .process_phase import process_kdp_leastsquare_single_window
from .process_phase import process_kdp_leastsquare_double_window
from .process_phase import process_phidp_kdp_lp, process_phidp_kdp_Maesaka
from .process_phase import process_attenuation

from .process_calib import process_correct_bias, process_correct_noise_rhohv
from .process_calib import process_rhohv_rain, process_zdr_rain
from .process_calib import process_estimate_phidp0
from .process_calib import process_selfconsistency_kdp_phidp
from .process_calib import process_selfconsistency_bias
from .process_calib import process_monitoring
from .process_calib import process_time_avg, process_weighted_time_avg
from .process_calib import process_time_avg_flag
from .process_calib import process_colocated_gates, process_intercomp
from .process_calib import process_intercomp_time_avg
from .process_calib import process_sun_hits

from .process_retrieve import process_signal_power, process_snr
from .process_retrieve import process_l, process_cdr
from .process_retrieve import process_rainrate, process_wind_vel
from .process_retrieve import process_windshear

from .process_cosmo import process_cosmo, process_cosmo_lookup_table
from .process_cosmo import process_cosmo_coord

__all__ = [s for s in dir() if not s.startswith('_')]
