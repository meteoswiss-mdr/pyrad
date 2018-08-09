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
    process_roi
    process_grid
    process_qvp
    process_rqvp
    process_svp
    process_evp
    process_time_height

Echo classification and filtering
=================================

.. autosummary::
    :toctree: generated/

    process_echo_id
    process_birds_id
    process_clt_to_echo_id
    process_echo_filter
    process_cdf
    process_filter_snr
    process_filter_visibility
    process_outlier_filter
    process_hydroclass
    process_melting_layer
    process_filter_vel_diff
    process_zdr_column

Phase processing and attenuation correction
===========================================

.. autosummary::
    :toctree: generated/

    process_correct_phidp0
    process_smooth_phidp_single_window
    process_smooth_phidp_double_window
    process_kdp_leastsquare_single_window
    process_kdp_leastsquare_double_window
    process_phidp_kdp_Vulpiani
    process_phidp_kdp_Kalman
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
    process_zdr_precip
    process_zdr_snow
    process_estimate_phidp0
    process_sun_hits
    process_selfconsistency_kdp_phidp
    process_selfconsistency_bias
    process_occurrence
    process_occurrence_period
    process_monitoring
    process_gc_monitoring
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
    process_vol_refl
    process_bird_density

Doppler processing
==================

.. autosummary::
    :toctree: generated/

    process_dealias_fourdd
    process_dealias_region_based
    process_dealias_unwrap_phase
    process_wind_vel
    process_windshear
    process_vad

Time series functions
====================

.. autosummary::
    :toctree: generated/

    process_point_measurement
    process_qvp
    process_rqvp
    process_svp
    process_evp
    process_time_height

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
    process_hzt
    process_hzt_lookup_table
    process_hzt_coord


"""

from .process_aux import get_process_func, process_raw, process_save_radar
from .process_aux import process_grid, process_roi

from .process_timeseries import process_point_measurement, process_qvp
from .process_timeseries import process_rqvp, process_evp, process_svp
from .process_timeseries import process_time_height

from .process_traj import process_trajectory, process_traj_atplane
from .process_traj import process_traj_antenna_pattern

from .process_echoclass import process_echo_id, process_birds_id
from .process_echoclass import process_echo_filter, process_clt_to_echo_id
from .process_echoclass import process_filter_snr, process_filter_visibility
from .process_echoclass import process_outlier_filter, process_hydroclass
from .process_echoclass import process_cdf, process_melting_layer
from .process_echoclass import process_filter_vel_diff, process_zdr_column

from .process_phase import process_correct_phidp0
from .process_phase import process_smooth_phidp_single_window
from .process_phase import process_smooth_phidp_double_window
from .process_phase import process_kdp_leastsquare_single_window
from .process_phase import process_kdp_leastsquare_double_window
from .process_phase import process_phidp_kdp_lp, process_phidp_kdp_Maesaka
from .process_phase import process_phidp_kdp_Vulpiani
from .process_phase import process_phidp_kdp_Kalman
from .process_phase import process_attenuation

from .process_calib import process_correct_bias, process_correct_noise_rhohv
from .process_calib import process_occurrence, process_occurrence_period
from .process_calib import process_gc_monitoring, process_sun_hits

from .process_intercomp import process_time_avg, process_weighted_time_avg
from .process_intercomp import process_time_avg_flag
from .process_intercomp import process_colocated_gates
from .process_intercomp import process_intercomp, process_intercomp_time_avg

from .process_monitoring import process_zdr_snow, process_zdr_precip
from .process_monitoring import process_estimate_phidp0, process_rhohv_rain
from .process_monitoring import process_selfconsistency_kdp_phidp
from .process_monitoring import process_selfconsistency_bias
from .process_monitoring import process_monitoring

from .process_retrieve import process_signal_power, process_snr
from .process_retrieve import process_l, process_cdr, process_bird_density
from .process_retrieve import process_rainrate, process_vol_refl

from .process_Doppler import process_wind_vel, process_windshear
from .process_Doppler import process_dealias_fourdd
from .process_Doppler import process_dealias_region_based
from .process_Doppler import process_dealias_unwrap_phase
from .process_Doppler import process_vad

from .process_cosmo import process_cosmo, process_cosmo_lookup_table
from .process_cosmo import process_cosmo_coord, process_hzt
from .process_cosmo import process_hzt_lookup_table, process_hzt_coord

__all__ = [s for s in dir() if not s.startswith('_')]
