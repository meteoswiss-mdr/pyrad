"""
==================================
Input and output (:mod:`pyrad.io`)
==================================

.. currentmodule:: pyrad.io

Functions to read and write data and configuration files.

Reading configuration files
===========================

.. autosummary::
    :toctree: generated/

    read_config

Reading radar data
==================

.. autosummary::
    :toctree: generated/

    get_data

Reading cosmo data
==================

.. autosummary::
    :toctree: generated/

    cosmo2radar_data
    cosmo2radar_coord
    hzt2radar_data
    hzt2radar_coord
    get_cosmo_fields
    get_iso0_field
    read_cosmo_data
    read_cosmo_coord
    read_hzt_data

Reading other data
==================

.. autosummary::
    :toctree: generated/

    read_last_state
    read_status
    read_rad4alp_cosmo
    read_rad4alp_vis
    read_excess_gates
    read_colocated_gates
    read_colocated_data
    read_timeseries
    read_ts_cum
    read_monitoring_ts
    read_intercomp_scores_ts
    get_sensor_data
    read_smn
    read_smn2
    read_disdro_scattering
    read_sun_hits
    read_sun_hits_multiple_days
    read_sun_retrieval
    read_solar_flux
    read_selfconsistency
    read_antenna_pattern
    read_lightning
    read_lightning_traj
    read_trt_data
    read_trt_traj_data
    read_rhi_profile
    read_histogram
    read_quantiles
    read_profile_ts
    read_histogram_ts
    read_quantiles_ts
    read_ml_ts

Writing data
==================

.. autosummary::
    :toctree: generated/

    send_msg
    write_alarm_msg
    write_last_state
    write_smn
    write_colocated_gates
    write_colocated_data
    write_colocated_data_time_avg
    write_timeseries
    write_ts_polar_data
    write_ts_cum
    write_monitoring_ts
    write_excess_gates
    write_intercomp_scores_ts
    write_sun_hits
    write_sun_retrieval
    write_cdf
    write_rhi_profile
    write_field_coverage
    write_histogram
    write_quantiles


Auxiliary functions
===================

.. autosummary::
    :toctree: generated/

    map_hydro
    map_Doppler
    get_save_dir
    make_filename
    get_datetime
    get_datasetfields
    get_file_list
    get_trtfile_list
    get_datatype_fields
    get_field_unit
    get_fieldname_pyart
    get_fieldname_cosmo
    generate_field_name_str
    find_raw_cosmo_file
    find_hzt_file
    add_field
    interpol_field
    get_new_rainbow_file_name
    _get_datetime

Trajectory
==========

.. autosummary::
    :toctree: generated/

    Trajectory

TimeSeries
==========

.. autosummary::
    :toctree: generated/

    TimeSeries

"""

from .config import read_config

from .read_data_radar import get_data, add_field, interpol_field

from .read_data_cosmo import read_cosmo_data, read_cosmo_coord
from .read_data_cosmo import cosmo2radar_data, cosmo2radar_coord
from .read_data_cosmo import get_cosmo_fields

from .read_data_hzt import read_hzt_data, hzt2radar_data, hzt2radar_coord
from .read_data_hzt import get_iso0_field

from .read_data_other import read_status, read_rad4alp_cosmo, read_rad4alp_vis
from .read_data_other import read_timeseries, read_monitoring_ts, read_ts_cum
from .read_data_other import read_intercomp_scores_ts, read_quantiles
from .read_data_other import read_selfconsistency, read_colocated_gates
from .read_data_other import read_colocated_data, read_antenna_pattern
from .read_data_other import read_last_state, read_rhi_profile
from .read_data_other import read_excess_gates, read_histogram
from .read_data_other import read_profile_ts, read_histogram_ts
from .read_data_other import read_quantiles_ts, read_ml_ts

from .read_data_sensor import read_lightning, read_lightning_traj
from .read_data_sensor import get_sensor_data, read_smn, read_smn2
from .read_data_sensor import read_disdro_scattering, read_trt_data
from .read_data_sensor import read_trt_traj_data

from .read_data_sun import read_sun_hits_multiple_days, read_sun_hits
from .read_data_sun import read_sun_retrieval, read_solar_flux

from .write_data import write_smn, write_ts_polar_data, write_ts_cum
from .write_data import write_monitoring_ts, write_intercomp_scores_ts
from .write_data import write_sun_hits, write_sun_retrieval
from .write_data import write_colocated_gates, write_colocated_data
from .write_data import write_colocated_data_time_avg, write_cdf
from .write_data import write_rhi_profile, write_field_coverage
from .write_data import write_last_state, write_alarm_msg, send_msg
from .write_data import write_excess_gates, write_trt_cell_data
from .write_data import write_histogram, write_quantiles

from .io_aux import get_save_dir, make_filename, get_new_rainbow_file_name
from .io_aux import get_datetime, get_dataset_fields, map_hydro, map_Doppler
from .io_aux import get_file_list, get_trtfile_list, get_datatype_fields
from .io_aux import get_fieldname_pyart, get_field_unit, get_fieldname_cosmo
from .io_aux import generate_field_name_str, find_raw_cosmo_file
from .io_aux import find_hzt_file, _get_datetime

from .trajectory import Trajectory

from .timeseries import TimeSeries

__all__ = [s for s in dir() if not s.startswith('_')]
