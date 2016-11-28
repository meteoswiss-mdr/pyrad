"""
pyrad.proc.process_calib
===========================

Functions for monitoring data quality and correct bias and noise effects

.. autosummary::
    :toctree: generated/

    process_correct_bias
    process_correct_noise_rhohv
    process_selfconsistency_kdp_phidp
    process_selfconsistency_bias
    process_estimate_phidp0
    process_rhohv_rain
    process_zdr_rain
    process_monitoring
    process_time_avg
    process_weighted_time_avg
    process_time_avg_flag
    process_colocated_gates
    process_intercomp
    process_sun_hits

"""

from copy import deepcopy
from warnings import warn
import datetime

import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename
from ..io.read_data_other import read_selfconsistency, read_colocated_gates
from ..io.read_data_other import read_sun_hits_multiple_days, read_solar_flux
from ..io.read_data_other import read_colocated_data
from ..io.read_data_radar import interpol_field

from ..util.radar_utils import get_closest_solar_flux, get_histogram_bins
from ..util.radar_utils import time_avg_range, get_range_bins_to_avg
from ..util.radar_utils import find_ray_index, find_rng_index


def process_correct_bias(procstatus, dscfg, radar_list=None):
    """
    Corrects a bias on the data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type to correct for bias
        bias : float. Dataset keyword
            The bias to be corrected [dB]. Default 0
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        break
    field_name = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if field_name not in radar.fields:
        warn('Unable to correct for bias field ' + field_name +
             '. Field not available')
        return None, None

    bias = 0.
    if 'bias' in dscfg:
        bias = dscfg['bias']

    corrected_field = pyart.correct.correct_bias(
        radar, bias=bias, field_name=field_name)

    if field_name.startswith('corrected_'):
        new_field_name = field_name
    else:
        new_field_name = 'corrected_'+field_name

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field(new_field_name, corrected_field)

    return new_dataset, ind_rad


def process_correct_noise_rhohv(procstatus, dscfg, radar_list=None):
    """
    identifies echoes as 0: No data, 1: Noise, 2: Clutter,
    3: Precipitation

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The data types used in the correction
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if datatype == 'uRhoHV':
            urhohv = 'uncorrected_cross_correlation_ratio'
        if datatype == 'SNRh':
            snr = 'signal_to_noise_ratio_hh'
        if datatype == 'ZDR':
            zdr = 'differential_reflectivity'
        if datatype == 'ZDRc':
            zdr = 'corrected_differential_reflectivity'
        if datatype == 'Nh':
            nh = 'noisedBZ_hh'
        if datatype == 'Nv':
            nv = 'noisedBZ_vv'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((urhohv not in radar.fields) or
            (snr not in radar.fields) or
            (zdr not in radar.fields) or
            (nh not in radar.fields) or
            (nv not in radar.fields)):
        warn('Unable to correct RhoHV field for noise. Missing fields')
        return None, None

    rhohv = pyart.correct.correct_noise_rhohv(
        radar, urhohv_field=urhohv, snr_field=snr, zdr_field=zdr,
        nh_field=nh, nv_field=nv, rhohv_field='cross_correlation_ratio')

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field('cross_correlation_ratio', rhohv)

    return new_dataset, ind_rad


def process_selfconsistency_kdp_phidp(procstatus, dscfg, radar_list=None):
    """
    Computes specific differential phase and differential phase in rain using
    the selfconsistency between Zdr, Zh and KDP

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of strings. Dataset keyword
            The input data types
        rsmooth : float. Dataset keyword
            length of the smoothing window [m]. Default 1000.
        min_rhohv : float. Dataset keyword
            minimum valid RhoHV. Default 0.92
        max_phidp : float. Dataset keyword
            maximum valid PhiDP [deg]. Default 20.
        ml_thickness : float. Dataset keyword
            assumed melting layer thickness [m]. Default 700.
        fzl : float. Dataset keyword
            The default freezing level height. It will be used if no
            temperature field name is specified or the temperature field is
            not in the radar object. Default 2000.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    temp = None
    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if datatype == 'dBZc':
            refl = 'corrected_reflectivity'
        if datatype == 'dBZ':
            refl = 'reflectivity'
        if datatype == 'ZDRc':
            zdr = 'corrected_differential_reflectivity'
        if datatype == 'ZDR':
            zdr = 'differential_reflectivity'
        if datatype == 'PhiDPc':
            phidp = 'corrected_differential_phase'
        if datatype == 'PhiDP':
            phidp = 'differential_phase'
        if datatype == 'TEMP':
            temp = 'temperature'
        if datatype == 'RhoHV':
            rhohv = 'cross_correlation_ratio'
        if datatype == 'RhoHVc':
            rhohv = 'corrected_cross_correlation_ratio'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl not in radar.fields) or
            (zdr not in radar.fields) or
            (phidp not in radar.fields) or
            (rhohv not in radar.fields)):
        warn('Unable to estimate reflectivity bias using selfconsistency. ' +
             'Missing data')
        return None, None

    if (temp is not None) and (temp not in radar.fields):
        warn('COSMO temperature field not available. ' +
             'Using fixed freezing level height')

    fzl = None
    if (temp is None) or (temp not in radar.fields):
        if 'fzl' in dscfg:
            fzl = dscfg['fzl']
        else:
            fzl = 2000.
            warn('Freezing level height not defined. Using default ' +
                 str(fzl)+' m')

    fname = (
        dscfg['configpath'] + 'selfconsistency/' +
        'selfconsistency_zdr_zhkdp_Xband_temp10_elev000_mu05.txt')
    zdr_kdpzh_table = read_selfconsistency(fname)

    # default values
    rsmooth = 1000.
    min_rhohv = 0.92
    max_phidp = 20.
    ml_thickness = 700.

    # get user defined values
    if 'rsmooth' in dscfg:
        rsmooth = dscfg['rsmooth']
    if 'min_rhohv' in dscfg:
        min_rhohv = dscfg['min_rhohv']
    if 'max_phidp' in dscfg:
        max_phidp = dscfg['max_phidp']
    if 'ml_thickness' in dscfg:
        ml_thickness = dscfg['ml_thickness']

    kdpsim_field = 'specific_differential_phase'
    phidpsim_field = 'differential_phase'
    r_res = radar.range['data'][1]-radar.range['data'][0]
    smooth_wind_len = int(rsmooth/r_res)

    kdpsim, phidpsim = pyart.correct.selfconsistency_kdp_phidp(
        radar, zdr_kdpzh_table, min_rhohv=min_rhohv, max_phidp=max_phidp,
        smooth_wind_len=smooth_wind_len, doc=15, fzl=fzl,
        thickness=ml_thickness, refl_field=refl, phidp_field=phidp,
        zdr_field=zdr, temp_field=temp, rhohv_field=rhohv,
        kdpsim_field=kdpsim_field, phidpsim_field=phidpsim_field)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    new_dataset.add_field(kdpsim_field, kdpsim)
    new_dataset.add_field(phidpsim_field, phidpsim)

    return new_dataset, ind_rad


def process_selfconsistency_bias(procstatus, dscfg, radar_list=None):
    """
    Estimates the reflectivity bias by means of the selfconsistency
    algorithm by Gourley

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rsmooth : float. Dataset keyword
            length of the smoothing window [m]. Default 1000.
        min_rhohv : float. Dataset keyword
            minimum valid RhoHV. Default 0.92
        max_phidp : float. Dataset keyword
            maximum valid PhiDP [deg]. Default 20.
        rcell : float. Dataset keyword
            length of continuous precipitation to consider the precipitation
            cell a valid phidp segment [m]. Default 1000.
        dphidp_min : float. Dataset keyword
            minimum phase shift [deg]. Default 2.
        dphidp_max : float. Dataset keyword
            maximum phase shift [deg]. Default 16.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    temp = None
    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if datatype == 'dBZc':
            refl = 'corrected_reflectivity'
        if datatype == 'dBZ':
            refl = 'reflectivity'
        if datatype == 'ZDRc':
            zdr = 'corrected_differential_reflectivity'
        if datatype == 'ZDR':
            zdr = 'differential_reflectivity'
        if datatype == 'PhiDPc':
            phidp = 'corrected_differential_phase'
        if datatype == 'PhiDP':
            phidp = 'differential_phase'
        if datatype == 'TEMP':
            temp = 'temperature'
        if datatype == 'RhoHV':
            rhohv = 'cross_correlation_ratio'
        if datatype == 'RhoHVc':
            rhohv = 'corrected_cross_correlation_ratio'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl not in radar.fields) or
            (zdr not in radar.fields) or
            (phidp not in radar.fields) or
            (rhohv not in radar.fields)):
        warn('Unable to estimate reflectivity bias using selfconsistency. ' +
             'Missing data')
        return None, None

    if (temp is not None) and (temp not in radar.fields):
        warn('COSMO temperature field not available. ' +
             'Using fixed freezing level height')

    fzl = None
    if (temp is None) or (temp not in radar.fields):
        if 'fzl' in dscfg:
            fzl = dscfg['fzl']
        else:
            fzl = 2000.
            warn('Freezing level height not defined. Using default ' +
                 str(fzl)+' m')

    # default values
    rsmooth = 1000.
    min_rhohv = 0.92
    max_phidp = 20.
    ml_thickness = 700.
    rcell = 1000.
    dphidp_min = 2.
    dphidp_max = 16.

    # get user defined values
    if 'rsmooth' in dscfg:
        rsmooth = dscfg['rsmooth']
    if 'min_rhohv' in dscfg:
        min_rhohv = dscfg['min_rhohv']
    if 'max_phidp' in dscfg:
        max_phidp = dscfg['max_phidp']
    if 'ml_thickness' in dscfg:
        ml_thickness = dscfg['ml_thickness']
    if 'rcell' in dscfg:
        rcell = dscfg['rcell']
    if 'dphidp_min' in dscfg:
        dphidp_min = dscfg['dphidp_min']
    if 'dphidp_max' in dscfg:
        dphidp_max = dscfg['dphidp_max']

    fname = (
        dscfg['configpath'] + 'selfconsistency/' +
        'selfconsistency_zdr_zhkdp_Xband_temp10_elev000_mu05.txt')
    zdr_kdpzh_table = read_selfconsistency(fname)

    r_res = radar.range['data'][1]-radar.range['data'][0]
    smooth_wind_len = int(rsmooth/r_res)
    min_rcons = int(rcell/r_res)

    step = None
    if 'step' in dscfg:
        step = dscfg['step']

    refl_bias = pyart.correct.selfconsistency_bias(
        radar, zdr_kdpzh_table, min_rhohv=min_rhohv, max_phidp=max_phidp,
        smooth_wind_len=smooth_wind_len, doc=15, fzl=fzl,
        thickness=ml_thickness, min_rcons=min_rcons, dphidp_min=dphidp_min,
        dphidp_max=dphidp_max, refl_field=refl, phidp_field=phidp,
        zdr_field=zdr, temp_field=temp, rhohv_field=rhohv)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    new_dataset.add_field('reflectivity_bias', refl_bias)

    return new_dataset, ind_rad


def process_estimate_phidp0(procstatus, dscfg, radar_list=None):
    """
    estimates the system differential phase offset at each ray

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rmin : float. Dataset keyword
            The minimum range where to look for valid data [m]
        rmax : float. Dataset keyword
            The maximum range where to look for valid data [m]
        rcell : float. Dataset keyword
            The length of a continuous cell to consider it valid precip [m]
        Zmin : float. Dataset keyword
            The minimum reflectivity [dBZ]
        Zmax : float. Dataset keyword
            The maximum reflectivity [dBZ]
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'PhiDP':
            psidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            psidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            psidp_field = 'uncorrected_differential_phase'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if (refl_field not in radar.fields) or (psidp_field not in radar.fields):
        warn('Unable to estimate PhiDP system offset. Missing data')
        return None, None

    ind_rmin = np.where(radar.range['data'] > dscfg['rmin'])[0][0]
    ind_rmax = np.where(radar.range['data'] < dscfg['rmax'])[0][-1]
    r_res = radar.range['data'][1]-radar.range['data'][0]
    min_rcons = int(dscfg['rcell']/r_res)

    step = None
    if 'step' in dscfg:
        step = dscfg['step']

    phidp0, first_gates = pyart.correct.det_sys_phase_ray(
        radar, ind_rmin=ind_rmin, ind_rmax=ind_rmax, min_rcons=min_rcons,
        zmin=dscfg['Zmin'], zmax=dscfg['Zmax'], phidp_field=psidp_field,
        refl_field=refl_field)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    new_dataset.add_field('system_differential_phase', phidp0)
    new_dataset.add_field('first_gate_differential_phase', first_gates)

    return new_dataset, ind_rad


def process_rhohv_rain(procstatus, dscfg, radar_list=None):
    """
    Keeps only suitable data to evaluate the 80 percentile of RhoHV in rain

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rmin : float. Dataset keyword
            minimum range where to look for rain [m]. Default 1000.
        rmax : float. Dataset keyword
            maximum range where to look for rain [m]. Default 50000.
        Zmin : float. Dataset keyword
            minimum reflectivity to consider the bin as precipitation [dBZ].
            Default 20.
        Zmax : float. Dataset keyword
            maximum reflectivity to consider the bin as precipitation [dBZ]
            Default 40.
        ml_thickness : float. Dataset keyword
            assumed thickness of the melting layer. Default 700.
        fzl : float. Dataset keyword
            The default freezing level height. It will be used if no
            temperature field name is specified or the temperature field is
            not in the radar object. Default 2000.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    temp_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if datatype == 'RhoHV':
            rhohv_field = 'cross_correlation_ratio'
        if datatype == 'RhoHVc':
            rhohv_field = 'corrected_cross_correlation_ratio'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'TEMP':
            temp_field = 'temperature'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl_field not in radar.fields) or
            (rhohv_field not in radar.fields)):
        warn('Unable to estimate RhoHV in rain. Missing data')
        return None, None

    if (temp_field is not None) and (temp_field not in radar.fields):
        warn('COSMO temperature field not available. ' +
             'Using fixed freezing level height')

    fzl = None
    if (temp_field is None) or (temp_field not in radar.fields):
        if 'fzl' in dscfg:
            fzl = dscfg['fzl']
        else:
            fzl = 2000.
            warn('Freezing level height not defined. Using default ' +
                 str(fzl)+' m')

    # default values
    rmin = 1000.
    rmax = 50000.
    zmin = 20.
    zmax = 40.
    thickness = 700.

    # user defined values
    if 'rmin' in dscfg:
        rmin = dscfg['rmin']
    if 'rmax' in dscfg:
        rmax = dscfg['rmax']
    if 'Zmin' in dscfg:
        zmin = dscfg['Zmin']
    if 'Zmax' in dscfg:
        zmax = dscfg['Zmax']
    if 'ml_thickness' in dscfg:
        thickness = dscfg['ml_thickness']

    ind_rmin = np.where(radar.range['data'] > rmin)[0][0]
    ind_rmax = np.where(radar.range['data'] < rmax)[0][-1]

    rhohv_rain = pyart.correct.est_rhohv_rain(
        radar, ind_rmin=ind_rmin, ind_rmax=ind_rmax, zmin=zmin,
        zmax=zmax, thickness=thickness, doc=15, fzl=fzl,
        rhohv_field=rhohv_field, temp_field=temp_field, refl_field=refl_field)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    new_dataset.add_field('cross_correlation_ratio_in_rain', rhohv_rain)

    return new_dataset, ind_rad


def process_zdr_rain(procstatus, dscfg, radar_list=None):
    """
    Keeps only suitable data to evaluate the differential reflectivity in
    moderate rain

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rmin : float. Dataset keyword
            minimum range where to look for rain [m]. Default 1000.
        rmax : float. Dataset keyword
            maximum range where to look for rain [m]. Default 50000.
        Zmin : float. Dataset keyword
            minimum reflectivity to consider the bin as precipitation [dBZ].
            Default 20.
        Zmax : float. Dataset keyword
            maximum reflectivity to consider the bin as precipitation [dBZ]
            Default 40.
        rhohvmin : float. Dataset keyword
            minimum RhoHV to consider the bin as precipitation
            Default 0.97
        phidpmax : float. Dataset keyword
            maximum PhiDP to consider the bin as precipitation [deg]
            Default 10.
        elmax : float. Dataset keyword
            maximum elevation angle where to look for precipitation [deg]
            Default 20.
        ml_thickness : float. Dataset keyword
            assumed thickness of the melting layer. Default 700.
        fzl : float. Dataset keyword
            The default freezing level height. It will be used if no
            temperature field name is specified or the temperature field is
            not in the radar object. Default 2000.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    temp_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if datatype == 'ZDR':
            zdr_field = 'differential_reflectivity'
        if datatype == 'ZDRc':
            zdr_field = 'corrected_differential_reflectivity'
        if datatype == 'PhiDP':
            phidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            phidp_field = 'corrected_differential_phase'
        if datatype == 'RhoHV':
            rhohv_field = 'cross_correlation_ratio'
        if datatype == 'RhoHVc':
            rhohv_field = 'corrected_cross_correlation_ratio'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'TEMP':
            temp_field = 'temperature'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl_field not in radar.fields) or
            (rhohv_field not in radar.fields)):
        warn('Unable to estimate ZDR in rain. Missing data')
        return None, None

    if (temp_field is not None) and (temp_field not in radar.fields):
        warn('COSMO temperature field not available. ' +
             'Using fixed freezing level height')

    fzl = None
    if (temp_field is None) or (temp_field not in radar.fields):
        if 'fzl' in dscfg:
            fzl = dscfg['fzl']
        else:
            fzl = 2000.
            warn('Freezing level height not defined. Using default ' +
                 str(fzl)+' m')

    # default values
    rmin = 1000.
    rmax = 50000.
    zmin = 20.
    zmax = 40.
    rhohvmin = 0.97
    phidpmax = 10.
    elmax = 20.
    thickness = 700.

    # user defined values
    if 'rmin' in dscfg:
        rmin = dscfg['rmin']
    if 'rmax' in dscfg:
        rmax = dscfg['rmax']
    if 'Zmin' in dscfg:
        zmin = dscfg['Zmin']
    if 'Zmax' in dscfg:
        zmax = dscfg['Zmax']
    if 'RhoHVmin' in dscfg:
        rhohvmin = dscfg['RhoHVmin']
    if 'PhiDPmax' in dscfg:
        phidpmax = dscfg['PhiDPmax']
    if 'elmax' in dscfg:
        elmax = dscfg['elmax']
    if 'ml_thickness' in dscfg:
        thickness = dscfg['ml_thickness']

    ind_rmin = np.where(radar.range['data'] > rmin)[0][0]
    ind_rmax = np.where(radar.range['data'] < rmax)[0][-1]

    zdr_rain = pyart.correct.est_zdr_rain(
        radar, ind_rmin=ind_rmin, ind_rmax=ind_rmax, zmin=zmin,
        zmax=zmax, rhohvmin=rhohvmin, phidpmax=phidpmax, elmax=elmax,
        thickness=thickness, doc=15, fzl=fzl, zdr_field=zdr_field,
        rhohv_field=rhohv_field, phidp_field=phidp_field,
        temp_field=temp_field, refl_field=refl_field)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    new_dataset.add_field('differential_reflectivity_in_rain', zdr_rain)

    return new_dataset, ind_rad


def process_monitoring(procstatus, dscfg, radar_list=None):
    """
    computes monitoring statistics

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        step : float. Dataset keyword
            The width of the histogram bin. Default is None. In that case the
            default step in function get_histogram_bins is used
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object containing histogram data
    ind_rad : int
        radar index

    """

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        for datatypedescr in dscfg['datatype']:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            field_name = get_fieldname_pyart(datatype)
            break
        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn(field_name+' not available.')
            return None, None

        step = None
        if 'step' in dscfg:
            step = dscfg['step']

        bins = get_histogram_bins(field_name, step=step)
        nbins = len(bins)-1

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.range['data'] = bins[0:-1]
        radar_aux.ngates = nbins

        field_dict = pyart.config.get_metadata(field_name)
        field_dict['data'] = np.ma.zeros((radar.nrays, nbins), dtype=int)

        for ray in range(radar.nrays):
            field_dict['data'][ray, :], bin_edges = np.histogram(
                radar.fields[field_name]['data'][ray, :].compressed(),
                bins=bins)

        radar_aux.add_field(field_name, field_dict)

        # keep histogram in Memory or add to existing histogram
        if dscfg['initialized'] == 0:

            dscfg['global_data'] = radar_aux
            dscfg['initialized'] = 1
        else:
            field_interp = interpol_field(
                dscfg['global_data'], radar_aux, field_name, fill_value=0)
            dscfg['global_data'].fields[field_name]['data'] += (
                field_interp['data'].filled(fill_value=0)).astype('int64')

        dataset = dict()
        dataset.update({'hist_obj': radar_aux})
        dataset.update({'hist_type': 'instant'})

        return dataset, ind_rad

    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None

        for datatypedescr in dscfg['datatype']:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            field_name = get_fieldname_pyart(datatype)
            break
        ind_rad = int(radarnr[5:8])-1

        dataset = dict()
        dataset.update({'hist_obj': dscfg['global_data']})
        dataset.update({'hist_type': 'cumulative'})

        return dataset, ind_rad


def process_time_avg(procstatus, dscfg, radar_list=None):
    """
    computes the temporal mean of a field

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        period : float. Dataset keyword
            the period to average [s]. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
        lin_trans: int. Dataset keyword
            If 1 apply linear transformation before averaging
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """
    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        break
    ind_rad = int(radarnr[5:8])-1

    lin_trans = 0
    if 'lin_trans' in dscfg:
        lin_trans = dscfg['lin_trans']

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn(field_name+' not available.')
            return None, None

        period = 3600.
        if 'period' in dscfg:
            period = dscfg['period']

        field = deepcopy(radar.fields[field_name])
        if lin_trans:
            field['data'] = np.ma.power(10., 0.1*field['data'])

        field['data'] = field['data'].filled(fill_value=0.)
        field['data'] = np.ma.asarray(field['data'])

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field(field_name, field)
        npoints_dict = pyart.config.get_metadata('number_of_samples')
        npoints_dict['data'] = np.ma.ones(
            (radar.nrays, radar.ngates), dtype=int)
        radar_aux.add_field('number_of_samples', npoints_dict)

        # first volume: initialize start and end time of averaging
        if dscfg['initialized'] == 0:
            start_average = 0.  # seconds from midnight
            if 'start_average' in dscfg:
                start_average = dscfg['start_average']

            date_00 = dscfg['timeinfo'].replace(
                hour=0, minute=0, second=0, microsecond=0)

            avg_par = dict()
            avg_par.update(
                {'starttime': date_00+datetime.timedelta(
                    seconds=start_average)})
            avg_par.update(
                {'endtime': avg_par['starttime']+datetime.timedelta(
                    seconds=period)})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        # no radar object in global data: create it
        if 'radar_obj' not in dscfg['global_data']:
            # get start and stop times of new radar object
            (dscfg['global_data']['starttime'],
             dscfg['global_data']['endtime']) = (
                time_avg_range(
                    dscfg['timeinfo'], dscfg['global_data']['starttime'],
                    dscfg['global_data']['endtime'], period))

            # check if volume time older than starttime
            if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
                dscfg['global_data'].update({'radar_obj': radar_aux})

            return None, None

        # still accumulating: add field to global field
        if dscfg['timeinfo'] < dscfg['global_data']['endtime']:
            field_interp = interpol_field(
                dscfg['global_data']['radar_obj'], radar_aux, field_name)
            npoints_interp = interpol_field(
                dscfg['global_data']['radar_obj'], radar_aux,
                'number_of_samples')
            dscfg['global_data']['radar_obj'].fields[field_name]['data'] += (
                field_interp['data'].filled(fill_value=0))
            dscfg['global_data']['radar_obj'].fields[
                'number_of_samples']['data'] += (
                npoints_interp['data'].filled(fill_value=0)).astype('int')

            return None, None

        # we have reached the end of the accumulation period: do the averaging
        # and start a new object
        dscfg['global_data']['radar_obj'].fields[field_name]['data'] /= (
            dscfg['global_data']['radar_obj'].fields[
                'number_of_samples']['data'])
        if lin_trans:
            dscfg['global_data']['radar_obj'].fields[field_name]['data'] = (
                10.*np.ma.log10(
                    dscfg['global_data']['radar_obj'].fields[
                        field_name]['data']))

        new_dataset = deepcopy(dscfg['global_data']['radar_obj'])

        dscfg['global_data']['starttime'] += datetime.timedelta(
            seconds=period)
        dscfg['global_data']['endtime'] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg['global_data'].pop('radar_obj', None)

        # get start and stop times of new radar object
        dscfg['global_data']['starttime'], dscfg['global_data']['endtime'] = (
            time_avg_range(
                dscfg['timeinfo'], dscfg['global_data']['starttime'],
                dscfg['global_data']['endtime'], period))

        # check if volume time older than starttime
        if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
            dscfg['global_data'].update({'radar_obj': radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'radar_obj' not in dscfg['global_data']:
            return None, None

        (dscfg['global_data']['radar_obj'].fields[field_name][
            'data']) /= (
            dscfg['global_data']['radar_obj'].fields[
                'number_of_samples']['data'])
        if lin_trans:
            dscfg['global_data']['radar_obj'].fields[field_name]['data'] = (
                10.*np.ma.log10(
                    dscfg['global_data']['radar_obj'].fields[
                        field_name]['data']))

        new_dataset = deepcopy(dscfg['global_data']['radar_obj'])

        return new_dataset, ind_rad


def process_weighted_time_avg(procstatus, dscfg, radar_list=None):
    """
    computes the temporal mean of a field weighted by the reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        period : float. Dataset keyword
            the period to average [s]. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """
    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if ((datatype == 'dBZ') or (datatype == 'dBZc') or
                (datatype == 'dBuZ') or (datatype == 'dBZv') or
                (datatype == 'dBZvc') or (datatype == 'dBuZv')):
            refl_name = get_fieldname_pyart(datatype)
        else:
            field_name = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]
        if field_name not in radar.fields or refl_name not in radar.fields:
            warn('Unable to compute weighted average. Missing data')
            return None, None

        period = 3600.
        if 'period' in dscfg:
            period = dscfg['period']

        field = deepcopy(radar.fields[field_name])
        field['data'] = field['data'].filled(fill_value=0.)
        field['data'] = np.ma.asarray(field['data'])

        refl_field = deepcopy(radar.fields[refl_name])
        refl_field['data'] = np.ma.power(10., 0.1*refl_field['data'])
        refl_field['data'] = refl_field['data'].filled(fill_value=0.)
        refl_field['data'] = np.ma.asarray(refl_field['data'])

        field['data'] *= refl_field['data']

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field(field_name, field)
        radar_aux.add_field(refl_name, refl_field)

        # first volume: initialize start and end time of averaging
        if dscfg['initialized'] == 0:
            start_average = 0.  # seconds from midnight
            if 'start_average' in dscfg:
                start_average = dscfg['start_average']

            date_00 = dscfg['timeinfo'].replace(
                hour=0, minute=0, second=0, microsecond=0)

            avg_par = dict()
            avg_par.update(
                {'starttime': date_00+datetime.timedelta(
                    seconds=start_average)})
            avg_par.update(
                {'endtime': avg_par['starttime']+datetime.timedelta(
                    seconds=period)})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        # no radar object in global data: create it
        if 'radar_obj' not in dscfg['global_data']:
            # get start and stop times of new radar object
            (dscfg['global_data']['starttime'],
             dscfg['global_data']['endtime']) = (
                time_avg_range(
                    dscfg['timeinfo'], dscfg['global_data']['starttime'],
                    dscfg['global_data']['endtime'], period))

            # check if volume time older than starttime
            if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
                dscfg['global_data'].update({'radar_obj': radar_aux})

            return None, None

        # still accumulating: add field to global field
        if dscfg['timeinfo'] < dscfg['global_data']['endtime']:
            field_interp = interpol_field(
                dscfg['global_data']['radar_obj'], radar_aux, field_name)
            dscfg['global_data']['radar_obj'].fields[field_name]['data'] += (
                field_interp['data'].filled(fill_value=0))

            refl_interp = interpol_field(
                dscfg['global_data']['radar_obj'], radar_aux, refl_name)
            dscfg['global_data']['radar_obj'].fields[refl_name]['data'] += (
                refl_interp['data'].filled(fill_value=0))

            return None, None

        # we have reached the end of the accumulation period: do the averaging
        # and start a new object
        dscfg['global_data']['radar_obj'].fields[field_name]['data'] /= (
            dscfg['global_data']['radar_obj'].fields[refl_name]['data'])

        new_dataset = deepcopy(dscfg['global_data']['radar_obj'])

        dscfg['global_data']['starttime'] += datetime.timedelta(
            seconds=period)
        dscfg['global_data']['endtime'] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg['global_data'].pop('radar_obj', None)

        # get start and stop times of new radar object
        dscfg['global_data']['starttime'], dscfg['global_data']['endtime'] = (
            time_avg_range(
                dscfg['timeinfo'], dscfg['global_data']['starttime'],
                dscfg['global_data']['endtime'], period))

        # check if volume time older than starttime
        if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
            dscfg['global_data'].update({'radar_obj': radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'radar_obj' not in dscfg['global_data']:
            return None, None

        dscfg['global_data']['radar_obj'].fields[field_name]['data'] /= (
            dscfg['global_data']['radar_obj'].fields[refl_name]['data'])

        new_dataset = deepcopy(dscfg['global_data']['radar_obj'])

        return new_dataset, ind_rad


def process_time_avg_flag(procstatus, dscfg, radar_list=None):
    """
    computes a flag field describing the conditions of the data used while
    averaging

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        period : float. Dataset keyword
            the period to average [s]. Default 3600.
        start_average : float. Dataset keyword
            when to start the average [s from midnight UTC]. Default 0.
        phidpmax: float. Dataset keyword
            maximum PhiDP
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """
    temp_name = None
    hydro_name = None
    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if ((datatype == 'PhiDP') or (datatype == 'PhiDPc')):
            phidp_name = get_fieldname_pyart(datatype)
        elif datatype == 'echoID':
            echo_name = get_fieldname_pyart(datatype)
        elif datatype == 'hydro':
            hydro_name = get_fieldname_pyart(datatype)
        elif datatype == 'TEMP':
            temp_name = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        phidpmax = 60.
        if 'phidpmax' in dscfg:
            phidpmax = dscfg['phidpmax']

        period = 3600.
        if 'period' in dscfg:
            period = dscfg['period']

        time_avg_flag = pyart.config.get_metadata('time_avg_flag')
        time_avg_flag['data'] = np.ma.zeros(
            (radar.nrays, radar.ngates), dtype=int)

        if phidp_name not in radar.fields:
            warn('Missing PhiDP data')
            time_avg_flag['data'] += 1
        else:
            phidp_field = radar.fields[phidp_name]
            time_avg_flag['data'][phidp_field['data'] > phidpmax] += 1

        if echo_name not in radar.fields:
            warn('Missing echo ID data')
            time_avg_flag['data'] += 100
        else:
            echo_field = radar.fields[echo_name]
            time_avg_flag['data'][echo_field['data'] == 2] += 100

        if hydro_name is not None:
            if ((hydro_name not in radar.fields) or
                    (echo_name not in radar.fields)):
                warn('Missing hydrometeor classification data')
                time_avg_flag['data'] += 10000
            else:
                hydro_field = radar.fields[hydro_name]
                # check where is no rain
                is_not_rain = np.logical_and(
                    hydro_field['data'] != 3, hydro_field['data'] != 5)
                # where is no rain should be precip
                is_not_rain = np.logical_and(
                    is_not_rain, echo_field['data'] == 3)
                time_avg_flag['data'][is_not_rain] += 10000
        elif temp_name is not None:
            if temp_name not in radar.fields:
                warn('Missing temperature data')
                time_avg_flag['data'] += 10000
            else:
                if 'radar_beam_width_h' in radar.instrument_parameters:
                    beamwidth = (
                        radar.instrument_parameters[
                            'radar_beam_width_h']['data'][0])
                else:
                    warn('Unknown radar antenna beamwidth.')
                    beamwidth = None

                mask_fzl, end_gate_arr = get_mask_fzl(
                    radar, fzl=None, doc=None, min_temp=0., thickness=700.,
                    beamwidth=beamwidth, temp_field=temp_name)
                time_avg_flag['data'][mask] += 10000

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.add_field('time_avg_flag', time_avg_flag)

        # first volume: initialize start and end time of averaging
        if dscfg['initialized'] == 0:
            start_average = 0.  # seconds from midnight
            if 'start_average' in dscfg:
                start_average = dscfg['start_average']

            date_00 = dscfg['timeinfo'].replace(
                hour=0, minute=0, second=0, microsecond=0)

            avg_par = dict()
            avg_par.update(
                {'starttime': date_00+datetime.timedelta(
                    seconds=start_average)})
            avg_par.update(
                {'endtime': avg_par['starttime']+datetime.timedelta(
                    seconds=period)})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        # no radar object in global data: create it
        if 'radar_obj' not in dscfg['global_data']:
            # get start and stop times of new radar object
            (dscfg['global_data']['starttime'],
             dscfg['global_data']['endtime']) = (
                time_avg_range(
                    dscfg['timeinfo'], dscfg['global_data']['starttime'],
                    dscfg['global_data']['endtime'], period))

            # check if volume time older than starttime
            if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
                dscfg['global_data'].update({'radar_obj': radar_aux})

            return None, None

        # still accumulating: add field to global field
        if dscfg['timeinfo'] < dscfg['global_data']['endtime']:
            flag_interp = interpol_field(
                dscfg['global_data']['radar_obj'], radar_aux, 'time_avg_flag')
            dscfg['global_data']['radar_obj'].fields[
                'time_avg_flag']['data'] += (
                    flag_interp['data'].filled(fill_value=0)).astype(int)

            return None, None

        # we have reached the end of the accumulation: start a new object
        new_dataset = deepcopy(dscfg['global_data']['radar_obj'])

        dscfg['global_data']['starttime'] += datetime.timedelta(
            seconds=period)
        dscfg['global_data']['endtime'] += datetime.timedelta(seconds=period)

        # remove old radar object from global_data dictionary
        dscfg['global_data'].pop('radar_obj', None)

        # get start and stop times of new radar object
        dscfg['global_data']['starttime'], dscfg['global_data']['endtime'] = (
            time_avg_range(
                dscfg['timeinfo'], dscfg['global_data']['starttime'],
                dscfg['global_data']['endtime'], period))

        # check if volume time older than starttime
        if dscfg['timeinfo'] > dscfg['global_data']['starttime']:
            dscfg['global_data'].update({'radar_obj': radar_aux})

        return new_dataset, ind_rad

    # no more files to process if there is global data pack it up
    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None
        if 'radar_obj' not in dscfg['global_data']:
            return None, None

        new_dataset = deepcopy(dscfg['global_data']['radar_obj'])

        return new_dataset, ind_rad


def process_colocated_gates(procstatus, dscfg, radar_list=None):
    """
    Find colocated gates within two radars

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        h_tol : float. Dataset keyword
            Tolerance in altitude difference between radar gates [m].
            Default 100.
        latlon_tol : float. Dataset keyword
            Tolerance in latitude and longitude position between radar gates
            [deg]. Default 0.0005
        vol_d_tol : float. Dataset keyword
            Tolerance in pulse volume diameter [m]. Default 100.
        vismin : float. Dataset keyword
            Minimum visibility [percent]. Default None.
        hmin : float. Dataset keyword
            Minimum altitude [m MSL]. Default None.
        hmax : float. Dataset keyword
            Maximum altitude [m MSL]. Default None.
        rmin : float. Dataset keyword
            Minimum range [m]. Default None.
        rmax : float. Dataset keyword
            Maximum range [m]. Default None.
        elmin : float. Dataset keyword
            Minimum elevation angle [deg]. Default None.
        elmax : float. Dataset keyword
            Maximum elevation angle [deg]. Default None.
        azmin : float. Dataset keyword
            Minimum azimuth angle [deg]. Default None.
        azmax : float. Dataset keyword
            Maximum azimuth angle [deg]. Default None.

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : radar object
        radar object containing the flag field
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # check how many radars are there
    radarnr_dict = dict()
    ind_radar_list = set()
    for datatypedescr in dscfg['datatype']:
        radarnr = datatypedescr.split(':')[0]
        radarnr_dict.update({radarnr: []})
        ind_radar_list.add(int(radarnr[5:8])-1)

    ind_radar_list = list(ind_radar_list)

    if (len(radarnr_dict) != 2) or (len(radar_list) < 2):
        warn('Intercomparison requires data from two different radars')
        return None, None

    # create the list of data types for each radar
    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if radarnr in radarnr_dict:
            radarnr_dict[radarnr].append(get_fieldname_pyart(datatype))

    radar1 = radar_list[ind_radar_list[0]]
    radar2 = radar_list[ind_radar_list[1]]

    coloc_gates_field = 'colocated_gates'

    h_tol = 100.
    if 'h_tol' in dscfg:
        h_tol = dscfg['h_tol']

    latlon_tol = 0.0005
    if 'latlon_tol' in dscfg:
        latlon_tol = dscfg['latlon_tol']

    vol_d_tol = 100.
    if 'vol_d_tol' in dscfg:
        vol_d_tol = dscfg['vol_d_tol']

    vismin = None
    if 'vismin' in dscfg:
        vismin = dscfg['vismin']

    hmin = None
    if 'hmin' in dscfg:
        hmin = dscfg['hmin']

    hmax = None
    if 'hmax' in dscfg:
        hmax = dscfg['hmax']

    rmin = None
    if 'rmin' in dscfg:
        rmin = dscfg['rmin']

    rmax = None
    if 'rmax' in dscfg:
        rmax = dscfg['rmax']

    elmin = None
    if 'elmin' in dscfg:
        elmin = dscfg['elmin']

    elmax = None
    if 'elmax' in dscfg:
        elmax = dscfg['elmax']

    azmin = None
    if 'azmin' in dscfg:
        azmin = dscfg['azmin']

    azmax = None
    if 'azmax' in dscfg:
        azmax = dscfg['azmax']

    visib_field = None
    if 'visibility' in radarnr_dict['RADAR'+'{:03d}'.format(
            ind_radar_list[0]+1)]:
        visib_field = 'visibility'
    gate_coloc_rad1_dict = pyart.util.intersection(
        radar1, radar2,
        h_tol=h_tol, latlon_tol=latlon_tol, vol_d_tol=vol_d_tol,
        vismin=vismin, hmin=hmin, hmax=hmax, rmin=rmin, rmax=rmax,
        elmin=elmin, elmax=elmax, azmin=azmin, azmax=azmax,
        visib_field=visib_field, intersec_field=coloc_gates_field)

    visib_field = None
    if 'visibility' in radarnr_dict['RADAR'+'{:03d}'.format(
            ind_radar_list[1]+1)]:
        visib_field = 'visibility'
    gate_coloc_rad2_dict = pyart.util.intersection(
        radar2, radar1,
        h_tol=h_tol, latlon_tol=latlon_tol, vol_d_tol=vol_d_tol,
        vismin=vismin, hmin=hmin, hmax=hmax, rmin=rmin, rmax=rmax,
        elmin=elmin, elmax=elmax, azmin=azmin, azmax=azmax,
        visib_field=visib_field, intersec_field=coloc_gates_field)

    new_rad1 = deepcopy(radar1)
    new_rad1.fields = dict()
    new_rad1.add_field('colocated_gates', gate_coloc_rad1_dict)

    new_rad2 = deepcopy(radar2)
    new_rad2.fields = dict()
    new_rad2.add_field('colocated_gates', gate_coloc_rad2_dict)

    coloc_rad1_dict, new_rad1.fields['colocated_gates'] = (
        pyart.util.colocated_gates(
            new_rad1, new_rad2, h_tol=h_tol,
            latlon_tol=latlon_tol, coloc_gates_field=coloc_gates_field))

    coloc_rad2_dict, new_rad2.fields['colocated_gates'] = (
        pyart.util.colocated_gates(
            new_rad2, new_rad1, h_tol=h_tol,
            latlon_tol=latlon_tol, coloc_gates_field=coloc_gates_field))

    # prepare output
    rad1_dict = {
        'coloc_dict': coloc_rad1_dict,
        'radar': new_rad1}

    rad2_dict = {
        'coloc_dict': coloc_rad2_dict,
        'radar': new_rad2}

    new_dataset = {
        'RADAR'+'{:03d}'.format(ind_radar_list[0]+1): rad1_dict,
        'RADAR'+'{:03d}'.format(ind_radar_list[1]+1): rad2_dict}

    return new_dataset, ind_radar_list


def process_intercomp(procstatus, dscfg, radar_list=None):
    """
    intercomparison between two radars

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    sun_hits_dict : dict
        dictionary containing a radar object, a sun_hits dict and a
        sun_retrieval dictionary
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        savedir = dscfg['colocgatespath']+dscfg['coloc_radars_name']+'/'

        fname = make_filename(
            'info', 'COLOCATED_GATES', dscfg['coloc_radars_name'], ['csv'],
            timeinfo=None)

        (rad1_ele, rad1_azi, rad1_rng,
         rad2_ele, rad2_azi, rad2_rng) = read_colocated_gates(savedir+fname[0])

        if rad1_ele is None:
            raise ValueError('Unable to intercompare radars. ' +
                             'Missing colocated gates file')

        dscfg['global_data'] = {
            'rad1_ele': rad1_ele,
            'rad1_azi': rad1_azi,
            'rad1_rng': rad1_rng,
            'rad2_ele': rad2_ele,
            'rad2_azi': rad2_azi,
            'rad2_rng': rad2_rng}

        return None, None

    if procstatus == 1:
        # check how many radars are there
        radarnr_dict = dict()
        ind_radar_list = set()
        for datatypedescr in dscfg['datatype']:
            radarnr = datatypedescr.split(':')[0]
            radarnr_dict.update({radarnr: []})
            ind_radar_list.add(int(radarnr[5:8])-1)

        ind_radar_list = list(ind_radar_list)

        if (len(radarnr_dict) != 2) or (len(radar_list) < 2):
            warn('Intercomparison requires data from two different radars')
            return None, None

        # create the list of data types for each radar
        for datatypedescr in dscfg['datatype']:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            field_name = get_fieldname_pyart(datatype)
            break

        radar1 = radar_list[ind_radar_list[0]]
        radar2 = radar_list[ind_radar_list[1]]

        if ((field_name not in radar1.fields) or
                (field_name not in radar2.fields)):
            warn('Unable to get values of field '+field_name +
                 ' at colocated range bins. ' +
                 'Field missing in one of the radars')
            return None, None

        if not dscfg['initialized']:
            dscfg['global_data'].update({'timeinfo': dscfg['timeinfo']})
            dscfg['initialized'] = 1

        rad1_field = radar1.fields[field_name]['data']
        rad2_field = radar2.fields[field_name]['data']

        azi_tol = 0.5
        ele_tol = 0.5
        rng_tol = 50.

        if 'azi_tol' in dscfg:
            azi_tol = dscfg['azi_tol']
        if 'ele_tol' in dscfg:
            ele_tol = dscfg['ele_tol']
        if 'rng_tol' in dscfg:
            rng_tol = dscfg['rng_tol']

        intercomp_dict = {
            'rad1_ele': [],
            'rad1_azi': [],
            'rad1_rng': [],
            'rad1_val': [],
            'rad2_ele': [],
            'rad2_azi': [],
            'rad2_rng': [],
            'rad2_val': []}

        # determine if radar data has to be averaged
        avg_rad1, avg_rad2, avg_rad_lim = get_range_bins_to_avg(
            radar1.range['data'], radar2.range['data'])

        for i in range(len(dscfg['global_data']['rad1_ele'])):
            ind_ray_rad1 = find_ray_index(
                radar1.elevation['data'], radar1.azimuth['data'],
                dscfg['global_data']['rad1_ele'][i],
                dscfg['global_data']['rad1_azi'][i],
                ele_tol=ele_tol, azi_tol=azi_tol)
            if ind_ray_rad1 is None:
                continue
            ind_rng_rad1 = find_rng_index(
                radar1.range['data'], dscfg['global_data']['rad1_rng'][i],
                rng_tol=rng_tol)
            if ind_rng_rad1 is None:
                continue

            ind_ray_rad2 = find_ray_index(
                radar2.elevation['data'], radar2.azimuth['data'],
                dscfg['global_data']['rad2_ele'][i],
                dscfg['global_data']['rad2_azi'][i],
                ele_tol=ele_tol, azi_tol=azi_tol)
            if ind_ray_rad2 is None:
                continue
            ind_rng_rad2 = find_rng_index(
                radar2.range['data'], dscfg['global_data']['rad2_rng'][i],
                rng_tol=rng_tol)
            if ind_rng_rad2 is None:
                continue

            val1 = np.ma.asarray(rad1_field[ind_ray_rad1[0], ind_rng_rad1[0]])
            val2 = np.ma.asarray(rad2_field[ind_ray_rad2[0], ind_rng_rad2[0]])
            if avg_rad1:
                if (ind_rng_rad1[0]+avg_rad_lim[1] >= radar1.ngates or
                        ind_rng_rad1[0]+avg_rad_lim[0] < 0):
                    continue
                ind_rng = list(range(
                    ind_rng_rad1[0]+avg_rad_lim[0],
                    ind_rng_rad1[0]+avg_rad_lim[1]+1))
                val1 = np.ma.asarray(np.ma.mean(
                    rad1_field[ind_ray_rad1[0], ind_rng]))
            elif avg_rad2:
                if (ind_rng_rad2[0]+avg_rad_lim[1] >= radar2.ngates or
                        ind_rng_rad2[0]+avg_rad_lim[0] < 0):
                    continue
                ind_rng = list(range(
                    ind_rng_rad2[0]+avg_rad_lim[0],
                    ind_rng_rad2[0]+avg_rad_lim[1]+1))
                val2 = np.ma.asarray(np.ma.mean(
                    rad2_field[ind_ray_rad2[0], ind_rng]))

            if val1.mask or val2.mask:
                continue

            intercomp_dict['rad1_ele'].append(
                radar1.elevation['data'][ind_ray_rad1[0]])
            intercomp_dict['rad1_azi'].append(
                radar1.azimuth['data'][ind_ray_rad1[0]])
            intercomp_dict['rad1_rng'].append(
                radar1.range['data'][ind_rng_rad1[0]])
            intercomp_dict['rad1_val'].append(val1)

            intercomp_dict['rad2_ele'].append(
                radar2.elevation['data'][ind_ray_rad2[0]])
            intercomp_dict['rad2_azi'].append(
                radar2.azimuth['data'][ind_ray_rad2[0]])
            intercomp_dict['rad2_rng'].append(
                radar2.range['data'][ind_rng_rad2[0]])
            intercomp_dict['rad2_val'].append(val2)

        new_dataset = {'intercomp_dict': intercomp_dict,
                       'final': False}
        return new_dataset, None

    if procstatus == 2:
        savedir = get_save_dir(
            dscfg['basepath'], dscfg['procname'], dscfg['dsname'],
            dscfg['coloc_data_dir'],
            timeinfo=dscfg['global_data']['timeinfo'], create_dir=False)

        fname = make_filename(
            'colocated_data', dscfg['type'], 'dBZc', ['csv'],
            timeinfo=dscfg['global_data']['timeinfo'], timeformat='%Y%m%d')

        fname = savedir+fname[0]

        coloc_data = read_colocated_data(fname)

        intercomp_dict = {
            'rad1_ele': coloc_data[0],
            'rad1_azi': coloc_data[1],
            'rad1_rng': coloc_data[2],
            'rad1_val': coloc_data[3],
            'rad2_ele': coloc_data[4],
            'rad2_azi': coloc_data[5],
            'rad2_rng': coloc_data[6],
            'rad2_val': coloc_data[7]}

        new_dataset = {'intercomp_dict': intercomp_dict,
                       'timeinfo': dscfg['global_data']['timeinfo'],
                       'final': True}

        return new_dataset, None


def process_sun_hits(procstatus, dscfg, radar_list=None):
    """
    monitoring of the radar using sun hits

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rmin : float. Dataset keyword
            minimum range where to look for a sun hit signal [m]. Default 20
        delev_max : float. Dataset keyword
            maximum elevation distance from nominal radar elevation where to
            look for a sun hit signal [deg]. Default 1.5
        dazim_max : float. Dataset keyword
            maximum azimuth distance from nominal radar elevation where to
            look for a sun hit signal [deg]. Default 1.5
        elmin : float. Dataset keyword
            minimum radar elevation where to look for sun hits [deg].
            Default 1.
        percent_bins : float. Dataset keyword.
            minimum percentage of range bins that have to contain signal to
            consider the ray a potential sun hit. Default 10.
        attg : float. Dataset keyword
            gaseous attenuation. Default None
        max_std : float. Dataset keyword
            maximum standard deviation to consider the data noise. Default 1.
        az_width_co : float. Dataset keyword
            co-polar antenna azimuth width (convoluted with sun width) [deg].
            Default None
        el_width_co : float. Dataset keyword
            co-polar antenna elevation width (convoluted with sun width)
            [deg]. Default None
        az_width_cross : float. Dataset keyword
            cross-polar antenna azimuth width (convoluted with sun width)
            [deg]. Default None
        el_width_cross : float. Dataset keyword
            cross-polar antenna elevation width (convoluted with sun width)
            [deg]. Default None
        ndays : int. Dataset keyword
            number of days used in sun retrieval. Default 1
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    sun_hits_dict : dict
        dictionary containing a radar object, a sun_hits dict and a
        sun_retrieval dictionary
    ind_rad : int
        radar index

    """

    if procstatus == 0:
        return None, None

    if procstatus == 1:
        for datatypedescr in dscfg['datatype']:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            if datatype == 'dBm':
                pwrh_field = 'signal_power_hh'
            if datatype == 'dBmv':
                pwrv_field = 'signal_power_vv'
            if datatype == 'ZDRu':
                zdr_field = 'unfiltered_differential_reflectivity'
            if datatype == 'ZDRuc':
                zdr_field = 'corrected_unfiltered_differential_reflectivity'

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if ((pwrh_field not in radar.fields) or
                (pwrv_field not in radar.fields) or
                (zdr_field not in radar.fields)):
            warn('Unable to get sun hits. Missing data')
            return None, None

        # initialize dataset
        if dscfg['initialized'] == 0:
            radar_par = dict()
            if 'frequency' in radar.instrument_parameters:
                radar_par.update(
                    {'wavelen': (
                        3e8 /
                        radar.instrument_parameters['frequency']['data'][0])})
            else:
                warn('Radar frequency unknown.')
            if 'radar_beam_width_h' in radar.instrument_parameters:
                radar_par.update(
                    {'beamwidth': (
                        radar.instrument_parameters[
                            'radar_beam_width_h']['data'][0])})
            elif 'radar_beam_width_v' in radar.instrument_parameters:
                radar_par.update(
                    {'beamwidth': (
                        radar.instrument_parameters[
                            'radar_beam_width_v']['data'][0])})
            else:
                warn('Antenna beam width unknown.')
            if 'pulse_width' in radar.instrument_parameters:
                radar_par.update(
                    {'pulse_width': (
                        radar.instrument_parameters[
                            'pulse_width']['data'][0])})
            else:
                warn('Pulse width unknown.')
            if radar.ray_angle_res is not None:
                radar_par.update(
                    {'angle_step': (radar.ray_angle_res['data'][0])})

            dscfg['global_data'] = radar_par
            dscfg['initialized'] = 1

        # default values
        rmin = 20.
        delev_max = 1.5
        dazim_max = 1.5
        elmin = 1.
        percent_bins = 10.
        attg = None
        max_std = 1.

        # user values
        if 'rmin' in dscfg:
            rmin = dscfg['rmin']
        if 'delev_max' in dscfg:
            delev_max = dscfg['delev_max']
        if 'dazim_max' in dscfg:
            dazim_max = dscfg['dazim_max']
        if 'elmin' in dscfg:
            elmin = dscfg['elmin']
        if 'percent_bins' in dscfg:
            percent_bins = dscfg['percent_bins']
        if 'attg' in dscfg:
            attg = dscfg['attg']
        if 'max_std' in dscfg:
            max_std = dscfg['max_std']

        ind_rmin = np.where(radar.range['data'] > rmin)[0][0]

        sun_hits, new_radar = pyart.correct.get_sun_hits(
            radar, delev_max=delev_max, dazim_max=dazim_max, elmin=elmin,
            ind_rmin=ind_rmin, percent_bins=percent_bins, max_std=max_std,
            attg=attg, pwrh_field=pwrh_field, pwrv_field=pwrv_field,
            zdr_field=zdr_field)

        if sun_hits is None:
            return None, None

        sun_hits_dataset = dict()
        sun_hits_dataset.update({'sun_hits': sun_hits})
        sun_hits_dataset.update({'radar': new_radar})

        return sun_hits_dataset, ind_rad

    if procstatus == 2:
        for datatypedescr in dscfg['datatype']:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            break

        ind_rad = int(radarnr[5:8])-1

        # default values
        az_width_co = None
        el_width_co = None
        az_width_cross = None
        el_width_cross = None
        nfiles = 1

        # user values
        if 'az_width_co' in dscfg:
            az_width_co = dscfg['az_width_co']
        if 'el_width_co' in dscfg:
            el_width_co = dscfg['el_width_co']
        if 'az_width_cross' in dscfg:
            az_width_cross = dscfg['az_width_cross']
        if 'el_width_cross' in dscfg:
            el_width_cross = dscfg['el_width_cross']
        if 'ndays' in dscfg:
            nfiles = dscfg['ndays']

        sun_hits = read_sun_hits_multiple_days(dscfg, nfiles=nfiles)

        if sun_hits[0] is None:
            return None, None

        sun_pwr_h = sun_hits[7]
        sun_pwr_v = sun_hits[11]

        sun_pwr_ref = np.ma.asarray(np.ma.masked)
        ref_time = None
        if (('pulse_width' in dscfg['global_data']) and
                ('wavelen' in dscfg['global_data']) and
                ('angle_step' in dscfg['global_data']) and
                ('beamwidth' in dscfg['global_data']) and
                (dscfg['AntennaGain'] is not None)):

            flx_dt, flx_val = read_solar_flux(
                dscfg['solarfluxpath']+'fluxtable.txt')

            if flx_dt is not None:
                flx_dt_closest, flx_val_closest = get_closest_solar_flux(
                    sun_hits[0], flx_dt, flx_val)

                sun_pwr_drao = pyart.correct.sun_power(
                    flx_val_closest, dscfg['global_data']['pulse_width'],
                    dscfg['global_data']['wavelen'], dscfg['AntennaGain'],
                    dscfg['global_data']['angle_step'],
                    dscfg['global_data']['beamwidth'])

                sun_pwr_ref = np.ma.asarray(sun_pwr_drao[-1])
                ref_time = flx_dt_closest[-1]

                # scaling of the power to account for solar flux variations.
                # The last sun hit is the reference
                scale_factor = sun_pwr_drao/sun_pwr_ref
                sun_pwr_h *= scale_factor
                sun_pwr_v *= scale_factor

        sun_retrieval_h = pyart.correct.sun_retrieval(
            sun_hits[4], sun_hits[6], sun_hits[3], sun_hits[5],
            sun_pwr_h, sun_hits[8],
            az_width_co=az_width_co, el_width_co=el_width_co,
            az_width_cross=az_width_cross, el_width_cross=el_width_cross,
            is_zdr=False)

        sun_retrieval_v = pyart.correct.sun_retrieval(
            sun_hits[4], sun_hits[6], sun_hits[3], sun_hits[5],
            sun_pwr_v, sun_hits[12],
            az_width_co=az_width_co, el_width_co=el_width_co,
            az_width_cross=az_width_cross, el_width_cross=el_width_cross,
            is_zdr=False)

        sun_retrieval_zdr = pyart.correct.sun_retrieval(
            sun_hits[4], sun_hits[6], sun_hits[3], sun_hits[5],
            sun_hits[15], sun_hits[16],
            az_width_co=az_width_co, el_width_co=el_width_co,
            az_width_cross=az_width_cross, el_width_cross=el_width_cross,
            is_zdr=True)

        sun_retrieval_dict = {
            'first_hit_time': sun_hits[0][0],
            'last_hit_time': sun_hits[0][-1],
            'dBm_sun_est': np.ma.asarray(np.ma.masked),
            'std(dBm_sun_est)': np.ma.asarray(np.ma.masked),
            'az_bias_h': np.ma.asarray(np.ma.masked),
            'el_bias_h': np.ma.asarray(np.ma.masked),
            'az_width_h': np.ma.asarray(np.ma.masked),
            'el_width_h': np.ma.asarray(np.ma.masked),
            'nhits_h': 0,
            'par_h': None,
            'dBmv_sun_est': np.ma.asarray(np.ma.masked),
            'std(dBmv_sun_est)': np.ma.asarray(np.ma.masked),
            'az_bias_v': np.ma.asarray(np.ma.masked),
            'el_bias_v': np.ma.asarray(np.ma.masked),
            'az_width_v': np.ma.asarray(np.ma.masked),
            'el_width_v': np.ma.asarray(np.ma.masked),
            'nhits_v': 0,
            'par_v': None,
            'ZDR_sun_est': np.ma.asarray(np.ma.masked),
            'std(ZDR_sun_est)': np.ma.asarray(np.ma.masked),
            'az_bias_zdr': np.ma.asarray(np.ma.masked),
            'el_bias_zdr': np.ma.asarray(np.ma.masked),
            'nhits_zdr': 0,
            'par_zdr': None,
            'dBm_sun_ref': sun_pwr_ref,
            'ref_time': ref_time}

        if sun_retrieval_h is not None:
            sun_retrieval_dict['dBm_sun_est'] = np.ma.asarray(
                sun_retrieval_h[0])
            sun_retrieval_dict['std(dBm_sun_est)'] = np.ma.asarray(
                sun_retrieval_h[1])
            sun_retrieval_dict['az_bias_h'] = np.ma.asarray(
                sun_retrieval_h[2])
            sun_retrieval_dict['el_bias_h'] = np.ma.asarray(
                sun_retrieval_h[3])
            sun_retrieval_dict['az_width_h'] = np.ma.asarray(
                sun_retrieval_h[4])
            sun_retrieval_dict['el_width_h'] = np.ma.asarray(
                sun_retrieval_h[5])
            sun_retrieval_dict['nhits_h'] = np.asarray(sun_retrieval_h[6])
            sun_retrieval_dict['par_h'] = np.ma.asarray(sun_retrieval_h[7])
        if sun_retrieval_v is not None:
            sun_retrieval_dict['dBmv_sun_est'] = np.ma.asarray(
                sun_retrieval_v[0])
            sun_retrieval_dict['std(dBmv_sun_est)'] = np.ma.asarray(
                sun_retrieval_v[1])
            sun_retrieval_dict['az_bias_v'] = np.ma.asarray(
                sun_retrieval_v[2])
            sun_retrieval_dict['el_bias_v'] = np.ma.asarray(
                sun_retrieval_v[3])
            sun_retrieval_dict['az_width_v'] = np.ma.asarray(
                sun_retrieval_v[4])
            sun_retrieval_dict['el_width_v'] = np.ma.asarray(
                sun_retrieval_v[5])
            sun_retrieval_dict['nhits_v'] = np.ma.asarray(sun_retrieval_v[6])
            sun_retrieval_dict['par_v'] = np.ma.asarray(sun_retrieval_v[7])
        if sun_retrieval_zdr is not None:
            sun_retrieval_dict['ZDR_sun_est'] = np.ma.asarray(
                sun_retrieval_zdr[0])
            sun_retrieval_dict['std(ZDR_sun_est)'] = np.ma.asarray(
                sun_retrieval_zdr[1])
            sun_retrieval_dict['az_bias_zdr'] = np.ma.asarray(
                sun_retrieval_zdr[2])
            sun_retrieval_dict['el_bias_dzr'] = np.ma.asarray(
                sun_retrieval_zdr[3])
            sun_retrieval_dict['nhits_zdr'] = np.asarray(sun_retrieval_zdr[6])
            sun_retrieval_dict['par_zdr'] = np.asarray(sun_retrieval_zdr[7])

        sun_hits_dict = {
            'time': sun_hits[0],
            'ray': sun_hits[1],
            'NPrng': sun_hits[2],
            'rad_el': sun_hits[3],
            'rad_az': sun_hits[4],
            'sun_el': sun_hits[5],
            'sun_az': sun_hits[6],
            'dBm_sun_hit': sun_hits[7],
            'std(dBm_sun_hit)': sun_hits[8],
            'NPh': sun_hits[9],
            'NPhval': sun_hits[10],
            'dBmv_sun_hit': sun_hits[11],
            'std(dBmv_sun_hit)': sun_hits[12],
            'NPv': sun_hits[13],
            'NPvval': sun_hits[14],
            'ZDR_sun_hit': sun_hits[15],
            'std(ZDR_sun_hit)': sun_hits[16],
            'NPzdr': sun_hits[17],
            'NPzdrval': sun_hits[18]}

        sun_hits_dataset = dict()
        sun_hits_dataset.update({'sun_hits_final': sun_hits_dict})
        sun_hits_dataset.update({'sun_retrieval': sun_retrieval_dict})

        return sun_hits_dataset, ind_rad
