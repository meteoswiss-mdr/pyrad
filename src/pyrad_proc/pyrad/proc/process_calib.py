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
    process_zdr_precip
    process_zdr_snow
    process_monitoring
    process_time_avg
    process_weighted_time_avg
    process_time_avg_flag
    process_colocated_gates
    process_intercomp
    process_intercomp_time_avg
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
from ..io.read_data_other import read_colocated_data_time_avg
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
    iso0 = None
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
        if datatype == 'H_ISO0':
            iso0 = 'height_over_iso0'
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

    # determine which freezing level reference
    temp_ref = 'temperature'
    if temp is None and iso0 is None:
        warn('Field to obtain the freezing level was not specified. ' +
             'Using fixed freezing level height')
        temp_ref = 'fixed_fzl'
    elif temp is not None:
        if temp not in radar.fields:
            warn('COSMO temperature field not available. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
    elif iso0 is not None:
        if iso0 not in radar.fields:
            warn('Height over iso0 field not available. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
        else:
            temp_ref = 'height_over_iso0'

    # determine freezing level height if necessary
    fzl = None
    if temp_ref == 'fixed_fzl':
        if 'fzl' in dscfg:
            fzl = dscfg['fzl']
        else:
            fzl = 2000.
            warn('Freezing level height not defined. Using default ' +
                 str(fzl)+' m')

    if dscfg['initialized'] == 0:
        # get frequency band
        freq_band = pyart.retrieve.get_freq_band(
            radar.instrument_parameters['frequency']['data'][0])

        # find unique elevations
        el_vec = np.unique(
            (10.*np.round(radar.elevation['data'], decimals=1)).astype(int))
        zdr_kdpzh_list = list()
        el_list = list()
        for i in range(len(el_vec)):
            fname = (
                dscfg['configpath'] + 'selfconsistency/' +
                'selfconsistency_zdr_zhkdp_'+freq_band+'band_temp10_elev' +
                '{:03d}'.format(el_vec[i])+'_mu05.txt')
            zdr_kdpzh_table = read_selfconsistency(fname)
            if zdr_kdpzh_table is not None:
                zdr_kdpzh_list.append(zdr_kdpzh_table)
                el_list.append((el_vec[i]/10.).astype(int))
        if not el_list:
            warn('Unable to retrieve PhiDP and KDP using self-consistency. ' +
                 'No selfconsistency files for the radar elevations.')

            return None

        zdr_kdpzh_dict = {'zdr_kdpzh': zdr_kdpzh_list,
                          'elev': el_list}
        dscfg['global_data'] = zdr_kdpzh_dict
        dscfg['initialized'] = 1

    if dscfg['initialized'] == 1:
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
            radar, dscfg['global_data'], min_rhohv=min_rhohv,
            max_phidp=max_phidp, smooth_wind_len=smooth_wind_len, doc=15,
            fzl=fzl, thickness=ml_thickness, refl_field=refl,
            phidp_field=phidp, zdr_field=zdr, temp_field=temp,
            iso0_field=iso0, rhohv_field=rhohv, kdpsim_field=kdpsim_field,
            phidpsim_field=phidpsim_field, temp_ref=temp_ref)

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
        fzl : float. Dataset keyword
            Default freezing level height. Default 2000.
        rsmooth : float. Dataset keyword
            length of the smoothing window [m]. Default 1000.
        min_rhohv : float. Dataset keyword
            minimum valid RhoHV. Default 0.92
        max_phidp : float. Dataset keyword
            maximum valid PhiDP [deg]. Default 20.
        ml_thickness : float. Dataset keyword
            Melting layer thickness [m]. Default 700.
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
    iso0 = None
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
        if datatype == 'H_ISO0':
            iso0 = 'height_over_iso0'
        if datatype == 'RhoHV':
            rhohv = 'cross_correlation_ratio'
        if datatype == 'RhoHVc':
            rhohv = 'corrected_cross_correlation_ratio'
        if datatype == 'uRhoHV':
            rhohv = 'uncorrected_cross_correlation_ratio'

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

    # determine which freezing level reference
    temp_ref = 'temperature'
    if temp is None and iso0 is None:
        warn('Field to obtain the freezing level was not specified. ' +
             'Using fixed freezing level height')
        temp_ref = 'fixed_fzl'
    elif temp is not None:
        if temp not in radar.fields:
            warn('COSMO temperature field not available. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
    elif iso0 is not None:
        if iso0 not in radar.fields:
            warn('Height over iso0 field not available. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
        else:
            temp_ref = 'height_over_iso0'

    # determine freezing level height if necessary
    fzl = None
    if temp_ref == 'fixed_fzl':
        if 'fzl' in dscfg:
            fzl = dscfg['fzl']
        else:
            fzl = 2000.
            warn('Freezing level height not defined. Using default ' +
                 str(fzl)+' m')

    if dscfg['initialized'] == 0:
        # get frequency band
        freq_band = pyart.retrieve.get_freq_band(
            radar.instrument_parameters['frequency']['data'][0])

        # find unique elevations
        el_vec = np.unique(
            (10.*np.round(radar.elevation['data'], decimals=1)).astype(int))
        zdr_kdpzh_list = list()
        el_list = list()
        for i in range(len(el_vec)):
            fname = (
                dscfg['configpath'] + 'selfconsistency/' +
                'selfconsistency_zdr_zhkdp_'+freq_band+'band_temp10_elev' +
                '{:03d}'.format(el_vec[i])+'_mu05.txt')
            zdr_kdpzh_table = read_selfconsistency(fname)
            if zdr_kdpzh_table is not None:
                zdr_kdpzh_list.append(zdr_kdpzh_table)
                el_list.append((el_vec[i]/10.).astype(int))
        if not el_list:
            warn('Unable to retrieve PhiDP and KDP using self-consistency. ' +
                 'No selfconsistency files for the radar elevations.')

            return None

        zdr_kdpzh_dict = {'zdr_kdpzh': zdr_kdpzh_list,
                          'elev': el_list}
        dscfg['global_data'] = zdr_kdpzh_dict
        dscfg['initialized'] = 1

    if dscfg['initialized'] == 1:
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

        r_res = radar.range['data'][1]-radar.range['data'][0]
        smooth_wind_len = int(rsmooth/r_res)
        min_rcons = int(rcell/r_res)

        refl_bias = pyart.correct.selfconsistency_bias(
            radar, dscfg['global_data'], min_rhohv=min_rhohv,
            max_phidp=max_phidp, smooth_wind_len=smooth_wind_len, doc=15,
            fzl=fzl, thickness=ml_thickness, min_rcons=min_rcons,
            dphidp_min=dphidp_min, dphidp_max=dphidp_max, refl_field=refl,
            phidp_field=phidp, zdr_field=zdr, temp_field=temp,
            iso0_field=iso0, rhohv_field=rhohv, temp_ref=temp_ref)

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
    iso0_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if datatype == 'RhoHV':
            rhohv_field = 'cross_correlation_ratio'
        if datatype == 'RhoHVc':
            rhohv_field = 'corrected_cross_correlation_ratio'
        if datatype == 'uRhoHV':
            rhohv_field = 'uncorrected_cross_correlation_ratio'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'TEMP':
            temp_field = 'temperature'
        if datatype == 'H_ISO0':
                iso0_field = 'height_over_iso0'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl_field not in radar.fields) or
            (rhohv_field not in radar.fields)):
        warn('Unable to estimate RhoHV in rain. Missing data')
        return None, None

    # determine which freezing level reference
    temp_ref = 'temperature'
    if temp_field is None and iso0_field is None:
        warn('Field to obtain the freezing level was not specified. ' +
             'Using fixed freezing level height')
        temp_ref = 'fixed_fzl'
    elif temp_field is not None:
        if temp_field not in radar.fields:
            warn('COSMO temperature field not available. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
    elif iso0_field is not None:
        if iso0_field not in radar.fields:
            warn('Height over iso0 field not available. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
        else:
            temp_ref = 'height_over_iso0'

    # determine freezing level height if necessary
    fzl = None
    if temp_ref == 'fixed_fzl':
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
        rhohv_field=rhohv_field, temp_field=temp_field, iso0_field=iso0_field,
        refl_field=refl_field, temp_ref=temp_ref)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    new_dataset.add_field('cross_correlation_ratio_in_rain', rhohv_rain)

    return new_dataset, ind_rad


def process_zdr_precip(procstatus, dscfg, radar_list=None):
    """
    Keeps only suitable data to evaluate the differential reflectivity in
    moderate rain or precipitation (for vertical scans)

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        ml_filter : boolean. Dataset keyword
            indicates if a filter on data in and above the melting layer is
            applied. Default True.
        rmin : float. Dataset keyword
            minimum range where to look for rain [m]. Default 1000.
        rmax : float. Dataset keyword
            maximum range where to look for rain [m]. Default 50000.
        Zmin : float. Dataset keyword
            minimum reflectivity to consider the bin as precipitation [dBZ].
            Default 20.
        Zmax : float. Dataset keyword
            maximum reflectivity to consider the bin as precipitation [dBZ]
            Default 22.
        RhoHVmin : float. Dataset keyword
            minimum RhoHV to consider the bin as precipitation
            Default 0.97
        PhiDPmax : float. Dataset keyword
            maximum PhiDP to consider the bin as precipitation [deg]
            Default 10.
        elmax : float. Dataset keyword
            maximum elevation angle where to look for precipitation [deg]
            Default None.
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
    iso0_field = None
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
        if datatype == 'uRhoHV':
            rhohv_field = 'uncorrected_cross_correlation_ratio'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'TEMP':
            temp_field = 'temperature'
        if datatype == 'H_ISO0':
            iso0_field = 'height_over_iso0'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl_field not in radar.fields) or
            (rhohv_field not in radar.fields) or
            (zdr_field not in radar.fields) or
            (phidp_field not in radar.fields)):
        warn('Unable to estimate ZDR in rain. Missing data')
        return None, None

    # if data in and above the melting layer has to be filtered determine the
    # field to use
    fzl = None
    ml_filter = True
    if 'ml_filter' in dscfg:
        ml_filter = dscfg['ml_filter']

    if ml_filter:
        # determine which freezing level reference
        temp_ref = 'temperature'
        if temp_field is None and iso0_field is None:
            warn('Field to obtain the freezing level was not specified. ' +
                 'Using fixed freezing level height')
            temp_ref = 'fixed_fzl'
        elif temp_field is not None:
            if temp_field not in radar.fields:
                warn('COSMO temperature field not available. ' +
                     'Using fixed freezing level height')
                temp_ref = 'fixed_fzl'
        elif iso0_field is not None:
            if iso0_field not in radar.fields:
                warn('Height over iso0 field not available. ' +
                     'Using fixed freezing level height')
                temp_ref = 'fixed_fzl'
            else:
                temp_ref = 'height_over_iso0'

        # determine freezing level height if necessary
        if temp_ref == 'fixed_fzl':
            if 'fzl' in dscfg:
                fzl = dscfg['fzl']
            else:
                fzl = 2000.
                warn('Freezing level height not defined. Using default ' +
                     str(fzl)+' m')
    else:
        temp_ref = None

    # default values
    rmin = 1000.
    rmax = 50000.
    zmin = 20.
    zmax = 22.
    rhohvmin = 0.97
    phidpmax = 10.
    elmax = None
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

    zdr_precip = pyart.correct.est_zdr_precip(
        radar, ind_rmin=ind_rmin, ind_rmax=ind_rmax, zmin=zmin,
        zmax=zmax, rhohvmin=rhohvmin, phidpmax=phidpmax, elmax=elmax,
        thickness=thickness, doc=15, fzl=fzl, zdr_field=zdr_field,
        rhohv_field=rhohv_field, phidp_field=phidp_field,
        temp_field=temp_field, iso0_field=iso0_field, refl_field=refl_field,
        temp_ref=temp_ref)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    new_dataset.add_field(
        'differential_reflectivity_in_precipitation', zdr_precip)

    return new_dataset, ind_rad


def process_zdr_snow(procstatus, dscfg, radar_list=None):
    """
    Keeps only suitable data to evaluate the differential reflectivity in
    snow

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
            minimum reflectivity to consider the bin as snow [dBZ].
            Default 0.
        Zmax : float. Dataset keyword
            maximum reflectivity to consider the bin as snow [dBZ]
            Default 30.
        SNRmin : float. Dataset keyword
            minimum SNR to consider the bin as snow [dB].
            Default 10.
        SNRmax : float. Dataset keyword
            maximum SNR to consider the bin as snow [dB]
            Default 50.
        RhoHVmin : float. Dataset keyword
            minimum RhoHV to consider the bin as snow
            Default 0.97
        PhiDPmax : float. Dataset keyword
            maximum PhiDP to consider the bin as snow [deg]
            Default 10.
        elmax : float. Dataset keyword
            maximum elevation angle where to look for snow [deg]
            Default None.
        KDPmax : float. Dataset keyword
            maximum KDP to consider the bin as snow [deg]
            Default None
        TEMPmin : float. Dataset keyword
            minimum temperature to consider the bin as snow [deg C].
            Default None
        TEMPmax : float. Dataset keyword
            maximum temperature to consider the bin as snow [deg C]
            Default None
        hydroclass : list of ints. Dataset keyword
            list of hydrometeor classes to keep for the analysis
            Default [1] (dry snow)
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
    kdp_field = None
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
        if datatype == 'uRhoHV':
            rhohv_field = 'uncorrected_cross_correlation_ratio'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'TEMP':
            temp_field = 'temperature'
        if datatype == 'PhiDP':
            kdp_field = 'specific_differential_phase'
        if datatype == 'PhiDPc':
            kdp_field = 'corrected_specific_differential_phase'
        if datatype == 'SNRh':
            snr_field = 'signal_to_noise_ratio_hh'
        if datatype == 'SNRv':
            snr_field = 'signal_to_noise_ratio_vv'
        if datatype == 'hydro':
            hydro_field = 'radar_echo_classification'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl_field not in radar.fields) or
            (rhohv_field not in radar.fields) or
            (zdr_field not in radar.fields) or
            (phidp_field not in radar.fields) or
            (snr_field not in radar.fields) or
            (hydro_field not in radar.fields)):
        warn('Unable to estimate ZDR in snow. Missing data')
        return None, None

    # default values
    rmin = 1000.
    rmax = 50000.
    zmin = 0.
    zmax = 30.
    snrmin = 10.
    snrmax = 50.
    rhohvmin = 0.97
    phidpmax = 10.
    elmax = None
    kdpmax = None
    tempmin = None
    tempmax = None
    hydroclass = [1]

    # user defined values
    if 'rmin' in dscfg:
        rmin = dscfg['rmin']
    if 'rmax' in dscfg:
        rmax = dscfg['rmax']
    if 'Zmin' in dscfg:
        zmin = dscfg['Zmin']
    if 'Zmax' in dscfg:
        zmax = dscfg['Zmax']
    if 'SNRmin' in dscfg:
        snrmin = dscfg['SNRmin']
    if 'SNRmax' in dscfg:
        snrmax = dscfg['SNRmax']
    if 'RhoHVmin' in dscfg:
        rhohvmin = dscfg['RhoHVmin']
    if 'PhiDPmax' in dscfg:
        phidpmax = dscfg['PhiDPmax']
    if 'KDPmax' in dscfg:
        kdpmax = dscfg['KDPmax']
    if 'TEMPmin' in dscfg:
        tempmin = dscfg['TEMPmin']
    if 'TEMPmax' in dscfg:
        tempmax = dscfg['TEMPmax']
    if 'elmax' in dscfg:
        elmax = dscfg['elmax']
    if 'hydroclass' in dscfg:
        hydroclass = dscfg['hydroclass']

    ind_rmin = np.where(radar.range['data'] > rmin)[0][0]
    ind_rmax = np.where(radar.range['data'] < rmax)[0][-1]

    zdr_snow = pyart.correct.est_zdr_snow(
        radar, ind_rmin=ind_rmin, ind_rmax=ind_rmax, zmin=zmin, zmax=zmax,
        snrmin=snrmin, snrmax=snrmax, rhohvmin=rhohvmin,
        kept_values=hydroclass, phidpmax=phidpmax, kdpmax=kdpmax,
        tempmin=tempmin, tempmax=tempmax, elmax=elmax, zdr_field=zdr_field,
        rhohv_field=rhohv_field, phidp_field=phidp_field,
        temp_field=temp_field, snr_field=snr_field, hydro_field=hydro_field,
        kdp_field=kdp_field, refl_field=refl_field)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    new_dataset.add_field(
        'differential_reflectivity_in_snow', zdr_snow)

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

        field = deepcopy(radar.fields[field_name]['data'])

        # put gates with values off limits to limit
        mask = np.ma.getmaskarray(field)
        ind = np.where(np.logical_and(mask == False, field < bins[0]))
        field[ind] = bins[0]

        ind = np.where(np.logical_and(mask == False, field > bins[-1]))
        field[ind] = bins[-1]

        for ray in range(radar.nrays):
            field_dict['data'][ray, :], bin_edges = np.histogram(
                field[ray, :].compressed(), bins=bins)

        radar_aux.add_field(field_name, field_dict)

        # keep histogram in Memory or add to existing histogram
        if dscfg['initialized'] == 0:
            dscfg['global_data'] = {'hist_obj': radar_aux,
                                    'timeinfo': dscfg['timeinfo']}
            dscfg['initialized'] = 1
        else:
            field_interp = interpol_field(
                dscfg['global_data']['hist_obj'], radar_aux, field_name,
                fill_value=0)
            dscfg['global_data']['hist_obj'].fields[field_name]['data'] += (
                field_interp['data'].filled(fill_value=0)).astype('int64')

            dscfg['global_data']['timeinfo'] = dscfg['timeinfo']

        dataset = dict()
        dataset.update({'hist_obj': radar_aux})
        dataset.update({'hist_type': 'instant'})
        dataset.update({'timeinfo': dscfg['timeinfo']})

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
        dataset.update({'hist_obj': dscfg['global_data']['hist_obj']})
        dataset.update({'hist_type': 'cumulative'})
        dataset.update({'timeinfo': dscfg['global_data']['timeinfo']})

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
            avg_par.update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        if dscfg['initialized'] == 0:
            return None, None

        dscfg['global_data']['timeinfo'] = dscfg['timeinfo']
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

        new_dataset = {
            'radar_obj': deepcopy(dscfg['global_data']['radar_obj']),
            'timeinfo': dscfg['global_data']['endtime']}

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

        new_dataset = {
            'radar_obj': deepcopy(dscfg['global_data']['radar_obj']),
            'timeinfo': dscfg['global_data']['endtime']}

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
            avg_par.update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        if dscfg['initialized'] == 0:
            return None, None

        dscfg['global_data']['timeinfo'] = dscfg['timeinfo']
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

        new_dataset = {
            'radar_obj': deepcopy(dscfg['global_data']['radar_obj']),
            'timeinfo': dscfg['global_data']['endtime']}

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

        new_dataset = {
            'radar_obj': deepcopy(dscfg['global_data']['radar_obj']),
            'timeinfo': dscfg['global_data']['endtime']}

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
    iso0_name = None
    echo_name = None
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
        elif datatype == 'H_ISO0':
            iso0_name = 'height_over_iso0'

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

        if echo_name is not None:
            if echo_name not in radar.fields:
                warn('Missing echo ID data')
                time_avg_flag['data'] += 100
            else:
                echo_field = radar.fields[echo_name]
                time_avg_flag['data'][echo_field['data'] == 2] += 100

        if hydro_name is not None and echo_name is not None:
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

                mask_fzl, end_gate_arr = pyart.correct.get_mask_fzl(
                    radar, fzl=None, doc=None, min_temp=0., max_h_iso0=0.,
                    thickness=700., beamwidth=beamwidth,
                    temp_field=temp_name, iso0_field=iso0_name,
                    temp_ref='temperature')
                time_avg_flag['data'][mask_fzl] += 10000
        elif iso0_name is not None:
            if iso0_name not in radar.fields:
                warn('Missing height relative to iso0 data')
                time_avg_flag['data'] += 10000
            else:
                if 'radar_beam_width_h' in radar.instrument_parameters:
                    beamwidth = (
                        radar.instrument_parameters[
                            'radar_beam_width_h']['data'][0])
                else:
                    warn('Unknown radar antenna beamwidth.')
                    beamwidth = None

                mask_fzl, end_gate_arr = pyart.correct.get_mask_fzl(
                    radar, fzl=None, doc=None, min_temp=0., max_h_iso0=0.,
                    thickness=700., beamwidth=beamwidth,
                    temp_field=temp_name, iso0_field=iso0_name,
                    temp_ref='height_over_iso0')
                time_avg_flag['data'][mask_fzl] += 10000

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
            avg_par.update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'] = avg_par
            dscfg['initialized'] = 1

        if dscfg['initialized'] == 0:
            return None, None

        dscfg['global_data']['timeinfo'] = dscfg['timeinfo']
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
        new_dataset = {
            'radar_obj': deepcopy(dscfg['global_data']['radar_obj']),
            'timeinfo': dscfg['global_data']['endtime']}

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

        new_dataset = {
            'radar_obj': deepcopy(dscfg['global_data']['radar_obj']),
            'timeinfo': dscfg['global_data']['endtime']}

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
        azrad1min : float. Dataset keyword
            Minimum azimuth angle [deg] for radar 1. Default None.
        azrad1max : float. Dataset keyword
            Maximum azimuth angle [deg] for radar 1. Default None.
        azrad2min : float. Dataset keyword
            Minimum azimuth angle [deg] for radar 2. Default None.
        azrad2max : float. Dataset keyword
            Maximum azimuth angle [deg] for radar 2. Default None.

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

    if radar1 is None or radar2 is None:
        warn('Unable to inter-compare radars. Missing radar')

    if 'instrument_name' in radar1.metadata:
        print('Radar 1: '+radar1.metadata['instrument_name'])
    if 'instrument_name' in radar2.metadata:
        print('Radar 2: '+radar2.metadata['instrument_name'])

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

    azrad1min = None
    if 'azrad1min' in dscfg:
        azrad1min = dscfg['azrad1min']

    azrad1max = None
    if 'azrad1max' in dscfg:
        azrad1max = dscfg['azrad1max']

    azrad2min = None
    if 'azrad2min' in dscfg:
        azrad2min = dscfg['azrad2min']

    azrad2max = None
    if 'azrad2max' in dscfg:
        azrad2max = dscfg['azrad2max']

    visib_field = None
    if 'visibility' in radarnr_dict['RADAR'+'{:03d}'.format(
            ind_radar_list[0]+1)]:
        visib_field = 'visibility'
    if vismin is not None and visib_field is None:
        warn('Unable to filter data according to visibility. ' +
             'Visibility field for RADAR'+'{:03d}'.format(
                ind_radar_list[0]+1)+' not available')

    gate_coloc_rad1_dict = pyart.util.intersection(
        radar1, radar2,
        h_tol=h_tol, latlon_tol=latlon_tol, vol_d_tol=vol_d_tol,
        vismin=vismin, hmin=hmin, hmax=hmax, rmin=rmin, rmax=rmax,
        elmin=elmin, elmax=elmax, azmin=azrad1min, azmax=azrad1max,
        visib_field=visib_field, intersec_field=coloc_gates_field)

    visib_field = None
    if 'visibility' in radarnr_dict['RADAR'+'{:03d}'.format(
            ind_radar_list[1]+1)]:
        visib_field = 'visibility'

    if vismin is not None and visib_field is None:
        warn('Unable to filter data according to visibility. ' +
             'Visibility field for RADAR'+'{:03d}'.format(
                ind_radar_list[1]+1)+' not available')

    gate_coloc_rad2_dict = pyart.util.intersection(
        radar2, radar1,
        h_tol=h_tol, latlon_tol=latlon_tol, vol_d_tol=vol_d_tol,
        vismin=vismin, hmin=hmin, hmax=hmax, rmin=rmin, rmax=rmax,
        elmin=elmin, elmax=elmax, azmin=azrad2min, azmax=azrad2max,
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
        coloc_data_dir : string. Dataset keyword
            name of the directory containing the csv file with colocated data
        coloc_radars_name : string. Dataset keyword
            string identifying the radar names
        azi_tol : float. Dataset keyword
            azimuth tolerance between the two radars. Default 0.5 deg
        ele_tol : float. Dataset keyword
            elevation tolerance between the two radars. Default 0.5 deg
        rng_tol : float. Dataset keyword
            range tolerance between the two radars. Default 50 m
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing a dictionary with intercomparison data and the
        key "final" which contains a boolean that is true when all volumes
        have been processed
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        savedir = dscfg['colocgatespath']+dscfg['coloc_radars_name']+'/'

        prdtype = 'info'
        if 'prdtype' in dscfg:
            prdtype = dscfg['prdtype']

        fname = make_filename(
            prdtype, 'COLOCATED_GATES', dscfg['coloc_radars_name'], ['csv'],
            timeinfo=None)[0]

        (rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
         rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng) = (
            read_colocated_gates(savedir+fname))

        if rad1_ele is None:
            raise ValueError('Unable to intercompare radars. ' +
                             'Missing colocated gates file')

        dscfg['global_data'] = {
            'rad1_ray_ind': rad1_ray_ind,
            'rad1_rng_ind': rad1_rng_ind,
            'rad1_ele': rad1_ele,
            'rad1_azi': rad1_azi,
            'rad1_rng': rad1_rng,
            'rad2_ray_ind': rad2_ray_ind,
            'rad2_rng_ind': rad2_rng_ind,
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

        if radar1 is None or radar2 is None:
            warn('Unable to inter-compare radars. Missing radar')
            return None, None

        if ((field_name not in radar1.fields) or
                (field_name not in radar2.fields)):
            warn('Unable to get values of field '+field_name +
                 ' at colocated range bins. ' +
                 'Field missing in one of the radars')
            return None, None

        if not dscfg['initialized']:
            dscfg['global_data'].update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'].update(
                {'rad1_name': dscfg['RadarName'][ind_radar_list[0]]})
            dscfg['global_data'].update(
                {'rad2_name': dscfg['RadarName'][ind_radar_list[1]]})
            dscfg['initialized'] = 1

        rad1_field = radar1.fields[field_name]['data']
        rad2_field = radar2.fields[field_name]['data']

        intercomp_dict = {
            'rad1_ray_ind': [],
            'rad1_rng_ind': [],
            'rad1_ele': [],
            'rad1_azi': [],
            'rad1_rng': [],
            'rad1_val': [],
            'rad2_ray_ind': [],
            'rad2_rng_ind': [],
            'rad2_ele': [],
            'rad2_azi': [],
            'rad2_rng': [],
            'rad2_val': []}

        # determine if radar data has to be averaged
        avg_rad1, avg_rad2, avg_rad_lim = get_range_bins_to_avg(
            radar1.range['data'], radar2.range['data'])

        # rays are indexed to regular grid
        rays_are_indexed = False
        if 'rays_are_indexed' in dscfg:
            rays_are_indexed = dscfg['rays_are_indexed']

        if not rays_are_indexed:
            azi_tol = 0.5
            ele_tol = 0.5
            rng_tol = 50.

            if 'azi_tol' in dscfg:
                azi_tol = dscfg['azi_tol']
            if 'ele_tol' in dscfg:
                ele_tol = dscfg['ele_tol']
            if 'rng_tol' in dscfg:
                rng_tol = dscfg['rng_tol']

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

                val1 = np.ma.asarray(rad1_field[ind_ray_rad1, ind_rng_rad1])
                val2 = np.ma.asarray(rad2_field[ind_ray_rad2, ind_rng_rad2])
                if avg_rad1:
                    if (ind_rng_rad1+avg_rad_lim[1] >= radar1.ngates or
                            ind_rng_rad1+avg_rad_lim[0] < 0):
                        continue
                    ind_rng = list(range(
                        ind_rng_rad1+avg_rad_lim[0],
                        ind_rng_rad1+avg_rad_lim[1]+1))
                    val1 = np.ma.asarray(np.ma.mean(
                        rad1_field[ind_ray_rad1, ind_rng]))
                elif avg_rad2:
                    if (ind_rng_rad2+avg_rad_lim[1] >= radar2.ngates or
                            ind_rng_rad2+avg_rad_lim[0] < 0):
                        continue
                    ind_rng = list(range(
                        ind_rng_rad2+avg_rad_lim[0],
                        ind_rng_rad2+avg_rad_lim[1]+1))
                    val2 = np.ma.asarray(np.ma.mean(
                        rad2_field[ind_ray_rad2, ind_rng]))

                if val1.mask or val2.mask:
                    continue

                intercomp_dict['rad1_ray_ind'].append(ind_ray_rad1)
                intercomp_dict['rad1_rng_ind'].append(ind_rng_rad1)
                intercomp_dict['rad1_ele'].append(
                    radar1.elevation['data'][ind_ray_rad1])
                intercomp_dict['rad1_azi'].append(
                    radar1.azimuth['data'][ind_ray_rad1])
                intercomp_dict['rad1_rng'].append(
                    radar1.range['data'][ind_rng_rad1])
                intercomp_dict['rad1_val'].append(val1)

                intercomp_dict['rad2_ray_ind'].append(ind_ray_rad2)
                intercomp_dict['rad2_rng_ind'].append(ind_rng_rad2)
                intercomp_dict['rad2_ele'].append(
                    radar2.elevation['data'][ind_ray_rad2])
                intercomp_dict['rad2_azi'].append(
                    radar2.azimuth['data'][ind_ray_rad2])
                intercomp_dict['rad2_rng'].append(
                    radar2.range['data'][ind_rng_rad2])
                intercomp_dict['rad2_val'].append(val2)
        else:
            rad1_ray_ind = deepcopy(dscfg['global_data']['rad1_ray_ind'])
            rad1_rng_ind = deepcopy(dscfg['global_data']['rad1_rng_ind'])
            rad2_ray_ind = deepcopy(dscfg['global_data']['rad2_ray_ind'])
            rad2_rng_ind = deepcopy(dscfg['global_data']['rad2_rng_ind'])

            val1_vec = rad1_field[rad1_ray_ind, rad1_rng_ind]
            val2_vec = rad2_field[rad1_ray_ind, rad1_rng_ind]

            mask_val1 = np.ma.getmaskarray(val1_vec)
            mask_val2 = np.ma.getmaskarray(val2_vec)

            isvalid = np.logical_not(np.logical_and(mask_val1, mask_val2))

            val1_vec = val1_vec[isvalid]
            val2_vec = val2_vec[isvalid]
            rad1_ray_ind = rad1_ray_ind[isvalid]
            rad1_rng_ind = rad1_rng_ind[isvalid]
            rad2_ray_ind = rad2_ray_ind[isvalid]
            rad2_rng_ind = rad2_rng_ind[isvalid]

            intercomp_dict['rad1_ray_ind'] = rad1_ray_ind
            intercomp_dict['rad1_rng_ind'] = rad1_rng_ind
            intercomp_dict['rad1_ele'] = radar1.elevation['data'][rad1_ray_ind]
            intercomp_dict['rad1_azi'] = radar1.azimuth['data'][rad1_ray_ind]
            intercomp_dict['rad1_rng'] = radar1.range['data'][rad1_rng_ind]
            intercomp_dict['rad1_dBZavg'] = refl1_vec
            intercomp_dict['rad1_PhiDPavg'] = phidp1_vec
            intercomp_dict['rad1_Flagavg'] = flag1_vec

            intercomp_dict['rad2_ray_ind'] = rad2_ray_ind
            intercomp_dict['rad2_rng_ind'] = rad2_rng_ind
            intercomp_dict['rad2_ele'] = radar2.elevation['data'][rad2_ray_ind]
            intercomp_dict['rad2_azi'] = radar2.azimuth['data'][rad2_ray_ind]
            intercomp_dict['rad2_rng'] = radar2.range['data'][rad2_rng_ind]
            intercomp_dict['rad2_dBZavg'] = refl2_vec
            intercomp_dict['rad2_PhiDPavg'] = phidp2_vec
            intercomp_dict['rad2_Flagavg'] = flag2_vec

        new_dataset = {'intercomp_dict': intercomp_dict,
                       'timeinfo': dscfg['global_data']['timeinfo'],
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
            'rad1_name': dscfg['global_data']['rad1_name'],
            'rad1_ray_ind': coloc_data[0],
            'rad1_rng_ind': coloc_data[1],
            'rad1_ele': coloc_data[2],
            'rad1_azi': coloc_data[3],
            'rad1_rng': coloc_data[4],
            'rad1_val': coloc_data[5],
            'rad2_name': dscfg['global_data']['rad2_name'],
            'rad2_ray_ind': coloc_data[6],
            'rad2_rng_ind': coloc_data[7],
            'rad2_ele': coloc_data[8],
            'rad2_azi': coloc_data[9],
            'rad2_rng': coloc_data[10],
            'rad2_val': coloc_data[11]}

        new_dataset = {'intercomp_dict': intercomp_dict,
                       'timeinfo': dscfg['global_data']['timeinfo'],
                       'final': True}

        return new_dataset, None


def process_intercomp_time_avg(procstatus, dscfg, radar_list=None):
    """
    intercomparison between the average reflectivity of two radars

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        coloc_data_dir : string. Dataset keyword
            name of the directory containing the csv file with colocated data
        coloc_radars_name : string. Dataset keyword
            string identifying the radar names
        azi_tol : float. Dataset keyword
            azimuth tolerance between the two radars. Default 0.5 deg
        ele_tol : float. Dataset keyword
            elevation tolerance between the two radars. Default 0.5 deg
        rng_tol : float. Dataset keyword
            range tolerance between the two radars. Default 50 m
        clt_max : int. Dataset keyword
            maximum number of samples that can be clutter contaminated.
            Default 100 i.e. all
        phi_excess_max : int. Dataset keyword
            maximum number of samples that can have excess instantaneous
            PhiDP. Default 100 i.e. all
        non_rain_max : int. Dataset keyword
            maximum number of samples that can be no rain. Default 100 i.e. all
        phi_avg_max : float. Dataset keyword
            maximum average PhiDP allowed. Default 600 deg i.e. any

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing a dictionary with intercomparison data and the
        key "final" which contains a boolean that is true when all volumes
        have been processed
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        savedir = dscfg['colocgatespath']+dscfg['coloc_radars_name']+'/'

        prdtype = 'info'
        if 'prdtype' in dscfg:
            prdtype = dscfg['prdtype']

        fname = make_filename(
            prdtype, 'COLOCATED_GATES', dscfg['coloc_radars_name'], ['csv'],
            timeinfo=None)[0]

        (rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
         rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng) = (
            read_colocated_gates(savedir+fname))

        if rad1_ele is None:
            raise ValueError('Unable to intercompare radars. ' +
                             'Missing colocated gates file')

        dscfg['global_data'] = {
            'rad1_ray_ind': rad1_ray_ind,
            'rad1_rng_ind': rad1_rng_ind,
            'rad1_ele': rad1_ele,
            'rad1_azi': rad1_azi,
            'rad1_rng': rad1_rng,
            'rad2_ray_ind': rad2_ray_ind,
            'rad2_rng_ind': rad2_rng_ind,
            'rad2_ele': rad2_ele,
            'rad2_azi': rad2_azi,
            'rad2_rng': rad2_rng}

        return None, None

    if procstatus == 1:
        # check how many radars are there
        ind_radar_list = set()
        for datatypedescr in dscfg['datatype']:
            radarnr = datatypedescr.split(':')[0]
            ind_radar_list.add(int(radarnr[5:8])-1)

        ind_radar_list = list(ind_radar_list)

        if (len(ind_radar_list) != 2) or (len(radar_list) < 2):
            warn('Intercomparison requires data from two different radars')
            return None, None

        radarnr_list = ['RADAR'+'{:03d}'.format(ind_radar_list[0]+1),
                        'RADAR'+'{:03d}'.format(ind_radar_list[1]+1)]

        # get field names
        for datatypedescr in dscfg['datatype']:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            if radarnr == radarnr_list[0]:
                if ((datatype == 'dBZ') or (datatype == 'dBZc') or
                        (datatype == 'dBuZ') or (datatype == 'dBZv') or
                        (datatype == 'dBZvc') or (datatype == 'dBuZv')):
                    rad1_refl_field = get_fieldname_pyart(datatype)
                elif (datatype == 'PhiDP') or (datatype == 'PhiDPc'):
                    rad1_phidp_field = get_fieldname_pyart(datatype)
                elif datatype == 'time_avg_flag':
                    rad1_flag_field = get_fieldname_pyart(datatype)
            elif radarnr == radarnr_list[1]:
                if ((datatype == 'dBZ') or (datatype == 'dBZc') or
                        (datatype == 'dBuZ') or (datatype == 'dBZv') or
                        (datatype == 'dBZvc') or (datatype == 'dBuZv')):
                    rad2_refl_field = get_fieldname_pyart(datatype)
                elif (datatype == 'PhiDP') or (datatype == 'PhiDPc'):
                    rad2_phidp_field = get_fieldname_pyart(datatype)
                elif datatype == 'time_avg_flag':
                    rad2_flag_field = get_fieldname_pyart(datatype)

        radar1 = radar_list[ind_radar_list[0]]
        radar2 = radar_list[ind_radar_list[1]]

        if radar1 is None or radar2 is None:
            warn('Unable to inter-compare radars. Missing radar')
            return None, None

        if ((rad1_refl_field not in radar1.fields) or
                (rad1_phidp_field not in radar1.fields) or
                (rad1_flag_field not in radar1.fields) or
                (rad2_refl_field not in radar2.fields) or
                (rad2_phidp_field not in radar2.fields) or
                (rad2_flag_field not in radar2.fields)):
            warn('Unable to compare radar time avg fields. ' +
                 'Fields missing')
            return None, None

        if not dscfg['initialized']:
            dscfg['global_data'].update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'].update(
                {'rad1_name': dscfg['RadarName'][ind_radar_list[0]]})
            dscfg['global_data'].update(
                {'rad2_name': dscfg['RadarName'][ind_radar_list[1]]})
            dscfg['initialized'] = 1

        refl1 = radar1.fields[rad1_refl_field]['data']
        refl2 = radar2.fields[rad2_refl_field]['data']

        phidp1 = radar1.fields[rad1_phidp_field]['data']
        phidp2 = radar2.fields[rad2_phidp_field]['data']

        flag1 = radar1.fields[rad1_flag_field]['data']
        flag2 = radar2.fields[rad2_flag_field]['data']

        intercomp_dict = {
            'rad1_ray_ind': [],
            'rad1_rng_ind': [],
            'rad1_ele': [],
            'rad1_azi': [],
            'rad1_rng': [],
            'rad1_dBZavg': [],
            'rad1_PhiDPavg': [],
            'rad1_Flagavg': [],
            'rad2_ray_ind': [],
            'rad2_rng_ind': [],
            'rad2_ele': [],
            'rad2_azi': [],
            'rad2_rng': [],
            'rad2_dBZavg': [],
            'rad2_PhiDPavg': [],
            'rad2_Flagavg': []}

        # determine if radar data has to be averaged
        avg_rad1, avg_rad2, avg_rad_lim = get_range_bins_to_avg(
            radar1.range['data'], radar2.range['data'])

        # rays are indexed to regular grid
        rays_are_indexed = False
        if 'rays_are_indexed' in dscfg:
            rays_are_indexed = dscfg['rays_are_indexed']

        if not rays_are_indexed:
            azi_tol = 0.5
            ele_tol = 0.5
            rng_tol = 50.

            if 'azi_tol' in dscfg:
                azi_tol = dscfg['azi_tol']
            if 'ele_tol' in dscfg:
                ele_tol = dscfg['ele_tol']
            if 'rng_tol' in dscfg:
                rng_tol = dscfg['rng_tol']

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

                refl1_val = np.ma.asarray(refl1[ind_ray_rad1, ind_rng_rad1])
                refl2_val = np.ma.asarray(refl2[ind_ray_rad2, ind_rng_rad2])

                phidp1_val = np.ma.asarray(
                    phidp1[ind_ray_rad1, ind_rng_rad1])
                phidp2_val = np.ma.asarray(
                    phidp2[ind_ray_rad2, ind_rng_rad2])

                flag1_val = flag1[ind_ray_rad1, ind_rng_rad1]
                flag2_val = flag2[ind_ray_rad2, ind_rng_rad2]

                if avg_rad1:
                    if (ind_rng_rad1+avg_rad_lim[1] >= radar1.ngates or
                            ind_rng_rad1+avg_rad_lim[0] < 0):
                        continue
                    ind_rng = list(range(
                        ind_rng_rad1+avg_rad_lim[0],
                        ind_rng_rad1+avg_rad_lim[1]+1))
                    refl1_val = np.ma.asarray(np.ma.mean(
                        refl1[ind_ray_rad1, ind_rng]))
                    phidp1_val = np.ma.asarray(np.ma.mean(
                        phidp1[ind_ray_rad1, ind_rng]))

                    rad1_flag = flag1[ind_ray_rad1, ind_rng_rad1]

                    rad1_excess_phi = rad1_flag % 100
                    rad1_clt = ((rad1_flag-rad1_excess_phi) % 10000) / 100
                    rad1_prec = (
                        ((rad1_flag-rad1_clt*100-rad1_excess_phi) % 1000000) /
                        10000)

                    flag1_val = int(
                        10000*np.max(rad1_prec)+100*np.max(rad1_clt) +
                        np.max(rad1_excess_phi))

                elif avg_rad2:
                    if (ind_rng_rad2+avg_rad_lim[1] >= radar2.ngates or
                            ind_rng_rad2+avg_rad_lim[0] < 0):
                        continue
                    ind_rng = list(range(
                        ind_rng_rad2+avg_rad_lim[0],
                        ind_rng_rad2+avg_rad_lim[1]+1))
                    refl2_val = np.ma.asarray(np.ma.mean(
                        refl2[ind_ray_rad2, ind_rng]))
                    phidp2_val = np.ma.asarray(np.ma.mean(
                        phidp2[ind_ray_rad1, ind_rng]))

                    rad2_flag = flag2[ind_ray_rad1, ind_rng_rad1]

                    rad2_excess_phi = rad2_flag % 100
                    rad2_clt = ((rad2_flag-rad2_excess_phi) % 10000) / 100
                    rad2_prec = (
                        ((rad2_flag-rad2_clt*100-rad2_excess_phi) % 1000000) /
                        10000)

                    flag2_val = int(
                        10000*np.max(rad2_prec)+100*np.max(rad2_clt) +
                        np.max(rad2_excess_phi))

                if (refl1_val.mask or refl2_val.mask or phidp1_val.mask or
                        phidp2_val.mask):
                    continue

                intercomp_dict['rad1_ray_ind'].append(ind_ray_rad1)
                intercomp_dict['rad1_rng_ind'].append(ind_rng_rad1)
                intercomp_dict['rad1_ele'].append(
                    radar1.elevation['data'][ind_ray_rad1])
                intercomp_dict['rad1_azi'].append(
                    radar1.azimuth['data'][ind_ray_rad1])
                intercomp_dict['rad1_rng'].append(
                    radar1.range['data'][ind_rng_rad1])
                intercomp_dict['rad1_dBZavg'].append(refl1_val)
                intercomp_dict['rad1_PhiDPavg'].append(phidp1_val)
                intercomp_dict['rad1_Flagavg'].append(flag1_val)

                intercomp_dict['rad2_ray_ind'].append(ind_ray_rad2)
                intercomp_dict['rad2_rng_ind'].append(ind_rng_rad2)
                intercomp_dict['rad2_ele'].append(
                    radar2.elevation['data'][ind_ray_rad2])
                intercomp_dict['rad2_azi'].append(
                    radar2.azimuth['data'][ind_ray_rad2])
                intercomp_dict['rad2_rng'].append(
                    radar2.range['data'][ind_rng_rad2])
                intercomp_dict['rad2_dBZavg'].append(refl2_val)
                intercomp_dict['rad2_PhiDPavg'].append(phidp2_val)
                intercomp_dict['rad2_Flagavg'].append(flag2_val)
        else:
            rad1_ray_ind = deepcopy(dscfg['global_data']['rad1_ray_ind'])
            rad1_rng_ind = deepcopy(dscfg['global_data']['rad1_rng_ind'])
            rad2_ray_ind = deepcopy(dscfg['global_data']['rad2_ray_ind'])
            rad2_rng_ind = deepcopy(dscfg['global_data']['rad2_rng_ind'])

            refl1_vec = refl1[rad1_ray_ind, rad1_rng_ind]
            phidp1_vec = phidp1[rad1_ray_ind, rad1_rng_ind]
            flag1_vec = flag1[rad1_ray_ind, rad1_rng_ind]

            refl2_vec = refl2[rad2_ray_ind, rad2_rng_ind]
            phidp2_vec = phidp2[rad2_ray_ind, rad2_rng_ind]
            flag2_vec = flag2[rad2_ray_ind, rad2_rng_ind]

            mask_refl1 = np.ma.getmaskarray(refl1_vec)
            mask_phidp1 = np.ma.getmaskarray(phidp1_vec)
            mask_refl2 = np.ma.getmaskarray(refl2_vec)
            mask_phidp2 = np.ma.getmaskarray(phidp2_vec)

            isvalid = np.logical_not(
                np.logical_and(np.logical_and(mask_refl1, mask_refl2),
                               np.logical_and(mask_phidp1, mask_phidp2)))

            refl1_vec = refl1_vec[isvalid]
            phidp1_vec = phidp1_vec[isvalid]
            flag1_vec = flag1_vec[isvalid]
            refl2_vec = refl2_vec[isvalid]
            phidp2_vec = phidp2_vec[isvalid]
            flag2_vec = flag2_vec[isvalid]
            rad1_ray_ind = rad1_ray_ind[isvalid]
            rad1_rng_ind = rad1_rng_ind[isvalid]
            rad2_ray_ind = rad2_ray_ind[isvalid]
            rad2_rng_ind = rad2_rng_ind[isvalid]

            intercomp_dict['rad1_ray_ind'] = rad1_ray_ind
            intercomp_dict['rad1_rng_ind'] = rad1_rng_ind
            intercomp_dict['rad1_ele'] = radar1.elevation['data'][rad1_ray_ind]
            intercomp_dict['rad1_azi'] = radar1.azimuth['data'][rad1_ray_ind]
            intercomp_dict['rad1_rng'] = radar1.range['data'][rad1_rng_ind]
            intercomp_dict['rad1_dBZavg'] = refl1_vec
            intercomp_dict['rad1_PhiDPavg'] = phidp1_vec
            intercomp_dict['rad1_Flagavg'] = flag1_vec

            intercomp_dict['rad2_ray_ind'] = rad2_ray_ind
            intercomp_dict['rad2_rng_ind'] = rad2_rng_ind
            intercomp_dict['rad2_ele'] = radar2.elevation['data'][rad2_ray_ind]
            intercomp_dict['rad2_azi'] = radar2.azimuth['data'][rad2_ray_ind]
            intercomp_dict['rad2_rng'] = radar2.range['data'][rad2_rng_ind]
            intercomp_dict['rad2_dBZavg'] = refl2_vec
            intercomp_dict['rad2_PhiDPavg'] = phidp2_vec
            intercomp_dict['rad2_Flagavg'] = flag2_vec

        new_dataset = {'intercomp_dict': intercomp_dict,
                       'timeinfo': dscfg['global_data']['timeinfo'],
                       'final': False}
        return new_dataset, None

    if procstatus == 2:
        # get field name
        refl_type = None
        for datatypedescr in dscfg['datatype']:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            if ((datatype == 'dBZ') or (datatype == 'dBZc') or
                    (datatype == 'dBuZ') or (datatype == 'dBZv') or
                    (datatype == 'dBZvc') or (datatype == 'dBuZv')):
                refl_type = datatype
                break

        if refl_type is None:
            warn('Unknown reflectivity type')
            return None, None

        savedir = get_save_dir(
            dscfg['basepath'], dscfg['procname'], dscfg['dsname'],
            dscfg['coloc_data_dir'],
            timeinfo=dscfg['global_data']['timeinfo'], create_dir=False)

        fname = make_filename(
            'colocated_data', dscfg['type'], refl_type, ['csv'],
            timeinfo=dscfg['global_data']['timeinfo'], timeformat='%Y%m%d')

        fname = savedir+fname[0]

        (rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng, rad1_dBZ,
         rad1_phi, rad1_flag, rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi,
         rad2_rng, rad2_dBZ, rad2_phi, rad2_flag) = (
            read_colocated_data_time_avg(fname))

        rad1_excess_phi = (rad1_flag % 100).astype(int)
        rad2_excess_phi = (rad2_flag % 100).astype(int)

        rad1_clt = (((rad1_flag-rad1_excess_phi) % 10000) / 100).astype(int)
        rad2_clt = (((rad2_flag-rad2_excess_phi) % 10000) / 100).astype(int)

        rad1_non_rain = (
            ((rad1_flag-rad1_clt*100-rad1_excess_phi) % 1000000) /
            10000).astype(int)
        rad2_non_rain = (
            ((rad2_flag-rad2_clt*100-rad2_excess_phi) % 1000000) /
            10000).astype(int)

        clt_max = 100
        phi_excess_max = 100
        non_rain_max = 100
        phi_avg_max = 600.

        if 'clt_max' in dscfg:
            clt_max = dscfg['clt_max']
        if 'phi_excess_max' in dscfg:
            phi_excess_max = dscfg['phi_excess_max']
        if 'non_rain_max' in dscfg:
            non_rain_max = dscfg['non_rain_max']
        if 'phi_avg_max' in dscfg:
            phi_avg_max = dscfg['phi_avg_max']

        # filter out invalid data
        ind_val = np.where(
            np.logical_and.reduce((
                rad1_clt <= clt_max, rad2_clt <= clt_max,
                rad1_excess_phi <= phi_excess_max,
                rad2_excess_phi <= phi_excess_max,
                rad1_non_rain <= non_rain_max, rad2_non_rain <= non_rain_max,
                rad1_phi <= phi_avg_max, rad2_phi <= phi_avg_max)))[0]

        intercomp_dict = {
            'rad1_name': dscfg['global_data']['rad1_name'],
            'rad1_ray_ind': rad1_ray_ind[ind_val],
            'rad1_rng_ind': rad1_rng_ind[ind_val],
            'rad1_ele': rad1_ele[ind_val],
            'rad1_azi': rad1_azi[ind_val],
            'rad1_rng': rad1_rng[ind_val],
            'rad1_val': rad1_dBZ[ind_val],
            'rad2_name': dscfg['global_data']['rad2_name'],
            'rad2_ray_ind': rad1_ray_ind[ind_val],
            'rad2_rng_ind': rad1_rng_ind[ind_val],
            'rad2_ele': rad2_ele[ind_val],
            'rad2_azi': rad2_azi[ind_val],
            'rad2_rng': rad2_rng[ind_val],
            'rad2_val': rad2_dBZ[ind_val]}

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
            minimum range where to look for a sun hit signal [m].
            Default 50000.
        hmin : float. Dataset keyword
            minimum altitude where to look for a sun hit signal [m MSL].
            Default 10000. The actual range from which a sun hit signal will
            be search will be the minimum between rmin and the range from
            which the altitude is higher than hmin.
        delev_max : float. Dataset keyword
            maximum elevation distance from nominal radar elevation where to
            look for a sun hit signal [deg]. Default 1.5
        dazim_max : float. Dataset keyword
            maximum azimuth distance from nominal radar elevation where to
            look for a sun hit signal [deg]. Default 1.5
        elmin : float. Dataset keyword
            minimum radar elevation where to look for sun hits [deg].
            Default 1.
        nbins_min : int. Dataset keyword.
            minimum number of range bins that have to contain signal to
            consider the ray a potential sun hit. Default 10.
        attg : float. Dataset keyword
            gaseous attenuation. Default None
        max_std_pwr : float. Dataset keyword
            maximum standard deviation of the signal power to consider the
            data a sun hit [dB]. Default 2.
        max_std_zdr : float. Dataset keyword
            maximum standard deviation of the ZDR to consider the
            data a sun hit [dB]. Default 2.
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
        coeff_band : float. Dataset keyword
            multiplicate coefficient to transform pulse width into receiver
            bandwidth
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
            if datatype == 'ZDR':
                zdr_field = 'differential_reflectivity'

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
            else:
                warn('Angular resolution unknown.')

            radar_par.update({'timeinfo': dscfg['timeinfo']})
            dscfg['global_data'] = radar_par
            dscfg['initialized'] = 1

        dscfg['global_data']['timeinfo'] = dscfg['timeinfo']

        # default values
        rmin = 50000.
        hmin = 10000.
        delev_max = 1.5
        dazim_max = 1.5
        elmin = 1.
        nbins_min = 20
        attg = None
        max_std_pwr = 2.
        max_std_zdr = 2.

        # user values
        if 'rmin' in dscfg:
            rmin = dscfg['rmin']
        if 'hmin' in dscfg:
            hmin = dscfg['hmin']
        if 'delev_max' in dscfg:
            delev_max = dscfg['delev_max']
        if 'dazim_max' in dscfg:
            dazim_max = dscfg['dazim_max']
        if 'elmin' in dscfg:
            elmin = dscfg['elmin']
        if 'nbins_min' in dscfg:
            nbins_min = dscfg['nbins_min']
        if 'attg' in dscfg:
            attg = dscfg['attg']
        if 'max_std_pwr' in dscfg:
            max_std_pwr = dscfg['max_std_pwr']
        if 'max_std_zdr' in dscfg:
            max_std_zdr = dscfg['max_std_zdr']

        ind_rmin = np.where(radar.range['data'] > rmin)[0][0]

        sun_hits, new_radar = pyart.correct.get_sun_hits(
            radar, delev_max=delev_max, dazim_max=dazim_max, elmin=elmin,
            rmin=rmin, hmin=hmin, nbins_min=nbins_min,
            max_std_pwr=max_std_pwr, max_std_zdr=max_std_zdr,
            attg=attg, pwrh_field=pwrh_field, pwrv_field=pwrv_field,
            zdr_field=zdr_field)

        if sun_hits is None:
            return None, None

        sun_hits_dataset = dict()
        sun_hits_dataset.update({'sun_hits': sun_hits})
        sun_hits_dataset.update({'radar': new_radar})
        sun_hits_dataset.update({'timeinfo': dscfg['timeinfo']})

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
        coeff_band = 1.2

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
        if 'coeff_band' in dscfg:
            coeff_band = dscfg['coeff_band']

        sun_hits = read_sun_hits_multiple_days(
            dscfg, dscfg['global_data']['timeinfo'], nfiles=nfiles)

        if sun_hits[0] is None:
            return None, None

        sun_pwr_h = sun_hits[7]
        sun_pwr_v = sun_hits[11]

        # get DRAO reference
        sf_ref = np.ma.asarray(np.ma.masked)
        ref_time = None
        if 'wavelen' in dscfg['global_data']:

            flx_dt, flx_val = read_solar_flux(
                dscfg['solarfluxpath']+'fluxtable.txt')

            if flx_dt is not None:
                flx_dt_closest, flx_val_closest = get_closest_solar_flux(
                    sun_hits[0], flx_dt, flx_val)

                # flux at radar wavelength
                sf_radar = pyart.correct.solar_flux_lookup(
                    flx_val_closest, dscfg['global_data']['wavelen'])

                sf_ref = np.ma.asarray(sf_radar[-1])
                ref_time = flx_dt_closest[-1]

                # scaling of the power to account for solar flux variations.
                # The last sun hit is the reference. The scale factor is in dB
                scale_factor = -10.*np.log10(sf_radar/sf_ref)
                sun_pwr_h += scale_factor
                sun_pwr_v += scale_factor
            else:
                warn('Unable to compute solar power reference. ' +
                     'Missing DRAO data')
        else:
            warn('Unable to compute solar power reference. ' +
                 'Missing radar wavelength')

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
            'sf_h': np.ma.asarray(np.ma.masked),
            'az_bias_h': np.ma.asarray(np.ma.masked),
            'el_bias_h': np.ma.asarray(np.ma.masked),
            'az_width_h': np.ma.asarray(np.ma.masked),
            'el_width_h': np.ma.asarray(np.ma.masked),
            'nhits_h': 0,
            'par_h': None,
            'dBmv_sun_est': np.ma.asarray(np.ma.masked),
            'std(dBmv_sun_est)': np.ma.asarray(np.ma.masked),
            'sf_v': np.ma.asarray(np.ma.masked),
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
            'sf_ref': np.ma.asarray(sf_ref),
            'ref_time': ref_time,
            'lant': np.ma.asarray(np.ma.masked)}

        if sun_retrieval_h is not None:
            # correct for scanning losses and the polarization of the antenna
            if (('angle_step' in dscfg['global_data']) and
                    ('beamwidth' in dscfg['global_data'])):
                lant = pyart.correct.scanning_losses(
                    dscfg['global_data']['angle_step'],
                    dscfg['global_data']['beamwidth'])

            else:
                warn('Unable to estimate scanning losses. ' +
                     'Missing radar parameters. ' +
                     'Antenna losses will be neglected')
                lant = 0.
            ptoa_h = sun_retrieval_h[0]+lant+3.

            # compute observed solar flux
            if (('pulse_width' in dscfg['global_data']) and
                    ('wavelen' in dscfg['global_data']) and
                    (dscfg['AntennaGain'] is not None)):
                sf_h = pyart.correct.ptoa_to_sf(
                    ptoa_h, dscfg['global_data']['pulse_width'],
                    dscfg['global_data']['wavelen'],
                    dscfg['AntennaGain'],
                    )
            else:
                warn('Unable to estimate observed solar flux. ' +
                     'Missing radar parameters')
                sf_h = np.ma.asarray(np.ma.masked)

            sun_retrieval_dict['dBm_sun_est'] = np.ma.asarray(ptoa_h)
            sun_retrieval_dict['std(dBm_sun_est)'] = np.ma.asarray(
                sun_retrieval_h[1])
            sun_retrieval_dict['sf_h'] = np.ma.asarray(sf_h)
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
            sun_retrieval_dict['lant'] = np.ma.asarray(lant)
        if sun_retrieval_v is not None:
            # correct for scanning losses and the polarization of the antenna
            if (('angle_step' in dscfg['global_data']) and
                    ('beamwidth' in dscfg['global_data'])):
                lant = pyart.correct.scanning_losses(
                    dscfg['global_data']['angle_step'],
                    dscfg['global_data']['beamwidth'])
            else:
                lant = 0.
                warn('Unable to estimate scanning losses. ' +
                     'Missing radar parameters. ' +
                     'Antenna losses will be neglected')
            ptoa_v = sun_retrieval_v[0]+lant+3.

            # compute observed solar flux
            if (('pulse_width' in dscfg['global_data']) and
                    ('wavelen' in dscfg['global_data']) and
                    (dscfg['AntennaGain'] is not None)):
                sf_v = pyart.correct.ptoa_to_sf(
                    ptoa_v, dscfg['global_data']['pulse_width'],
                    dscfg['global_data']['wavelen'],
                    dscfg['AntennaGain'],
                    )
            else:
                warn('Unable to estimate observed solar flux. ' +
                     'Missing radar parameters')
                sf_v = np.ma.asarray(np.ma.masked)

            sun_retrieval_dict['dBmv_sun_est'] = np.ma.asarray(ptoa_v)
            sun_retrieval_dict['std(dBmv_sun_est)'] = np.ma.asarray(
                sun_retrieval_v[1])
            sun_retrieval_dict['sf_v'] = np.ma.asarray(sf_v)
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
            sun_retrieval_dict['lant'] = np.ma.asarray(lant)
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
        sun_hits_dataset.update(
            {'timeinfo': dscfg['global_data']['timeinfo']})

        return sun_hits_dataset, ind_rad
