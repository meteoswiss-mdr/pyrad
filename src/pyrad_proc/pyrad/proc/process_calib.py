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
    process_sun_hits

"""

from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.read_data_other import read_selfconsistency
from ..io.read_data_other import read_sun_hits_multiple_days, read_solar_flux
from ..io.read_data_radar import interpol_field

from ..util.radar_utils import get_closest_solar_flux, get_histogram_bins


def process_correct_bias(procstatus, dscfg, radar=None):
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

    radar : Radar
        Optional. Radar object

    Returns
    -------
    new_dataset : Radar
        radar object

    """

    if procstatus != 1:
        return None

    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        break
    field_name = get_fieldname_pyart(datatype)

    if field_name not in radar.fields:
        warn('Unable to correct for bias field ' + field_name +
             '. Field not available')
        return None

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

    return new_dataset


def process_correct_noise_rhohv(procstatus, dscfg, radar=None):
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

    radar : Radar
        Optional. Radar object

    Returns
    -------
    new_dataset : Radar
        radar object

    """

    if procstatus != 1:
        return None

    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatype_fields(
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

    if ((urhohv not in radar.fields) or
            (snr not in radar.fields) or
            (zdr not in radar.fields) or
            (nh not in radar.fields) or
            (nv not in radar.fields)):
        warn('Unable to correct RhoHV field for noise. Missing fields')
        return None

    rhohv = pyart.correct.correct_noise_rhohv(
        radar, urhohv_field=urhohv, snr_field=snr, zdr_field=zdr,
        nh_field=nh, nv_field=nv, rhohv_field='cross_correlation_ratio')

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field('cross_correlation_ratio', rhohv)

    return new_dataset


def process_selfconsistency_kdp_phidp(procstatus, dscfg, radar=None):
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

    radar : Radar
        Optional. Radar object

    Returns
    -------
    radar : Radar
        radar object

    """

    if procstatus != 1:
        return None

    temp = None
    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatype_fields(
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

    if ((refl not in radar.fields) or
            (zdr not in radar.fields) or
            (phidp not in radar.fields) or
            (rhohv not in radar.fields)):
        warn('Unable to estimate reflectivity bias using selfconsistency. ' +
             'Missing data')
        return None

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

    return new_dataset


def process_selfconsistency_bias(procstatus, dscfg, radar=None):
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

    radar : Radar
        Optional. Radar object

    Returns
    -------
    radar : Radar
        radar object

    """

    if procstatus != 1:
        return None

    temp = None
    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatype_fields(
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

    if ((refl not in radar.fields) or
            (zdr not in radar.fields) or
            (phidp not in radar.fields) or
            (rhohv not in radar.fields)):
        warn('Unable to estimate reflectivity bias using selfconsistency. ' +
             'Missing data')
        return None

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

    return new_dataset


def process_estimate_phidp0(procstatus, dscfg, radar=None):
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

    radar : Radar
        Optional. Radar object

    Returns
    -------
    new_dataset : Radar
        radar object

    """

    if procstatus != 1:
        return None

    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatype_fields(
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

    if (refl_field not in radar.fields) or (psidp_field not in radar.fields):
        warn('Unable to estimate PhiDP system offset. Missing data')
        return None

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

    return new_dataset


def process_rhohv_rain(procstatus, dscfg, radar=None):
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

    radar : Radar
        Optional. Radar object

    Returns
    -------
    radar : Radar
        radar object

    """

    if procstatus != 1:
        return None

    temp_field = None
    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatype_fields(
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

    if ((refl_field not in radar.fields) or
            (rhohv_field not in radar.fields)):
        warn('Unable to estimate RhoHV in rain. Missing data')
        return None

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

    return new_dataset


def process_zdr_rain(procstatus, dscfg, radar=None):
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

    radar : Radar
        Optional. Radar object

    Returns
    -------
    radar : Radar
        radar object

    """

    if procstatus != 1:
        return None

    temp_field = None
    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatype_fields(
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

    if ((refl_field not in radar.fields) or
            (rhohv_field not in radar.fields)):
        warn('Unable to estimate ZDR in rain. Missing data')
        return None

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

    return new_dataset


def process_monitoring(procstatus, dscfg, radar=None):
    """
    computes statistics

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types

    radar : Radar
        Optional. Radar object

    Returns
    -------
    radar : Radar
        radar object

    """

    if procstatus == 0:
        return None

    if procstatus == 1:
        for datatypedescr in dscfg['datatype']:
            datagroup, datatype, dataset, product = get_datatype_fields(
                datatypedescr)
            field_name = get_fieldname_pyart(datatype)
            break

        if field_name not in radar.fields:
            warn(field_name+' not available.')
            return None

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

        return dataset

    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None

        dataset = dict()
        dataset.update({'hist_obj': dscfg['global_data']})
        dataset.update({'hist_type': 'cumulative'})

        return dataset


def process_sun_hits(procstatus, dscfg, radar=None):
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

    radar : Radar
        Optional. Radar object

    Returns
    -------
    sun_hits_dict : dict
        dictionary containing a radar object, a sun_hits dict and a
        sun_retrieval dictionary

    """

    if procstatus == 0:
        return None

    if procstatus == 1:
        for datatypedescr in dscfg['datatype']:
            datagroup, datatype, dataset, product = get_datatype_fields(
                datatypedescr)
            if datatype == 'dBm':
                pwrh_field = 'signal_power_hh'
            if datatype == 'dBmv':
                pwrv_field = 'signal_power_vv'
            if datatype == 'ZDRu':
                zdr_field = 'unfiltered_differential_reflectivity'
            if datatype == 'ZDRuc':
                zdr_field = 'corrected_unfiltered_differential_reflectivity'

        if ((pwrh_field not in radar.fields) or
                (pwrv_field not in radar.fields) or
                (zdr_field not in radar.fields)):
            warn('Unable to get sun hits. Missing data')
            return None

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
            return None

        sun_hits_dataset = dict()
        sun_hits_dataset.update({'sun_hits': sun_hits})
        sun_hits_dataset.update({'radar': new_radar})

        return sun_hits_dataset

    if procstatus == 2:
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
            return None

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

        return sun_hits_dataset
