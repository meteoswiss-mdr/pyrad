"""
pyrad.proc.process_phase
======================================

Functions for PhiDP and KDP processing and attenuation correction

.. autosummary::
    :toctree: generated/

    process_correct_phidp0
    process_smooth_phidp_single_window
    process_smooth_phidp_double_window
    process_kdp_leastsquare_single_window
    process_kdp_leastsquare_double_window
    process_phidp_kdp_Maesaka
    process_phidp_kdp_lp
    process_selfconsistency_kdp_phidp
    process_selfconsistency_bias
    process_attenuation

"""

from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields


def process_correct_phidp0(procstatus, dscfg, radar_list=None):
    """
    corrects phidp of the system phase

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
        warn('Unable to correct PhiDP system offset. Missing data')
        return None, None

    ind_rmin = np.where(radar.range['data'] > dscfg['rmin'])[0][0]
    ind_rmax = np.where(radar.range['data'] < dscfg['rmax'])[0][-1]
    r_res = radar.range['data'][1]-radar.range['data'][0]
    min_rcons = int(dscfg['rcell']/r_res)

    if psidp_field.startswith('corrected_'):
        phidp_field = psidp_field
    elif psidp_field.startswith('uncorrected_'):
        phidp_field = psidp_field.replace('uncorrected_', 'corrected_', 1)
    else:
        phidp_field = 'corrected_'+psidp_field

    phidp = pyart.correct.correct_sys_phase(
        radar, ind_rmin=ind_rmin, ind_rmax=ind_rmax, min_rcons=min_rcons,
        zmin=dscfg['Zmin'], zmax=dscfg['Zmax'], psidp_field=psidp_field,
        refl_field=refl_field, phidp_field=phidp_field)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field(phidp_field, phidp)

    return new_dataset, ind_rad


def process_smooth_phidp_single_window(procstatus, dscfg, radar_list=None):
    """
    corrects phidp of the system phase and smoothes it using one window

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
        rwind : float. Dataset keyword
            The length of the smoothing window [m]
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
        warn('Unable to smooth PhiDP. Missing data')
        return None, None

    ind_rmin = np.where(radar.range['data'] > dscfg['rmin'])[0][0]
    ind_rmax = np.where(radar.range['data'] < dscfg['rmax'])[0][-1]
    r_res = radar.range['data'][1]-radar.range['data'][0]
    min_rcons = int(dscfg['rcell']/r_res)
    wind_len = int(dscfg['rwind']/r_res)
    min_valid = int(wind_len/2+1)

    if psidp_field.startswith('corrected_'):
        phidp_field = psidp_field
    elif psidp_field.startswith('uncorrected_'):
        phidp_field = psidp_field.replace('uncorrected_', 'corrected_', 1)
    else:
        phidp_field = 'corrected_'+psidp_field

    phidp = pyart.correct.smooth_phidp_single_window(
        radar, ind_rmin=ind_rmin, ind_rmax=ind_rmax, min_rcons=min_rcons,
        zmin=dscfg['Zmin'], zmax=dscfg['Zmax'], wind_len=wind_len,
        min_valid=min_valid, psidp_field=psidp_field, refl_field=refl_field,
        phidp_field=phidp_field)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field(phidp_field, phidp)

    return new_dataset, ind_rad


def process_smooth_phidp_double_window(procstatus, dscfg, radar_list=None):
    """
    corrects phidp of the system phase and smoothes it using one window

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
        rwinds : float. Dataset keyword
            The length of the short smoothing window [m]
        rwindl : float. Dataset keyword
            The length of the long smoothing window [m]
        Zmin : float. Dataset keyword
            The minimum reflectivity [dBZ]
        Zmax : float. Dataset keyword
            The maximum reflectivity [dBZ]
        Zthr : float. Dataset keyword
            The threshold defining wich smoothed data to used [dBZ]
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
        warn('Unable to smooth PhiDP. Missing data')
        return None, None

    ind_rmin = np.where(radar.range['data'] > dscfg['rmin'])[0][0]
    ind_rmax = np.where(radar.range['data'] < dscfg['rmax'])[0][-1]
    r_res = radar.range['data'][1]-radar.range['data'][0]
    min_rcons = int(dscfg['rcell']/r_res)
    swind_len = int(dscfg['rwinds']/r_res)
    smin_valid = int(swind_len/2+1)
    lwind_len = int(dscfg['rwindl']/r_res)
    lmin_valid = int(lwind_len/2+1)

    if psidp_field.startswith('corrected_'):
        phidp_field = psidp_field
    elif psidp_field.startswith('uncorrected_'):
        phidp_field = psidp_field.replace('uncorrected_', 'corrected_', 1)
    else:
        phidp_field = 'corrected_'+psidp_field

    phidp = pyart.correct.smooth_phidp_double_window(
        radar, ind_rmin=ind_rmin, ind_rmax=ind_rmax, min_rcons=min_rcons,
        zmin=dscfg['Zmin'], zmax=dscfg['Zmax'], swind_len=swind_len,
        smin_valid=smin_valid, lwind_len=lwind_len, lmin_valid=lmin_valid,
        zthr=dscfg['Zthr'], psidp_field=psidp_field, refl_field=refl_field,
        phidp_field=phidp_field)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field(phidp_field, phidp)

    return new_dataset, ind_rad


def process_phidp_kdp_Maesaka(procstatus, dscfg, radar_list=None):
    """
    Estimates PhiDP and KDP using the method by Maesaka

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
        if datatype == 'PhiDP':
            psidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            psidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            psidp_field = 'uncorrected_differential_phase'
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'TEMP':
            temp_field = 'temperature'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl_field not in radar.fields) or
            (psidp_field not in radar.fields) or
            (temp_field not in radar.fields)):
        warn('Unable to retrieve PhiDP KDP using the Maesaka approach. ' +
             'Missing data')
        return None, None

    phidp_field = 'corrected_differential_phase'
    kdp_field = 'corrected_specific_differential_phase'

    radar_aux = deepcopy(radar)

    # correct PhiDP0
    ind_rmin = np.where(radar.range['data'] > dscfg['rmin'])[0][0]
    ind_rmax = np.where(radar.range['data'] < dscfg['rmax'])[0][-1]
    r_res = radar.range['data'][1]-radar.range['data'][0]
    min_rcons = int(dscfg['rcell']/r_res)

    phidp = pyart.correct.correct_sys_phase(
        radar_aux, ind_rmin=ind_rmin, ind_rmax=ind_rmax, min_rcons=min_rcons,
        zmin=dscfg['Zmin'], zmax=dscfg['Zmax'], psidp_field=psidp_field,
        refl_field=refl_field, phidp_field=phidp_field)

    mask = np.ma.getmaskarray(phidp['data'])
    # filter out data in an above the melting layer
    if 'radar_beam_width_h' in radar.instrument_parameters:
        beamwidth = (
            radar.instrument_parameters['radar_beam_width_h']['data'][0])
    else:
        warn('Unknown radar antenna beamwidth.')
        beamwidth = None

    mask_fzl, end_gate_arr = get_mask_fzl(
        radar_aux, fzl=None, doc=None, min_temp=0., thickness=400.,
        beamwidth=beamwidth, temp_field=temp_field)
    mask = np.logical_or(mask, mask_refl)

    phidp['data'] = np.ma.masked_where(mask, phidp['data'])
    fill_value = pyart.config.get_fillvalue()
    radar_aux.add_field(phidp_field, phidp, replace_existing=True)

    # the return data is not a masked array
    kdp, phidpf, phidpr = pyart.retrieve.kdp_proc.kdp_maesaka(
        radar_aux, gatefilter=None, method='cg', backscatter=None, Clpf=1.,
        length_scale=None, first_guess=0.01, finite_order='low',
        fill_value=fill_value, psidp_field=phidp_field, kdp_field=kdp_field,
        phidp_field=phidp_field)

    kdp['data'] = np.ma.masked_where(mask, kdp['data'])
    phidpf['data'] = np.ma.masked_where(mask, phidpf['data'])

    # prepare for exit
    new_dataset = deepcopy(radar_aux)
    new_dataset.fields = dict()
    new_dataset.add_field(phidp_field, phidpf)
    new_dataset.add_field(kdp_field, kdp)

    return new_dataset, ind_rad


def process_phidp_kdp_lp(procstatus, dscfg, radar_list=None):
    """
    Estimates PhiDP and KDP using a linear programming algorithm

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
        if datatype == 'PhiDP':
            psidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            psidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            psidp_field = 'uncorrected_differential_phase'
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'
        if datatype == 'RhoHVc':
            rhv_field = 'corrected_cross_correlation_ratio'
        if datatype == 'RhoHV':
            rhv_field = 'cross_correlation_ratio'
        if datatype == 'SNRh':
            snr_field = 'signal_to_noise_ratio_hh'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl_field not in radar.fields) or
            (psidp_field not in radar.fields) or
            (rhv_field not in radar.fields) or
            (snr_field not in radar.fields)):
        warn('Unable to retrieve PhiDP and KDP using the LP approach. ' +
             'Missing data')
        return None, None

    mask = np.ma.getmaskarray(radar.fields[psidp_field]['data'])

    phidp_field = 'corrected_differential_phase'
    kdp_field = 'corrected_specific_differential_phase'

    phidp, kdp = pyart.correct.phase_proc_lp(
        radar, 0, debug=False, self_const=60000.0,
        low_z=10.0, high_z=53.0, min_phidp=0.01, min_ncp=10.,
        min_rhv=0.6, fzl=4000.0, sys_phase=0.0,
        overide_sys_phase=False, nowrap=None, really_verbose=False,
        LP_solver='cvxopt', refl_field=refl_field, ncp_field=snr_field,
        rhv_field=rhv_field, phidp_field=psidp_field, kdp_field=kdp_field,
        unf_field=phidp_field, window_len=35, proc=1)

    kdp['data'] = np.ma.masked_where(mask, kdp['data'])
    phidp['data'] = np.ma.masked_where(mask, phidp['data'])

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field(phidp_field, phidp)
    new_dataset.add_field(kdp_field, kdp)

    return new_dataset, ind_rad


def process_kdp_leastsquare_single_window(procstatus, dscfg, radar_list=None):
    """
    Computes specific differential phase using a piecewise least square method

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rwind : float. Dataset keyword
            The length of the segment for the least square method [m]
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
        if datatype == 'PhiDP':
            phidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            phidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            phidp_field = 'uncorrected_differential_phase'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if phidp_field not in radar.fields:
        warn('Unable to retrieve KDP from PhiDP using least square. ' +
             'Missing data')
        return None, None

    r_res = radar.range['data'][1]-radar.range['data'][0]
    wind_len = int(dscfg['rwind']/r_res)
    min_valid = int(wind_len/2+1)
    kdp_field = 'corrected_specific_differential_phase'

    kdp = pyart.retrieve.kdp_leastsquare_single_window(
        radar, wind_len=wind_len, min_valid=min_valid, phidp_field=phidp_field,
        kdp_field=kdp_field)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field(kdp_field, kdp)

    return new_dataset, ind_rad


def process_kdp_leastsquare_double_window(procstatus, dscfg, radar_list=None):
    """
    Computes specific differential phase using a piecewise least square method

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        rwinds : float. Dataset keyword
            The length of the short segment for the least square method [m]
        rwindl : float. Dataset keyword
            The length of the long segment for the least square method [m]
        Zthr : float. Dataset keyword
            The threshold defining which estimated data to use [dBZ]
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
        if datatype == 'PhiDP':
            phidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            phidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            phidp_field = 'uncorrected_differential_phase'
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBZc':
            refl_field = 'corrected_reflectivity'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if (phidp_field not in radar.fields) or (refl_field not in radar.fields):
        warn('Unable to retrieve KDP from PhiDP using least square. ' +
             'Missing data')
        return None, None

    r_res = radar.range['data'][1]-radar.range['data'][0]
    swind_len = int(dscfg['rwinds']/r_res)
    smin_valid = int(swind_len/2+1)
    lwind_len = int(dscfg['rwindl']/r_res)
    lmin_valid = int(lwind_len/2+1)

    kdp_field = 'corrected_specific_differential_phase'

    kdp = pyart.retrieve.kdp_leastsquare_double_window(
        radar, swind_len=swind_len, smin_valid=smin_valid,
        lwind_len=lwind_len, lmin_valid=lmin_valid, zthr=dscfg['Zthr'],
        phidp_field=phidp_field, refl_field=refl_field, kdp_field=kdp_field)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field(kdp_field, kdp)

    return new_dataset, ind_rad


def process_attenuation(procstatus, dscfg, radar_list=None):
    """
    Computes specific attenuation and specific differential attenuation using
    the Z-Phi method and corrects reflectivity and differential reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        ATT_METHOD : float. Dataset keyword
            The attenuation estimation method used. One of the following:
            ZPhi, Philin
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
        if datatype == 'PhiDPc':
            phidp = 'corrected_differential_phase'
        if datatype == 'ZDRc':
            zdr = 'corrected_differential_reflectivity'
        if datatype == 'dBZ':
            refl = 'reflectivity'
        if datatype == 'PhiDP':
            phidp = 'differential_phase'
        if datatype == 'ZDR':
            zdr = 'differential_reflectivity'
        if datatype == 'TEMP':
            temp = 'temperature'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((phidp not in radar.fields) or
            (refl not in radar.fields) or
            (zdr not in radar.fields)):
        warn('Unable to retrieve KDP from PhiDP using least square. ' +
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

    att_method = 'ZPhi'
    if 'ATT_METHOD' in dscfg:
        att_method = dscfg['ATT_METHOD']

    if (att_method != 'ZPhi') and (att_method != 'Philin'):
        raise ValueError(
            'Unknown attenuation correction method. ' +
            'Must be one of the following: [ZPhi, Philin]')

    if att_method == 'ZPhi':
        spec_at, cor_z, spec_diff_at, cor_zdr = (
            pyart.correct.calculate_attenuation_zphi(
                radar, doc=15, fzl=fzl, smooth_window_len=0, a_coef=None,
                beta=None, c=None, d=None, refl_field=refl, phidp_field=phidp,
                zdr_field=zdr, temp_field=temp, spec_at_field=None,
                corr_refl_field=None, spec_diff_at_field=None,
                corr_zdr_field=None))
    elif att_method == 'Philin':
        spec_at, cor_z, spec_diff_at, cor_zdr = (
            pyart.correct.calculate_attenuation_philinear(
                radar, doc=15, fzl=fzl, pia_coef=None, pida_coef=None,
                refl_field=refl, phidp_field=phidp, zdr_field=zdr,
                temp_field=temp, spec_at_field=None, corr_refl_field=None,
                spec_diff_at_field=None, corr_zdr_field=None))

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    new_dataset.add_field('specific_attenuation', spec_at)
    new_dataset.add_field('corrected_reflectivity', cor_z)

    if (spec_diff_at is not None) and (cor_zdr is not None):
        new_dataset.add_field(
            'specific_differential_attenuation', spec_diff_at)
        new_dataset.add_field('corrected_differential_reflectivity', cor_zdr)
    else:
        warn(
            ' Specific differential attenuation and attenuation ' +
            'corrected differential reflectivity not available')

    return new_dataset, ind_rad
