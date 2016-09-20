"""
pyrad.proc.process_dataset
==========================

Functions for processing Pyrad datasets

.. autosummary::
    :toctree: generated/

    get_process_type
    process_raw
    process_save_radar
    process_snr
    process_l
    process_cdr
    process_correct_noise_rhohv
    process_correct_bias
    process_echo_id
    process_echo_filter
    process_filter_snr
    process_attenuation
    process_rainrate
    process_hydroclass
    process_correct_phidp0
    process_smooth_phidp_single_window
    process_smooth_phidp_double_window
    process_kdp_leastsquare_single_window
    process_kdp_leastsquare_double_window
    process_phidp_kdp_Maesaka
    process_phidp_kdp_lp

"""

from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from ..io.read_data import get_datatypefields, get_fieldname_rainbow
from netCDF4 import num2date


def get_process_type(dataset_type):
    """
    maps the dataset type into its processing function and data set format

    Parameters
    ----------
    dataset_type : str
        data set type, i.e. 'RAW', 'SAN', etc.

    Returns
    -------
    func_name : str
        pyrad function used to process the data set type

    dsformat : str
        data set format, i.e.: 'VOL', etc.

    """

    dsformat = 'VOL'
    if dataset_type == 'RAW':
        func_name = 'process_raw'
    elif dataset_type == 'NCVOL':
        func_name = 'process_save_radar'
    elif dataset_type == 'SNR':
        func_name = 'process_snr'
    elif dataset_type == 'RHOHV_CORRECTION':
        func_name = 'process_correct_noise_rhohv'
    elif dataset_type == 'BIAS_CORRECTION':
        func_name = 'process_correct_bias'
    elif dataset_type == 'L':
        func_name = 'process_l'
    elif dataset_type == 'CDR':
        func_name = 'process_cdr'
    elif dataset_type == 'SAN':
        func_name = 'process_echo_id'
    elif dataset_type == 'ECHO_FILTER':
        func_name = 'process_echo_filter'
    elif dataset_type == 'SNR_FILTER':
        func_name = 'process_filter_snr'
    elif dataset_type == 'PHIDP0_CORRECTION':
        func_name = 'process_correct_phidp0'
    elif dataset_type == 'PHIDP_SMOOTH_1W':
        func_name = 'process_smooth_phidp_single_window'
    elif dataset_type == 'PHIDP_SMOOTH_2W':
        func_name = 'process_smooth_phidp_double_window'
    elif dataset_type == 'PHIDP_KDP_MAESAKA':
        func_name = 'process_phidp_kdp_Maesaka'
    elif dataset_type == 'PHIDP_KDP_LP':
        func_name = 'process_phidp_kdp_lp'
    elif dataset_type == 'KDP_LEASTSQUARE_1W':
        func_name = 'process_kdp_leastsquare_single_window'
    elif dataset_type == 'KDP_LEASTSQUARE_2W':
        func_name = 'process_kdp_leastsquare_double_window'
    elif dataset_type == 'ATTENUATION':
        func_name = 'process_attenuation'
    elif dataset_type == 'RAINRATE':
        func_name = 'process_rainrate'
    elif dataset_type == 'HYDROCLASS':
        func_name = 'process_hydroclass'
    elif dataset_type == 'POINT_MEASUREMENT':
        func_name = 'process_point_measurement'
        dsformat = 'TIMESERIES'
    else:
        raise ValueError('ERROR: Unknown dataset type')

    return func_name, dsformat


def process_raw(procstatus, dscfg, radar=None):
    """
    dummy function that returns the initial input data set

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

    radar : Radar
        Optional. Radar object

    Returns
    -------
    new_dataset : Radar
        radar object

    """

    if procstatus != 1:
        return None

    new_dataset = deepcopy(radar)
    return new_dataset


def process_save_radar(procstatus, dscfg, radar=None):
    """
    dummy function that allows to save the entire radar object

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

    radar : Radar
        Optional. Radar object

    Returns
    -------
    new_dataset : Radar
        radar object

    """

    if procstatus != 1:
        return None

    new_dataset = deepcopy(radar)
    return new_dataset


def process_correct_bias(procstatus, dscfg, radar=None):
    """
    Corrects a bias on the data

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

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
        datagroup, datatype, dataset, product = get_datatypefields(
            datatypedescr)
        break
    field_name = get_fieldname_rainbow(datatype)

    corrected_field = pyart.correct.correct_bias(
        radar, bias=dscfg['bias'], field_name=field_name)

    if field_name.startswith('corrected_'):
        new_field_name = field_name
    else:
        new_field_name = 'corrected_'+field_name

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field(new_field_name, corrected_field)

    return new_dataset


def process_snr(procstatus, dscfg, radar=None):
    """
    Computes SNR

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

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
        datagroup, datatype, dataset, product = get_datatypefields(
            datatypedescr)
        if datatype == 'dBZ':
            refl = 'reflectivity'
        if datatype == 'dBuZ':
            refl = 'unfiltered_reflectivity'
        if datatype == 'dBZv':
            refl = 'reflectivity_vv'
        if datatype == 'dBuZv':
            refl = 'unfiltered_reflectivity_vv'
        if datatype == 'Nh':
            noise = 'noisedBZ_hh'
        if datatype == 'Nv':
            noise = 'noisedBZ_vv'

    if 'output_type' in dscfg:
        if dscfg['output_type'] == 'SNRh':
            snr_field = 'signal_to_noise_ratio_hh'
        elif dscfg['output_type'] == 'SNRv':
            snr_field = 'signal_to_noise_ratio_vv'
    else:
        snr_field = 'signal_to_noise_ratio_hh'

    snr = pyart.retrieve.compute_snr(
        radar, refl_field=refl, noise_field=noise,
        snr_field=snr_field)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field(snr_field, snr)

    return new_dataset


def process_l(procstatus, dscfg, radar=None):
    """
    Computes L parameter

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

    radar : Radar
        Optional. Radar object

    Returns
    -------
    new_dataset : Radar
        radar object

    """

    if procstatus != 1:
        return None

    datagroup, datatype, dataset, product = get_datatypefields(
        dscfg['datatype'])
    rhohv = get_fieldname_rainbow(datatype)

    l = pyart.retrieve.compute_l(
        radar, rhohv_field=rhohv,
        l_field='logarithmic_cross_correlation_ratio')

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field('logarithmic_cross_correlation_ratio', l)

    return new_dataset


def process_cdr(procstatus, dscfg, radar=None):
    """
    Computes Circular Depolarization Ratio

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

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
        datagroup, datatype, dataset, product = get_datatypefields(
            datatypedescr)
        if datatype == 'RhoHV':
            rhohv = 'cross_correlation_ratio'
        if datatype == 'uRhoHV':
            rhohv = 'uncorrected_cross_correlation_ratio'
        if datatype == 'RhoHVu':
            rhohv = 'unfiltered_cross_correlation_ratio'
        if datatype == 'ZDR':
            zdr = 'differential_reflectivity'
        if datatype == 'ZDRc':
            zdr = 'corrected_differential_reflectivity'

    cdr = pyart.retrieve.compute_cdr(
        radar, rhohv_field=rhohv, zdr_field=zdr,
        cdr_field='circular_depolarization_ratio')

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field('circular_depolarization_ratio', cdr)

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
        data set configuration

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
        datagroup, datatype, dataset, product = get_datatypefields(
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

    rhohv = pyart.correct.correct_noise_rhohv(
        radar, urhohv_field=urhohv, snr_field=snr, zdr_field=zdr,
        nh_field=nh, nv_field=nv, rhohv_field='cross_correlation_ratio')

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field('cross_correlation_ratio', rhohv)

    return new_dataset


def process_echo_id(procstatus, dscfg, radar=None):
    """
    identifies echoes as 0: No data, 1: Noise, 2: Clutter,
    3: Precipitation

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

    radar : Radar
        Optional. Radar object

    Returns
    -------
    new_dataset : Radar
        radar object

    """

    if procstatus != 1:
        return None

    id = np.zeros((radar.nrays, radar.ngates), dtype='int32')+3

    # look for clutter
    gatefilter = pyart.filters.moment_and_texture_based_gate_filter(
        radar, zdr_field=None, rhv_field=None, phi_field=None,
        refl_field=None, textzdr_field=None, textrhv_field=None,
        textphi_field=None, textrefl_field=None, wind_size=7,
        max_textphi=20., max_textrhv=0.3, max_textzdr=2.85,
        max_textrefl=8., min_rhv=0.6)

    is_clutter = gatefilter.gate_excluded == 1
    id[is_clutter] = 2

    # look for noise
    is_noise = radar.fields['reflectivity']['data'].data == (
        pyart.config.get_fillvalue())
    id[is_noise] = 1

    id_field = pyart.config.get_metadata('radar_echo_id')
    id_field['data'] = id

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field('radar_echo_id', id_field)

    return new_dataset


def process_echo_filter(procstatus, dscfg, radar=None):
    """
    filters out undesired echo types.
    TODO: make the selection of the echo types to filter and to keep user
    defined

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

    radar : Radar
        Optional. Radar object

    Returns
    -------
    new_dataset : Radar
        radar object

    """

    if procstatus != 1:
        return None

    is_not_precip = radar.fields['radar_echo_id']['data'] != 3

    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatypefields(
            datatypedescr)
        field_name = get_fieldname_rainbow(datatype)

        radar_field = deepcopy(radar.fields[field_name])
        radar_field['data'] = np.ma.masked_where(
            is_not_precip, radar_field['data'])

        if field_name.startswith('corrected_'):
            new_field_name = field_name
        elif field_name.startswith('uncorrected_'):
            new_field_name = field_name.replace(
                'uncorrected_', 'corrected_', 1)
        else:
            new_field_name = 'corrected_'+field_name
        new_dataset.add_field(new_field_name, radar_field)

    return new_dataset


def process_filter_snr(procstatus, dscfg, radar=None):
    """
    filters out low SNR echoes

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

    radar : Radar
        Optional. Radar object

    Returns
    -------
    new_dataset : Radar
        radar object

    """

    if procstatus != 1:
        return None

    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatypefields(
            datatypedescr)
        if (datatype == 'SNRh') or (datatype == 'SNRv'):
            snr_field = get_fieldname_rainbow(datatype)
            break

    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatypefields(
            datatypedescr)

        if (datatype != 'SNRh') and (datatype != 'SNRv'):
            field_name = get_fieldname_rainbow(datatype)

            gatefilter = pyart.filters.snr_based_gate_filter(
                radar, snr_field=snr_field, min_snr=dscfg['SNRmin'])
            is_lowSNR = gatefilter.gate_excluded == 1

            radar_field = deepcopy(radar.fields[field_name])
            radar_field['data'] = np.ma.masked_where(
                is_lowSNR, radar_field['data'])

            if field_name.startswith('corrected_'):
                new_field_name = field_name
            elif field_name.startswith('uncorrected_'):
                new_field_name = field_name.replace(
                    'uncorrected_', 'corrected_', 1)
            else:
                new_field_name = 'corrected_'+field_name
            new_dataset.add_field(new_field_name, radar_field)

    return new_dataset


def process_correct_phidp0(procstatus, dscfg, radar=None):
    """
    corrects phidp of the system phase

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

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
        datagroup, datatype, dataset, product = get_datatypefields(
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

    return new_dataset


def process_smooth_phidp_single_window(procstatus, dscfg, radar=None):
    """
    corrects phidp of the system phase and smoothes it using one window

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

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
        datagroup, datatype, dataset, product = get_datatypefields(
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

    return new_dataset


def process_smooth_phidp_double_window(procstatus, dscfg, radar=None):
    """
    corrects phidp of the system phase and smoothes it using one window

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

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
        datagroup, datatype, dataset, product = get_datatypefields(
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

    return new_dataset


def process_phidp_kdp_Maesaka(procstatus, dscfg, radar=None):
    """
    Estimates PhiDP and KDP using the method by Maesaka

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

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
        datagroup, datatype, dataset, product = get_datatypefields(
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

    # filter out data in an above the melting layer
    gatefilter = pyart.filters.temp_based_gate_filter(
        radar_aux, temp_field=temp_field, min_temp=0., thickness=400.)
    is_notrain = gatefilter.gate_excluded == 1

    phidp['data'][is_notrain] = np.ma.masked
    mask = np.ma.getmaskarray(phidp['data'])
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

    return new_dataset


def process_phidp_kdp_lp(procstatus, dscfg, radar=None):
    """
    Estimates PhiDP and KDP using a linear programming algorithm

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

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
        datagroup, datatype, dataset, product = get_datatypefields(
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

    return new_dataset


def process_kdp_leastsquare_single_window(procstatus, dscfg, radar=None):
    """
    Computes specific differential phase using a piecewise least square method

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

    radar : Radar
        Optional. Radar object

    Returns
    -------
    radar : Radar
        radar object

    """
    if procstatus != 1:
        return None

    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatypefields(
            datatypedescr)
        if datatype == 'PhiDP':
            phidp_field = 'differential_phase'
        if datatype == 'PhiDPc':
            phidp_field = 'corrected_differential_phase'
        if datatype == 'uPhiDP':
            phidp_field = 'uncorrected_differential_phase'

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

    return new_dataset


def process_kdp_leastsquare_double_window(procstatus, dscfg, radar=None):
    """
    Computes specific differential phase using a piecewise least square method

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

    radar : Radar
        Optional. Radar object

    Returns
    -------
    radar : Radar
        radar object

    """
    if procstatus != 1:
        return None

    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatypefields(
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

    return new_dataset


def process_attenuation(procstatus, dscfg, radar=None):
    """
    Computes specific attenuation and specific differential attenuation using
    the Z-Phi method and corrects reflectivity and differential reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

    radar : Radar
        Optional. Radar object

    Returns
    -------
    radar : Radar
        radar object

    """

    if procstatus != 1:
        return None

    for datatypedescr in dscfg['datatype']:
        datagroup, datatype, dataset, product = get_datatypefields(
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

    spec_at, cor_z, spec_diff_at, cor_zdr = (
        pyart.correct.calculate_attenuation(
            radar, doc=15, fzl=None, smooth_window_len=0, a_coef=None,
            beta=None, c=None, d=None, refl_field=refl, phidp_field=phidp,
            zdr_field=zdr, temp_field=temp, spec_at_field=None,
            corr_refl_field=None, spec_diff_at_field=None,
            corr_zdr_field=None))

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

    return new_dataset


def process_rainrate(procstatus, dscfg, radar=None):
    """
    Estimates rainfall rate from polarimetric moments

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

    radar : Radar
        Optional. Radar object

    Returns
    -------
    radar : Radar
        radar object

    """
    if procstatus != 1:
        return None

    if dscfg['RR_METHOD'] == 'Z':
        datagroup, datatype, dataset, product = get_datatypefields(
            dscfg['datatype'][0])
        refl_field = get_fieldname_rainbow(datatype)
        rain = pyart.retrieve.est_rain_rate_z(
            radar, alpha=0.0376, beta=0.6112, refl_field=refl_field,
            rr_field=None)

    elif dscfg['RR_METHOD'] == 'ZPoly':
        datagroup, datatype, dataset, product = get_datatypefields(
            dscfg['datatype'][0])
        refl_field = get_fieldname_rainbow(datatype)
        rain = pyart.retrieve.est_rain_rate_zpoly(
            radar, refl_field=refl_field, rr_field=None)

    elif dscfg['RR_METHOD'] == 'KDP':
        datagroup, datatype, dataset, product = get_datatypefields(
            dscfg['datatype'][0])
        kdp_field = get_fieldname_rainbow(datatype)
        rain = pyart.retrieve.est_rain_rate_kdp(
            radar, alpha=None, beta=None, kdp_field=kdp_field, rr_field=None)

    elif dscfg['RR_METHOD'] == 'A':
        datagroup, datatype, dataset, product = get_datatypefields(
            dscfg['datatype'][0])
        a_field = get_fieldname_rainbow(datatype)
        rain = pyart.retrieve.est_rain_rate_a(
            radar, alpha=None, beta=None, a_field=a_field, rr_field=None)

    elif dscfg['RR_METHOD'] == 'ZKDP':
        for datatypedescr in dscfg['datatype']:
            datagroup, datatype, dataset, product = get_datatypefields(
                datatypedescr)
            if datatype == 'dBZc':
                refl_field = 'corrected_reflectivity'
            if datatype == 'KDPc':
                kdp_field = 'corrected_specific_differential_phase'
            if datatype == 'dBZ':
                refl_field = 'reflectivity'
            if datatype == 'KDP':
                kdp_field = 'specific_differential_phase'

        rain = pyart.retrieve.est_rain_rate_zkdp(
            radar, alphaz=0.0376, betaz=0.6112, alphakdp=None, betakdp=None,
            refl_field=refl_field, kdp_field=kdp_field, rr_field=None,
            master_field=refl_field, thresh=10., thresh_max=True)

    elif dscfg['RR_METHOD'] == 'ZA':
        for datatypedescr in dscfg['datatype']:
            datagroup, datatype, dataset, product = get_datatypefields(
                datatypedescr)
            if datatype == 'dBZc':
                refl_field = 'corrected_reflectivity'
            if datatype == 'Ahc':
                a_field = 'corrected_specific_attenuation'
            if datatype == 'dBZ':
                refl_field = 'reflectivity'
            if datatype == 'Ah':
                a_field = 'specific_attenuation'

        rain = pyart.retrieve.est_rain_rate_za(
            radar, alphaz=0.0376, betaz=0.6112, alphaa=None, betaa=None,
            refl_field=refl_field, a_field=a_field, rr_field=None,
            master_field=refl_field, thresh=0.04, thresh_max=False)

    elif dscfg['RR_METHOD'] == 'hydro':
        for datatypedescr in dscfg['datatype']:
            datagroup, datatype, dataset, product = get_datatypefields(
                datatypedescr)
            if datatype == 'dBZc':
                refl_field = 'corrected_reflectivity'
            if datatype == 'Ahc':
                a_field = 'corrected_specific_attenuation'
            if datatype == 'hydro':
                hydro_field = 'radar_echo_classification'
            if datatype == 'dBZ':
                refl_field = 'reflectivity'
            if datatype == 'Ah':
                a_field = 'specific_attenuation'

        rain = pyart.retrieve.est_rain_rate_hydro(
            radar, alphazr=0.0376, betazr=0.6112, alphazs=0.1, betazs=0.5,
            alphaa=None, betaa=None, mp_factor=0.6, refl_field=refl_field,
            a_field=a_field, hydro_field=hydro_field, rr_field=None,
            master_field=refl_field, thresh=0.04, thresh_max=False)

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field('radar_estimated_rain_rate', rain)

    return new_dataset


def process_hydroclass(procstatus, dscfg, radar=None):
    """
    Classifies precipitation echoes

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

    radar : Radar
        Optional. Radar object

    Returns
    -------
    radar : Radar
        radar object

    """
    if procstatus != 1:
        return None

    if dscfg['HYDRO_METHOD'] == 'SEMISUPERVISED':
        for datatypedescr in dscfg['datatype']:
            datagroup, datatype, dataset, product = get_datatypefields(
                datatypedescr)
            if datatype == 'dBZc':
                refl_field = 'corrected_reflectivity'
            if datatype == 'ZDRc':
                zdr_field = 'corrected_differential_reflectivity'
            if datatype == 'RhoHVc':
                rhv_field = 'corrected_cross_correlation_ratio'
            if datatype == 'KDPc':
                kdp_field = 'corrected_specific_differential_phase'
            if datatype == 'TEMP':
                temp_field = 'temperature'
            if datatype == 'dBZ':
                refl_field = 'reflectivity'
            if datatype == 'ZDR':
                zdr_field = 'differential_reflectivity'
            if datatype == 'RhoHV':
                rhv_field = 'cross_correlation_ratio'
            if datatype == 'KDP':
                kdp_field = 'specific_differential_phase'

        mass_centers = np.zeros((9, 5))
        if dscfg['RADARCENTROIDS'] == 'A':
            #      Zh      ZDR     kdp   RhoHV   delta_Z
            mass_centers[0, :] = [
                13.5829,  0.4063, 0.0497, 0.9868,  1330.3]  # DS
            mass_centers[1, :] = [
                02.8453,  0.2457, 0.0000, 0.9798,  0653.8]  # CR
            mass_centers[2, :] = [
                07.6597,  0.2180, 0.0019, 0.9799, -1426.5]  # LR
            mass_centers[3, :] = [
                31.6815,  0.3926, 0.0828, 0.9978,  0535.3]  # GR
            mass_centers[4, :] = [
                39.4703,  1.0734, 0.4919, 0.9876, -1036.3]  # RN
            mass_centers[5, :] = [
                04.8267, -0.5690, 0.0000, 0.9691,  0869.8]  # VI
            mass_centers[6, :] = [
                30.8613,  0.9819, 0.1998, 0.9845, -0066.1]  # WS
            mass_centers[7, :] = [
                52.3969,  2.1094, 2.4675, 0.9730, -1550.2]  # MH
            mass_centers[8, :] = [
                50.6186, -0.0649, 0.0946, 0.9904,  1179.9]  # IH/HDG
        elif dscfg['RADARCENTROIDS'] == 'L':
            #       Zh      ZDR     kdp   RhoHV   delta_Z
            mass_centers[0, :] = [
                13.8231,  0.2514, 0.0644, 0.9861,  1380.6]  # DS
            mass_centers[1, :] = [
                03.0239,  0.1971, 0.0000, 0.9661,  1464.1]  # CR
            mass_centers[2, :] = [
                04.9447,  0.1142, 0.0000, 0.9787, -0974.7]  # LR
            mass_centers[3, :] = [
                34.2450,  0.5540, 0.1459, 0.9937,  0945.3]  # GR
            mass_centers[4, :] = [
                40.9432,  1.0110, 0.5141, 0.9928, -0993.5]  # RN
            mass_centers[5, :] = [
                03.5202, -0.3498, 0.0000, 0.9746,  0843.2]  # VI
            mass_centers[6, :] = [
                32.5287,  0.9751, 0.2640, 0.9804, -0055.5]  # WS
            mass_centers[7, :] = [
                52.6547,  2.7054, 2.5101, 0.9765, -1114.6]  # MH
            mass_centers[8, :] = [
                46.4998,  0.1978, 0.6431, 0.9845,  1010.1]  # IH/HDG
        elif dscfg['RADARCENTROIDS'] == 'P':
            #       Zh      ZDR     kdp   RhoHV   delta_Z
            mass_centers[0, :] = [
                13.9882,  0.2470, 0.0690, 0.9939,  1418.1]  # DS
            mass_centers[1, :] = [
                00.9834,  0.4830, 0.0043, 0.9834,  0950.6]  # CR
            mass_centers[2, :] = [
                05.3962,  0.2689, 0.0000, 0.9831, -0479.5]  # LR
            mass_centers[3, :] = [
                35.3411,  0.1502, 0.0940, 0.9974,  0920.9]  # GR
            mass_centers[4, :] = [
                35.0114,  0.9681, 0.1106, 0.9785, -0374.0]  # RN
            mass_centers[5, :] = [
                02.5897, -0.3879, 0.0282, 0.9876,  0985.5]  # VI
            mass_centers[6, :] = [
                32.2914,  0.7789, 0.1443, 0.9075, -0153.5]  # WS
            mass_centers[7, :] = [
                53.2413,  1.8723, 0.3857, 0.9454, -0470.8]  # MH
            mass_centers[8, :] = [
                44.7896,  0.0015, 0.1349, 0.9968,  1116.7]  # IH/HDG
        elif dscfg['RADARCENTROIDS'] == 'DX50':
            #       Zh      ZDR     kdp   RhoHV   delta_Z
            mass_centers[0, :] = [
                19.0770,  0.4139, 0.0099, 0.9841,  1061.7]  # DS
            mass_centers[1, :] = [
                03.9877,  0.5040, 0.0000, 0.9642,  0856.6]  # CR
            mass_centers[2, :] = [
                20.7982,  0.3177, 0.0004, 0.9858, -1375.1]  # LR
            mass_centers[3, :] = [
                34.7124, -0.3748, 0.0988, 0.9828,  1224.2]  # GR
            mass_centers[4, :] = [
                33.0134,  0.6614, 0.0819, 0.9802, -1169.8]  # RN
            mass_centers[5, :] = [
                08.2610, -0.4681, 0.0000, 0.9722,  1100.7]  # VI
            mass_centers[6, :] = [
                35.1801,  1.2830, 0.1322, 0.9162, -0159.8]  # WS
            mass_centers[7, :] = [
                52.4539,  2.3714, 1.1120, 0.9382, -1618.5]  # MH
            mass_centers[8, :] = [
                44.2216, -0.3419, 0.0687, 0.9683,  1272.7]  # IH/HDG
        else:
            warn(
                ' Unknown radar. ' +
                'Default centroids will be used in classification.')
            mass_centers = None

        hydro = pyart.retrieve.hydroclass_semisupervised(
            radar, mass_centers=mass_centers,
            weights=np.array([1., 1., 1., 0.75, 0.5]), refl_field=refl_field,
            zdr_field=zdr_field, rhv_field=rhv_field, kdp_field=kdp_field,
            temp_field=temp_field, hydro_field=None)

        # prepare for exit
        new_dataset = deepcopy(radar)
        new_dataset.fields = dict()
        new_dataset.add_field('radar_echo_classification', hydro)

        return new_dataset


def process_point_measurement(procstatus, dscfg, radar=None):
    """
    Obtains the radar data at a point measurement

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing

    dscfg : dictionary of dictionaries
        data set configuration

    radar : Radar
        Optional. Radar object

    Returns
    -------
    radar : Radar
        radar object

    """
    if procstatus != 1:
        return None

    projparams = dict()
    projparams.update({'proj': 'pyart_aeqd'})
    projparams.update({'lon_0': radar.longitude['data']})
    projparams.update({'lat_0': radar.latitude['data']})

    if dscfg['latlon']:
        lon = dscfg['lon']
        lat = dscfg['lat']
        alt = dscfg['alt']
        x, y = pyart.core.geographic_to_cartesian(lon, lat, projparams)

        if not dscfg['truealt']:
            ke = 4./3.  # constant for effective radius
            a = 6378100.  # earth radius
            re = a * ke  # effective radius

            elrad = dscfg['ele'] * np.pi / 180.
            r_ground = np.sqrt(x ** 2. + y ** 2.)
            r = r_ground / np.cos(elrad)
            alt_radar = radar.altitude['data']+np.sqrt(
                r ** 2. + re ** 2. + 2. * r * re * np.sin(elrad)) - re
            alt_radar = alt_radar[0]
        else:
            alt_radar = dscfg['alt']

        r, az, el = pyart.core.cartesian_to_antenna(
            x, y, alt_radar-radar.altitude['data'])
        r = r[0]
        az = az[0]
        el = el[0]
    else:
        r = dscfg['rng']
        az = dscfg['azi']
        el = dscfg['ele']

        x, y, alt = antenna_to_cartesian(r, az, el)
        lon, lat = cartesian_to_geographic(x, y, projparams)

    d_az = np.min(np.abs(radar.azimuth['data'] - az))
    if d_az > dscfg['AziTol']:
        warn(' No radar bin found for point (az, el, r):(' +
             str(az)+', '+str(el)+', '+str(r) +
             '). Minimum distance to radar azimuth '+str(d_az) +
             ' larger than tolerance')
        return None

    d_el = np.min(np.abs(radar.elevation['data'] - el))
    if d_el > dscfg['EleTol']:
        warn(' No radar bin found for point (az, el, r):(' +
             str(az)+', '+str(el)+', '+str(r) +
             '). Minimum distance to radar elevation '+str(d_el) +
             ' larger than tolerance')
        return None

    d_r = np.min(np.abs(radar.range['data'] - r))
    if d_r > dscfg['RngTol']:
        warn(' No radar bin found for point (az, el, r):(' +
             str(az)+', '+str(el)+', '+str(r) +
             '). Minimum distance to radar range bin '+str(d_r) +
             ' larger than tolerance')
        return None

    datagroup, datatype, dataset, product = get_datatypefields(
            dscfg['datatype'][0])
    field_name = get_fieldname_rainbow(datatype)

    ind_ray = np.argmin(np.abs(radar.azimuth['data'] - az) +
                        np.abs(radar.elevation['data'] - el))
    ind_r = np.argmin(np.abs(radar.range['data'] - r))

    val = radar.fields[field_name]['data'].data[ind_ray, ind_r]
    time = num2date(radar.time['data'][ind_ray], radar.time['units'],
                    radar.time['calendar'])

    # prepare for exit
    new_dataset = dict()
    new_dataset.update({'value': val})
    new_dataset.update({'datatype': datatype})
    new_dataset.update({'time': time})
    new_dataset.update(
        {'point_coordinates_WGS84_lon_lat_alt': [lon, lat, alt]})
    new_dataset.update({'antenna_coordinates_az_el_r': [az, el, r]})
    new_dataset.update(
        {'used_antenna_coordinates_az_el_r': [radar.azimuth['data'][ind_ray],
         radar.elevation['data'][ind_ray],
         radar.range['data'][ind_r]]})

    return new_dataset
