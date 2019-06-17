"""
pyrad.proc.process_spectra
==========================

Functions to processes spectral data.

.. autosummary::
    :toctree: generated/

    process_raw_spectra
    process_filter_0Doppler
    process_filter_srhohv
    process_filter_spectra_noise
    process_spectral_power
    process_spectral_phase
    process_spectral_reflectivity
    process_spectral_differential_reflectivity
    process_spectral_differential_phase
    process_spectral_rhohv
    process_pol_variables
    process_reflectivity
    process_differential_reflectivity
    process_differential_phase
    process_rhohv
    process_Doppler_velocity
    process_Doppler_width

"""

from copy import deepcopy
from warnings import warn
import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart


def process_raw_spectra(procstatus, dscfg, radar_list=None):
    """
    Dummy function that returns the initial input data set

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break
    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    new_dataset = {'radar_out': deepcopy(radar_list[ind_rad])}

    return new_dataset, ind_rad


def process_filter_0Doppler(procstatus, dscfg, radar_list=None):
    """
    Function to filter the 0-Doppler line bin and neighbours of the
    Doppler spectra

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
        filter_width : float
            The Doppler filter width. Default 0.
        filter_units : str
            Can be 'm/s' or 'Hz'. Default 'm/s'
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    field_name_list = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name_list.append(get_fieldname_pyart(datatype))

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    filter_width = dscfg.get('filter_width', 0.)
    filter_units = dscfg.get('filter_units', 'm/s')

    if filter_units == 'm/s':
        axis = psr.Doppler_velocity['data']
    else:
        axis = psr.Doppler_frequency['data']

    fields = dict()
    field_name_list_aux = []
    for field_name in field_name_list:
        if field_name not in psr.fields:
            warn('Unable to filter 0-Doppler. Missing field '+field_name)
            continue

        field = deepcopy(psr.fields[field_name])
        for ray in range(psr.nrays):
            ind = np.ma.where(np.logical_and(
                axis[ray, :] >= -filter_width/2.,
                axis[ray, :] <= filter_width/2.))
            field['data'][ray, :, ind] = np.ma.masked
        fields.update({field_name: field})
        field_name_list_aux.append(field_name)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(psr)}
    new_dataset['radar_out'].fields = dict()
    for field_name in field_name_list_aux:
        new_dataset['radar_out'].add_field(field_name, fields[field_name])

    return new_dataset, ind_rad


def process_filter_srhohv(procstatus, dscfg, radar_list=None):
    """
    Filter Doppler spectra as a function of spectral RhoHV

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
        sRhoHV_threshold : float
            Data with sRhoHV module above this threshold will be filtered.
            Default 1.
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    field_name_list = []
    sRhoHV_found = False
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'sRhoHV' and not sRhoHV_found:
            sRhoHV_field = get_fieldname_pyart(datatype)
            sRhoHV_found = True
        else:
            field_name_list.append(get_fieldname_pyart(datatype))

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if sRhoHV_field not in psr.fields:
        warn('Unable to obtain apply sRhoHV filter. Missing field ' +
             sRhoHV_field)
        return None, None

    sRhoHV_threshold = dscfg.get('sRhoHV_threshold', 1.)
    sRhoHV = psr.fields[sRhoHV_field]['data']

    fields = dict()
    field_name_list_aux = []
    for field_name in field_name_list:
        if field_name not in psr.fields:
            warn('Unable to filter 0-Doppler. Missing field '+field_name)
            continue

        field = deepcopy(psr.fields[field_name])
        field['data'][np.ma.abs(sRhoHV) >= sRhoHV_threshold] = np.ma.masked
        fields.update({field_name: field})
        field_name_list_aux.append(field_name)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(psr)}
    new_dataset['radar_out'].fields = dict()
    for field_name in field_name_list_aux:
        new_dataset['radar_out'].add_field(field_name, fields[field_name])

    return new_dataset, ind_rad


def process_filter_spectra_noise(procstatus, dscfg, radar_list=None):
    """
    Filter the noise of the Doppler spectra by clipping any data below
    the noise level plus a margin

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
        clipping_level : float
            The clipping level [dB above noise level]. Default 10.
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ('ShhADU', 'SvvADU'):
            signal_field = get_fieldname_pyart(datatype)
        elif datatype in ('NADUh', 'NADUv'):
            noise_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if signal_field not in psr.fields or noise_field not in psr.fields:
        warn('Unable to obtain apply spectral noise filter. Missing fields')
        return None, None

    clipping_level = dscfg.get('clipping_level', 10.)

    clip_pwr = (
        psr.fields[noise_field]['data']*np.power(10., 0.1*clipping_level))

    s_pwr = pyart.retrieve.compute_spectral_power(
        psr, units='ADU', signal_field=signal_field,
        noise_field=noise_field)

    mask = np.ma.less_equal(s_pwr['data'], clip_pwr)

    field = deepcopy(psr.fields[signal_field])
    field['data'][mask] = np.ma.masked

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(psr)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(signal_field, field)

    return new_dataset, ind_rad


def process_spectral_power(procstatus, dscfg, radar_list=None):
    """
    Computes the spectral power

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
        units : str
            The units of the returned signal. Can be 'ADU', 'dBADU' or 'dBm'
        subtract_noise : Bool
            If True noise will be subtracted from the signal
        smooth_window : int or None
            Size of the moving Gaussian smoothing window. If none no smoothing
            will be applied
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    noise_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ('ShhADU', 'SvvADU'):
            signal_field = get_fieldname_pyart(datatype)
        elif datatype in ('NADUh', 'NADUv'):
            noise_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if signal_field not in psr.fields:
        warn('Unable to obtain spectral signal power. Missing field ' +
             signal_field)
        return None, None

    units = dscfg.get('units', 'dBADU')
    subtract_noise = dscfg.get('subtract_noise', False)
    smooth_window = dscfg.get('smooth_window', None)

    s_pwr = pyart.retrieve.compute_spectral_power(
        psr, units=units, subtract_noise=subtract_noise,
        smooth_window=smooth_window, signal_field=signal_field,
        noise_field=noise_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(psr)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(s_pwr['standard_name'], s_pwr)

    return new_dataset, ind_rad


def process_spectral_phase(procstatus, dscfg, radar_list=None):
    """
    Computes the spectral phase

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ('ShhADU', 'SvvADU'):
            signal_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if signal_field not in psr.fields:
        warn('Unable to obtain spectral phase. Missing field ' +
             signal_field)
        return None, None

    s_phase = pyart.retrieve.compute_spectral_phase(
        psr, signal_field=signal_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(psr)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(s_phase['standard_name'], s_phase)

    return new_dataset, ind_rad


def process_spectral_reflectivity(procstatus, dscfg, radar_list=None):
    """
    Computes spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
        subtract_noise : Bool
            If True noise will be subtracted from the signal
        smooth_window : int or None
            Size of the moving Gaussian smoothing window. If none no smoothing
            will be applied
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    noise_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ('ShhADU', 'SvvADU'):
            signal_field = get_fieldname_pyart(datatype)
        elif datatype in ('NADUh', 'NADUv'):
            noise_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if signal_field not in psr.fields:
        warn('Unable to obtain spectral reflectivity. Missing field ' +
             signal_field)
        return None, None

    subtract_noise = dscfg.get('subtract_noise', False)
    smooth_window = dscfg.get('smooth_window', None)

    sdBZ = pyart.retrieve.compute_spectral_reflectivity(
        psr, subtract_noise=subtract_noise, smooth_window=smooth_window,
        signal_field=signal_field, noise_field=noise_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(psr)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(sdBZ['standard_name'], sdBZ)

    return new_dataset, ind_rad


def process_spectral_differential_reflectivity(procstatus, dscfg, radar_list=None):
    """
    Computes spectral differential reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
        subtract_noise : Bool
            If True noise will be subtracted from the signal
        smooth_window : int or None
            Size of the moving Gaussian smoothing window. If none no smoothing
            will be applied
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    noise_h_field = None
    noise_v_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'ShhADU':
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype == 'SvvADU':
            signal_v_field = get_fieldname_pyart(datatype)
        elif datatype == 'NADUh':
            noise_h_field = get_fieldname_pyart(datatype)
        elif datatype == 'NADUv':
            noise_v_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if signal_h_field not in psr.fields or signal_v_field not in psr.fields:
        warn('Unable to obtain spectral differential reflectivity. ' +
             'Missing fields')
        return None, None

    subtract_noise = dscfg.get('subtract_noise', False)
    smooth_window = dscfg.get('smooth_window', None)

    sZDR = pyart.retrieve.compute_spectral_differential_reflectivity(
        psr, subtract_noise=subtract_noise, smooth_window=smooth_window,
        signal_h_field=signal_h_field, signal_v_field=signal_v_field,
        noise_h_field=noise_h_field, noise_v_field=noise_v_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(psr)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(sZDR['standard_name'], sZDR)

    return new_dataset, ind_rad


def process_spectral_differential_phase(procstatus, dscfg, radar_list=None):
    """
    Computes the spectral differential phase

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'ShhADU':
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype == 'SvvADU':
            signal_v_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if signal_h_field not in psr.fields or signal_v_field not in psr.fields:
        warn('Unable to obtain spectral signal differential phase. ' +
             'Missing fields')
        return None, None

    sPhiDP = pyart.retrieve.compute_spectral_differential_phase(
        psr, signal_h_field=signal_h_field, signal_v_field=signal_v_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(psr)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(sPhiDP['standard_name'], sPhiDP)

    return new_dataset, ind_rad


def process_spectral_rhohv(procstatus, dscfg, radar_list=None):
    """
    Computes the spectral RhoHV

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
        subtract_noise : Bool
            If True noise will be subtracted from the signal
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    noise_h_field = None
    noise_v_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'ShhADU':
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype == 'SvvADU':
            signal_v_field = get_fieldname_pyart(datatype)
        elif datatype == 'NADUh':
            noise_h_field = get_fieldname_pyart(datatype)
        elif datatype == 'NADUv':
            noise_v_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if signal_h_field not in psr.fields or signal_v_field not in psr.fields:
        warn('Unable to obtain spectral RhoHV. ' +
             'Missing fields')
        return None, None

    subtract_noise = dscfg.get('subtract_noise', False)

    sRhoHV = pyart.retrieve.compute_spectral_rhohv(
        psr, subtract_noise=subtract_noise, signal_h_field=signal_h_field,
        signal_v_field=signal_v_field, noise_h_field=noise_h_field,
        noise_v_field=noise_v_field)

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(psr)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(sRhoHV['standard_name'], sRhoHV)

    return new_dataset, ind_rad


def process_pol_variables(procstatus, dscfg, radar_list=None):
    """
    Computes the polarimetric variables from the complex spectra

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
        subtract_noise : Bool
            If True noise will be subtracted from the signal
        smooth_window : int or None
            Size of the moving Gaussian smoothing window. If none no smoothing
            will be applied
        variables : list of str
            list of variables to compute. Default dBZ
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    noise_h_field = None
    noise_v_field = None
    signal_h_field = None
    signal_v_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'ShhADU':
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype == 'SvvADU':
            signal_v_field = get_fieldname_pyart(datatype)
        elif datatype == 'NADUh':
            noise_h_field = get_fieldname_pyart(datatype)
        elif datatype == 'NADUv':
            noise_v_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if signal_h_field not in psr.fields and signal_v_field not in psr.fields:
        warn('Unable to obtain polarimetric variables. ' +
             'Missing fields')
        return None, None

    subtract_noise = dscfg.get('subtract_noise', False)
    smooth_window = dscfg.get('smooth_window', None)
    variables = dscfg.get('variables', ['dBZ'])

    fields_list = []
    for variable in variables:
        fields_list.append(get_fieldname_pyart(variable))

    radar = pyart.retrieve.compute_pol_variables(
        psr, fields_list, subtract_noise=subtract_noise,
        smooth_window=smooth_window, signal_h_field=signal_h_field,
        signal_v_field=signal_v_field, noise_h_field=noise_h_field,
        noise_v_field=noise_v_field)

    # prepare for exit
    new_dataset = {'radar_out': radar}

    return new_dataset, ind_rad


def process_reflectivity(procstatus, dscfg, radar_list=None):
    """
    Computes reflectivity from the spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ('sdBZ', 'sdBZv'):
            sdBZ_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if sdBZ_field not in psr.fields:
        warn('Unable to obtain reflectivity. ' +
             'Missing field '+sdBZ_field)
        return None, None

    dBZ = pyart.retrieve.compute_reflectivity(
        psr, sdBZ_field=sdBZ_field)

    reflectivity_field = 'reflectivity'
    if datatype == 'sdBZv':
        reflectivity_field += 'vv'

    # prepare for exit
    new_dataset = {'radar_out': pyart.util.radar_from_spectra(psr)}
    new_dataset['radar_out'].add_field(reflectivity_field, dBZ)

    return new_dataset, ind_rad


def process_differential_reflectivity(procstatus, dscfg, radar_list=None):
    """
    Computes differential reflectivity from the horizontal and vertical
    spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'sdBZ':
            sdBZ_field = get_fieldname_pyart(datatype)
        elif datatype == 'sdBZv':
            sdBZv_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if sdBZ_field not in psr.fields or sdBZv_field not in psr.fields:
        warn('Unable to obtain differential reflectivity. ' +
             'Missing fields.')
        return None, None

    ZDR = pyart.retrieve.compute_differential_reflectivity(
        psr, sdBZ_field=sdBZ_field, sdBZv_field=sdBZv_field)

    # prepare for exit
    new_dataset = {'radar_out': pyart.util.radar_from_spectra(psr)}
    new_dataset['radar_out'].add_field('differential_reflectivity', ZDR)

    return new_dataset, ind_rad


def process_differential_phase(procstatus, dscfg, radar_list=None):
    """
    Computes the differential phase from the spectral differential phase and
    the spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ('sdBZ', 'sdBZv'):
            sdBZ_field = get_fieldname_pyart(datatype)
        elif datatype == 'sPhiDP':
            sPhiDP_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if sdBZ_field not in psr.fields or sPhiDP_field not in psr.fields:
        warn('Unable to obtain PhiDP. Missing fields')
        return None, None

    PhiDP = pyart.retrieve.compute_differential_phase(
        psr, sdBZ_field=sdBZ_field, sPhiDP_field=sPhiDP_field)

    # prepare for exit
    new_dataset = {'radar_out': pyart.util.radar_from_spectra(psr)}
    new_dataset['radar_out'].add_field('differential_phase', PhiDP)

    return new_dataset, ind_rad


def process_rhohv(procstatus, dscfg, radar_list=None):
    """
    Computes RhoHV from the complex spectras

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
        subtract_noise : Bool
            If True noise will be subtracted from the signal
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    noise_h_field = None
    noise_v_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'ShhADU':
            signal_h_field = get_fieldname_pyart(datatype)
        elif datatype == 'SvvADU':
            signal_v_field = get_fieldname_pyart(datatype)
        elif datatype == 'NADUh':
            noise_h_field = get_fieldname_pyart(datatype)
        elif datatype == 'NADUv':
            noise_v_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if signal_h_field not in psr.fields or signal_v_field not in psr.fields:
        warn('Unable to obtain RhoHV. ' +
             'Missing fields')
        return None, None

    subtract_noise = dscfg.get('subtract_noise', False)

    RhoHV = pyart.retrieve.compute_rhohv(
        psr, subtract_noise=subtract_noise, signal_h_field=signal_h_field,
        signal_v_field=signal_v_field, noise_h_field=noise_h_field,
        noise_v_field=noise_v_field)

    # prepare for exit
    new_dataset = {'radar_out': pyart.util.radar_from_spectra(psr)}
    new_dataset['radar_out'].add_field('cross_correlation_ratio', RhoHV)

    return new_dataset, ind_rad


def process_Doppler_velocity(procstatus, dscfg, radar_list=None):
    """
    Compute the Doppler velocity from the spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ('sdBZ', 'sdBZv'):
            sdBZ_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if sdBZ_field not in psr.fields:
        warn('Unable to obtain Doppler velocity. ' +
             'Missing field '+sdBZ_field)
        return None, None

    vel = pyart.retrieve.compute_Doppler_velocity(
        psr, sdBZ_field=sdBZ_field)

    # prepare for exit
    new_dataset = {'radar_out': pyart.util.radar_from_spectra(psr)}
    new_dataset['radar_out'].add_field('velocity', vel)

    return new_dataset, ind_rad


def process_Doppler_width(procstatus, dscfg, radar_list=None):
    """
    Compute the Doppler spectrum width from the spectral reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted configuration keywords::

        datatype : list of string. Dataset keyword
            The input data types
    radar_list : list of spectra objects
        Optional. list of spectra objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype in ('sdBZ', 'sdBZv'):
            sdBZ_field = get_fieldname_pyart(datatype)

    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    psr = radar_list[ind_rad]

    if sdBZ_field not in psr.fields:
        warn('Unable to obtain Doppler spectrum width. ' +
             'Missing field '+sdBZ_field)
        return None, None

    width = pyart.retrieve.compute_Doppler_width(
        psr, sdBZ_field=sdBZ_field)

    # prepare for exit
    new_dataset = {'radar_out': pyart.util.radar_from_spectra(psr)}
    new_dataset['radar_out'].add_field('spectrum_width', width)

    return new_dataset, ind_rad
