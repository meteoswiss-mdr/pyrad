"""
pyrad.proc.process_monitoring
=============================

Functions for monitoring of the polarimetric variables

.. autosummary::
    :toctree: generated/

    process_selfconsistency_kdp_phidp
    process_selfconsistency_bias
    process_estimate_phidp0
    process_rhohv_rain
    process_zdr_precip
    process_zdr_snow
    process_monitoring

"""

from copy import deepcopy
from warnings import warn
import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.read_data_other import read_selfconsistency
from ..io.read_data_radar import interpol_field

from ..util.radar_utils import get_histogram_bins


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
        frequency : float. Dataset keyword
            the radar frequency [Hz]. If None that of the key
            frequency in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            the selfconsistency will not be computed
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    temp = None
    iso0 = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
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
        freq = dscfg.get('frequency', None)
        if freq is None:
            if (radar.instrument_parameters is not None and
                    'frequency' in radar.instrument_parameters):
                freq = radar.instrument_parameters['frequency']['data'][0]
            else:
                warn('Unable to retrieve PhiDP and KDP using ' +
                     'self-consistency. Unknown radar frequency')
                return None, None

        # get frequency band
        freq_band = pyart.retrieve.get_freq_band(
            radar.instrument_parameters['frequency']['data'][0])

        # find unique elevations
        el_vec = np.unique(
            (10.*np.round(radar.elevation['data'], decimals=1)).astype(int))
        zdr_kdpzh_list = list()
        el_list = list()
        for el in el_vec:
            fname = (
                dscfg['configpath'] + 'selfconsistency/' +
                'selfconsistency_zdr_zhkdp_'+freq_band+'band_temp10_elev' +
                '{:03d}'.format(el)+'_mu05.txt')
            zdr_kdpzh_table = read_selfconsistency(fname)
            if zdr_kdpzh_table is not None:
                zdr_kdpzh_list.append(zdr_kdpzh_table)
                el_list.append((el/10.).astype(int))
        if not el_list:
            warn('Unable to retrieve PhiDP and KDP using self-consistency. ' +
                 'No selfconsistency files for the radar elevations.')

            return None, None

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
        new_dataset = {'radar_out': deepcopy(radar)}
        new_dataset['radar_out'].fields = dict()

        new_dataset['radar_out'].add_field(kdpsim_field, kdpsim)
        new_dataset['radar_out'].add_field(phidpsim_field, phidpsim)

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
        frequency : float. Dataset keyword
            the radar frequency [Hz]. If None that of the key
            frequency in attribute instrument_parameters of the radar
            object will be used. If the key or the attribute are not present
            the selfconsistency will not be computed
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    temp = None
    iso0 = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
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
        freq = dscfg.get('frequency', None)
        if freq is None:
            if (radar.instrument_parameters is not None and
                    'frequency' in radar.instrument_parameters):
                freq = radar.instrument_parameters['frequency']['data'][0]
            else:
                warn('Unable to retrieve PhiDP and KDP using ' +
                     'self-consistency. Unknown radar frequency')
                return None, None

        # get frequency band
        freq_band = pyart.retrieve.get_freq_band(
            radar.instrument_parameters['frequency']['data'][0])

        # find unique elevations
        el_vec = np.unique(
            (10.*np.round(radar.elevation['data'], decimals=1)).astype(int))
        zdr_kdpzh_list = list()
        el_list = list()
        for el in el_vec:
            fname = (
                dscfg['configpath'] + 'selfconsistency/' +
                'selfconsistency_zdr_zhkdp_'+freq_band+'band_temp10_elev' +
                '{:03d}'.format(el)+'_mu05.txt')
            zdr_kdpzh_table = read_selfconsistency(fname)
            if zdr_kdpzh_table is not None:
                zdr_kdpzh_list.append(zdr_kdpzh_table)
                el_list.append((el/10.).astype(int))
        if not el_list:
            warn('Unable to retrieve PhiDP and KDP using self-consistency. ' +
                 'No selfconsistency files for the radar elevations.')

            return None, None

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
        new_dataset = {'radar_out': deepcopy(radar)}
        new_dataset['radar_out'].fields = dict()

        new_dataset['radar_out'].add_field('reflectivity_bias', refl_bias)

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
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
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
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()

    new_dataset['radar_out'].add_field('system_differential_phase', phidp0)
    new_dataset['radar_out'].add_field(
        'first_gate_differential_phase', first_gates)

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
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    temp_field = None
    iso0_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
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
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(
        'cross_correlation_ratio_in_rain', rhohv_rain)

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
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    temp_field = None
    iso0_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
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
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()

    new_dataset['radar_out'].add_field(
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
            Default [2] (dry snow)
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index
    """

    if procstatus != 1:
        return None, None

    temp_field = None
    kdp_field = None
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
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
            (hydro_field not in radar.fields)):
        warn('Unable to estimate ZDR in snow. Missing data')
        return None, None

    # User defined values
    rmin = dscfg.get('rmin', 1000.)
    rmax = dscfg.get('rmax', 50000.)
    zmin = dscfg.get('Zmin', 0.)
    zmax = dscfg.get('Zmax', 30.)
    snrmin = dscfg.get('SNRmin', 10.)
    snrmax = dscfg.get('SNRmax', 50.)
    rhohvmin = dscfg.get('RhoHVmin', 0.97)
    phidpmax = dscfg.get('PhiDPmax', 10.)
    elmax = dscfg.get('elmax', None)
    kdpmax = dscfg.get('KDPmax', None)
    tempmin = dscfg.get('TEMPmin', None)
    tempmax = dscfg.get('TEMPmax', None)
    hydroclass = dscfg.get('hydroclass', [2])

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
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()

    new_dataset['radar_out'].add_field(
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
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
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

        step = dscfg.get('step', None)

        bin_edges = get_histogram_bins(field_name, step=step)
        nbins = len(bin_edges)-1
        step = bin_edges[1]-bin_edges[0]
        bin_centers = bin_edges[:-1]+step/2.

        radar_aux = deepcopy(radar)
        radar_aux.fields = dict()
        radar_aux.range['data'] = bin_centers
        radar_aux.ngates = nbins

        field_dict = pyart.config.get_metadata(field_name)
        field_dict['data'] = np.ma.zeros((radar.nrays, nbins), dtype=int)

        field = deepcopy(radar.fields[field_name]['data'])

        # put gates with values off limits to limit
        mask = np.ma.getmaskarray(field)
        ind = np.where(np.logical_and(mask == False, field < bin_centers[0]))
        field[ind] = bin_centers[0]

        ind = np.where(np.logical_and(mask == False, field > bin_centers[-1]))
        field[ind] = bin_centers[-1]

        for ray in range(radar.nrays):
            field_dict['data'][ray, :], bin_edges = np.histogram(
                field[ray, :].compressed(), bins=bin_edges)

        radar_aux.add_field(field_name, field_dict)
        start_time = pyart.graph.common.generate_radar_time_begin(radar_aux)

        # keep histogram in Memory or add to existing histogram
        if dscfg['initialized'] == 0:
            dscfg['global_data'] = {'hist_obj': radar_aux,
                                    'timeinfo': start_time}
            dscfg['initialized'] = 1
        else:
            field_interp = interpol_field(
                dscfg['global_data']['hist_obj'], radar_aux, field_name,
                fill_value=0)
            dscfg['global_data']['hist_obj'].fields[field_name]['data'] += (
                field_interp['data'].filled(fill_value=0)).astype('int64')

        #    dscfg['global_data']['timeinfo'] = dscfg['timeinfo']

        dataset = dict()
        dataset.update({'hist_obj': radar_aux})
        dataset.update({'hist_type': 'instant'})
        dataset.update({'timeinfo': start_time})

        return dataset, ind_rad

    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None

        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_name = get_fieldname_pyart(datatype)
            break
        ind_rad = int(radarnr[5:8])-1

        dataset = dict()
        dataset.update({'hist_obj': dscfg['global_data']['hist_obj']})
        dataset.update({'hist_type': 'cumulative'})
        dataset.update({'timeinfo': dscfg['global_data']['timeinfo']})

        return dataset, ind_rad
