"""
pyrad.proc.process_echoclass
===============================

Functions for echo classification and filtering

.. autosummary::
    :toctree: generated/

    process_echo_id
    process_echo_filter
    process_filter_snr
    process_filter_visibility
    process_outlier_filter
    process_hydroclass

"""

from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart


def process_echo_id(procstatus, dscfg, radar_list=None):
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
        if datatype == 'dBZ':
            refl_field = 'reflectivity'
        if datatype == 'dBuZ':
            refl_field = 'unfiltered_reflectivity'
        if datatype == 'ZDR':
            zdr_field = 'differential_reflectivity'
        if datatype == 'ZDRu':
            zdr_field = 'unfiltered_differential_reflectivity'
        if datatype == 'RhoHV':
            rhv_field = 'cross_correlation_ratio'
        if datatype == 'uPhiDP':
            phi_field = 'uncorrected_differential_phase'

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if ((refl_field not in radar.fields) or
            (zdr_field not in radar.fields) or
            (rhv_field not in radar.fields) or
            (phi_field not in radar.fields)):
        warn('Unable to create radar_echo_id dataset. Missing data')
        return None, None

    id = np.zeros((radar.nrays, radar.ngates), dtype='int32')+3

    # look for clutter
    gatefilter = pyart.filters.moment_and_texture_based_gate_filter(
        radar, zdr_field=zdr_field, rhv_field=rhv_field, phi_field=phi_field,
        refl_field=refl_field, textzdr_field=None, textrhv_field=None,
        textphi_field=None, textrefl_field=None, wind_size=7,
        max_textphi=20., max_textrhv=0.3, max_textzdr=2.85,
        max_textrefl=8., min_rhv=0.6)

    is_clutter = gatefilter.gate_excluded == 1
    id[is_clutter] = 2

    # look for noise
    is_noise = radar.fields[refl_field]['data'].data == (
        pyart.config.get_fillvalue())
    id[is_noise] = 1

    id_field = pyart.config.get_metadata('radar_echo_id')
    id_field['data'] = id

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field('radar_echo_id', id_field)

    return new_dataset, ind_rad


def process_echo_filter(procstatus, dscfg, radar_list=None):
    """
    Masks all echo types that are not of the class specified in
    keyword echo_type

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        echo_type : int
            The type of echo to keep: 1 noise, 2 clutter, 3 precipitation.
            Default 3
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
        if (datatype == 'echoID'):
            echoid_field = get_fieldname_pyart(datatype)
            break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if echoid_field not in radar.fields:
        warn('Unable to filter data. Missing echo ID field')
        return None, None

    echo_type = 3
    if 'echo_type' in dscfg:
        echo_type = dscfg['echo_type']

    mask = radar.fields[echoid_field]['data'] != echo_type

    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        if field_name not in echoid_field:
            if field_name in radar.fields:
                radar_field = deepcopy(radar.fields[field_name])
                radar_field['data'] = np.ma.masked_where(
                    mask, radar_field['data'])

                if field_name.startswith('corrected_'):
                    new_field_name = field_name
                elif field_name.startswith('uncorrected_'):
                    new_field_name = field_name.replace(
                        'uncorrected_', 'corrected_', 1)
                else:
                    new_field_name = 'corrected_'+field_name
                new_dataset.add_field(new_field_name, radar_field)
            else:
                warn('Unable to filter '+field_name+' according to echo ID. ' +
                     'No valid input fields')

    if not new_dataset.fields:
        return None, None

    return new_dataset, ind_rad


def process_filter_snr(procstatus, dscfg, radar_list=None):
    """
    filters out low SNR echoes

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        SNRmin : float. Dataset keyword
            The minimum SNR to keep the data.
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
        if (datatype == 'SNRh') or (datatype == 'SNRv'):
            snr_field = get_fieldname_pyart(datatype)
            break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    if snr_field not in radar.fields:
        warn('Unable to filter dataset according to SNR. Missing SNR field')
        return None, None

    gatefilter = pyart.filters.snr_based_gate_filter(
        radar, snr_field=snr_field, min_snr=dscfg['SNRmin'])
    is_lowSNR = gatefilter.gate_excluded == 1

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)

        if (datatype != 'SNRh') and (datatype != 'SNRv'):
            field_name = get_fieldname_pyart(datatype)
            if field_name in radar.fields:
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
            else:
                warn('Unable to filter '+field_name +
                     ' according to SNR. '+'No valid input fields')

    if not new_dataset.fields:
        return None, None

    return new_dataset, ind_rad


def process_filter_visibility(procstatus, dscfg, radar_list=None):
    """
    filters out rays gates with low visibility and corrects the reflectivity

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        VISmin : float. Dataset keyword
            The minimum visibility to keep the data.
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
        if datatype == 'VIS':
            vis_field = get_fieldname_pyart(datatype)
            break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()

    if vis_field not in radar.fields:
        warn('Unable to filter dataset according to visibility. ' +
             'Missing visibility field')
        return None, None

    gatefilter = pyart.filters.visibility_based_gate_filter(
        radar, vis_field=vis_field, min_vis=dscfg['VISmin'])
    is_lowVIS = gatefilter.gate_excluded == 1

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)

        if datatype != 'VIS':
            field_name = get_fieldname_pyart(datatype)
            if field_name in radar.fields:
                radar_aux = deepcopy(radar)
                radar_aux.fields[field_name]['data'] = np.ma.masked_where(
                    is_lowVIS, radar_aux.fields[field_name]['data'])

                if ((datatype == 'dBZ') or (datatype == 'dBZc') or
                        (datatype == 'dBuZ') or (datatype == 'dBZv') or
                        (datatype == 'dBZvc') or (datatype == 'dBuZv')):
                    radar_field = pyart.correct.correct_visibility(
                        radar_aux, vis_field=vis_field, field_name=field_name)
                else:
                    radar_field = radar_aux.fields[field_name]

                if field_name.startswith('corrected_'):
                    new_field_name = field_name
                elif field_name.startswith('uncorrected_'):
                    new_field_name = field_name.replace(
                        'uncorrected_', 'corrected_', 1)
                else:
                    new_field_name = 'corrected_'+field_name
                new_dataset.add_field(new_field_name, radar_field)
            else:
                warn('Unable to filter '+field_name +
                     ' according to visibility. No valid input fields')

    if not new_dataset.fields:
        return None, None

    return new_dataset, ind_rad


def process_outlier_filter(procstatus, dscfg, radar_list=None):
    """
    filters out gates which are outliers respect to the surrounding

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        threshold : float. Dataset keyword
            The distance between the value of the examined range gate and the
            median of the surrounding gates to consider the gate an outlier
        nb : int. Dataset keyword
            The number of neighbours (to one side) to analyse. i.e. 2 would
            correspond to 24 gates
        nb_min : int. Dataset keyword
            Minimum number of neighbouring gates to consider the examined gate
            valid
        percentile_min, percentile_max : float. Dataset keyword
            gates below (above) these percentiles (computed over the sweep) are
            considered potential outliers and further examined
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

    radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
        dscfg['datatype'][0])

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    field_name = get_fieldname_pyart(datatype)
    if field_name not in radar.fields:
        warn('Unable to perform outlier removal. No valid data')
        return None, None

    threshold = 10.
    if 'threshold' in dscfg:
        threshold = dscfg['threshold']
    nb = 2
    if 'nb' in dscfg:
        nb = dscfg['nb']
    nb_min = 3
    if 'nb_min' in dscfg:
        nb_min = dscfg['nb_min']
    percentile_min = 5.
    if 'percentile_min' in dscfg:
        percentile_min = dscfg['percentile_min']
    percentile_max = 95.
    if 'percentile_max' in dscfg:
        percentile_max = dscfg['percentile_max']

    field = radar.fields[field_name]
    field_out = deepcopy(field)
    for sweep in range(radar.nsweeps):
        # find gates suspected to be outliers
        sweep_start = radar.sweep_start_ray_index['data'][sweep]
        sweep_end = radar.sweep_end_ray_index['data'][sweep]
        nrays_sweep = radar.rays_per_sweep['data'][sweep]
        data_sweep = field['data'][sweep_start:sweep_end+1, :]
        percent_vals = np.nanpercentile(
            data_sweep.filled(fill_value=np.nan),
            (percentile_min, percentile_max))
        ind_ray, ind_rng = np.ma.where(
            np.ma.logical_or(
                data_sweep < percent_vals[0], data_sweep > percent_vals[1]))

        for gate in range(len(ind_ray)):
            # find neighbours of suspected outlier gate
            data_cube = []
            for ray_nb in range(-nb, nb+1):
                for rng_nb in range(-nb, nb+1):
                    if ray_nb == 0 and rng_nb == 0:
                        continue
                    if ((ind_ray[gate]+ray_nb >= 0) and
                            (ind_ray[gate]+ray_nb < nrays_sweep) and
                            (ind_rng[gate]+rng_nb >= 0) and
                            (ind_rng[gate]+rng_nb < radar.ngates)):
                        if (data_sweep[ind_ray[gate]+ray_nb,
                                       ind_rng[gate]+rng_nb] is not
                                np.ma.masked):
                            data_cube.append(
                                data_sweep[ind_ray[gate]+ray_nb,
                                           ind_rng[gate]+rng_nb])

            # remove data far from median of neighbours or with not enough
            # valid neighbours
            if len(data_cube) < nb_min:
                field_out['data'][
                    sweep_start+ind_ray[gate], ind_rng[gate]] = np.ma.masked
            elif abs(np.ma.median(data_cube)-data_sweep[ind_ray[gate],
                     ind_rng[gate]]) > threshold:
                field_out['data'][
                    sweep_start+ind_ray[gate], ind_rng[gate]] = np.ma.masked

    if field_name.startswith('corrected_'):
        new_field_name = field_name
    elif field_name.startswith('uncorrected_'):
        new_field_name = field_name.replace(
            'uncorrected_', 'corrected_', 1)
    else:
        new_field_name = 'corrected_'+field_name

    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field(new_field_name, field_out)

    return new_dataset, ind_rad


def process_hydroclass(procstatus, dscfg, radar_list=None):
    """
    Classifies precipitation echoes

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types
        HYDRO_METHOD : string. Dataset keyword
            The hydrometeor classification method. One of the following:
            SEMISUPERVISED
        RADARCENTROIDS : string. Datset keyword
            Used with HYDRO_METHOD SEMISUPERVISED. The name of the radar of
            which the derived centroids will be used. One of the following: A
            Albis, L Lema, P Plaine Morte, DX50
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

    if dscfg['HYDRO_METHOD'] == 'SEMISUPERVISED':
        for datatypedescr in dscfg['datatype']:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            if datatype == 'dBZ':
                refl_field = 'reflectivity'
            if datatype == 'dBZc':
                refl_field = 'corrected_reflectivity'
            if datatype == 'ZDR':
                zdr_field = 'differential_reflectivity'
            if datatype == 'ZDRc':
                zdr_field = 'corrected_differential_reflectivity'
            if datatype == 'RhoHV':
                rhv_field = 'cross_correlation_ratio'
            if datatype == 'RhoHVc':
                rhv_field = 'corrected_cross_correlation_ratio'
            if datatype == 'KDP':
                kdp_field = 'specific_differential_phase'
            if datatype == 'KDPc':
                kdp_field = 'corrected_specific_differential_phase'
            if datatype == 'TEMP':
                temp_field = 'temperature'

        ind_rad = int(radarnr[5:8])-1
        if radar_list[ind_rad] is None:
            warn('No valid radar')
            return None, None
        radar = radar_list[ind_rad]

        if ((refl_field not in radar.fields) or
                (zdr_field not in radar.fields) or
                (rhv_field not in radar.fields) or
                (kdp_field not in radar.fields) or
                (temp_field not in radar.fields)):
            warn('Unable to create hydrometeor classification field. ' +
                 'Missing data')
            return None, None

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

        return new_dataset, ind_rad
