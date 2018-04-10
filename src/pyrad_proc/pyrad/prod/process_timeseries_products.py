"""
pyrad.prod.process_timeseries_products
======================================

Functions for obtaining Pyrad products from a time series datasets

.. autosummary::
    :toctree: generated/

    generate_timeseries_products

"""

from copy import deepcopy
from warnings import warn

import numpy as np

from ..io.io_aux import get_save_dir, make_filename
from ..io.io_aux import generate_field_name_str

from ..io.read_data_sensor import get_sensor_data
from ..io.read_data_other import read_timeseries

from ..io.write_data import write_ts_polar_data, write_ts_cum

from ..graph.plots_timeseries import plot_timeseries, plot_timeseries_comp
from ..graph.plots import plot_scatter_comp

from ..util.radar_utils import rainfall_accumulation


def generate_timeseries_products(dataset, prdcfg):
    """
    Generates time series products

    Parameters
    ----------
    dataset : dictionary
        radar object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    no return

    """

    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'PLOT_AND_WRITE_POINT':
        if dataset['final']:
            return None

        dpi = 72
        if 'dpi' in prdcfg:
            dpi = prdcfg['dpi']

        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')[0]

        csvfname = savedir+csvfname

        write_ts_polar_data(dataset, csvfname)
        print('saved CSV file: '+csvfname)

        date, value = read_timeseries(csvfname)

        if date is None:
            warn(
                'Unable to plot time series. No valid data')
            return None

        figfname_list = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=date[0], timeformat='%Y%m%d')

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
        titl = ('Time Series '+date[0].strftime('%Y-%m-%d'))

        labely = generate_field_name_str(dataset['datatype'])

        plot_timeseries(
            date, [value], figfname_list, labelx='Time UTC',
            labely=labely, labels=[label1], title=titl, dpi=dpi)
        print('----- save to '+' '.join(figfname_list))

        return figfname_list

    elif prdcfg['type'] == 'PLOT_CUMULATIVE_POINT':
        if dataset['final']:
            return None

        dpi = 72
        if 'dpi' in prdcfg:
            dpi = prdcfg['dpi']

        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')[0]

        csvfname = savedir+csvfname

        date, value = read_timeseries(csvfname)

        if date is None:
            warn(
                'Unable to plot accumulation time series. No valid data')
            return None

        figfname_list = make_filename(
            'ts_cum', prdcfg['dstype'], dataset['datatype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=date[0], timeformat='%Y%m%d')

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
        titl = ('Time Series Acc. '+date[0].strftime('%Y-%m-%d'))

        labely = 'Radar estimated rainfall accumulation (mm)'

        plot_timeseries(
            date, [value], figfname_list, labelx='Time UTC',
            labely=labely, labels=[label1], title=titl,
            period=prdcfg['ScanPeriod']*60., dpi=dpi)
        print('----- save to '+' '.join(figfname_list))

        return figfname_list

    elif prdcfg['type'] == 'COMPARE_POINT':
        if dataset['final']:
            return None

        dpi = 72
        if 'dpi' in prdcfg:
            dpi = prdcfg['dpi']

        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir_ts = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')

        csvfname = savedir_ts+csvfname

        radardate, radarvalue = read_timeseries(csvfname)
        if radardate is None:
            warn(
                'Unable to plot sensor comparison at point of interest. ' +
                'No valid radar data')
            return None

        sensordate, sensorvalue, sensortype, period = get_sensor_data(
            radardate[0], dataset['datatype'], prdcfg)
        if sensordate is None:
            warn(
                'Unable to plot sensor comparison at point of interest. ' +
                'No valid sensor data')
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=radardate[0])

        figfname_list = make_filename(
            'ts_comp', prdcfg['dstype'], dataset['datatype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=radardate[0], timeformat='%Y%m%d')

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
        label2 = sensortype+' '+prdcfg['sensorid']
        titl = 'Time Series Comp. '+radardate[0].strftime('%Y-%m-%d')
        labely = generate_field_name_str(dataset['datatype'])

        plot_timeseries_comp(
            radardate, radarvalue, sensordate, sensorvalue, figfname_list,
            labelx='Time UTC', labely=labely, label1=label1, label2=label2,
            titl=titl, dpi=dpi)
        print('----- save to '+' '.join(figfname_list))

        return figfname_list

    elif prdcfg['type'] == 'COMPARE_CUMULATIVE_POINT':
        if dataset['final']:
            return None

        dpi = 72
        if 'dpi' in prdcfg:
            dpi = prdcfg['dpi']

        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir_ts = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')[0]

        csvfname = savedir_ts+csvfname

        radardate, radarvalue = read_timeseries(csvfname)
        if radardate is None:
            warn(
                'Unable to plot sensor comparison at point of interest. ' +
                'No valid radar data')
            return None

        sensordate, sensorvalue, sensortype, period2 = get_sensor_data(
            radardate[0], dataset['datatype'], prdcfg)
        if sensordate is None:
            warn(
                'Unable to plot sensor comparison at point of interest. ' +
                'No valid sensor data')
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=radardate[0])

        figfname_list = make_filename(
            'ts_cumcomp', prdcfg['dstype'], dataset['datatype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=radardate[0], timeformat='%Y%m%d')

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
        label2 = sensortype+' '+prdcfg['sensorid']
        titl = ('Time Series Acc. Comp. ' +
                radardate[0].strftime('%Y-%m-%d'))
        labely = 'Rainfall accumulation (mm)'

        plot_timeseries_comp(
            radardate, radarvalue, sensordate, sensorvalue,
            figfname_list, labelx='Time UTC', labely=labely,
            label1=label1, label2=label2, titl=titl,
            period1=prdcfg['ScanPeriod']*60., period2=period2, dpi=dpi)
        print('----- save to '+' '.join(figfname_list))

        return figfname_list

    elif prdcfg['type'] == 'COMPARE_TIME_AVG':
        if not dataset['final']:
            return None

        dpi = 72
        if 'dpi' in prdcfg:
            dpi = prdcfg['dpi']

        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir_ts = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=dataset['time'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=dataset['time'],
            timeformat='%Y%m%d')[0]

        radardate, radarvalue = read_timeseries(savedir_ts+csvfname)
        if radardate is None:
            warn(
                'Unable to compared time averaged data at POI. ' +
                'No valid radar data')
            return None

        sensordate, sensorvalue, sensortype, period2 = get_sensor_data(
            radardate[0], dataset['datatype'], prdcfg)
        if sensordate is None:
            warn(
                'Unable to compared time averaged data at POI. ' +
                'No valid sensor data')
            return None

        cum_time = 3600
        if 'cum_time' in prdcfg:
            cum_time = prdcfg['cum_time']

        base_time = 0
        if 'base_time' in prdcfg:
            base_time = prdcfg['base_time']

        sensordate_cum, sensorvalue_cum, np_sensor_cum = rainfall_accumulation(
            sensordate, sensorvalue, cum_time=cum_time, base_time=base_time,
            dropnan=False)

        radardate_cum, radarvalue_cum, np_radar_cum = rainfall_accumulation(
            radardate, radarvalue, cum_time=cum_time, base_time=base_time,
            dropnan=False)

        # find common time stamps
        ind = np.where(np.in1d(radardate_cum, sensordate_cum))[0]
        if not ind:
            warn('No sensor data for radar data time stamps')
        radardate_cum2 = radardate_cum[ind]
        radarvalue_cum2 = radarvalue_cum[ind]
        np_radar_cum2 = np_radar_cum[ind]

        ind = np.where(np.in1d(sensordate_cum, radardate_cum2))[0]
        sensorvalue_cum2 = sensorvalue_cum[ind]
        np_sensor_cum2 = np_sensor_cum[ind]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=radardate[0])

        fname = make_filename(
            str(cum_time)+'s_acc_ts_comp', prdcfg['dstype'],
            dataset['datatype'], ['csv'], prdcfginfo=gateinfo,
            timeinfo=radardate[0], timeformat='%Y%m%d')[0]

        fname = savedir+fname

        new_dataset = deepcopy(dataset)
        new_dataset.update({
            'time': radardate_cum2,
            'sensor_value': sensorvalue_cum2,
            'np_sensor': np_sensor_cum2,
            'radar_value': radarvalue_cum2,
            'np_radar': np_radar_cum2,
            'cum_time': cum_time})
        new_dataset.update(prdcfg)

        write_ts_cum(new_dataset, fname)

        print('saved CSV file: '+fname)

        figfname_list = make_filename(
            str(cum_time)+'s_acc_ts_comp', prdcfg['dstype'],
            dataset['datatype'], prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=radardate[0], timeformat='%Y%m%d')

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        labelx = sensortype+' '+prdcfg['sensorid']+' (mm)'
        labely = 'Radar (az, el, r): ('+az+', '+el+', '+r+') (mm)'
        titl = (str(cum_time)+'s Acc. Comp. ' +
                radardate_cum[0].strftime('%Y-%m-%d'))

        plot_scatter_comp(
            sensorvalue_cum2, radarvalue_cum2, figfname_list, labelx=labelx,
            labely=labely, titl=titl, axis='equal', dpi=dpi)

        print('----- save to '+' '.join(figfname_list))

        return figfname_list

    # ================================================================
    elif prdcfg['type'] == 'PLOT_AND_WRITE':

        timeinfo = dataset.time_vector[0]

        savedir = get_save_dir(prdcfg['basepath'], prdcfg['procname'],
                               dssavedir, prdcfg['prdname'],
                               timeinfo=timeinfo)

        dstype_str = prdcfg['dstype'].lower().replace('_', '')
        fname = make_filename('ts', dstype_str, dataset.datatype,
                              ['csv'],
                              prdcfginfo=None, timeinfo=timeinfo,
                              timeformat='%Y%m%d%H%M%S',
                              runinfo=prdcfg['runinfo'])

        dataset.write(savedir + fname[0])

        fname = make_filename('ts', dstype_str, dataset.datatype,
                              prdcfg['imgformat'],
                              prdcfginfo=None, timeinfo=timeinfo,
                              timeformat='%Y%m%d%H%M%S',
                              runinfo=prdcfg['runinfo'])

        ymin = None
        ymax = None
        if 'ymin' in prdcfg:
            ymin = prdcfg['ymin']
        if 'ymax' in prdcfg:
            ymax = prdcfg['ymax']

        dataset.plot(savedir + fname[0], ymin=ymin, ymax=ymax)

        return None

    elif prdcfg['type'] == 'PLOT_HIST':

        timeinfo = dataset.time_vector[0]

        savedir = get_save_dir(prdcfg['basepath'], prdcfg['procname'],
                               dssavedir, prdcfg['prdname'],
                               timeinfo=timeinfo)

        dstype_str = prdcfg['dstype'].lower().replace('_', '')

        fname = make_filename('hist', dstype_str, dataset.datatype,
                              prdcfg['imgformat'],
                              prdcfginfo=None, timeinfo=timeinfo,
                              timeformat='%Y%m%d%H%M%S',
                              runinfo=prdcfg['runinfo'])

        step = None
        if 'step' in prdcfg:
            step = prdcfg['step']

        dataset.plot_hist(savedir + fname[0], step=step)

        return None

    # ================================================================
    else:
        raise Exception("ERROR: Unsupported product type: '%s' of dataset '%s'"
                        % (prdcfg['type'], prdcfg['dsname']))
