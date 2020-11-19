#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

This program processes and post-processes data over a time span

To run the processing framework type:
    python main_process_data.py \
[config_file] --starttime [process_start_time] --endtime [process_end_time] \
--postproc_cfgfile [postproc_config_file] --cfgpath [cfgpath]

If startime and endtime are not specified the program determines them from
the trajectory file or the last processed volume.
postproc_cfgfile is an optional argument with default: None
cfgpath is an optional argument with default: \
'$HOME/pyrad/config/processing/'
The trajectory file can be of type plane, lightning or proc_periods. If it is \
of type lightning the flash number can be specified

Example:
    python main_process_data.py 'paradiso_fvj_vol.txt' --starttime \
'20140523000000' --endtime '20140523001000' --postproc_cfgfile \
'paradiso_fvj_vol_postproc.txt' --cfgpath '$HOME/pyrad/config/processing/'

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import os
import subprocess
from warnings import warn
import glob

import pandas as pd
import numpy as np

print(__doc__)


def main():
    """
    """

    print("====== PYRAD data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD data processing finished: ")

    basepath = '/store/msrad/radar/pyrad_products/rad4alp_POH/'

    # read the 3 CSV files containing POH statistics
    df_2018 = pd.read_csv(
        basepath+'2018_POH_above_90_percent.csv',
        usecols=['nominal time', 'day/night', 'npixels with POH > 90%'])
    df_2018['nominal time'] = pd.to_datetime(df_2018['nominal time'])
    print('Number of time steps with POH > 90% in 2018:', df_2018.shape[0])

    df_2019 = pd.read_csv(
        basepath+'2019_POH_above_90_percent.csv',
        usecols=['nominal time', 'day/night', 'npixels with POH > 90%'])
    df_2019['nominal time'] = pd.to_datetime(df_2019['nominal time'])
    print('Number of time steps with POH > 90% in 2019:', df_2019.shape[0])

    df_2020 = pd.read_csv(
        basepath+'2020_POH_above_90_percent.csv',
        usecols=['nominal time', 'day/night', 'npixels with POH > 90%'])
    df_2020['nominal time'] = pd.to_datetime(df_2020['nominal time'])
    print('Number of time steps with POH > 90% in 2020:', df_2020.shape[0])

    df = pd.concat([df_2018, df_2019, df_2020])
    print('Total number of time steps with POH > 90%:', df.shape[0])

    # separate the data frame into day, night and transition. Check the number of images
    df_day = df[df['day/night'] == 'day']
    df_night = df[df['day/night'] == 'night']
    df_trans = df[df['day/night'] == 'transition']
    print('day light images:', df_day.shape[0])
    print('night time images:', df_night.shape[0])
    print('transition images:', df_trans.shape[0])
    
    df_day['just_date'] = df_day['nominal time'].dt.date
    
    days = np.unique(df_day['just_date'])
    
    print('number of different days:', days.size)
    
    for day in days:
        print('POH data for day '+day.strftime('%Y%m%d'))
        # cp = subprocess.run(
        #     ["./get_rad4alp_cartesian_data_CSCS_2.sh", '-d '+day.strftime('%Y%m%d'), '--prod BZC', '-p /store/msrad/radar/rad4alp/POH/'],
        #     shell=False, check=True, capture_output=True)
        # print(cp)
        os.system('./get_rad4alp_cartesian_data_CSCS_2.sh -d '+day.strftime('%Y%m%d')+' --prod BZC -p /store/msrad/radar/rad4alp/POH/')
        
   

def _print_end_msg(text):
    """
    prints end message

    Parameters
    ----------
    text : str
        the text to be printed

    Returns
    -------
    Nothing

    """
    print(text + datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))


# ---------------------------------------------------------
# Start main:
# ---------------------------------------------------------
if __name__ == "__main__":
    main()   