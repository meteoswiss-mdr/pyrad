#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
generate_tasks_file
================================================

This script generates a tasks file for parallel computing at CSCS

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import os
import pandas as pd

print(__doc__)


def main():
    """
    """
    ftasks_name = os.path.expanduser('~')+'/pyrad/tools/processData/tasks'
    task = os.path.expanduser('~')+'/pyrad/tools/processData/get_and_process_rad4alp_data_cscs.sh'

    # generate list of dates
    start_date = datetime.datetime(2018, 1, 1)
    ndays = 2
    datelist = pd.date_range(start_date, periods=ndays).tolist()

    datelist.append(datetime.datetime(2018, 3, 7))

    start_time = '000001'
    end_time = '001000'

    cfgfiles = 'cscs_rad4alp_PLA_postproc.txt'
    radar = 'A'
    res = 'L'

    get_data = 1
    rm_data = 0
    data_destbase = '/store/msrad/radar/rad4alp/tmp/'

    get_cosmo = 0
    rm_cosmo = 0
    cosmo_destbase = '/store/msrad/cosmo/tmp/TEMP/raw1/'

    get_hzt = 0
    rm_hzt = 0
    hzt_destbase = '/store/msrad/cosmo/tmp/HZT/'

    print("====== tasks file generation started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== tasks file generation finished: ")

    with open(ftasks_name, 'w', newline='') as txtfile:
        for day in datelist:
            txtfile.write(
                task+' -d '+day.strftime("%Y%m%d")+' --start_time '+start_time +
                ' --end_time '+end_time+' -c '+cfgfiles +' -r '+radar+' -e '+res +
                ' --get_data '+str(get_data)+' --rm_data '+str(rm_data) +
                ' --data_destbase '+data_destbase +
                ' --get_cosmo '+str(get_cosmo)+' --rm_cosmo '+str(rm_cosmo) +
                ' --cosmo_destbase '+cosmo_destbase +
                ' --get_hzt '+str(get_hzt)+' --rm_hzt '+str(rm_hzt) +
                ' --hzt_destbase '+hzt_destbase+'\n')

        txtfile.close()


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
