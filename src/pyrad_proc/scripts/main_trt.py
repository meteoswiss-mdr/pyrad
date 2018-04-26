#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_trt
================================================

This program processes TRT data

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import os

import numpy as np

from pyrad.io import get_trtfile_list, read_trt_data, write_trt_cell_data

print(__doc__)


def main():
    """
    """
    database_path = '/data/TRT/'

    start_time_list = [
        '20170629000000', '20170630000000', '20170710000000', '20170718000000']
    end_time_list = [
        '20170629235900', '20170630235900', '20170710235900', '20170718235900']
    nsteps_min = 3

    print("====== TRT cell processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== comparison finished: ")

    for i, start_time_str in enumerate(start_time_list):
        end_time_str = end_time_list[i]

        starttime = datetime.datetime.strptime(start_time_str, '%Y%m%d%H%M%S')
        endtime = datetime.datetime.strptime(end_time_str, '%Y%m%d%H%M%S')

        data_input_path = (
            database_path+starttime.strftime('%Y-%m-%d')+'/TRTC/')
        data_output_path = (
            database_path+starttime.strftime('%Y-%m-%d')+'/TRTC_cell/')
        if not os.path.isdir(data_output_path):
            os.makedirs(data_output_path)

        flist = get_trtfile_list(data_input_path, starttime, endtime)

        traj_ID = np.array([], dtype=int)
        yyyymmddHHMM = np.array([], dtype=datetime.datetime)
        lon = np.array([], dtype=float)
        lat = np.array([], dtype=float)
        ell_L = np.array([], dtype=float)
        ell_S = np.array([], dtype=float)
        ell_or = np.array([], dtype=float)
        area = np.array([], dtype=float)
        vel_x = np.ma.array([], dtype=float)
        vel_y = np.ma.array([], dtype=float)
        det = np.ma.array([], dtype=float)
        RANKr = np.array([], dtype=int)
        CG_n = np.array([], dtype=int)
        CG_p = np.array([], dtype=int)
        CG = np.array([], dtype=int)
        CG_percent_p = np.ma.array([], dtype=float)
        ET45 = np.ma.array([], dtype=float)
        ET45m = np.ma.array([], dtype=float)
        ET15 = np.ma.array([], dtype=float)
        ET15m = np.ma.array([], dtype=float)
        VIL = np.ma.array([], dtype=float)
        maxH = np.ma.array([], dtype=float)
        maxHm = np.ma.array([], dtype=float)
        POH = np.ma.array([], dtype=float)
        RANK = np.ma.array([], dtype=float)
        Dvel_x = np.ma.array([], dtype=float)
        Dvel_y = np.ma.array([], dtype=float)
        cell_contour = []

        for fname in flist:
            print('Reading TRT file '+fname)

            (traj_ID_aux, yyyymmddHHMM_aux, lon_aux, lat_aux, ell_L_aux,
             ell_S_aux, ell_or_aux, area_aux, vel_x_aux, vel_y_aux, det_aux,
             RANKr_aux, CG_n_aux, CG_p_aux, CG_aux, CG_percent_p_aux,
             ET45_aux, ET45m_aux, ET15_aux, ET15m_aux, VIL_aux, maxH_aux,
             maxHm_aux, POH_aux, RANK_aux, Dvel_x_aux, Dvel_y_aux,
             cell_contour_aux) = read_trt_data(fname)

            traj_ID = np.append(traj_ID, traj_ID_aux)
            yyyymmddHHMM = np.append(yyyymmddHHMM, yyyymmddHHMM_aux)
            lon = np.append(lon, lon_aux)
            lat = np.append(lat, lat_aux)
            ell_L = np.append(ell_L, ell_L_aux)
            ell_S = np.append(ell_S, ell_S_aux)
            ell_or = np.append(ell_or, ell_or_aux)
            area = np.append(area, area_aux)
            vel_x = np.append(vel_x, vel_x_aux)
            vel_y = np.append(vel_y, vel_y_aux)
            det = np.append(det, det_aux)
            RANKr = np.append(RANKr, RANKr_aux)
            CG_n = np.append(CG_n, CG_n_aux)
            CG_p = np.append(CG_p, CG_p_aux)
            CG = np.append(CG, CG_aux)
            CG_percent_p = np.append(CG_percent_p, CG_percent_p_aux)
            ET45 = np.append(ET45, ET45_aux)
            ET45m = np.append(ET45m, ET45m_aux)
            ET15 = np.append(ET15, ET15_aux)
            ET15m = np.append(ET15m, ET15m_aux)
            VIL = np.append(VIL, VIL_aux)
            maxH = np.append(maxH, maxH_aux)
            maxHm = np.append(maxHm, maxHm_aux)
            POH = np.append(POH, POH_aux)
            RANK = np.append(RANK, RANK_aux)
            Dvel_x = np.append(Dvel_x, Dvel_x_aux)
            Dvel_y = np.append(Dvel_y, Dvel_y_aux)
            cell_contour.append(cell_contour_aux)

        traj_ID_unique_list = np.unique(traj_ID)

        print('Total Number of cells: '+str(traj_ID_unique_list.size))

        ncells = 0
        for traj_ID_unique in traj_ID_unique_list:
            ind = np.where(traj_ID == traj_ID_unique)[0]

            traj_ID_cell = traj_ID[ind]
            yyyymmddHHMM_cell = yyyymmddHHMM[ind]
            lon_cell = lon[ind]
            lat_cell = lat[ind]
            ell_L_cell = ell_L[ind]
            ell_S_cell = ell_S[ind]
            ell_or_cell = ell_or[ind]
            area_cell = area[ind]
            vel_x_cell = vel_x[ind]
            vel_y_cell = vel_y[ind]
            det_cell = det[ind]
            RANKr_cell = RANKr[ind]
            CG_n_cell = CG_n[ind]
            CG_p_cell = CG_p[ind]
            CG_cell = CG[ind]
            CG_percent_p_cell = CG_percent_p[ind]
            ET45_cell = ET45[ind]
            ET45m_cell = ET45m[ind]
            ET15_cell = ET15[ind]
            ET15m_cell = ET15m[ind]
            VIL_cell = VIL[ind]
            maxH_cell = maxH[ind]
            maxHm_cell = maxHm[ind]
            POH_cell = POH[ind]
            RANK_cell = RANK[ind]
            Dvel_x_cell = Dvel_x[ind]
            Dvel_y_cell = Dvel_y[ind]
            cell_contour_cell = None

            if traj_ID_cell.size < nsteps_min:
                continue

            fname = data_output_path+str(traj_ID_unique)+'.trt'
            fname = write_trt_cell_data(
                traj_ID_cell, yyyymmddHHMM_cell, lon_cell, lat_cell,
                ell_L_cell, ell_S_cell, ell_or_cell, area_cell, vel_x_cell,
                vel_y_cell, det_cell, RANKr_cell, CG_n_cell, CG_p_cell,
                CG_cell, CG_percent_p_cell, ET45_cell, ET45m_cell, ET15_cell,
                ET15m_cell, VIL_cell, maxH_cell, maxHm_cell, POH_cell,
                RANK_cell, Dvel_x_cell,
                Dvel_y_cell, cell_contour_cell, fname)

            print('Written individual TRT cell file '+fname)
            ncells += 1

        print('Number of cells with '+str(nsteps_min) +
              ' or more time steps: '+str(ncells))


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
