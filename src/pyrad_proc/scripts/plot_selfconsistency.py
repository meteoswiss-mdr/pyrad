#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
plot_antenna_pattern
================================================

This program plots an antenna pattern provided in a .csv file

"""

# Author: fvj
# License: BSD 3 clause
import numpy as np

from pyrad.io import read_selfconsistency
from pyrad.graph import plot_selfconsitency

print(__doc__)


def main():
    """
    """
    
    zdr = np.arange(0., 6., 0.05)
    file_path = '/users/jfigui/pyrad/config/selfconsistency/'
    img_ext = '.png'
    
    # Vaccarono
    kdpzh = 1.77e-4*np.power(np.power(10., 0.1*zdr), -2.09)    
    zdrkdp_table = [zdr, kdpzh]
    base_name = 'selfconsistency_zdr_zh091kdp085_Cband_temp00_elev000_mu05_Vaccarono'
    plot_selfconsitency(
            zdrkdp_table, [file_path+base_name+img_ext], labelx='ZDR [dB]',
            title=base_name, ymin=None, ymax=None, save_fig=True)
            
    return
    
    # Gorgucci expression, probably Zh and Zdr in linear units
    kdpzh = 1.82e-4*np.power(np.power(10., 0.1*zdr), -1.28)    
    zdrkdp_table = [zdr, kdpzh]
    base_name = 'selfconsistency_zdr_zh095kdp_Cband_temp00_elev000_mu05_Gorgucci'
    plot_selfconsitency(
            zdrkdp_table, [file_path+base_name+img_ext], labelx='ZDR [dB]',
            title=base_name, ymin=None, ymax=None, save_fig=True)
            
    return

    file_path = '/users/jfigui/pyrad/config/selfconsistency/'
    base_name_vec = ['selfconsistency_zdr_zhkdp_Cband_temp10_elev000_mu05']
    # base_name_vec = ['selfconsistency_zdr_zhkdp_Xband_temp10_elev000_mu05']
    img_ext = '.png'

    for base_name in base_name_vec:
        zdrkdp_table = read_selfconsistency(file_path+base_name+'.txt')
        fig, ax = plot_selfconsitency(
            zdrkdp_table, [file_path+base_name+img_ext], labelx='ZDR [dB]',
            title=base_name, ymin=0, ymax=0.00008, save_fig=False)

    zdr = np.arange(0., 6., 0.05)
    
    # Gorgucci expression, probably Zh and Zdr in linear units
    kdpzh = 1.82e-4*np.power(np.power(10., 0.1*zdr), -1.28)    
    zdrkdp_table = [zdr, kdpzh]
    base_name = 'selfconsistency_zdr_zh095kdp_Cband_temp00_elev000_mu05_Gorgucci'
    plot_selfconsitency(
            zdrkdp_table, [file_path+base_name+img_ext], labelx='ZDR [dB]',
            title=base_name, ymin=0, ymax=0.00008, save_fig=True)
            
     
    
    # Louf
    kdpzh = 1e-5*(6.607-4.577*zdr+1.577*zdr*zdr-0.23*zdr*zdr*zdr)
    zdrkdp_table = [zdr, kdpzh]

    base_name = 'selfconsistency_zdr_zhkdp_Cband_tempXX_elev000_Louf'
    plot_selfconsitency(
            zdrkdp_table, [file_path+base_name+img_ext], labelx='ZDR [dB]',
            title=base_name, ymin=0, ymax=0.00008, save_fig=True, fig=fig,
            ax=ax)
            
    return
    
    
    # Wolfensberger C-band
    kdpzh = 0.000053*np.power(zdr, -0.08054)*np.exp(-0.247351*zdr)
    zdrkdp_table = [zdr, kdpzh]

    base_name = 'selfconsistency_zdr_zhkdp_Cband_tempXX_elev016_Wolfensberger'
    plot_selfconsitency(
            zdrkdp_table, [file_path+base_name+img_ext], labelx='ZDR [dB]',
            title=base_name, ymin=0, ymax=0.00008, save_fig=True, fig=fig,
            ax=ax)
            
    return
    
    # Gourley X-band
    kdpzh = 1e-5*(11.74-4.020*zdr-0.140*zdr*zdr+0.130*zdr*zdr*zdr)
    zdrkdp_table = [zdr, kdpzh]

    base_name = 'selfconsistency_zdr_zhkdp_Xband_temp00_elev000_mu05_Gourley'
    plot_selfconsitency(
            zdrkdp_table, [file_path+base_name+img_ext], labelx='ZDR [dB]',
            title=base_name, ymin=0, ymax=0.00015, save_fig=True, fig=fig,
            ax=ax)
            
    return
            

    # Gourley C-band
    kdpzh = 1e-5*(6.746-2.970*zdr+0.711*zdr*zdr-0.079*zdr*zdr*zdr)
    zdrkdp_table = [zdr, kdpzh]

    base_name = 'selfconsistency_zdr_zhkdp_Cband_temp00_elev000_mu05_Gourley'
    plot_selfconsitency(
            zdrkdp_table, [file_path+base_name+img_ext], labelx='ZDR [dB]',
            title=base_name, ymin=0, ymax=0.00008, save_fig=True, fig=fig,
            ax=ax)

    # Gourley S-band
    kdpzh = 1e-5*(3.696-1.963*zdr+0.504*zdr*zdr-0.051*zdr*zdr*zdr)
    zdrkdp_table = [zdr, kdpzh]

    base_name = 'selfconsistency_zdr_zhkdp_Sband_temp00_elev000_mu05_Gourley'
    plot_selfconsitency(
            zdrkdp_table, [file_path+base_name+img_ext], labelx='ZDR [dB]',
            title=base_name, ymin=0, ymax=0.00008, save_fig=True)

    # Gourley X-band
    kdpzh = 1e-5*(11.74-4.020*zdr-0.140*zdr*zdr+0.130*zdr*zdr*zdr)
    zdrkdp_table = [zdr, kdpzh]

    base_name = 'selfconsistency_zdr_zhkdp_Xband_temp00_elev000_mu05_Gourley'
    plot_selfconsitency(
            zdrkdp_table, [file_path+base_name+img_ext], labelx='ZDR [dB]',
            title=base_name, ymin=0, ymax=0.00008, save_fig=True)

    # Gorgucci expression, probably Zh and Zdr in linear units
    kdpzh = 1.82e-4*np.power(np.power(10., 0.1*zdr), -1.28)    
    zdrkdp_table = [zdr, kdpzh]
    base_name = 'selfconsistency_zdr_zh095kdp_Cband_temp00_elev000_mu05_Gorgucci'
    plot_selfconsitency(
            zdrkdp_table, [file_path+base_name+img_ext], labelx='ZDR [dB]',
            title=base_name, ymin=0, ymax=0.00008, save_fig=True)


# ---------------------------------------------------------
# Start main:
# ---------------------------------------------------------
if __name__ == "__main__":
    main()
