# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 10:48:31 2016

@author: fvanden

Configuration file for mxpol pyart.core.Radar class. Some information may be 
redundant because this file is a copy from the ProfileLab toolkit.

Functions to retrieve data from this file may be found in pyrad.io.read_data_mxpol
under the utilities section

"""

# radar information 

MCH_elev=[-0.2,0.4,1,1.6,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,11,13,16,20,25,30,35,40]
NYQUIST_VEL=[8.3, 9.6,8.3,12.4,11.0,12.4,13.8,12.4,13.8,16.5,16.5,16.5,20.6,20.6,20.6,20.6,20.6,20.6,20.6,20.6]
 
RADAR_INFO = {
    'coordinates' : {
        'ALB' : [47.284,8.512],
        'DOL' : [46.425,6.099],
        'PPM' : [46.371,7.487],
        'MLE' : [46.041,8.833],
        'DX50' : [46.8425,6.9184],
        'MXPOL' : [46.8133,6.9428]
        },
    'altitude': {
        'ALB' : 938,
        'DOL' : 1682,
        'PPM' : 2937,
        'MLE' : 1626,
        'DX50': 451,
        'MXPOL': 489
        },
    'searchkey' : {
        'ALB' : 'PHA*hdf*',
        'DOL' : 'PHD*hdf*',
        'PPM' : 'PHP*hdf*',
        'MLE' : 'PHL*hdf*',
        'DX50' : None,
        'MXPOL' : None        
        },
    'radarID' : {
        'ALB' : 'ALB',
        'A':'ALB',
        'DOL' : 'DOL',
        'D':'DOL',
        'PPM' : 'PPM',
        'P':'PPM',
        'MLE' : 'MLE',
        'M':'MLE',
        'DX50' : 'DX50',
        'MXPOL' : 'MXPOL'
        },
    'dbbeam' : {
        'ALB' : 1.,
        'DOL' : 1.,
        'PPM' : 1.,
        'MLE' : 1.,
        'MXPOL' : 1.4,
        'DX50' : 1.27
        },
    'elevations': {
        'ALB': MCH_elev,
        'DOL': MCH_elev,
        'PPM': MCH_elev,
        'MLE' : MCH_elev,
        'DX50': None,
        'MXPOL': None
        }
    }

MY_METADATA = {

        'nyq_vel' : NYQUIST_VEL,

        # Metadata for instrument tables

        'Radar_info' : {
            'searchkey' : None,
            'coordinates' : None,
            'altitude' : None, 
            'dbbeam' : None,
            'radarID' : None},
            
        'Parsivel_info' : {
            'parsID' : None, 
            'coordinates' : None,
            'altitude' : None,
            'filepath' : None,
            'long_name': None},
    
        # Metadata for VerticalProfile attributes
        
        'Time' : {
            'units': 'Python time struct',
            'standard_name': 'time',
            'long_name': 'time of measurement'},                    
        
        'Profile_Type' : {
            'standard_name' : 'Profile type',
            'long_name' : 'time_average, spatial_average or instantaneous profile'},
            
        'Polvar' : {
            'units' : None,
            'standard_name' : None,
            'short_name' : None,
            'long_name' : None,
            'valid_min': None, 
            'valid_max': None,
            'plot_interval' : None},
            
        'Time_average_over' : {
            'units' : 'minutes',
            'standard_name' : 'time_average_over',
            'long_name' : 'time period over which average vertical profile has been calculated'},
            
        'Spatial_average_over' : {
            'units' : 'metres',
            'standard_name' : 'spatial_average_over',
            'long_name' : 'spatial extent over which average vertical profile has been calculated'},
            
        'Scan_list' : {
            'standard_name' : 'scan_list',
            'long_name' : 'list of databasenames for scans used to create vertical profile'},
            
        'Vp_calculation' : {
            'standard_name' : 'vp_calculation',
            'long_name' : 'dictionary with information on how vp calculation was performed'},
            
        'Profiles' : {
            'standard_name' : 'profiles',
            'long_name' : 'extracted vertical profiles',
            'cartesian_locations_on_grid': None,
            'cartesian_locations_distance_from_radar' : None,
            'latitude_locations' : None,
            'longitude_locations' : None,
            'type' : None},
            
        'Grid_info' : {
            'standard_name' : 'grid_info',
            'long_name' : 'information concerning grid projection on from which profiles were extracted',
            'latcoords' : None,
            'loncoords' : None },
            
        # metadata for grid
            
        'projection' : {
            'standard_name' : 'projection',
            'long_name' : 'information on the algorithm used for the projection of the instrument',
            'proj': None,
            'script': None,
            'date_script': None },
            
        # Metadata for MetInfo attributes
            
        'timeperiod' : {
            'units' : 'Python time struct',
            'standard_name' : 'timeperiod',
            'long_name' : 'time periods for start/end of measurement'},
            
        'location' : {
            'standard_name' : 'location',
            'long_name' : 'geographical location and name of the meteorological station of the meteorological data',
            'ID' : None, 
            'station': None,
            'canton' : None,            
            'lat' : None,
            'lon' : None,
            'altitude' : None,
            'Ind_OMM' : None,
            'Ind_Nat' : None,
            'distance' : None,
            'variables' : None },
            
        'met_var' : {
            'standard_name' : None,
            'long_name' : None, 
            'units' : None,
            'valid_min' : None,
            'valid_max' : None },
            
        # Metadata for extra field names     
            
        'radar_echo_id': {
            'units': '-',
            'standard_name': 'radar_echo_id',
            'long_name': 'Radar Echo Identification',
            'coordinates': 'elevation azimuth range'},
            
            }
            
MY_POLARNAMES = {

        # Metadata for polarimetric short and long names

    'Zh' : ['reflectivity','reflectivity','dBZ', 0., 55.,1.],
    'Zdr' : ['differential_reflectivity','Differential reflectivity', 'dB', -1., 5.,0.1],
    'Kdp' : ['specific_differential_phase','Specific differential phase','deg/km',-2., 7., 0.1],
    'Phidp' : ['differential_phase','Differential phase', 'deg',0., 150.,1.],
    'Rhohv' : ['cross_correlation_ratio','Copolar correlation coefficient', '-',0.57, 1., 0.05],
    'ZhCorr' : ['corrected_unfiltered_reflectivity','Attenuation corrected reflectivity', 'dBZ', 0., 55.,1.], 
    'ZdrCorr' : ['corrected_unfiltered_reflectivity_vv','Attenuation corrected reflectivity','dB', 0., 3., 0.1],
    'RVel' : ['velocity','Mean doppler velocity','m/s', -15., 15.,0.5],
    'Sw' : ['spectrum_width','Spectral Width','m2/s2', 0., 3., 0.1],
    'Zv' : ['reflectivity_vv', 'Vertical reflectivity','dBZ', 0., 45., 1.],
    'Clut' : ['clutter', 'Output clutter algorithm','-',0.,100.,10.], 
    'corrected_Z' : ['corrected_reflectivity', 'Clutter filtered reflectivity', 'dBZ', 0., 55., 1.]        

    }