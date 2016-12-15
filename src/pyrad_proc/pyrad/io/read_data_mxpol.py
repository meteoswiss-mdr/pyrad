"""
pyrad.io.read_data_mxpol
========================

Functions for reading radar mxpol data files

.. autosummary::
    :toctree: generated/

    pyrad_MXPOL
    findTimes
    row_stack
    int2float_radar
    

"""
import pyart
import imp
import numpy as np
import re
import os
import datetime
import netCDF4


class pyrad_MXPOL(pyart.core.Radar):
    def __init__(self,filename, field_names, max_range=np.Inf,min_range=10000):
        
        print(field_names)
        
        # find information based on filename
        
        all_files = [filename]
        
        fname_basename=os.path.basename(filename)
        
        if 'PPI' in fname_basename:
            scan_type='ppi'
        elif 'RHI' in fname_basename:
            scan_type='rhi'
            
        strdate=re.findall(r"([0-9]{8}-[0-9]{6})",fname_basename)[0]
        date=datetime.datetime.strptime(strdate,'%Y%m%d-%H%M%S')
        
        # convert fieldname if necessary
        
        varnames = []
        for fieldname in field_names:
            newname = convert_polvar_name('LTE', fieldname)
            varnames.append(newname)
            
        print(varnames)
            
        # get labels, units etc
            
        long_names = []
        standard_names = []
        units = []
        vmin = []
        vmax = []
        
        for varname in varnames:
            metadata = generate_polvar_metadata(varname)
            standard_names.append( metadata['standard_name'] )
            long_names.append( metadata['long_name'] )
            units.append( metadata['units'] )
            vmin.append( metadata['valid_min'] )
            vmax.append( metadata['valid_max'] )
            
        # initiate empty vectors
            
        N_sweeps=len(all_files)
        fields={}
        fixed_angle={}
        fixed_angle['data']=np.zeros(N_sweeps,)
        
        sweep_start_ray_index={}
        sweep_start_ray_index['data']=[]
        sweep_stop_ray_index={}
        sweep_stop_ray_index['data']=[]
        
        for i,k in enumerate(varnames):
            fields[k]={}
            fields[k]['data']=[]
            fields[k]['long_name']=long_names[i]
            fields[k]['standard_name'] = standard_names[i]
            fields[k]['units']=units[i]
            fields[k]['valid_min']=vmin[i]
            fields[k]['valid_max']=vmax[i]
            
        idx_start=0
        idx_stop=0
        elevations=[]
        azimuths=[]
        nyquist=[]
        
        # read data and create dictionaries

        for i in range(N_sweeps):
            metadata, data = readMXPOLRadData(all_files[i],varnames,max_range)
            if scan_type == 'rhi':
                fixed_angle['data'] = data['azimuth']
            elif scan_type == 'ppi':
                fixed_angle['data'] = data['elevation']
                
            [N_az,N_ranges]=data[varnames[0]].shape
            idx_stop=idx_start+N_az-1
            sweep_start_ray_index['data'].append(idx_start)
            sweep_stop_ray_index['data'].append(idx_stop)
            idx_start=idx_stop+1
            elevations.extend(list(data['elevation']))
            nyquist.extend([data['nyquist_vel']]*N_az)
            azimuths.extend(list(data['azimuth']))
            
            for j,v in enumerate(varnames):
                if v in data.keys():
                    if not len(fields[v]['data']):
                        fields[v]['data']=data[v]
                    else:
                        fields[v]['data']=row_stack(fields[v]['data'],data[v])
                else:
                    print('Variable '+v+' was not found in file!')
        
        # mask NaNs
            
        for v in varnames:
            fields[v]['data']=np.ma.masked_equal(fields[v]['data'], -99900.0)
            
        [a,N_ranges]=fields[varnames[0]]['data'].shape
        
        # create dictionaries according to pyART standard  
        
        latitude={'data' : np.array([data['latitude']]), 'units' : data['lat_units']}
        longitude={'data' :np.array([data['longitude']]), 'units' : data['lon_units']}    
        altitude={'data' : np.array([data['altitude']]), 'units' : data['alt_units']}     
        sweep_number={'data' : np.arange(0,len(all_files))}     
        sweep_mode={'data' : [scan_type]*N_sweeps}   
        instrument_parameters={'nyquist_velocity': {'data':np.array(nyquist)}}
        azimuth={'data' : np.array(azimuths), 'units' : data['azim_units']}
        rrange={'data':np.arange(N_ranges)*data['resolution'], 'units' : data['range_units']}
        elevation={'data' :np.array(elevations), 'units' : data['elev_units']}
        
        time_units='seconds since '+str(date)
        time={'data' : data['time'],'units': time_units}
            
        # change keys to match pyART metranet keys
        for keys in fields:
            newkey = fields[keys]['standard_name']
            fields[newkey] = fields.pop(keys)
            
        # Create PyART instance
        pyart.core.Radar.__init__(self,time,rrange,fields,metadata,scan_type,latitude,longitude,altitude,sweep_number,sweep_mode,fixed_angle,\
        sweep_start_ray_index,sweep_stop_ray_index,azimuth,elevation,instrument_parameters=instrument_parameters)
            
        
        
########################## utilities - read ##################################
        
def row_stack(a1,a2):
    """
    """
    [N1,M1]=a1.shape
    [N2,M2]=a2.shape

    if M1>M2:       # changed M1>M2 instead of N1>N2
        a2=np.pad(a2,((0,0),(0,M1-M2)),mode='constant',constant_values=-9999999)
    elif M2<M1:     # changed M1<M2 instead of N1<N2
        a1=np.pad(a2,((0,0),(0,M2-M1)),mode='constant',constant_values=-9999999)
    
    out = np.vstack((a1,a2))
    out[out == -9999999] = np.nan
    
    return out
        
def readMXPOLRadData(filename, variableList, max_range=np.Inf,min_range=0):
    """
    Reads a netcdf containing processed radar data in polar coordinates
    
    Parameters
    ----------
    filename: str
        complete path of the file
        
    variableList: list
        list of variables to be read
        
    Returns
    -------
    varPol: dict
        dictionary containing the variables, the azimuth and the range
        
    metadata: dict
        dictionary containing the metadata of the file
        
    """

    varPol={}
    metadata={}
    ncid=netCDF4.Dataset(filename)
    
    time=ncid.variables['Time']
    time-=time[0] # To get time in seconds from beginning of scan

    rrange= ncid.variables['Range'][:]
        
    # Get indexes between min_range and max_range
    idx2keep=np.where(np.logical_and(rrange<max_range,rrange>min_range))[0]
    rrange=rrange[idx2keep]
    
    # Get variables in polar coordinates
    for varname in variableList:
        try:
            varPol[varname]=ncid.variables[varname][:].T
        except:
            pass
        
    varPol['resolution']=ncid.__dict__['RangeResolution-value']
    varPol['range']=rrange
    varPol['range_units'] = ncid.__dict__['RangeResolution-unit']
    varPol['azimuth']=ncid.variables['Azimuth'][:]
    varPol['azim_units']=ncid.__dict__['Azimuth-unit'] 
    varPol['elevation']=ncid.variables['Elevation'][:]
    varPol['elev_units']=ncid.__dict__['Elevation-unit']
    varPol['nyquist_vel']=ncid.__dict__['NyquistVelocity-value']
    varPol['longitude']=ncid.__dict__['Longitude-value']
    varPol['lon_units'] = ncid.__dict__['Longitude-unit']
    varPol['latitude']=ncid.__dict__['Latitude-value']   
    varPol['lat_units']=ncid.__dict__['Latitude-unit']
    varPol['altitude']=ncid.__dict__['Altitude-value']
    varPol['alt_units']=ncid.__dict__['Altitude-unit']
    varPol['time']=time
    
    metadata['Source'] = ncid.__dict__['Source']
    metadata['Institution'] = ncid.__dict__['Institution']
    metadata['History'] = ncid.__dict__['History']
    metadata['ContactInformation'] = ncid.__dict__['ContactInformation']
    
    # Close netcdf
    ncid.close()

    return metadata, varPol


######################### utilities - config #################################

_dirname = os.path.dirname(__file__)
_DEFAULT_CONFIG_FILE = os.path.join(_dirname, 'mxpol_config.py')

def load_myconfig(filename = None):
    """
    Load configuration from a config file.
    
    Parameters
    ----------
    filename: str
        Filename of the configuration file. If None the default configuration
        file is loaded from the directory.
        
    Returns
    -------
    _DEFAULT_METADATA: dict
        Dictionary with metadata
        
    """
    
    if filename is None:
        filename = _DEFAULT_CONFIG_FILE
        
    # private:
    
    global cfile
    global _DEFAULT_POLARNAMES
    global _DEFAULT_METADATA
    
    cfile = imp.load_source('metadata_config', filename)
    _DEFAULT_METADATA = cfile.MY_METADATA
    _DEFAULT_POLARNAMES = cfile.MY_POLARNAMES

    return _DEFAULT_METADATA
    
def get_mymetadata(p, filename = None):
    """
    Return a dictionary of metadata for a given parameter, p.

    An empty dictionary will be returned if no metadata dictionary exists for
    parameter p.
    
    Parameters
    ----------
    p: str
        parameter name (i.e. Polvar) for which to return metadata
    
    filename: str
        Filename of the configuration file. If None the default configuration
        file is loaded from the directory.
        
    Returns
    -------
    _DEFAULT_METADATA[p].copy(): dict
        a copy of the parameter of interest from the metadata dictionary
        
    """
    load_myconfig(filename = filename)    
    
    if p in _DEFAULT_METADATA:
        return _DEFAULT_METADATA[p].copy()
    else:
        return {}

def generate_polvar_metadata(polvar, filename = None):
    """
    Generates a dictionary with metadata for a polarimetric variable
    
    Parameters
    ----------
    polvar: str
        polatimetric variable of interest
        
    filename: str
        Filename of the configuration file. If None the default configuration
        file is loaded from the directory.
        
    Returns
    -------
    polvar_metadata: dict
        dictionary with metatdata for polarimetric variable of interest
        
    """
    load_myconfig(filename = filename)   
    polvar = convert_polvar_name('LTE', polvar)
    
    if polvar in _DEFAULT_POLARNAMES:
        standard_name, long_name, units, valid_min, valid_max, plot_interval =  _DEFAULT_POLARNAMES[polvar]
    else:
        standard_name, long_name, units, valid_min, valid_max, plot_interval = None, None, None, None, None, None
        
    polvar_metadata = get_mymetadata('Polvar', filename)
    polvar_metadata['units'] = units
    polvar_metadata['standard_name'] = standard_name
    polvar_metadata['short_name'] = convert_polvar_name('MCH', polvar)
    polvar_metadata['long_name'] = long_name
    polvar_metadata['valid_min'] = valid_min
    polvar_metadata['valid_max'] = valid_max
    polvar_metadata['plot_interval'] = plot_interval
    
    return polvar_metadata

def convert_polvar_name(convention, polvar):
    """
    Finds the correct variable name for a given convention (MXPOL, MCH) and 
    a given variable name which was spelled with a different case or 
    according to a different convention. For example, MXPOL convention uses 
    'Z' for the reflectivity variable, but if a user inserted 'Zh' this 
    function will convert it to 'Z'. 
    
    Parameters
    ----------
    convention : str, destination convention; either MCH or LTE 
            
    polvar : str, key of polarimetric variable to be converted
        
    Returns
    -------
    mykey : str, polarimertric variable key as used within the ProfileLab
        toolbox context
            
    """
    # Generate dictionary for the conversion
    
    MCH = ['ZH','ZV','ZDR','PHI','VEL','WID', 'RHO','CLUT', 'MPH','STA1', 'STA2', 'WBN']
    LTE = ['Zh','Zv','Zdr','Phidp','RVel','Sw','Rhohv','Clut', 'MPH','STA1', 'STA2', 'WBN']

    
    convertkeys = {}
    convertkeys['MCH'] = {}
    for i in range(0,len(MCH)):
        convertkeys['MCH'][MCH[i]] = LTE[i]
        
    convertkeys['LTE'] = {}
    for i in range(0,len(LTE)):
        convertkeys['LTE'][LTE[i]] = MCH[i]
    
    # translate between conventions
    mykey = polvar
    
    for key, value in convertkeys[convention].items():
        if polvar in value:
            mykey = key
            break
                
    return mykey
    