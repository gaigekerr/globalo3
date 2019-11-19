#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create GEOSChem emission inventory for passive species

This module uses an existing emission inventory for GEOSChem CO2 emissions
and mimicks this file structure and its attributes to create a time-invarient 
emission inventory for passive tracers as well as adapting CCMI CO50 tracer
for use in GEOSChem
    
Revision History
----------------
    20092019 -- initial version created
    01102019 -- added tracer (GLOBAL_50) with global emissions
    18112019 -- changed bands of tracer emissions such that instead of just 
                having three bands (i.e., 0-30N/S, 30-60N/S, and 60-90N/S), 
                there are distinct tracers every 10 deg 
    19112019 -- adapt RETRO emissions for CO (year 2000) for use as a passive
                tracer in GEOSChem
"""
def create_passive10degtracer(): 
    """Data variables in output netCDF file are for a passive species with 
    a 50 d lifetime with annual, time-invariant emissions for 10 degree 
    latitude bands (i.e., 0-10 deg N/S, 10-20 deg N/S, etc.) and a tracer 
    with global emissions. 

    Parameters
    ----------
    None
    
    Returns
    -------
    None
    """
    import datetime
    import xarray as xr
    import numpy as np
    INVENTORYPATH = '/Users/ghkerr/phd/globalo3/GEOSChem/'
    inventoryname = 'fossilCO2.annual.geos.2x25.nc'
    # Open existing inventory
    inventory = xr.open_dataset(INVENTORYPATH+inventoryname)
    # Extract coordinate data
    lng = inventory.lon.data
    lat = inventory.lat.data
    time = inventory.time.data
    # Extract emissions data for CO2
    emiss = inventory.CO2.data
    # Based on a visual inspection of the CO2 emissions inventory, a 
    # reasonable CO2 emission rate (kg/m2/s) in the midlatitudes is ~5e-8
    #import matplotlib.pyplot as plt
    #ax = plt.subplot2grid((1, 1), (0, 0))
    #mb = ax.contourf(np.nanmean(emiss, axis=0), np.linspace(0, 5e-8, 10), 
    # extend='max')
    #plt.colorbar(mb, extend='max')
    emiss_rate = 5e-9
    # 0-10 deg N/S
    trac50_0_10 = np.zeros(emiss.shape)
    trac50_0_10[:, (np.abs(lat-0.)).argmin():
        (np.abs(lat-10.)).argmin(), :] = emiss_rate
    trac50_0_10[:, (np.abs(lat-(-10.))).argmin():
        (np.abs(lat-(-0.))).argmin(), :] = emiss_rate
    # 10-20 deg N/S
    trac50_10_20 = np.zeros(emiss.shape)
    trac50_10_20[:, (np.abs(lat-10.)).argmin():
        (np.abs(lat-20.)).argmin(), :] = emiss_rate
    trac50_10_20[:, (np.abs(lat-(-20.))).argmin():
        (np.abs(lat-(-10.))).argmin(), :] = emiss_rate
    # 20-30 deg N/S
    trac50_20_30 = np.zeros(emiss.shape)
    trac50_20_30[:, (np.abs(lat-20.)).argmin():
        (np.abs(lat-30.)).argmin(), :] = emiss_rate
    trac50_20_30[:, (np.abs(lat-(-30.))).argmin():
        (np.abs(lat-(-20.))).argmin(), :] = emiss_rate
    # 30-40 deg N/S
    trac50_30_40 = np.zeros(emiss.shape)
    trac50_30_40[:, (np.abs(lat-30.)).argmin():
        (np.abs(lat-40.)).argmin(), :] = emiss_rate
    trac50_30_40[:, (np.abs(lat-(-40.))).argmin():
        (np.abs(lat-(-30.))).argmin(), :] = emiss_rate
    # 40-50 deg N/S    
    trac50_40_50 = np.zeros(emiss.shape)
    trac50_40_50[:, (np.abs(lat-40.)).argmin():
        (np.abs(lat-50.)).argmin(), :] = emiss_rate
    trac50_40_50[:, (np.abs(lat-(-50.))).argmin():
        (np.abs(lat-(-40.))).argmin(), :] = emiss_rate
    # 50-60 deg N/S    
    trac50_50_60 = np.zeros(emiss.shape)
    trac50_50_60[:, (np.abs(lat-50.)).argmin():
        (np.abs(lat-60.)).argmin(), :] = emiss_rate
    trac50_50_60[:, (np.abs(lat-(-60.))).argmin():
        (np.abs(lat-(-50.))).argmin(), :] = emiss_rate    
    # 60-70 deg N/S    
    trac50_60_70 = np.zeros(emiss.shape)
    trac50_60_70[:, (np.abs(lat-60.)).argmin():
        (np.abs(lat-70.)).argmin(), :] = emiss_rate
    trac50_60_70[:, (np.abs(lat-(-70.))).argmin():
        (np.abs(lat-(-60.))).argmin(), :] = emiss_rate   
    # 70-80 deg N/S    
    trac50_70_80 = np.zeros(emiss.shape)
    trac50_70_80[:, (np.abs(lat-70.)).argmin():
        (np.abs(lat-80.)).argmin(), :] = emiss_rate
    trac50_70_80[:, (np.abs(lat-(-80.))).argmin():
        (np.abs(lat-(-70.))).argmin(), :] = emiss_rate       
    # 80-90 deg N/S    
    trac50_80_90 = np.zeros(emiss.shape)
    trac50_80_90[:, (np.abs(lat-80.)).argmin():
        (np.abs(lat-90.)).argmin()+1, :] = emiss_rate
    trac50_80_90[:, (np.abs(lat-(-90.))).argmin():
        (np.abs(lat-(-80.))).argmin(), :] = emiss_rate       
    # Fill global region 
    globalr = np.zeros(emiss.shape)
    globalr[:] = emiss_rate
    # History for new tracer emission inventory
    date = datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    history = ('%s: created by Gaige Hunter Kerr (gaige.kerr@jhu.edu). '
        'Emissions are for three passive species with yearly, time-invariant '
        'emissions fields. All emission rates are constant'
        ', %s kg/m2/s. File contains these emissions for ' 
        '10 degree latitude bands (i.e., 0-10 deg N/S, 10-20 deg N/S) and '
        'a tracer with global emissions (GLOBAL_50).') %(date, str(emiss_rate))
    # Write output file 
    inventory_new = xr.Dataset(
        {'TRAC50_0_10': (['time', 'lat', 'lon'],  trac50_0_10),
         'TRAC50_10_20': (['time', 'lat', 'lon'],  trac50_10_20),
         'TRAC50_20_30': (['time', 'lat', 'lon'],  trac50_20_30),
         'TRAC50_30_40': (['time', 'lat', 'lon'],  trac50_30_40),
         'TRAC50_40_50': (['time', 'lat', 'lon'],  trac50_40_50),
         'TRAC50_50_60': (['time', 'lat', 'lon'],  trac50_50_60),
         'TRAC50_60_70': (['time', 'lat', 'lon'],  trac50_60_70),     
         'TRAC50_70_80': (['time', 'lat', 'lon'],  trac50_70_80),
         'TRAC50_80_90': (['time', 'lat', 'lon'],  trac50_80_90),     
         'GLOBAL_50': (['time', 'lat', 'lon'],  globalr)},        
        coords={'lat': (['lat'], lat),
                'lon': (['lon'], lng),            
                'time': (['time'], time)},
        attrs={'CDI': inventory.CDI,
               'Conventions': inventory.Conventions,
               'history': history,
               'Title': inventory.Title,
               'Model': inventory.Model,
               'Delta_Lon': inventory.Delta_Lon,
               'Delta_Lat': inventory.Delta_Lat,
               'NLayers': inventory.NLayers,
               'Start_Date': inventory.Start_Date,
               'Start_Time': inventory.Start_Time,
               'End_Date': inventory.End_Date,
               'End_Time': inventory.End_Time,
               'Delta_Time': inventory.Delta_Time,
               'CDO': inventory.CDO})
    # Fill in missing attribute data for variables 
    # For lon
    inventory_new.lon.attrs['long_name'] = inventory.lon.long_name
    inventory_new.lon.attrs['units'] = inventory.lon.units
    # For lat
    inventory_new.lat.attrs['standard_name'] = inventory.lat.standard_name
    inventory_new.lat.attrs['long_name'] = inventory.lat.long_name
    inventory_new.lat.attrs['units'] = inventory.lat.units
    inventory_new.lat.attrs['axis'] = inventory.lat.axis
    # For time
    inventory_new.time.attrs['standard_name'] = inventory.time.standard_name
    # For tracers
    inventory_new.TRAC50_0_10.attrs['long_name'] = 'TRAC50_0_10'
    inventory_new.TRAC50_10_20.attrs['long_name'] = 'TRAC50_10_20'
    inventory_new.TRAC50_20_30.attrs['long_name'] = 'TRAC50_20_30'
    inventory_new.TRAC50_30_40.attrs['long_name'] = 'TRAC50_30_40'
    inventory_new.TRAC50_40_50.attrs['long_name'] = 'TRAC50_40_50'
    inventory_new.TRAC50_50_60.attrs['long_name'] = 'TRAC50_50_60'
    inventory_new.TRAC50_60_70.attrs['long_name'] = 'TRAC50_60_70'
    inventory_new.TRAC50_70_80.attrs['long_name'] = 'TRAC50_70_80'
    inventory_new.TRAC50_80_90.attrs['long_name'] = 'TRAC50_80_90'
    inventory_new.GLOBAL_50.attrs['long_name'] = 'GLOBAL_50'
    inventory_new.TRAC50_0_10.attrs['units'] = 'kg/m2/s'
    inventory_new.TRAC50_10_20.attrs['units'] = 'kg/m2/s'
    inventory_new.TRAC50_20_30.attrs['units'] = 'kg/m2/s'
    inventory_new.TRAC50_30_40.attrs['units'] = 'kg/m2/s'
    inventory_new.TRAC50_40_50.attrs['units'] = 'kg/m2/s'
    inventory_new.TRAC50_50_60.attrs['units'] = 'kg/m2/s'
    inventory_new.TRAC50_60_70.attrs['units'] = 'kg/m2/s'
    inventory_new.TRAC50_70_80.attrs['units'] = 'kg/m2/s'
    inventory_new.TRAC50_80_90.attrs['units'] = 'kg/m2/s'
    inventory_new.GLOBAL_50.attrs['units'] = 'kg/m2/s'
    inventory_new.to_netcdf(INVENTORYPATH+
        'tracers_globalsource_ghk.annual.geos.2x25.nc')
    return

def adapt_co_50_forgeoschem(): 
    """adapt CO_50 emissions inventory from CCMI for use in GEOSChem by 
    adding the correct file and variable attributes following existing GEOSChem
    emissions inventories.

    Parameters
    ----------
    None
    
    Returns
    -------
    None    
    
    """
    import datetime
    import xarray as xr
    import numpy as np
    INVENTORYPATH = '/Users/ghkerr/phd/globalo3/GEOSChem/'
    inventoryname = 'RETRO_ANTHRO_V2_2000_CO_aggregated.nc'
    # Open existing inventory
    inventory = xr.open_dataset(INVENTORYPATH+inventoryname)
    # Extract coordinate data
    lng = inventory.lon.data
    lat = inventory.lat.data
    time = inventory.time.data
    # Extract emissions data for CO_50
    emiss = inventory.emission_flux.data
    # Attributes for GEOSChem-compatible CO_50 emission inventory
    date = datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
    history = ('%s: created by Gaige Hunter Kerr (gaige.kerr@jhu.edu). '
        'This inventory is an adaptation of RETRO emissions for '
        'species CO, year 2000 from Martin Schultz at MPI-H. '
        'This file is compatiable with HEMCO in GEOSChem. '%(date))
    history = history+inventory.history
    Title = inventory.title
    # The following is kludgey and copy and pasted from the GEOSChem 
    # fossilCO2.annual.geos.2x25.nc emissions inventory
    Model = 'GEOS4'
    Delta_Lon = np.float32(np.diff(lng[1:-1]).mean())
    Delta_Lat = np.float32(np.diff(lat[1:-1]).mean())
    NLayers = np.int32(55)
    Start_Date = np.int32(20000101)
    Start_Time = np.int32(0)
    End_Date = np.int32(20000101)
    End_Time = np.int32(0)
    Delta_Time = np.int32(0)
    CDO = 'Climate Data Operators version 1.5.5 (http://code.zmaw.de/'+\
        'projects/cdo)'
    CDI = 'Climate Data Interface version 1.5.5 (http://code.zmaw.de/'+\
        'projects/cdi)'
    Conventions = 'COARDS'
    # Write output file 
    inventory_new = xr.Dataset(
        {'CO_50': (['time', 'lat', 'lon'],  emiss)},        
        coords={'lat': (['lat'], lat),
                'lon': (['lon'], lng),            
                'time': (['time'], time)},
        attrs={'CDI': CDI,
               'Conventions': Conventions, 
               'history': history,
               'Title': Title,
               'Model': Model, 
               'Delta_Lon': Delta_Lon, 
               'Delta_Lat': Delta_Lat, 
               'NLayers': NLayers, 
               'Start_Date': Start_Date, 
               'Start_Time': Start_Time,
               'End_Date': End_Date,
               'End_Time': End_Time,
               'Delta_Time': Delta_Time,
               'CDO': CDO})
    # Convert longitude coordinates from 0-359 to -180-179:
    inventory_new.assign_coords(lon=(((inventory_new.lon+180)%360)-180))
    # Fill in missing attribute data for variables 
    # For lon
    inventory_new.lon.attrs['standard_name'] = 'longitude'
    inventory_new.lon.attrs['long_name'] = 'Longitude'
    inventory_new.lon.attrs['units'] = 'degrees_east'
    # For lat
    inventory_new.lat.attrs['standard_name'] = 'latitude'
    inventory_new.lat.attrs['long_name'] = 'Latitude'
    inventory_new.lat.attrs['units'] = 'degrees_north'
    inventory_new.lat.attrs['axis'] = 'Y'
    # For time
    inventory_new.time.attrs['standard_name'] = 'time'
    # For tracers
    inventory_new.CO_50.attrs['long_name'] = 'CO_50'
    inventory_new.CO_50.attrs['units'] = 'kg/m2/s'
    inventory_new.to_netcdf(INVENTORYPATH+
        'RETRO_ANTHRO_V2_2000_CO_aggregated.geos.nc')
    return
