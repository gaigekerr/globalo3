#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create GEOSChem emission inventory for passive species

This module uses an existing emission inventory for GEOSChem CO2 emissions
and mimicks this file structure and its attributes to create a time-invarient 
emission inventory for passive tracers

Notes
-----
    There are three data variables in the output netCDF files: POLAR_50 for 
    a passive species with a 50 d lifetime with annual, time-invariant 
    emissions from 60-90˚N/S, MIDLAT_50 for emissions from 30-60˚N/S, and 
    TROPIC_50 from 0-30˚N/S. One could, in the future create separate variables
    for each hemisphere to investigate exchanges from each hemisphere.
    
Revision History
----------------
    20092019 -- initial version created
    01102019 -- added tracer (GLOBAL_50) with global emissions
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
emiss_rate = 5e-8
# Find locations of latitudinal bands (in both hemispheres) for each 
# tracer
# 90              /------\
# 75            /  POLAR   \
# 60          /TRACER (60-90)\
# 45        /  MID-LAT TRACER  \
# 30      /       (30-60)        \
# 15    /  TROPICS TRACER (0-30)   \
#  0  /______________________________\  
# Fill polar region
polar = np.zeros(emiss.shape)
polar[:, (np.abs(lat-60.)).argmin():
    (np.abs(lat-90.)).argmin(), :] = emiss_rate
polar[:, (np.abs(lat-(-90.))).argmin():
    (np.abs(lat-(-60.))).argmin(), :] = emiss_rate
# Fill mid-latitude regions 
midlat = np.zeros(emiss.shape)
midlat[:, (np.abs(lat-30.)).argmin():
    (np.abs(lat-60.)).argmin(), :] = emiss_rate
midlat[:, (np.abs(lat-(-60.))).argmin():
    (np.abs(lat-(-30.))).argmin(), :] = emiss_rate    
# Fill tropics regions 
tropics = np.zeros(emiss.shape)
tropics[:, (np.abs(lat-0.)).argmin():
    (np.abs(lat-30.)).argmin(), :] = emiss_rate
tropics[:, (np.abs(lat-(-30.))).argmin():
    (np.abs(lat-(-0.))).argmin(), :] = emiss_rate    
# Fill global region 
globalr = np.zeros(emiss.shape)
globalr[:] = emiss_rate
# History for new tracer emission inventory
date = datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
history = ('%s: created by Gaige Hunter Kerr (gaige.kerr@jhu.edu). '
    'Emissions are for three passive species with yearly, time-invariant '
    'emissions fields. All emission rates are constant'
    ', %s kg/m2/s. File contains these emissions for ' 
    'three different regions: 0-30˚N/S (data variable TROPIC_50),'
    '30-60˚N/S (MIDLAT_50), and 60-90˚N/S (POLAR_50) as well as a tracer '
    'with global emissions (GLOBAL_50).') %(
    date, str(emiss_rate))
Title = 'GEOSChem passive species emissions'
Model = inventory.Model
Delta_Lon = inventory.Delta_Lon
Delta_Lat = inventory.Delta_Lat
NLayers = inventory.NLayers
Start_Date = inventory.Start_Date
Start_Time = inventory.Start_Time
End_Date = inventory.End_Date
End_Time = inventory.End_Time
Delta_Time = inventory.Delta_Time
CDO = inventory.CDO
# Write output file 
inventory_new = xr.Dataset(
    {'TROPIC_50': (['time', 'lat', 'lon'],  tropics),
     'MIDLAT_50': (['time', 'lat', 'lon'],  midlat),
     'POLAR_50': (['time', 'lat', 'lon'],  polar),
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
inventory_new.TROPIC_50.attrs['long_name'] = 'TROPIC_50'
inventory_new.MIDLAT_50.attrs['long_name'] = 'MIDLAT_50'
inventory_new.POLAR_50.attrs['long_name'] = 'POLAR_50'
inventory_new.GLOBAL_50.attrs['long_name'] = 'POLAR_50'
inventory_new.TROPIC_50.attrs['units'] = 'kg/m2/s'
inventory_new.MIDLAT_50.attrs['units'] = 'kg/m2/s'
inventory_new.POLAR_50.attrs['units'] = 'kg/m2/s'
inventory_new.GLOBAL_50.attrs['units'] = 'kg/m2/s'
inventory_new.to_netcdf(INVENTORYPATH+
    'tracers_globalsource_ghk.annual.geos.2x25.nc')