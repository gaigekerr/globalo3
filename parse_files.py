#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Open MERRA-2 and GMI CTM files and parse out needed information into 
smaller files just containing JJA 2008-2010 fields (as well as JJA 2016-2017 


fields to compare with Chinese data) for the Northern Hemisphere.
If the "full" files needed to create these compact files are no longer stored
locally to run this, consult "kerr_surface_2020_README.pdf" to see where the 
full files are remotely stored and where they should be locally to run this 
script (n.b., some of the code to open full files is rather kludgey with 
regard to the path names and directories).

Revision History
    03032020 -- initial version created
"""
import numpy as np
import xarray as xr
from datetime import date
today = date.today().strftime('%d/%m/%Y')
import sys
sys.path.append('/Users/ghkerr/phd/globalo3/')
import globalo3_open, globalo3_calculate, observations_open
sys.path.append('/Users/ghkerr/phd/transporto3/')
import transporto3_open

# Directory where parsed files will be created
path_parsed_files = '/Users/ghkerr/phd/globalo3/data/parsed/'

# # # # Load data 
# GMI CTM control O3 
lat_gmi, lng_gmi, times_gmi, o3_gmi = \
    globalo3_open.open_overpass2_specifieddomain([2008, 2009, 2010], 
    ['jun', 'jul', 'aug'], -1., 90., 0., 360., 'O3', 'HindcastMR2')
o3_gmi = o3_gmi*1e9
lat_gmi_china, lng_gmi_china, times_gmi_china, o3_gmi_china = \
    globalo3_open.open_overpass2_specifieddomain([2016, 2017], 
    ['jun', 'jul', 'aug'], -1., 90., 0., 360., 'O3', 'HindcastMR2')
o3_gmi_china = o3_gmi_china*1e9

# GMI CTM transport-only O3
lat_gmi, lng_gmi, times_gmi, o3_transport_gmi = \
    globalo3_open.open_overpass2_specifieddomain([2008, 2009, 2010], 
    ['jun', 'jul', 'aug'], -1., 90., 0., 360., 'O3', 'HindcastMR2-DiurnalAvgTQ')
o3_transport_gmi = o3_transport_gmi*1e9 

# MERRA-2 2-meter temperature 
lat_merra, lng_merra, t2m_merra = \
    globalo3_open.open_merra2t2m_specifieddomain([2008, 2009, 2010], 
    ['jun', 'jul', 'aug'], -1., 90., 0., 360.)
t2m_merra = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
    lng_gmi, lat_merra, lng_merra, t2m_merra)
lat_merra_china, lng_merra_china, t2m_merra_china = \
    globalo3_open.open_merra2t2m_specifieddomain([2016, 2017], 
    ['jun', 'jul', 'aug'], -1., 90., 0., 360.)        
t2m_merra_china = globalo3_open.interpolate_merra_to_ctmresolution(
    lat_gmi, lng_gmi, lat_merra_china, lng_merra_china, 
    t2m_merra_china)

# MERRA-2 2-meter specific humidity
qv2m_merra, mtime, lat_merra, lng_merra = transporto3_open.open_merra2(
    [2008, 2009, 2010], [0, 3, 9, 12, 15, 18, 21], 'QV2M', 'tavg1_2d_slv_Nx', 
    'JJA_rh.nc', 0., 90., 360., -1., dailyavg='yes')
qv2m_merra = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
    lng_gmi, lat_merra, lng_merra, qv2m_merra)   
qv2m_merra = qv2m_merra*1000. # Convert from kg kg-1 to g kg-1
qv2m_merra_china, mtime_china, lat_merra_china, lng_merra_china = \
    transporto3_open.open_merra2([2016, 2017], 
    [0, 3, 9, 12, 15, 18, 21], 'QV2M', 'tavg1_2d_slv_Nx', 'JJA_rh.nc', 
    0., 90., 360., -1., dailyavg='yes')
qv2m_merra_china = globalo3_open.interpolate_merra_to_ctmresolution(
    lat_gmi, lng_gmi, lat_merra_china, lng_merra_china, 
    qv2m_merra_china)
qv2m_merra_china = qv2m_merra_china*1000. 

# MERRA-2 10-m northward (V) wind
u10m, mtime, lat_merra, lng_merra = transporto3_open.open_merra2(
    [2008, 2009, 2010], np.arange(0, 24, 1), 'U10M', 'tavg1_2d_slv_Nx', 
    'JJA_wind.nc', 0., 90., 360., -1., dailyavg='yes')
u10m = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, lng_gmi, 
    lat_merra, lng_merra, u10m)

# MERRA-2 10-m northward (V) wind
v10m, mtime, lat_merra, lng_merra = transporto3_open.open_merra2(
    [2008, 2009, 2010], np.arange(0, 24, 1), 'V10M', 'tavg1_2d_slv_Nx', 
    'JJA_wind.nc', 0., 90., 360., -1., dailyavg='yes')
v10m = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, lng_gmi, 
    lat_merra, lng_merra, v10m)

# MERRA-2 PBL height
mtime, lat_merra, lng_merra, pblh_merra = \
    globalo3_open.open_merra2pblh_specifieddomain([2008, 2009, 2010], 
    np.arange(0, 24, 1), -1., 90., 0., 360., dailyavg='yes')
pblh_merra = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
    lng_gmi, lat_merra, lng_merra, pblh_merra, checkplot='yes')

# Jet stream latitude, derived from MERRA-2 500-hPa eastward (U) wind
U500, mtime, lat_merra, lng_merra = \
    globalo3_open.open_merra2_specifieddomain([2008, 2009, 2010], 
    ['jun', 'jul', 'aug'], [0, 3, 6, 9, 12, 15, 18, 21], 'U', 
    'inst3_3d_asm_Np_500hPa', 0., 90., 360., -1., dailyavg='yes')
U500 = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
    lng_gmi, lat_merra, lng_merra, U500, checkplot='yes')
# Subset fields in mid-latitudes
U500_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(U500, lat_gmi, 
    lng_gmi, 0., 360., 20., 70.)
# Find latitude of eddy-driven jet and fields at jet
lat_jet_ml, scratch = globalo3_calculate.find_field_atjet(U500_ml, U500_ml, 
    lat_ml, lng_ml, 20, anom = True)
# Smooth jet latitude with window size of ~10 deg
lat_jet_ml = globalo3_calculate.convolve_jet(lat_jet_ml, 9)    
lat_grid = np.tile(lat_gmi, (lng_gmi.shape[0], 1)).T
# Grid will be populated with the distance from the eddy-driven jet (in 
# degrees) at each lat, lon coordinate
edj_dist = np.empty(shape=(len(lat_jet_ml),lat_grid.shape[0],lat_grid.shape[1]))
edj_dist[:] = np.nan
# Loop through days in measuring period; and, at each day, subtract the 
# latitude of the jet from the local latitude  
for day in np.arange(0,len(lat_jet_ml),1):
    edj_dist[day] = lat_jet_ml[day]-lat_grid
    
# Cyclones
cyclones = globalo3_open.open_merra2_cyclones('06/01/2008', 
    '08/31/2010', ['jun', 'jul', 'aug'])
# Bin cyclones (~4.125˚ lat x ~4.045˚ lng) and determine frequency
cyclones_binned, lat_cyclones_binned, lng_cyclones_binned = \
    globalo3_calculate.field_binner(np.linspace(lat_gmi[0], 
    lat_gmi[-1], 23), np.linspace(lng_gmi[0], lng_gmi[-1], 90), 
    times_gmi, cyclones['Latitude'].values, 
    cyclones['Longitude'].values, 
    np.ones(shape=cyclones['Latitude'].values.shape), 
    cyclones['Date'].values, 'sum')
cyclones_binned = np.nan_to_num(cyclones_binned)
# Same as above but temperature on days with equatorward and poleward 
# jet 
(eqjet_lat_cyclones, eqjet_lng_cyclones, eqjet_time_cyclones,
    pwjet_lat_cyclones, pwjet_lng_cyclones, pwjet_time_cyclones, 
    eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_o3, 
    eqjet_o3) = globalo3_calculate.segregate_field_bylat(o3_gmi, lng_gmi, 
    lat_jet_ml, times_gmi, cyclones=cyclones)
# Same as above but humidity on days with equatorward and poleward 
# jet              
(eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_qv2m, 
    eqjet_qv2m) = globalo3_calculate.segregate_field_bylat(qv2m_merra, lng_gmi, 
    lat_jet_ml, times_gmi)             
# Same as above but find O3 on days with equatorward and poleward jet 
(eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_t2m, 
    eqjet_t2m) = globalo3_calculate.segregate_field_bylat(t2m_merra, lng_gmi, 
    lat_jet_ml, times_gmi)   
# Bin cyclones (~4.125˚ lat x ~4.045˚ lng) on days with equatorward jet
# and determine frequency
eqjet_cyclones_binned, eqjet_lat_cyclones_binned, eqjet_lng_cyclones_binned = \
    globalo3_calculate.field_binner(np.linspace(lat_gmi[0], 
    lat_gmi[-1], 23), np.linspace(lng_gmi[0], lng_gmi[-1], 90), 
    times_gmi, eqjet_lat_cyclones, eqjet_lng_cyclones, 
    np.ones(shape=eqjet_lat_cyclones.shape[0]), eqjet_time_cyclones,
    'sum')
eqjet_cyclones_binned = np.nan_to_num(eqjet_cyclones_binned)  
# Same as above but for days with poleward jet
pwjet_cyclones_binned, pwjet_lat_cyclones_binned, pwjet_lng_cyclones_binned = \
    globalo3_calculate.field_binner(np.linspace(lat_gmi[0], 
    lat_gmi[-1], 23), np.linspace(lng_gmi[0], lng_gmi[-1], 90), 
    times_gmi, pwjet_lat_cyclones, pwjet_lng_cyclones, 
    np.ones(shape=pwjet_lat_cyclones.shape[0]), pwjet_time_cyclones, 
    'sum')
pwjet_cyclones_binned = np.nan_to_num(pwjet_cyclones_binned)
# Add cyclic point to longitude coordinates so that it wraps around the 
# Prime Meridian when plotting
lng_cyclones_binned[0] = 0.
lng_cyclones_binned[-1] = 360.
eqjet_lng_cyclones_binned[0] = 0.
eqjet_lng_cyclones_binned[-1] = 360.
pwjet_lng_cyclones_binned[0] = 0.
pwjet_lng_cyclones_binned[-1] = 360. 

# # # # Write files
# GMI CTM ontrol O3 
ds = xr.Dataset({'O3_control': (('time', 'lat', 'lng'), o3_gmi)},
    coords={'time' : mtime, 'lat': lat_gmi, 'lng': lng_gmi})
ds.attrs['title'] = 'GMI CTM surface level O3'
ds.attrs['history'] = 'O3 from the control simulation of the GMI CTM '+\
    'for JJA 2008-2010 for the Northern Hemisphere. Model output is from '+\
    'the lowest model level (eta mid is 992 hPa; eta top is 985 hPa). '+\
    'This file is derived from monthly files in '+\
    'kerr_surface_2020/data/full/o3/YYYY/'
ds.attrs['date'] =  today
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.O3_control.attrs['units'] = 'ppbv'
ds.to_netcdf(path_parsed_files+'gmi_O3control_JJA2008-2010.nc')
del ds

ds = xr.Dataset({'O3_control': (('time', 'lat', 'lng'), o3_gmi_china)},
    coords={'time' : mtime_china, 'lat': lat_gmi, 'lng': lng_gmi})
ds.attrs['title'] = 'GMI CTM surface level O3'
ds.attrs['history'] = 'O3 from the control simulation of the GMI CTM '+\
    'for JJA 2016-2017 for the Northern Hemisphere. Model output is from '+\
    'the lowest model level (eta mid is 992 hPa; eta top is 985 hPa). '+\
    'This file is derived from monthly files in '+\
    'kerr_surface_2020/data/full/o3/YYYY/'
ds.attrs['date'] =  today
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.O3_control.attrs['units'] = 'ppbv' 
ds.to_netcdf(path_parsed_files+'gmi_O3control_JJA2016-2017.nc')
del ds

# GMI CTM transport-only O3
ds = xr.Dataset({'O3_transportonly': (('time', 'lat', 'lng'), 
    o3_transport_gmi)}, coords={'time' : mtime, 'lat': lat_gmi, 
    'lng': lng_gmi})
ds.attrs['title'] = 'GMI CTM transport-only surface level O3'
ds.attrs['history'] = 'O3 from the transport-only simulation of the GMI CTM '+\
    'for JJA 2008-2010 for the Northern Hemisphere. Model output is from '+\
    'the lowest model level (eta mid is 992 hPa; eta top is 985 hPa). '+\
    'This file is derived from monthly files in '+\
    'kerr_surface_2020/data/full/o3/YYYY/'
ds.attrs['date'] =  today
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.O3_transportonly.attrs['units'] = 'ppbv'
ds.to_netcdf(path_parsed_files+'gmi_O3transportonly_JJA2008-2010.nc')
del ds

# MERRA-2 2-meter temperature 
ds = xr.Dataset({'T2M': (('time', 'lat', 'lng'), 
    t2m_merra)}, coords={'time' : mtime, 'lat': lat_gmi, 
    'lng': lng_gmi})
ds.attrs['title'] = 'MERRA-2 daily maximum temperature'
ds.attrs['history'] = 'Daily maximum MERRA-2 2-m temperatures derived '+\
    'from hourly data. Temperatures are for JJA 2008-2010 for the Northern '+\
    'Hemisphere and are interpolated to the resolution of the GMI CTM. '+\
    'This file is derived from monthly files in '+\
    'kerr_surface_2020/data/full/t2m/YYYY/'
ds.attrs['date'] =  today
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.T2M.attrs['units'] = 'K'
ds.to_netcdf(path_parsed_files+'merra2_T2M_JJA2008-2010.nc')
del ds

ds = xr.Dataset({'T2M': (('time', 'lat', 'lng'), 
    t2m_merra_china)}, coords={'time' : mtime_china, 'lat': lat_gmi, 
    'lng': lng_gmi})
ds.attrs['title'] = 'MERRA-2 daily maximum temperature'
ds.attrs['history'] = 'Daily maximum MERRA-2 2-m temperatures derived '+\
    'from hourly data. Temperatures are for JJA 2016-2017 for the Northern '+\
    'Hemisphere and are interpolated to the resolution of the GMI CTM. '+\
    'This file is derived from monthly files in '+\
    'kerr_surface_2020/data/full/t2m/YYYY/'
ds.attrs['date'] =  today
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.T2M.attrs['units'] = 'K'
ds.to_netcdf(path_parsed_files+'merra2_T2M_JJA2016-2017.nc')
del ds

# MERRA-2 2-meter specific humidity
ds = xr.Dataset({'QV2M': (('time', 'lat', 'lng'), 
    qv2m_merra)}, coords={'time' : mtime, 'lat': lat_gmi, 
    'lng': lng_gmi})
ds.attrs['title'] = 'MERRA-2 daily mean specific humidity'
ds.attrs['history'] = 'Daily mean MERRA-2 2-m specific humidity derived '+\
    'from three-hourly data. Specific humidity is for JJA 2008-2010 for '+\
    'the Northern Hemisphere and is interpolated to the resolution of '+\
    'the GMI CTM. This file is derived from monthly files in '+\
    'kerr_surface_2020/data/full/qv2m/YYYY/'
ds.attrs['date'] =  today
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.QV2M.attrs['units'] = 'g kg-1'
ds.to_netcdf(path_parsed_files+'merra2_QV2M_JJA2008-2010.nc')
del ds

ds = xr.Dataset({'QV2M': (('time', 'lat', 'lng'), 
    qv2m_merra_china)}, coords={'time' : mtime_china, 'lat': lat_gmi, 
    'lng': lng_gmi})
ds.attrs['title'] = 'MERRA-2 daily mean specific humidity'
ds.attrs['history'] = 'Daily mean MERRA-2 2-m specific humidity derived '+\
    'from three-hourly data. Specific humidity is for JJA 2016-2017 for '+\
    'the Northern Hemisphere and is interpolated to the resolution of '+\
    'the GMI CTM. This file is derived from monthly files in '+\
    'kerr_surface_2020/data/full/qv2m/YYYY/'
ds.attrs['date'] =  today
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.QV2M.attrs['units'] = 'g kg-1'
ds.to_netcdf(path_parsed_files+'merra2_QV2M_JJA2016-2017.nc')
del ds

# MERRA-2 10-m eastward (U) wind
ds = xr.Dataset({'U10M': (('time', 'lat', 'lng'), 
    u10m)}, coords={'time' : mtime, 'lat': lat_gmi, 
    'lng': lng_gmi})
ds.attrs['title'] = 'MERRA-2 daily mean northward wind'
ds.attrs['history'] = 'Daily mean MERRA-2 10-m northward wind derived '+\
    'from hourly data. Wind is for JJA 2008-2010 for the Northern '+\
    'Hemisphere and is interpolated to the resolution of the GMI CTM. '+\
    'This file is derived from monthly files in '+\
    'kerr_surface_2020/data/full/wind10m/YYYY/'
ds.attrs['date'] =  today
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.U10M.attrs['units'] = 'm s-1'
ds.to_netcdf(path_parsed_files+'merra2_U10M_JJA2008-2010.nc')
del ds

# MERRA-2 10-m northward (V) wind
ds = xr.Dataset({'V10M': (('time', 'lat', 'lng'), 
    v10m)}, coords={'time' : mtime, 'lat': lat_gmi, 
    'lng': lng_gmi})
ds.attrs['title'] = 'MERRA-2 daily mean eastward wind'
ds.attrs['history'] = 'Daily mean MERRA-2 10-m eastward wind derived '+\
    'from hourly data. Wind is for JJA 2008-2010 for the Northern '+\
    'Hemisphere and is interpolated to the resolution of the GMI CTM. '+\
    'This file is derived from monthly files in '+\
    'kerr_surface_2020/data/full/wind10m/YYYY/'
ds.attrs['date'] =  today
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.V10M.attrs['units'] = 'm s-1'
ds.to_netcdf(path_parsed_files+'merra2_V10M_JJA2008-2010.nc')
del ds

# MERRA-2 PBL height
ds = xr.Dataset({'PBLH': (('time', 'lat', 'lng'), 
    pblh_merra)}, coords={'time' : mtime, 'lat': lat_gmi, 
    'lng': lng_gmi})
ds.attrs['title'] = 'MERRA-2 daily mean planetary boundary layer height'
ds.attrs['history'] = 'Daily mean MERRA-2 planetary boundary layer height '+\
    '(PBLH) from hourly data. PBLH is for JJA 2008-2010 for the Northern '+\
    'Hemisphere and is interpolated to the resolution of the GMI CTM. '+\
    'This file is derived from monthly files in '+\
    'kerr_surface_2020/data/full/pblh/YYYY/'
ds.attrs['date'] =  today
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.PBLH.attrs['units'] = 'm'
ds.to_netcdf(path_parsed_files+'merra2_PBLH_JJA2008-2010.nc')
del ds

# Jet stream latitude, derived from MERRA-2 500-hPa eastward (U) wind
ds = xr.Dataset({'jetdistance': (('time', 'lat', 'lng'), 
    edj_dist), 'jetlatitude' : (('time', 'lng'), lat_jet_ml)}, 
    coords={'time' : mtime, 'lat': lat_gmi, 'lng': lng_gmi})
ds.attrs['title'] = 'Jet stream'
ds.attrs['history'] = 'Jet stream latitude and distance from the jet '+\
    'stream. The distance from the jet stream is defined, in degrees, '+\
    'as the latitude of the jet stream minus the local latitude (i.e., '+\
    'the latitude of the grid cell), such that a '+\
    'positive distance means that the jet stream is north. The jet '+\
    'stream position is defined as the latitude of maximum daily mean '+\
    '(based on 3-hourly data) 500 hPa eastward wind interpolated to the '+\
    'resolution of the GMI CTM. The latitude of the jet stream is '+\
    'smoothed with simple moving average ("boxcar smoothing") and '+\
    'replaces each data value with the average value of the ~9 '+\
    'neighboring values. This file is derived from daily files in '+\
    'kerr_surface_2020/data/full/jetdistance/YYYY/'
ds.attrs['date'] =  today
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.jetdistance.attrs['units'] = 'degrees'
ds.jetlatitude.attrs['units'] = 'degrees'
ds.to_netcdf(path_parsed_files+'merra2_JET_JJA2008-2010.nc')
del ds

# Cyclone frequency
ds = xr.Dataset({'cyclones': (('lat', 'lng'), cyclones_binned), 
    'pwjet_cyclones': (('lat', 'lng'), pwjet_cyclones_binned), 
    'ewjet_cyclones': (('lat', 'lng'), eqjet_cyclones_binned)},
    coords={'lat': lat_cyclones_binned, 'lng': lng_cyclones_binned})
ds.attrs['title'] = 'Binned cyclone frequency'
ds.attrs['history'] = 'Cyclone frequency binned to a ~4.125 deg lat x '+\
    '~4.045 deg lng grid for all days in JJA 2008-2010 and days with '+\
    'poleward and equatorward jets. Cyclone locations come from the '+\
    'MAP Climatology for Midlatitude Storminess. This file is derived '+\
    'from files corresponding to individual cyclones and is located at '+\
    'kerr_surface_2020/data/full/cyclones/merra2frontsYYYY/'
ds.attrs['date'] =  today
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.cyclones.attrs['units'] = ''
ds.pwjet_cyclones.attrs['units'] = ''
ds.ewjet_cyclones.attrs['units'] = ''
ds.to_netcdf(path_parsed_files+'mcms_CYCLONEFREQ_JJA2008-2010.nc')
del ds

# Cyclone dataframe (save as a pickle)
import pickle
cyclones.to_pickle(path_parsed_files+'cyclones.pickle')