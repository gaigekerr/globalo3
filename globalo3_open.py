#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Open GMI CTM trace gas concentrations, MERRA-2 2-meter temperatures, 
observations, and emissions inventories across the global domain.

Revision History
    06122018 -- initial version created
    15012019 -- function 'open_edgar_specificdomian' added
    24012019 -- function 'open_schnello3' added; edited 
                function 'interpolate_merra_to_ctmresolution' to interpolate
                any gridded dataset to the resolution of the CTM
    30012019 -- function 'open_overpass2_specifieddomain' modified to handle 
                any named GMI simulation's overpass2 output
    07022019 -- function 'open_gmideposition_specifieddomain' added
    08022019 -- np.roll issue with longitude corrected in functions 
                'open_gmideposition_specifieddomain' and 
                'open_merra2t2m_specifieddomain'
"""

def open_overpass2_specifieddomain(years, months, latmin, latmax, lngmin,
    lngmax, gas, exper):
    """load daily overpass2 files for GMI CTM simulation "HindcastMR2" for 
    specified trace gas, time period, and spatial domain.
    
    Parameters
    ----------
    years : list
        Year or range of years in measuring period
    months : list
        Three letter abbreviations (lowercase) for months in measuring period
    latmin : float    
        Latitude (degrees north) of bottom edge of bounding box for focus 
        region. For this parameter and others defining the bounding box, 
        function finds the closest index to bounding box edges
    latmax : float
        Latitude (degrees north) of upper edge of bounding box for focus region
    lngmin : float
        Longitude (degrees east, 0-360) of left edge of bounding box for focus 
        region        
    lngmax : float
        Longitude (degrees east, 0-360) of right edge of bounding box for focus 
        region
    gas : str
        Chemical formula of desired trace gas; e.g., O3, OH, NO, NO2, CO, CH2O
    exper : str
        GMI experiment's simulation name; n.b. as of 30 Jan 2019, downloaded 
        simulations include 'HindcastMR2', 'HindcastFFIgac2', and 
        'Hindcast3Igac2'
        
    Returns
    -------
    lat : numpy.ndarray
        GMI latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        GMI longitude coordinates, units of degrees east, [lng,]
    times : numpy.ndarray
        datetime.date objects corresponding to every day in measuring period. 
        n.b., information about the hour is given although the daily arrays 
        correspond to 1300-1400 hours local time, [time,]
    gas_all : numpy.ndarray
        Daily gridded trace gas concentrations at 1300-1400 hours local time, 
        units of volume mixing ratio, [time,lat,lng]
    """
    import numpy as np
    from netCDF4 import Dataset
    import datetime, calendar
    import time
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    start_time = time.time()
    gas_all = []
    times = []
    # Loop through measuring period    
    for year in years:
        PATH_GMI = '/Users/ghkerr/phd/globalo3/data/GMI/%s/%s/'%(exper,year)        
        for month in months:
            # Open monthly overpass2 file 
            infile = Dataset(PATH_GMI+'gmic_%s_%s.%s_const.overpass2.nc'
                             %(exper,year,month),'r') 
            # On first iteration, extract dimensions and find indicies 
            # corresponding to (closest to) desired domain
            if (year == years[0]) and (month == months[0]):
                lat = infile.variables['latitude_dim'][:]
                lng = infile.variables['longitude_dim'][:]
                latmin = geo_idx(latmin, lat)
                latmax = geo_idx(latmax, lat)
                lngmin = geo_idx(lngmin, lng)
                lngmax = geo_idx(lngmax, lng)
                # Restrict coordinates over focus region 
                lat = lat[latmin:latmax+1]
                lng = lng[lngmin:lngmax+1]
            # Load chemical species names    
            species_raw = infile.variables['overpass_labels'][:]
            species = []
            for row in np.arange(0, np.shape(species_raw)[0]):
                temp = []
                for element in species_raw[row]:
                    temp.append(element.decode('UTF-8'))
                si = ''.join(temp)
                si = si.rstrip()
                species.append(si[:])      
            del species_raw
            # Find position of desired trace gas 
            gas_where = np.where(np.array(species) == gas)[0][0]
            # Extract trace gas concentrations for month 
            gas_month = infile.variables['const_overpass'][:,gas_where,0,
                                        latmin:latmax+1,lngmin:lngmax+1]
            # Drop leap days from analysis
            if (calendar.isleap(year)==True) and (month == 'feb'):
                gas_month = gas_month[:-1]
            # GMI simulation 'HindcastFFIgac2' only has 30 days in 
            # December 2008 (wtf?) unlike all over Hindcast simulations. 
            # For this simulation, repeat 30 December values for 31 December
            if (exper == 'HindcastFFIgac2') and (month == 'dec') and \
                (year == 2008):                    
                    gas_month = np.append(gas_month, gas_month[-1:],0)            
            # Append daily values to multimonth, multiyear list
            gas_all.append(gas_month)
            # Create time dataset for month; conversions from month 
            # abbreviation to integer from 
            # https://stackoverflow.com/questions/1549641/
            # how-to-capitalize-the-first-letter-of-each-word-in-a-string-python
            num_days = calendar.monthrange(year, 
                list(calendar.month_abbr).index(month.title()))[1]
            days = [datetime.date(year, 
                    list(calendar.month_abbr).index(month.title()), 
                    day) for day in range(1,num_days+1)]
            times.append(days)
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
          'GMI %s output for %d-%d loaded in %.2f seconds' %(gas,years[0],
          years[-1],(time.time()-start_time)))
    # Plotting in Cartopy has a bug that doesn't "wrap" around the Prime 
    # Meridian. To counter this for global/hemispheric output, I manually 
    # change the longitude entry for 358.75 deg to 360 deg.
    if (lng[-1] == 358.75) or (lng[-1] == 357.5):
        lng[-1] = 360.    
    return lat, lng, np.hstack(times), np.vstack(gas_all)

def open_gmideposition_specifieddomain(years, latmin, latmax, lngmin, lngmax, 
    gas, process):
    """load daily deposition or scavenging data from MERRA-2 GMI simulation
    
    Parameters
    ----------
    years : list
        Year or range of years in measuring period
    latmin : float    
        Latitude (degrees north) of bottom edge of bounding box for focus 
        region. For this parameter and others defining the bounding box, 
        function finds the closest index to bounding box edges
    latmax : float
        Latitude (degrees north) of upper edge of bounding box for focus region
    lngmin : float
        Longitude (degrees east, 0-360) of left edge of bounding box for focus 
        region        
    lngmax : float
        Longitude (degrees east, 0-360) of right edge of bounding box for focus 
        region
    gas : str
        Chemical formula of desired trace gas; e.g., O3, OH, NO, NO2, CO, CH2O
    process : str
        DD for dry deposition, WD for wet deposition, and SCAV for scavenging
        
    Returns
    -------
    lat : numpy.ndarray
        GMI latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        GMI longitude coordinates, units of degrees east, [lng,]
    dep : numpy.ndarray
        MERRA-2 GMI daily deposition data, for example, if process is DD
        (SCAV) and gas is O3, then the data are the dry deposition of Ox
        (savenging of Ox), units of kg m-3 s-1 (SCAV) or kg m-2 s-1 (DD/WD), 
        [years, days in year, lat, lng] (DD/WD) or [years, days in year, 
        level, lat, lng] (SCAV)
    """
    import time
    import numpy as np
    import xarray as xr
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx    
    start_time = time.time()
    PATH_DEPOSITION = '/Users/ghkerr/phd/globalo3/data/GMI/MERRA2_GMI/'
    dep = []
    for year in years:
        ds = xr.open_dataset(PATH_DEPOSITION+
                             'MERRA2_GMI.tavg24_3d_dep_Nv.%dJJA_%s.nc4'
                             %(year,gas))
        if (year==years[0]):
            lat = ds['lat'].data
            lng = ds['lon'].data
            # Convert longitude from (-180-180) to (0-360)
            lng = lng % 360       
            # Shift this grid such that it spans (0-360) rather than 
            # (180-360, 0-180)
            lng = np.roll(lng,int(lng.shape[0]/2)-1)
            latmin = geo_idx(latmin,lat)
            latmax = geo_idx(latmax,lat)
            lngmin = geo_idx(lngmin,lng)
            lngmax = geo_idx(lngmax,lng)
        # Extract desired data, roll deposition grid similar to longitude grid
        # Restrict data to focus region
        dep_yr = ds['%s_%s'%(process,gas)].data
        dep_yr = np.roll(dep_yr, int(dep_yr.shape[-1]/2)-1, axis=2)    
        dep_yr = dep_yr[:,latmin:latmax+1,lngmin:lngmax+1]
        lat = lat[latmin:latmax+1]
        lng = lng[lngmin:lngmax+1]       
        dep.append(dep_yr)
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
          'MERRA-2 GMI %s %s data for %d-%d loaded in %.2f seconds' 
          %(gas,process,years[0],years[-1],(time.time()-start_time)))
    return lat, lng, np.array(dep)

def open_merra2t2m_specifieddomain(years, months, latmin, latmax, lngmin, 
    lngmax):
    """load daily maximum 2-meter temperatures from MERRA-2 over the specific 
    spatial and temporal domains. 
    
    Parameters
    ----------
    years : list
        Year or range of years in measuring period
    months : list
        Three letter abbreviations (lowercase) for months in measuring period
    latmin : float    
        Latitude (degrees north) of bottom edge of bounding box for focus 
        region. For this parameter and others defining the bounding box, 
        function finds the closest index to bounding box edges
    latmax : float
        Latitude (degrees north) of upper edge of bounding box for focus region
    lngmin : float
        Longitude (degrees east, 0-360) of left edge of bounding box for focus 
        region        
    lngmax : float
        Longitude (degrees east, 0-360) of right edge of bounding box for focus 
        region
        
    Returns
    -------
    lat : numpy.ndarray
        MERRA-2 latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        MERRA-2 longitude coordinates, units of degrees east, [lng,]
    t2m_all : numpy.ndarray
        Daily maximum 2-meter temperatures, [time, lat, lng]
    """
    import numpy as np
    from netCDF4 import Dataset
    import calendar
    import time
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    start_time = time.time()
    months_int = []
    t2m_all = []
    # Convert month abbreviations to integers
    for month in months:
        months_int.append(list(calendar.month_abbr).index(month.title()))
    # Loop through measuring period    
    for year in years:
        PATH_MERRA = '/Users/ghkerr/phd/globalo3/data/MERRA-2/%d/'%year        
        for month in months_int:
            # Open monthly overpass2 file 
            infile = Dataset(PATH_MERRA+'MERRA2_300.inst1_2d_asm_Nx.%d%.2d.SUB.nc'
                             %(year,month),'r')        
            if (month==months_int[0]) and (year==years[0]):
                lat = infile.variables['lat'][:]
                lng = infile.variables['lon'][:]
                # Convert longitude from (-180-180) to (0-360)
                lng = lng % 360
                # Shift this grid such that it spans (0-360) rather than 
                # (180-360, 0-180)
                lng = np.roll(lng,int(lng.shape[0]/2)-1)
                latmin = geo_idx(latmin,lat)
                latmax = geo_idx(latmax,lat)
                lngmin = geo_idx(lngmin,lng)
                lngmax = geo_idx(lngmax,lng)
                # Restrict coordinates over focus region 
                lat = lat[latmin:latmax+1]
                lng = lng[lngmin:lngmax+1]           
            # Extract 2-meter temperatures for the month
            t2m_month = infile.variables['T2M'][:]
            # Drop leap days from analysis
            if (calendar.isleap(year)==True) and (month == 2):
                t2m_month = t2m_month[:-1]
            # Roll grid similar to longitude grid
            t2m_month = np.roll(t2m_month, int(t2m_month.shape[-1]/2)-1, axis = 2)
            t2m_month = t2m_month[:,latmin:latmax+1,lngmin:lngmax+1]
            t2m_all.append(t2m_month)
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
          'MERRA-2 for %d-%d loaded in %.2f seconds' %(years[0],years[-1],
          (time.time()-start_time)))
    return lat, lng, np.vstack(t2m_all)

def open_edgar_specifieddomain(years, latmin, latmax, lngmin, lngmax, 
    substance):
    """load annual mean emissions from EDGAR files for the specified 
    substance, time period, and spatial domain.
    
    Parameters
    ----------
    years : list
        Year or range of years in measuring period
    latmin : float    
        Latitude (degrees north) of bottom edge of bounding box for focus 
        region. For this parameter and others defining the bounding box, 
        function finds the closest index to bounding box edges
    latmax : float
        Latitude (degrees north) of upper edge of bounding box for focus region
    lngmin : float
        Longitude (degrees east, 0-360) of left edge of bounding box for focus 
        region        
    lngmax : float
        Longitude (degrees east, 0-360) of right edge of bounding box for focus 
        region
    substance : str
        Chemical formula of desired trace gas; e.g., NOx, CO
        
    Returns
    -------
    lat : numpy.ndarray
        EDGAR latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        EDGAR longitude coordinates, units of degrees east, [lng,]
    gas_all : numpy.ndarray
        Annual mean emissions of substance from EDGAR, units of kg substance
        m-2 s-1, [years, lat, lng]
    """
    import numpy as np
    from netCDF4 import Dataset
    import time
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    start_time = time.time()
    PATH_EMISS = '/Users/ghkerr/phd/globalo3/data/EDGAR/'
    edgar_all = []
    # Loop through measuring period    
    for year in years:
        # Open annual mean emissions of substance from EDGAR 
        infile = Dataset(PATH_EMISS+'v431_v2_REFERENCE_%s_%d.0.1x0.1.nc'
                         %(substance, year), 'r')        
        # On first iteration, extract dimensions and find indicies 
        # corresponding to (closest to) desired domain
        if (year == years[0]):
            lat = infile.variables['lat'][:]
            lng = infile.variables['lon'][:]
            latmin = geo_idx(latmin,lat)
            latmax = geo_idx(latmax,lat)
            lngmin = geo_idx(lngmin,lng)
            lngmax = geo_idx(lngmax,lng)
            # Restrict coordinates over focus region 
            lat = lat[latmin:latmax+1]
            lng = lng[lngmin:lngmax+1]
        # Extract constituent 
        edgar_year = infile.variables['emi_%s' %(substance.lower())][:]
        edgar_year = edgar_year[latmin:latmax+1,lngmin:lngmax+1]
        edgar_all.append(edgar_year)
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
          'EDGAR emissions of %s for %d-%d loaded in %.2f seconds' 
          %(substance, years[0], years[-1], (time.time()-start_time)))
    return lat, lng, np.stack(edgar_all)

def interpolate_merra_to_ctmresolution(lat_gmi, lng_gmi, lat_merra, lng_merra,
    t2m, checkplot='yes'): 
    """interpolate MERRA-2 field (0.5 deg latitude x 0.625 deg longitude) to 
    resolution of GMI CTM (1 deg latitude x 1.25 deg longitude) using xESMF. 
    
    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI longitude coordinates, units of degrees east, [lng,]
    lat_merra : numpy.ndarray
        MERRA-2 latitude coordinates, units of degrees north, [lat,]
    lng_merra : numpy.ndarray
        MERRA-2 longitude coordinates, units of degrees east, [lng,]        
    t2m : numpy.ndarray
        Daily maximum 2-meter temperatures, [time, lat, lng]
    checkplot : str
        If 'yes' the mean fields before/after interpolation are plotted
        
    Returns
    -------
    t2m_interp : numpy.ndarray 
        Interpolated daily maximum 2-meter temperatures, [time, lat, lng]
    """
    import numpy as np
    import xesmf as xe
    import matplotlib.pyplot as plt
    import time
    start_time = time.time()
    # Interpolate finer grid (MERRA-2) to the resolution of the CTM
    grid_in = {'lon':lng_merra,'lat':lat_merra}
    # Output grid has a coarser resolution
    grid_out = {'lon':lng_gmi,'lat':lat_gmi}
    # Use xESMF with numpy array
    regridder = xe.Regridder(grid_in,grid_out,'bilinear')
    regridder.clean_weight_file()
    t2m_interp = regridder(t2m)
    # Check to ensure that interpolation worked
    fig = plt.figure()
    ax1 = plt.subplot2grid((2,2),(0,0),colspan=2)
    ax2 = plt.subplot2grid((2,2),(1,0),colspan=2)
    clevs = np.linspace(t2m.min(), t2m.max(), 10)
    if checkplot == 'yes':
        a = t2m.mean(axis=tuple(range(0, t2m.ndim-2)))
        b = t2m_interp.mean(axis=tuple(range(0, t2m_interp.ndim-2)))        
        ax1.contourf(a, clevs)
        ax1.set_title('Original')
        m2 = ax2.contourf(b, clevs)
        ax2.set_title('Interpolated')
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(m2, cax=cbar_ax)
        plt.subplots_adjust(hspace=0.4)
        plt.show()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
          'MERRA-2 interpolated to CTM resolution in %.2f seconds' 
          %((time.time()-start_time)))
    return t2m_interp

def interpolate_edgar_to_ctmresolution(lat_gmi, lng_gmi, lat_edgar, lng_edgar,
    edgar): 
    """interpolate EDGAR field (0.1 deg latitude x 0.1 deg longitude) to 
    resolution of GMI CTM (1 deg latitude x 1.25 deg longitude) using xESMF. 
    
    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI longitude coordinates, units of degrees east, [lng,]
    lat_edgar : numpy.ndarray
        EDGAR latitude coordinates, units of degrees north, [lat,]
    lng_edgar : numpy.ndarray
        EDGAR longitude coordinates, units of degrees east, [lng,]        
    edgar : numpy.ndarray
        Annual mean emissions of substance from EDGAR, units of kg substance
        m-2 s-1, [years, lat, lng]
        
    Returns
    -------
    edgar_interp : numpy.ndarray 
        Interpolated annual mean emissions of substance from EDGAR, [years, 
        lat, lng]
    """
    import numpy as np
    import xesmf as xe
    import matplotlib.pyplot as plt
    import time
    start_time = time.time()
    # Interpolate EDGAR emissions to the resolution of the CTM
    grid_in = {'lon':lng_edgar, 'lat':lat_edgar}
    # Output grid has a coarser resolution
    grid_out = {'lon':lng_gmi, 'lat':lat_gmi}
    # Use xESMF with numpy array
    regridder = xe.Regridder(grid_in,grid_out,'bilinear')
    regridder.clean_weight_file()
    edgar_interp = regridder(edgar)
    # Check to ensure that interpolation worked
    fig = plt.figure()
    ax1 = plt.subplot2grid((2,2),(0,0),colspan=2)
    ax2 = plt.subplot2grid((2,2),(1,0),colspan=2)
    clevs = np.linspace(0, edgar.max()/200, 20)
    ax1.contourf(np.mean(edgar, axis=0), clevs)
    ax1.set_title('Original')
    m2 = ax2.contourf(np.mean(edgar_interp, axis=0), clevs)
    ax2.set_title('Interpolated')
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(m2, cax=cbar_ax)
    plt.subplots_adjust(hspace=0.4)
    plt.show()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
          'EDGAR interpolated to CTM resolution in %.2f seconds' 
          %((time.time()-start_time)))
    return edgar_interp

def open_schnello3(years, months, domain):
    """open O3 observations from Schnell et al. (2014) interpolated onto 
    1˚ x 1˚ gridded surface for United States or European domain for specified 
    time period.  
    
    Parameters
    ----------
    years : list
        Year or range of years in measuring period
    months : list
        Three letter abbreviations (lowercase) for months in measuring period
    domain : str
        Either 'US' or 'EU'
        
    Returns
    -------
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]
    o3 : numpy.ndarray
        Maximum daily 8-hour average (MDA8) surface-level O3, [time, lat, lng]    
    
    References
    ----------
    [1] J. L. Schnell, C. D. Holmes, A. Jangam, and M. J. Prather, "Skill in 
    forecasting extreme ozone pollution episodes with a global atmospheric 
    chemistry model," Atmos. Chem. Phys., 14, 7721-7739, 2014. 
    """
    import time
    import numpy as np
    import calendar
    import xarray as xr
    start_time = time.time()
    PATH_OBS = '/Users/ghkerr/phd/globalo3/data/'
    ds = xr.open_dataset(PATH_OBS+'mda8.surfO3.%s.2000.2009.nc'%domain)
    whereyear = np.where(np.in1d(ds.time.dt.year, np.array(years)) == 
                         True)[0]
    # Convert month abbreviations to integers
    months_int = []
    for month in months:
        months_int.append(list(calendar.month_abbr).index(month.title()))
    wheremonth = np.where(np.in1d(ds.time.dt.month, np.array(months_int)) == 
                          True)[0]
    # Find years and months of interest
    wheretime = np.where(np.in1d(whereyear, wheremonth) == True)[0]
    # Extract data
    o3 = ds['MDA8_SurfO3'][wheretime]
    lat = ds['lat']
    lng = ds['lon']
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
          'Schnell et al. (2014) gridded O3 for %d-%d loaded in %.2f seconds' 
          %(years[0], years[-1],(time.time()-start_time)))
    return lat.data, lng.data, o3.data