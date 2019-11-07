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
    07042019 -- function 'open_toar' added
    11062019 -- function 'open_geos_c1sd' added to open output from GSFC
                GEOS-C1SD tracer simulations
    18062019 -- 'open_geos_c1sd' modified to handle datasets with no vertical
                coordinates (i.e., MSLP) and transport output arrays such that 
                time is the first dimension
    16062019 -- function 'open_merra2_specifieddomain' moved from 
                transporto3_open and edited 
    20072019 -- function 'open_merra2_cyclones' added
    22072019 -- edited function 'interpolate_merra_to_ctmresolution' to 
                indicate grid cells outside of domain as NaNs
    23072019 -- function 'open_napso3' added 
    13082019 -- function 'open_geos_c1sd' removed and replaced with 
                'open_m2g_c90'
    20082019 -- cyclones' longitude converted from (-180˚-180˚) to (0˚-360˚) in 
                function 'open_merra2_cyclones'
    18092019 -- function 'open_merra2pblh_specifieddomain' added
    29092019 -- function 'open_merra2u_specifieddomain' added
    17102019 -- function 'open_geoschem_merra2_2x25_RnPbBe' added
    31122019 -- function 'open_merra2pblh_specifieddomain' edited to handle 
                hourly PBLH rather than daily-averaged heights; function can 
                now produce daily averages or averages over selected hours
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
    sys.path.append('/Users/ghkerr/phd/utils/')
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
#            # Drop leap days from analysis
#            if (calendar.isleap(year)==True) and (month == 'feb'):
#                gas_month = gas_month[:-1]
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
            lat = lat[latmin:latmax+1]
            lng = lng[lngmin:lngmax+1]                   
        # Extract desired data, roll deposition grid similar to longitude grid
        # Restrict data to focus region
        dep_yr = ds['%s_%s'%(process,gas)].data
        dep_yr = np.roll(dep_yr, int(dep_yr.shape[-1]/2)-1, axis=2)    
        dep_yr = dep_yr[:,latmin:latmax+1,lngmin:lngmax+1]
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
        Daily maximum 2-meter temperatures, units of K, [time, lat, lng]
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
        PATH_MERRA = '/Users/ghkerr/phd/globalo3/data/MERRA-2/%d/inst1_2d_asm_Nx/'%year        
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
    import time
    start_time = time.time()
    # Using xESMF yields no problems if the input and output domains cover the
    # exact same range, however, if one domain is smaller than other (or vice
    # versa), all values outside this range are set to 0.0. This is a problem 
    # when plotting: I cannot simply set all values that are 0.0 to NaNs as 
    # some values might be real. Therefore, I got around this with the 
    # following: https://github.com/JiaweiZhuang/xESMF/issues/15
    C = np.pi*np.e**5
    t2m = t2m + C
    # Interpolate finer grid (MERRA-2) to the resolution of the CTM
    grid_in = {'lon':lng_merra,'lat':lat_merra}
    # Output grid has a coarser resolution
    grid_out = {'lon':lng_gmi,'lat':lat_gmi}
    # Use xESMF with numpy array
    regridder = xe.Regridder(grid_in,grid_out,'bilinear')
    regridder.clean_weight_file()
    t2m_interp = regridder(t2m)
    # To get around the aforementioned issue of 0.0 being placed for regions
    # outside the domain
    t2m_interp[t2m_interp==0.0] = np.nan
    t2m_interp = t2m_interp - C
    t2m = t2m - C
    if checkplot == 'yes':
        import matplotlib.pyplot as plt        
        # Check to ensure that interpolation worked
        fig = plt.figure()
        ax1 = plt.subplot2grid((2,2),(0,0),colspan=2)
        ax2 = plt.subplot2grid((2,2),(1,0),colspan=2)
        clevs = np.linspace(np.nanmin(t2m), np.nanmax(t2m), 10)
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
    PATH_OBS = '/Users/ghkerr/phd/globalo3/data/observations/schnell/'
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
    return lat.data, lng.data, np.swapaxes(o3.data, 2, 1)

def open_merra2_rh2mvars(years, hours, lngmin, latmax, lngmax, latmin):
    """Using function "open_merra2," function opens the components needed 
    to calculated relative humidity (specific humidity, surface pressure, 
    temperature) at 2 meters and produces a daily average over the hours 
    specified in input parameter. 

    Parameters
    ----------
    years : list
        Year or range of years in measuring period, [years,]
    hours : list 
        Hours (in Zulu time) during which reanalysis data, n.b., for 
        inst3_3d_asm_Np collection, output is 3-hourly (0, 3, 9, etc.), 
        [hours,]
    lngmin : float
        Longitude coordinate of the left side (minimum) of the bounding box 
        containing the focus region, units of degrees east        
    latmax : float 
        Latitude coordinate of the top side (maximum) of the bounding box 
        containing the focus region, units of degrees north    
    lngmax : float 
        Longitude coordinate of the right side (maximum) of the bounding box 
        containing the focus region, units of degrees east        
    latmin : float
        Latitude coordinate of the bottom side (minimum) of the bounding box 
        containing the focus region, units of degrees north           

    Returns
    -------    
    QV2M : numpy.ndarray
        MERRA-2 2-meter specific humidity regridded to GMI resolution, units
        of kg kg-1, [time, lat, lng]
    T2M : numpy.ndarray
        MERRA-2 2-meter temperature regridded to GMI resolution, units of K, 
        [time, lat, lng]    
    PS : numpy.ndarray
        MERRA-2 surface pressure regridded to GMI resolution, units of Pa, 
        [time, lat, lng]
    lat_gmi_n : numpy.ndarray
        GMI latitude coordinates, units of degrees north, [lat,]
    lng_gmi_n : numpy.ndarray
        GMI longitude coordinates, units of degrees east, [lng,]
    """
    import sys
    sys.path.append('/Users/ghkerr/phd/transporto3/')
    import transporto3_open
    sys.path.append('/Users/ghkerr/phd/globalo3/')
    import globalo3_open    
    # Open GMI output to obtain lat/lng
    lat_gmi_n, lng_gmi_n, times_n, o3_n = \
        globalo3_open.open_overpass2_specifieddomain([2008], 
        ['jun','jul','aug'], latmin, latmax, lngmin, lngmax, 'O3', 
        'HindcastMR2')
    # Load 2-meter specific humidity, temperature and surface pressure and 
    # interpolate
    QV2M, mtime, lat_merra, lng_merra = transporto3_open.open_merra2(years, hours, 
        'QV2M', 'tavg1_2d_slv_Nx', 'JJA_rh.nc', lngmin, latmax, lngmax, latmin, 
        dailyavg='yes')
    QV2M = interpolate_merra_to_ctmresolution(lat_gmi_n, lng_gmi_n, 
        lat_merra, lng_merra, QV2M)
    T2M, mtime, lat_merra, lng_merra = transporto3_open.open_merra2(years, hours, 
        'T2M', 'tavg1_2d_slv_Nx', 'JJA_rh.nc', lngmin, latmax, lngmax, latmin, 
        dailyavg='yes')
    T2M = interpolate_merra_to_ctmresolution(lat_gmi_n, lng_gmi_n, 
        lat_merra, lng_merra, T2M)
    PS, mtime, lat_merra, lng_merra = transporto3_open.open_merra2(years, hours, 
        'PS', 'tavg1_2d_slv_Nx', 'JJA_rh.nc', lngmin, latmax, lngmax, latmin, 
        dailyavg='yes')
    PS = interpolate_merra_to_ctmresolution(lat_gmi_n, lng_gmi_n, 
        lat_merra, lng_merra, PS)
    return QV2M, T2M, PS, lat_gmi_n, lng_gmi_n

def open_toar(years, months, varname, res):
    """For the variable/metric of interest function opens monthly mean 
    daily maximum 8-hour averaged gridded ozone for the specified resolution. 
    
    Parameters
    ----------
    years : list
        Year or range of years in measuring period, [years,]
    months : list
        Three letter abbreviations (lowercase) for months in measuring period
    varname : str
        Variable of interest (e.g., 'rural_mean', 'rural_median', 'urban_mean',
        urban_median')
    res : int 
        Resolution of TOAR gridded set (resolutions of 10°×10°, 5°×5°, and 
        2°×2° are available). Recommended use of the 5° longitude ×5° latitude 
        products is encouraged as they provide a reasonable compromise between 
        global coverage and regional differentiation.
    
    Returns
    -------
    var : numpy.ndarray
         TOAR output for specified metric, units of ppbv, [time, lat, lng]
    ttime : pandas.core.indexes.datetimes.DatetimeIndex
        Timestamps of months for which TOAR data is desired, [time,]
    lat : numpy.ndarray
        TOAR latitude coordinates for resolution of interest, units of degrees 
        north, [lat,]
    lng : numpy.ndarray
        TOAR longitude coordinates for resolution of interest, units of degrees 
        east, [lng,]    
    """
    import time
    start_time = time.time()
    import numpy as np
    import pandas as pd
    import calendar
    import xarray as xr
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading TOAR %s O3...' %varname)
    PATH_TOAR = '/Users/ghkerr/phd/globalo3/data/TOAR/'
    # n.b., this will open the monthly timeseries TOAR data at specified 
    # resolution
    ds = xr.open_dataset(PATH_TOAR+'%dx%d_degrees/'%(res,res)+
                         'timeseries_1990-2014/'+
                         'toar_monthly_dma8epax_1990-2014_2x2.nc')
    whereyear = np.where(np.in1d(ds.time.dt.year, np.array(years)) == 
                         True)[0]
    # Convert month abbreviations to integers
    months_int = []
    for month in months:
        months_int.append(list(calendar.month_abbr).index(month.title()))
    wheremonth = np.where(np.in1d(ds.time.dt.month, np.array(months_int)) == 
                          True)[0]
    # Find years and months of interest; n.b., the first wheretime finds the 
    # indices intersecting years and months with respect to whereyear, so
    # whereyear must be indexed to find the overall position of months/years of
    # interest in the array
    wheretime = np.where(np.in1d(whereyear, wheremonth) == True)[0]
    wheretime = whereyear[wheretime]
    #dt = pd.to_datetime(ds.time.data)
    #dt = dt[wheretime]
    # Extract relevant variable 
    ds = ds[[varname]]
    ds = ds.assign_coords(lon=(ds.lon % 360)).roll(
            lon=(ds.dims['lon']//2), roll_coords=True) # ghk: removed -1 from 
            # ds.dims['lon']//2-1
    ds = ds.isel(time=wheretime)  
    var = ds[varname].values
    ttime = pd.to_datetime(ds.time.values)
    lat = ds.lat.values
    lng = ds.lon.values
    print('TOAR %s for %d-%d loaded in %.2f seconds'%(varname, years[0],
             years[-1], time.time() - start_time))
    return var, ttime, lat, lng

def open_m2g_c90(years, varname, levmax, levmin, lngmin, latmax, lngmax, 
    latmin, columnmean=False):
    """function opens daily output for the specified months and years for the
    CCMI tracer of interest from MERRA2-GMI (Replay) M2G_c90 simulations for 
    the region and pressure level(s) of interest. If specified, function can 
    also compute the column-averaged value.

    Parameters
    ----------
    years : list
        Year or range of years in measuring period, [years,]
    varname : str
        Variable of interest from simulations. Options include
        ps (surface pressure), u (wind), t (temperature), co50 (CO emissions
        over mainly Asia with a 50 day life time), co25 (same as co50 but with 
        a 25 day lifetime), nh50 (fixed mixing ratio at the surface between 
        30 and 50˚N with a 50 day lifetime), st80_25 (fixed mixing ratio in the 
        stratosphere, < 80 hPa, and 25 day lifetime in the troposphere).
    levmax : int/float
        Desired pressure level closest to the surface; n.b. function finds 
        pressure level closest to value in coordinates
    levmin : int/float
        Desired pressure level aloft 
    lngmin : float
        Longitude coordinate of the left side (minimum) of the bounding box 
        containing the focus region, units of degrees east        
    latmax : float 
        Latitude coordinate of the top side (maximum) of the bounding box 
        containing the focus region, units of degrees north    
    lngmax : float 
        Longitude coordinate of the right side (maximum) of the bounding box 
        containing the focus region, units of degrees east        
    latmin : float
        Latitude coordinate of the bottom side (minimum) of the bounding box 
        containing the focus region, units of degrees north      
    columnmean : bool
        If True, the field is averaged over all pressure levels bounded by 
        levmax and levmin
        
    Returns
    -------
    var : numpy.ndarray
        Model output for specified variable, units of ppbv, if columnsum = True 
        then shape is [time, lat, lng] if False then shape is [time, lev, 
        lat, lng]
    lat : numpy.ndarray
        Model latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Model numpy.ndarray coordinates, units of degrees east, [lng,]    
    lev : numpy.ndarray
        Model pressure levels, units of hPa, [lev]
    """
    import time
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading %s from MERRA2-GMI M2G_c90 simulation...' %varname)
    start_time = time.time()
    import os
    import numpy as np
    import xarray as xr
    # Search for appropriate files given variable of interested specified in 
    # varname
    infiles = []
    for year in years:
        PATH_TRACER = '/Users/ghkerr/phd/globalo3/data/GMI/M2G_c90/%d/'%year
        infiles_ty = [PATH_TRACER+fn for fn in os.listdir(PATH_TRACER) if 
            any(ext in fn for ext in [str(y) for y in years])]
        infiles_ty.sort()
        infiles.append(infiles_ty)
    infiles = list(np.hstack(infiles))
    # Open multiple NetCDF files and store in Dataset
    ds = xr.open_mfdataset(infiles, concat_dim='time')
    # Extract relevant variable 
    ds = ds[[varname]]
    ds = ds.assign_coords(lon=(ds.lon % 360)).roll(
        lon = (ds.dims['lon']//2-1), roll_coords = True)
    # Subset Dataset with respect to spatial coordinates 
    ds = ds.sel(lat = slice(latmin, latmax), lon = slice(lngmin, lngmax))
    # Sample appropriate level if pressure coordinates exist 
    ds = ds.sel(lev = slice(levmax,levmin))
    # Extract
    var = ds[varname].values
    lat = ds.lat.values
    lng = ds.lon.values
    lev = ds.lev.values
    # If specified, average variable over all pressure levels (i.e., if 
    # column-averaged mixing ratios for a tracer is desired)
    if (columnmean == True) and (len(lev) > 1):
        lev_dim = var.shape.index(lev.shape[0])
        var = np.nanmean(var, axis = lev_dim)
    elif (columnmean == True) and (len(lev) == 1):
        var = var[:,:,0]
    var = var*1e9
    print('%d-%d hPa %s for %d-%d loaded in '%(lev[0], lev[-1], 
          varname, years[0], years[-1])+'%.2f seconds!'
          %(time.time() - start_time))
    return var, lat, lng, lev     
    
def open_merra2_specifieddomain(years, months, hours, var, collection, lngmin, 
    latmax, lngmax, latmin, dailyavg='yes'): 
    """load MERRA-2 reanalysis over the specified measuring period and spatial 
    domain for desired variable (n.b. the function works assuming that the 
    input .nc file contains variables for a particular level).
    
    Parameters
    ----------
    years : list
        Year or range of years in measuring period, [years,]
    months : list
        Three letter abbreviations (lowercase) for months in measuring period
    hours : list 
        Hours (in Zulu time) during which reanalysis data, n.b., for 
        inst3_3d_asm_Np collection, output is 3-hourly (0, 3, 9, etc.), 
        [hours,]
    var : str
        MERRA-2 met field variable in collection, n.b. for inst3_3d_asm_Np
        collection stored locally, variables H, OMEGA, QV, SLP, T, U, V
        available
    collection : str
        The name of the MERRA-2 collection (e.g., inst3_3d_asm_Np)
    lngmin : float
        Longitude coordinate of the left side (minimum) of the bounding box 
        containing the focus region, units of degrees east        
    latmax : float 
        Latitude coordinate of the top side (maximum) of the bounding box 
        containing the focus region, units of degrees north    
    lngmax : float 
        Longitude coordinate of the right side (maximum) of the bounding box 
        containing the focus region, units of degrees east        
    latmin : float
        Latitude coordinate of the bottom side (minimum) of the bounding box 
        containing the focus region, units of degrees north           
    dailyavg : str
        If 'yes' a daily average over the hours specified in variable 'hours'
        is calculated. 

    Returns
    -------
    var : numpy.ndarray
         MERRA-2 output for specified variable, if argument 'dailyavg' is 'yes' 
         shape is [days, lat, lng] and if 'dailyavg' is 'no' shape is 
         [days*hours, lat, lng]
    mtime : pandas.core.indexes.datetimes.DatetimeIndex
        Timestamps of dates for which MERRA-2 data is desired, [days,]
    lat : numpy.ndarray
        MERRA-2 latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        MERRA-2 longitude coordinates, units of degrees east, [lng,]
    """
    import time
    start_time = time.time()
    import os
    import numpy as np
    import pandas as pd
    import xarray as xr
    import fnmatch
    from time import strptime
    # Convert month abbrevations to integers
    months_str = []
    for m in months:
        months_str.append(str(strptime(m,'%b').tm_mon).zfill(2))
    # Search files for years, months of interest
    infiles = []
    for year in years: 
        PATH_METFIELDS = '/Users/ghkerr/phd/globalo3/data/MERRA-2/%d/%s/'%(
                year, collection)
        infilesty = [PATH_METFIELDS+fn for fn in os.listdir(PATH_METFIELDS) if
                     any(fn.__contains__(str(year)+m) for m in months_str)]
        infilesty.sort()
        infiles+=infilesty
    # Open multiple NetCDF files and store in Dataset
    ds = xr.open_mfdataset(infiles, mask_and_scale=False, combine='nested', 
        concat_dim='time')
    varname = ds[var].long_name
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading MERRA-2 %s...' %varname)
    # Extract relevant variable 
    ds = ds[[var]]
    # Subset Dataset with respect to spatial coordinates
    ds = ds.assign_coords(lon=(ds.lon % 360)).roll(
            lon=(ds.dims['lon']//2-1), roll_coords=True)
    ds = ds[var].sel(lat=slice(latmin, latmax), 
                lon=slice(lngmin, lngmax))
    # Convert np.datetime64 to DatetimeIndex
    dt = pd.to_datetime(ds.time.data)
    # Find relevant hours and subset Dataset
    dt_where = np.where(np.in1d(dt.hour, hours)==True)[0]
    dt = dt[dt_where]
    # Select reanalysis at time of interest and sample appropriate level (if)
    # appropriate)
    try:
        ds.lev
    except AttributeError:
        ds = ds.isel(time=dt_where)       
    else: 
        ds = ds.isel(time=dt_where, lev=0)
    # Sort by time
    ds = ds.sortby('time')
    # Produce daily average
    if dailyavg == 'yes':
        ds = ds.resample({'time': '1D'}).mean(dim='time', skipna=True)
        # Resampling adds time dimensions to smooth out the discontinuous time
        # series (i.e., examining JJA 2008-2010 will add in September through 
        # May values for years); find indices of days not in JAA
        months_int = []
        for m in months:
            months_int.append(strptime(m,'%b').tm_mon)        
        season = np.where(np.in1d(ds.time.dt.month.values, 
            np.array(months_int)) == True)[0]
        ds = ds.isel(time=season)
    # Extract
    var = ds.values
    mtime = pd.to_datetime(ds.time.values)
    lat = ds.lat.values
    lng = ds.lon.values
    print('MERRA-2 %s for %d-%d loaded in %.2f seconds!'%(varname, years[0],
          years[-1], time.time() - start_time))
    return var, mtime, lat, lng 

def open_merra2_cyclones(sday, eday, months_str):
    """for a given date or range of dates, function opens netCDF files 
    corresponding to extratropical cyclones derived from the MERRA-2 database. 
    Although there are several variables within netCDF files, only the 
    latitude, longitude, sea level pressure (SLP), storm ID, and fraction of 
    land cover occupied by the cyclone are extracted. 
    
    From the cyclone database website: frontal boundaries are obtained using 
    the second version of the Modern Era Retrospective Analysis for Research 
    and Applications (MERRA-2, Gelaro et al., 2017) temperature, wind and 
    geopotential height output. The period covered is September 2006 to 
    December 2016, and each year of data files is contained in one TGZ archive. 
    Each file corresponds to one 6-hourly cyclone occurrence, and the file 
    names contain information on the cyclone, as described in the README.
    
    Parameters
    ----------
    sday : str
        Starting day of measuring period in MM/DD/YYYY format 
    eday : str
        Ending day of measuring period in MM/DD/YYYY format 
    months_str : list
        Three letter abbreviations (lowercase) for months in measuring period
    
    Returns
    -------
    cyclones : pandas.core.frame.DataFrame
        DataFrame containing the date, fraction of land cover, latitude, 
        longitude, SLP, and storm ID
    
    Notes
    -----
    Database of Cyclones with Warm and Cold Front Locations is found at 
    https://data.giss.nasa.gov/storms/obs-etc/
    
    References
    ----------
    Gelaro, R., W. McCarty, M.J. Suárez, R. Todling, A. Molod, L. Takacs, 
    C.A. Randles, A. Darmenov, M.G. Bosilovich, R. Reichle, K. Wargan, 
    L. Coy, R. Cullather, C. Draper, S. Akella, V. Buchard, A. Conaty, A.M. da 
    Silva, W. Gu, G. Kim, R. Koster, R. Lucchesi, D. Merkova, J.E. Nielsen, 
    G. Partyka, S. Pawson, W. Putman, M. Rienecker, S.D. Schubert, M. 
    Sienkiewicz, and B. Zhao, 2017: The Modern-Era Retrospective Analysis for 
    research and Applications, version 2 (MERRA-2). J. Climate, 30, 5419-5454, 
    doi:10.1175/JCLI-D-16-0758.1.
    """
    import time
    start_time = time.time()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading MERRA-2 cyclone database for %s-%s...' %(sday, eday))
    from time import strptime
    import os
    import numpy as np
    import xarray as xr
    import pandas as pd
    PATH_CYCLONES = '/Users/ghkerr/phd/globalo3/data/cyclones/'
    # Convert month abbrevations to integers
    months = []
    for m in months_str:
        months.append(strptime(m,'%b').tm_mon) 
    # Formulate days in measuring period
    days_mp = pd.date_range(start = sday, end = eday)
    # Filter out days in range that are not in the season of interest (i.e., 
    # if the start date is 06/01/2008 and the end date is 08/31/2009 and only 
    # JJA cyclones are desired, then this could remove September 2008 - May 
    # 2009 dates from the range) 
    inrange = np.where(np.in1d(np.array(days_mp.month), 
        np.array(months)) == True)[0]
    days_mp = days_mp[inrange]
    # Lists to be filled with values for entire measuring period
    lng_all, lat_all, slp_all, fracland_all, idoftrack_all, datestr_all = \
        [],[],[],[],[], []
    # Find files of cyclones on particular day 
    for day in days_mp:
        year = str(day.year)
        month = str(day.month).zfill(2)
        directory = PATH_CYCLONES+'/merra2fronts%s/'%(year)+\
            '%s%s/'%(year, month)
        infiles = [directory+fn for fn in os.listdir(directory) if
                   fn.__contains__(day.strftime('%Y%m%d'))]
        # Open all cyclones on a particular day 
        dropvars = ['longitude', 'latitude', 'MERRA2SLP', 'CF_hewson1km', 
                    'CF_simmonds850', 'CF_combined', 'WF_Hewson850', 
                    'WF_Hewson1km', 'WF_HewsonWB']
        ds = xr.open_mfdataset(infiles, data_vars=['storm_info'],
            combine='nested', concat_dim='time', coords='minimal', 
            drop_variables=dropvars)
        # Extract latitude, longitude, SLP, and fraction of land cover from 
        # variable 'storm_info' and create date entry for each cyclone
        storm_info = ds.storm_info.values
        lng = storm_info[:, 0]
        lat = storm_info[:, 1]
        slp = storm_info[:, 2]
        fracland = storm_info[:, 3]
        datestr = np.repeat(day.strftime('%Y-%m-%d'), len(lng))
        # From https://portal.nccs.nasa.gov/datashare/Obs-ETC/Fronts-ETC/
        # Readme_MERRA2fronts.pdf, the filename convention means that each 
        # cyclone has an ID specified in the filename: e.g., 
        # MERRA2fronts_YYYYMMDD_UT_Latstorm_Lonstorm_surfacetype_IDofTrack
        idoftrack = [fn.split('_')[6][:-5] for fn in infiles]
        # Append daily values to lists for entire measuring period
        datestr_all.append(datestr)    
        lng_all.append(lng)
        lat_all.append(lat)
        slp_all.append(slp)
        fracland_all.append(fracland)
        idoftrack_all.append(idoftrack)
    datestr_all = np.hstack(datestr_all)
    # Convert longitude from (-180-180) to (0-360)
    lng_all = np.hstack(lng_all)
    lng_all = lng_all % 360
    lat_all = np.hstack(lat_all)
    slp_all = np.hstack(slp_all)
    fracland_all = np.hstack(fracland_all)
    idoftrack_all = np.hstack(idoftrack_all)
    # Create DataFrame and fill 
    cyclones = pd.DataFrame({'Date' : datestr_all, 
        'Longitude' : lng_all,
        'Latitude' : lat_all, 
        'SLP' : slp_all,
        'Land Cover' : fracland_all,
        'Storm ID' : idoftrack_all})
    print('Cyclone coordinates loaded in %.2f seconds!' %(time.time() - 
                                                          start_time))
    return cyclones

def open_merra2pblh_specifieddomain(years, hours, latmin, latmax, lngmin, 
    lngmax, dailyavg='no'):
    """Load hourly PBL height from MERRA-2 over the specified spatial and
    temporal domains and, if desired, produce daily-averaged values. 
    
    Parameters
    ----------
    years : list
        Year or range of years in measuring period
    hours : numpy.ndarray
        Hours (in Zulu time) during which reanalysis data is desired
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
    dailyavg : str
        If 'yes' a daily average over the hours specified in variable 'hours'
        is calculated.         

    Returns
    -------
    mtime : pandas.core.indexes.datetimes.DatetimeIndex
        Timestamps of MERRA-2 data, [time,]    
    lat : numpy.ndarray
        MERRA-2 latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        MERRA-2 longitude coordinates, units of degrees east, [lng,]
    var : numpy.ndarray
        Daily mean or hourly PBL height, units of m, [time, lat, lng]
    """
    import time
    start_time = time.time()
    import os
    import numpy as np
    import pandas as pd
    import xarray as xr
    # Find relevant files in years of interest
    infiles = []
    for year in years:
        PATH_PBLH = '/Users/ghkerr/phd/meteorology/data/'+\
            'tavg1_2d_flx_Nx/Y%s/'%year    
        infiles_ty = [PATH_PBLH+fn for fn in os.listdir(PATH_PBLH) if 
            any(fn.__contains__(ext) for ext in [str(y) for y in years])]
        infiles_ty.sort()
        infiles.append(infiles_ty)
    infiles = list(np.hstack(infiles))
    # Open multiple NetCDF files and store in Dataset
    ds = xr.open_mfdataset(infiles, combine='nested', concat_dim='time')
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading MERRA-2 PBL height...')
    # Extract PBL height
    ds = ds[['PBLH']]
    # Subset Dataset with respect to spatial coordinates
    ds = ds.assign_coords(lon=(ds.lon % 360)).roll(
        lon=(ds.dims['lon']//2-1), roll_coords=True)
    ds = ds['PBLH'].sel(lat=slice(latmin, latmax), 
        lon=slice(lngmin, lngmax))
    # Convert np.datetime64 to DatetimeIndex
    dt = pd.to_datetime(ds.time.data)
    # Find relevant hours and subset Dataset
    dt_where = np.where(np.in1d(dt.hour, hours)==True)[0]
    dt = dt[dt_where]
    ds = ds.isel(time=dt_where)       
    # Produce daily average
    if dailyavg == 'yes':
        ds = ds.resample({'time': '1D'}).mean(dim='time')
        # Resampling adds time dimensions to smooth out the discontinuous time
        # series (i.e., examining JJA 2008-2010 will add in September through 
        # May values for years); find indices of days not in JAA
        notJJA = np.where((ds.time.dt.month.values == 6) |
                          (ds.time.dt.month.values == 7) |
                          (ds.time.dt.month.values == 8))[0]
        ds = ds.isel(time=notJJA)
    # Extract
    var = ds.values
    mtime = pd.to_datetime(ds.time.values)
    lat = ds.lat.values
    lng = ds.lon.values
    print('MERRA-2 PBLH for %d-%d loaded in %.2f seconds!'%(years[0],years[-1], 
        time.time() - start_time))
    return mtime, lat, lng, var

def open_merra2u_specifieddomain(years, varname, levmax, levmin, lngmin, 
    latmax, lngmax, latmin, columnmean=False):
    """function opens daily mean MERRA2 inst3_3d_asm_Np fields (3d, 3-Hourly,
    Instantaneous, Pressure-Level, Assimilation, Assimilated Meteorological 
    Fields) for the specified months and years. The variable of interest 
    (either U or V) is extracted for the region and pressure level(s) of 
    interest. If specified, function can also compute the column-averaged 
    value.

    Parameters
    ----------
    years : list
        Year or range of years in measuring period, [years,]
    varname : str
        Variable of interest from GSFC/Held-Suarez runs. Options include
        ps (surface pressure), u (wind), t (temperature), co50 (CO emissions
        over mainly Asia with a 50 day life time), co25 (same as co50 but with 
        a 25 day lifetime), nh50 (fixed mixing ratio at the surface between 
        30 and 50˚N with a 50 day lifetime), st80_25 (fixed mixing ratio in the 
        stratosphere, < 80 hPa, and 25 day lifetime in the troposphere).
    levmax : int/float
        Desired pressure level closest to the surface; n.b. function finds 
        pressure level closest to value in coordinates
    levmin : int/float
        Desired pressure level aloft 
    lngmin : float
        Longitude coordinate of the left side (minimum) of the bounding box 
        containing the focus region, units of degrees east        
    latmax : float 
        Latitude coordinate of the top side (maximum) of the bounding box 
        containing the focus region, units of degrees north    
    lngmax : float 
        Longitude coordinate of the right side (maximum) of the bounding box 
        containing the focus region, units of degrees east        
    latmin : float
        Latitude coordinate of the bottom side (minimum) of the bounding box 
        containing the focus region, units of degrees north      
    columnmean : bool
        If True, the field is averaged over all pressure levels bounded by 
        levmax and levmin
        
    Returns
    -------
    var : numpy.ndarray
        Model output for specified variable, units of ppbv, if columnsum = True 
        then shape is [time, lat, lng] if False then shape is [time, lev, 
        lat, lng]
    lat : numpy.ndarray
        Model latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Model numpy.ndarray coordinates, units of degrees east, [lng,]    
    lev : numpy.ndarray
        Model pressure levels, units of hPa, [lev]
    """
    import time
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading %s from MERRA2 simulation...' %varname)
    start_time = time.time()
    import os
    import numpy as np
    import xarray as xr
    # Search for appropriate files given variable of interested specified in 
    # varname
    infiles = []
    for year in years:
        PATH_MERRA = '/Users/ghkerr/phd/meteorology/data/inst6_3d_ana_Np/%d/'%year
        infiles_ty = [PATH_MERRA+fn for fn in os.listdir(PATH_MERRA) if 
            any(ext in fn for ext in [str(y) for y in years])]
        infiles_ty.sort()
        infiles.append(infiles_ty)
    infiles = list(np.hstack(infiles))
    # Open multiple NetCDF files and store in Dataset
    ds = xr.open_mfdataset(infiles, concat_dim='time')
    # Extract relevant variable 
    ds = ds[[varname]]
    ds = ds.assign_coords(lon=(ds.lon % 360)).roll(
        lon = (ds.dims['lon']//2-1), roll_coords = True)
    # Subset Dataset with respect to spatial coordinates 
    ds = ds.sel(lat = slice(latmin, latmax), lon = slice(lngmin, lngmax))
    # Sample appropriate level if pressure coordinates exist 
    ds = ds.sel(lev = slice(levmax,levmin))
    # Extract
    var = ds[varname].values
    lat = ds.lat.values
    lng = ds.lon.values
    lev = ds.lev.values
    # If specified, average variable over all pressure levels (i.e., if 
    # column-averaged mixing ratios for a tracer is desired)
    if (columnmean == True) and (len(lev) > 1):
        lev_dim = var.shape.index(lev.shape[0])
        var = np.nanmean(var, axis = lev_dim)
    elif (columnmean == True) and (len(lev) == 1):
        var = var[:,:,0]
    print('%d-%d hPa %s for %d-%d loaded in '%(lev[0], lev[-1], 
          varname, years[0], years[-1])+'%.2f seconds!'
          %(time.time() - start_time))
    return var, lat, lng, lev

def open_geoschem_merra2_2x25_RnPbBe(years, months, varname, latmin, latmax,
    lngmin, lngmax, pmin, pmax, operation):
    """function opens daily mean GEOSChem RnPbBe output for tracer of interest
    and extracts columned diagnostics for the latitudes, longitudes, and levels
    of interest. A sum, average, standard deviation is performed over the 
    column, if desired.
    
    Parameters
    ----------
    years : list
        Year or range of years in measuring period
    months : list
        Three letter abbreviations (lowercase) for months in measuring period
    varname : str
        Variable name (options: SpeciesConc_GLOBAL_50, SpeciesConc_MIDLAT_50, 
        SpeciesConc_POLAR_50, SpeciesConc_TROPIC_50)
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
    pmin : float
        The upper pressure level of interest, units of hPa
    pmax : float
        The lower pressure level of interest, units of hPa    
    operation : str
        This operation will be conducted on the column of tracer diagnostics
        (options: mean, sum, standard deviation, none)
        
    Returns
    -------
    var : numpy.ndarray
        GEOSChem tracer diagnostics for the region of interest, units of dry 
        mixing ratio (mol mol-1), [time, lat, lng] or [time, lev, lat, lng]
    lat : numpy.ndarray
        GEOSChem latitude coordinates, units of degrees north, [lat,]    
    lng : numpy.ndarray
        GEOSChem longitude coordinates, units of degrees east, [lat,]    
    lev : numpy.ndarray
        GEOSChem pressure level coordinates, units of hPa, [lev,]
    """
    import time
    start_time = time.time()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading GEOSChem %s...' %varname)        
    import numpy as np
    from netCDF4 import Dataset
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    # Convert month abbrevations to integers
    months_int = []
    for m in months:
        months_int.append(str(time.strptime(m,'%b').tm_mon).zfill(2))
    # List will filled with monthly GEOSChem output for variable of interest
    var = []    
    # Loop through years, months of interest    
    for year in years:
        PATH_GEOSCHEM='/Users/ghkerr/phd/globalo3/data/GEOSChem/%d/'%(year)
        for month in months_int:
            infile = Dataset(PATH_GEOSCHEM+'GEOSChem.SpeciesConc.%d%s.nc4'
                %(year,month),'r') 
            # On first iteration, extract dimensions and find indicies 
            # corresponding to (closest to) desired domain
            if (year==years[0]) and (month==months_int[0]):
                lat = infile.variables['lat'][:]
                lng = infile.variables['lon'][:] 
                lev = infile.variables['lev'][:]*1000. # Convert to hPa
                # Convert longitude from (-180-180) to (0-360)
                lng = lng%360          
                # Shift grid such that it spans (0-360) rather than (180-360, w
                # 0-180)
                lng = np.roll(lng,int(lng.shape[0]/2))
                latmin = geo_idx(latmin,lat)
                latmax = geo_idx(latmax,lat)
                lngmin = geo_idx(lngmin,lng)
                lngmax = geo_idx(lngmax,lng)
                pmax = (np.abs(lev-pmax)).argmin()
                pmin = (np.abs(lev-pmin)).argmin()
                lnglen = lng.shape[0]
                # Restrict coordinates over focus region 
                lat = lat[latmin:latmax+1]
                lng = lng[lngmin:lngmax+1] 
                lev = lev[pmax:pmin+1]
            # Extract variable for the month
            var_month = infile.variables[varname][:]
            # Roll grid similar to longitude grid
            var_month = np.roll(var_month, int(var_month.shape[-1]/2), 
                axis=np.where(np.array(var_month.shape)==lnglen)[0][0])
            # Extract output for horizontal and verticle domain of interest;
            # for this to work, dimensions must be (time, lev, lat, lng)
            var_month = var_month[:, pmax:pmin+1, latmin:latmax+1, 
                lngmin:lngmax+1]
            if operation == 'mean':
                var_month = np.nanmean(var_month, axis=1)
                var.append(var_month)
            elif operation == 'sum':
                var_month = np.nansum(var_month, axis=1)            
                var.append(var_month)
            elif operation == 'std':
                var_month = np.nanstd(var_month, axis=1)
                var.append(var_month)
            elif operation == 'none':
                var.append(var_month)
    # Stack 
    var = np.vstack(var)
    print('%s for %d-%d loaded in %.2f seconds!' %(varname, years[0], 
        years[-1], (time.time()-start_time)))
    return var, lat, lng, lev