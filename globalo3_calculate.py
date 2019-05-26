#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module calculates statistical analysis of GMI outpu and MERRA-2 and emissions
data across the global domain.

Revision History
    23012019 -- initial version created
    28012019 -- changed functions to calculate trends to (1) calculate the 
                least-squares trend using spicy.stats.linregress (in order to 
                determine statistical significance) instead of np.polyfit and 
                (2) renamed these functions to specify that the least-squares
                (ls) trends were being calculated
    29012019 -- function 'calculate_mktrend' added
    30012019 -- trend-calculating functions consolidated into single function
    04042019 -- functions 'find_grid_in_bb' and 'find_grid_overland' added
    12042019 -- function 'find_field_atjet' added
    15042019 -- function 'calculate_o3jet_relationship' added
    18042019 -- function 'convert_uv_tocardinal' added
    14052019 -- function 'identify_SLPcenter' added
    17052019 -- function 'filter_center_bylat' added
    21052019 -- function 'filter_center_byjet' added
"""

def calculate_do3dt(t2m, o3, lat_gmi, lng_gmi):
    """calculate dO3/dT (Ordinary Least Squares linear regression of O3 against 
    temperature). 

    Parameters
    ----------
    t2m : numpy.ndarray
        Interpolated daily maximum 2-meter temperatures, [time, lat, lng]
    o3 : numpy.ndarray
        O3 concentrations at 1300-1400 hours local time, units of ppbv, [time,
        lat, lng]
    lat_gmi : numpy.ndarray
        GMI latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI longitude coordinates, units of degrees east, [lng,]        
        
    Returns
    -------
    do3dt : numpy.ndarray     
        dO3/dT, units of ppbv K-1, [lat, lng]
    """
    import numpy as np
    import time
    start_time = time.time()
    do3dt = np.empty(o3.shape[1:])
    for i,lat in enumerate(lat_gmi):
        for j,lng in enumerate(lng_gmi):
            try:
                do3dt[i,j] = np.polyfit(t2m[:,i,j], o3[:,i,j], deg=1)[0]
            except ValueError:
                do3dt[i,j] = np.nan
    print('dO3/dT in calculated in %.2f seconds' %((time.time()-start_time)))
    return do3dt

def calculate_r(x, y, lat, lng):
    """calculate the Pearson correlation coefficient between two gridded fields 
    with the same resolution.

    Parameters
    ----------
    x : numpy.ndarray
        Independent variable [time, lat, lng]
    y : numpy.ndarray
        Dependent variable [time, lat, lng]
    lat : numpy.ndarray
        Gridded latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Gridded longitude coordinates, units of degrees east, [lng,]        
        
    Returns
    -------
    r : numpy.ndarray     
        Pearson correlation coefficient between x and y, [lat, lng]
    """
    import numpy as np
    import time
    start_time = time.time()
    r = np.empty(y.shape[1:])
    for i,ilat in enumerate(lat):
        for j,ilng in enumerate(lng):
            r[i,j] = np.corrcoef(x[:,i,j], y[:,i,j])[0,1]
    print('correlation coefficient in calculated in %.2f seconds' 
          %((time.time()-start_time)))
    return r

def separate_years(x, lat_gmi, lng_gmi, years):
    """function separates continuous timeseries into equal-sized chunks with 
    the number of chunks corresponding to the number of years in the measuring 
    period; for example, the gridded O3 dataset might have dimensions [time,
    lat, lng] where time is (number of days in a year) x (no. years).
    This time dimension would be split into two: [years, days in year]
    
    Parameters
    ----------
    x : numpy.ndarray
        Continuous timeseries of gridded data, [time,lat,lng]
    lat_gmi : numpy.ndarray
        GMI latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI longitude coordinates, units of degrees east, [lng,]      
    years : list
        Year or range of years in measuring period
        
    Returns
    -------
    x_byyr : numpy.ndarray     
        Input array with added dimension corresponding to the year, 
        [no. years,no. days in year,lat,lng]
    """
    import numpy as np
    # Separate continuous timeseries into yearly chunks
    x_byyr = np.reshape(x, (len(years), int(x.shape[0]/len(years)), 
                              lat_gmi.shape[0], lng_gmi.shape[0]))
    return x_byyr

def calculate_trends(x, lat_gmi, lng_gmi, years, alpha=0.05):
    """function calculates yearly/seasonal averages and thereafter calculates 
    the (1) least-squares trend and significance and (2) the Mann-Kendall test 
    statistics and levels of significance to detect whether increasing or 
    decreasing monotonic trends exist. 
    
    Parameters
    ----------
    x : numpy.ndarray
        Timeseries of gridded data, [time,lat,lng] OR [no. years,no. days in 
        year,lat,lng]
    lat_gmi : numpy.ndarray
        GMI latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI longitude coordinates, units of degrees east, [lng,]      
    years : list
        Year or range of years in measuring period
    alpha : float
        The significance level
        
    Returns
    -------
    x : numpy.ndarray     
        Yearly/seasonally-averaged data, [years,lat,lng]
    trend : numpy.ndarray
        The slope of the regression line representing the trend in the yearly/
        seasonally-averaged data, units are data units yr-1, [lat,lng]
    p : numpy.ndarray
        Two-sided p-values for a hypotheis test whose null hypothesis is that
        the slope is 0, [lat,lng]
    mkz : numpy.ndarray
        Normalized  statistic (Z) based on the Mann-Kendall test statistic S; 
        the presence of a statistically significant trend is given by the Z 
        value. A positive (negative) value of Z indicates an upward (downward)
        trend, [lat,lng]
    mkp : numpy.ndarray
        To test for either an upward or downward monotone trend (a two-tailed 
        test) at α level of significance, H0 (which is, the slope is 0) is 
        rejected if the |Z| > Z1-α/2, where Z1-α/2 is obtained from the 
        standard normal cumulative distribution tables (i.e., for α = 0.05, |Z| 
        must be greater than 1.96), [lat,lng]         
    """
    import numpy as np
    from scipy.stats import linregress
    from mk_test import mk_test    
    import time
    start_time = time.time()
    # If data does not have a dimension for year, separate continuous 
    # timeseries into yearly chunks
    if x.shape[0] != len(years):
        x = separate_years(x, lat_gmi, lng_gmi, years)
        # Find yearly/seasonally-averaged values
        daydim = np.where((np.array(x.shape) != len(years)) &
                          (np.array(x.shape) != len(lat_gmi)) & 
                          (np.array(x.shape) != len(lng_gmi)))[0][0]
        x = np.mean(x, axis = daydim)
    # Calculate trend through yearly/seasonal values
    trend = np.empty(x.shape[1:])
    p = np.empty(x.shape[1:])    
    tl = np.arange(0, len(years), 1)
    for i,lat in enumerate(lat_gmi):
        for j,lng in enumerate(lng_gmi):
            trend[i,j] = linregress(tl, x[:,i,j]).slope
            p[i,j] = linregress(tl, x[:,i,j]).pvalue
    # Calculate Mann-Kendall trend through yearly/seasonal values
    mkp = np.empty(x.shape[1:])
    mkz = np.empty(x.shape[1:])
    for i,lat in enumerate(lat_gmi):
        for j,lng in enumerate(lng_gmi):
            pij, z = mk_test(x[:,i,j], alpha=alpha)
            mkp[i,j] = pij
            mkz[i,j] = z                    
    print('LS/MK trends in calculated in %.2f seconds' 
          %((time.time()-start_time)))
    return x, trend, p, mkz, mkp
    
def calculate_trends_do3dt(t2m, o3, lat_gmi, lng_gmi, years, alpha=0.05):
    """function first calculates dO3/dT for each year (in this
    case for each summer season) and thereafter fits (1) least-squares trends
    and significance and (2) the Mann-Kendall test statistics and levels of 
    significance to detect whether increasing or decreasing monotonic trends 
    exist. 
    
    Parameters
    ----------
    t2m : numpy.ndarray
        Interpolated daily maximum 2-meter temperatures, [time,lat,lng]
    o3 : numpy.ndarray
        O3 concentrations at 1300-1400 hours local time, units of ppbv, [time,
        lat,lng]
    lat_gmi : numpy.ndarray
        GMI latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI longitude coordinates, units of degrees east, [lng,]      
    years : list
        Year or range of years in measuring period
    alpha : float
        The significance level        
        
    Returns
    -------
    do3dt_byyr : numpy.ndarray
        Yearly/seasonally-averaged dO3/dT, units of ppbv K-1, [years,lat,lng]
    trend : numpy.ndarray
        The slope of the regression line representing the trend in the seasonal 
        dO3/dT, units of ppbv K-1 yr-1, [lat, lng]
    p : numpy.ndarray
        Two-sided p-values for a hypotheis test whose null hypothesis is that
        the slope is 0, [lat,lng]        
    mkz : numpy.ndarray
        Normalized  statistic (Z) based on the Mann-Kendall test statistic S; 
        the presence of a statistically significant trend is given by the Z 
        value. A positive (negative) value of Z indicates an upward (downward)
        trend, [lat,lng]
    mkp : numpy.ndarray
        To test for either an upward or downward monotone trend (a two-tailed 
        test) at α level of significance, H0 (which is, the slope is 0) is 
        rejected if the |Z| > Z1-α/2, where Z1-α/2 is obtained from the 
        standard normal cumulative distribution tables (i.e., for α = 0.05, |Z| 
        must be greater than 1.96), [lat,lng]          
    """
    import numpy as np
    from scipy.stats import linregress    
    import time
    start_time = time.time()
    # Separate continuous timeseries into yearly chunks
    o3_byyr = separate_years(o3, lat_gmi, lng_gmi, years)
    t2m_byyr = separate_years(t2m, lat_gmi, lng_gmi, years)
    # Calculate yearly dO3/dT
    do3dt_byyr = np.empty(shape=(len(years), lat_gmi.shape[0], 
                                 lng_gmi.shape[0]))
    for y,year in enumerate(years):
        for i,lat in enumerate(lat_gmi):
            for j,lng in enumerate(lng_gmi):
                do3dt_byyr[y,i,j] = linregress(t2m_byyr[y,:,i,j],
                          o3_byyr[y,:,i,j]).slope
    # Calculate trend through seasonal values of dO3/dT
    x, trend, p, mkz, mkp = calculate_trends(do3dt_byyr, lat_gmi, lng_gmi, 
                                             years, alpha)
    print('dO3/dT LS/MK trend in calculated in %.2f seconds' 
          %((time.time()-start_time)))
    return do3dt_byyr, trend, p, mkz, mkp

def find_grid_in_bb(ingrid, lat, lng, left, right, down, up): 
    """Given a bounding box (i.e., coordinates of minimum and maximum latitudes
    and longitudes), reduce a given grid and the dimensional coordinates to 
    only that focus region.
    
    Parameters
    ----------
    ingrid : numpy.ndarray
        Input grid
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]
    left : float
        Longitude coordinate of the left side (minimum) of the bounding box 
        containing the focus region, units of degrees east        
    right : float 
        Longitude coordinate of the right side (maximum) of the bounding box 
        containing the focus region, units of degrees east        
    down : float
        Latitude coordinate of the bottom side (minimum) of the bounding box 
        containing the focus region, units of degrees north            
    up : float 
        Latitude coordinate of the top side (maximum) of the bounding box 
        containing the focus region, units of degrees north    

    Returns
    -------
    outgrid : numpy.ndarray
        Output grid
    lat : numpy.ndarray
        Focus region latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Focus region longitude coordinates, units of degrees east, [lng,]
    """
    import numpy as np
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    # Find spines of focus region 
    left = geo_idx(left, lng)
    right = geo_idx(right, lng)
    up = geo_idx(up, lat)
    down = geo_idx(down, lat)
    # Reduce grid to only focus region
    lng_dim = np.where(np.array(ingrid.shape) == lng.shape[0])[0][0]
    lat_dim = np.where(np.array(ingrid.shape) == lat.shape[0])[0][0]
    # Eventually fix this section? Depending on whether the grid has a time
    # dimension or if lat/lng are reversed, indexing the grid will be handled
    # differently
    if (lat_dim == 2) and (lng_dim == 3):
        outgrid = ingrid[:, :, down:up+1, left:right+1]        
    elif (lat_dim == 1) and (lng_dim == 2):
        outgrid = ingrid[:, down:up+1, left:right+1]    
    elif (lat_dim == 0) and (lng_dim == 1):
        outgrid = ingrid[down:up+1, left:right+1]    
    else: 
        print('Dimensional information does not match!'+
              ' Cannot reduce grid.')    
    # Reduce lat/lng to focus region 
    lat = lat[down:up+1]
    lng = lng[left:right+1]
    return outgrid, lat, lng

def find_grid_overland(lat, lng): 
    """Determine which grid cells reside over land. From basemap documentation,
    the definition of land is based upon the GSHHS coastline polygons 
    associated with the class instance. Points over lakes inside land regions 
    are not counted as land points.

    Parameters
    ----------
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]

    Returns
    -------
    island : numpy.ndarray
        NaN corresponds to point over water, 1 over land [lat, lng]
    """
    import numpy as np
    from mpl_toolkits.basemap import Basemap
    # Default is projection is 'cyl'
    bm = Basemap(llcrnrlon=lng.min(),llcrnrlat=lat.min(),
                 urcrnrlon=lng.max(),urcrnrlat=lat.max())
    island = np.empty(shape=(lat.shape[0],lng.shape[0]), dtype='float')
    island[:] = np.nan
    # Loop through all unique lat/lng coordinate pairs
    for i, ilat in enumerate(lat):
        for j, ilng in enumerate(lng): 
            result = bm.is_land(ilng,ilat)
            # Set 'True' (land) to 1.; 'False' (not land) remain NaN.
            if result == True:
                island[i,j] = 1.
    return island

def calculate_rh_from_q(years, hours, lngmin, latmax, lngmax, latmin):
    """Calculate the relative humidity from specific humidity, temperature, and 
    pressure using the ratio of vapor pressure to saturation vapor pressures.
    Code adapted from earthscience.stackexchange.com/questions/2360/
    how-do-i-convert-specific-humidity-to-relative-humidity
    
    Parameters
    ----------    
    T : numpy.ndarray
        Atmospheric temperature, units of K, [time, lat, lng]
    Q : numpy.ndarray
        Specific humidity of air, units of kg kg-1, [time, lat, lng]
    PS : numpy.ndarray
         Atmospheric pressure, units of Pa, [time, lat, lng]

    Returns
    -------
    rh : numpy.ndarray
        Relative humidity as a unitless ratio, [time, lat, lng]
    """
    import time
    start_time = time.time()    
    import numpy as np
    import sys
    sys.path.append('/Users/ghkerr/phd/globalo3/')
    import globalo3_open    
    # Load fields needed for calculation
    Q, T, PS, lat_gmi, lng_gmi = globalo3_open.open_merra2_rh2mvars(years, 
        hours, lngmin, latmax, lngmax, latmin)
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
          'Calculating relative humidity using the Clausius-Clapeyron '+
          'equation...')   
    es = 6.112 * np.exp((17.67 * (T-273.15))/((T-273.15) + 243.5))
    e = Q * (PS/100.) / (0.378 * Q + 0.622)
    rh = e / es
    rh[rh > 1] = 1
    rh[rh < 0] = 0
    print('Relative humidity for %d-%d loaded in %.2f seconds!'%(years[0],
          years[-1], time.time() - start_time))    
    return rh

def find_field_atjet(field, U500_fr, lat_fr, lng_fr, jetdistance, anom=False): 
    """function identifies the latitude of maximum zonal (U) winds at 
    500 hPa (proxy for jet stream) and finds the values of the field of 
    interest at the jet and within +/- jetdistance of that latitude. If the 
    jet latitude is too close to the top or bottom of the focus region, then 
    the values of the field are found up until the focus region boundary and 
    remaining values are filled with NaNs. If anom=True, then the value of 
    the field at the jet's center is found and is subtracted from every 
    longitude band such that at the jet center the field=0 and values north/
    south represent departures from the value at the jet's center. 
    
    Parameters
    ----------
    field : numpy.ndarray
        Field of interest, [lat, lng] or [time, lat, lng]  
    U500_fr : numpy.ndarray
        Zonal (U) wind at 500 hPa in region, units of m s-1, [time, lat, lng]
    lat_fr : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng_fr : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]
    jetdistance : int
        Number of grid cells (with respect to latitude) on each side (north 
        and south) the jet (maximum zonal wind at 500 hPa) over which the 
        field will be evaluated
    anom : bool
        If True, field_jet will be calculated relative to values at the jet's 
        center

    Returns
    -------
    lat_jet : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in region, units of degrees north[time, lng]
    field_jet : numpy.ndarray
        The value of the field at the jet and within +/- jetdistance of jet, 
        [lat, lng] or [time, lat, lng]
    """
    import numpy as np
    # Array will be filled the daily latitude of the jet (defined as the 
    # location with maximum U wind at 500 mb) in the focus region
    lat_jet = np.empty(shape=(len(U500_fr), len(lng_fr)))
    lat_jet[:] = np.nan
    # Find latitude of jet location 
    for day in np.arange(0, len(U500_fr), 1):
        U500_day = U500_fr[day]
        # Loop through longitude
        for i in np.arange(0, len(lng_fr), 1):    
            # U wind at 500 hPa for longitude/day of interest 
            U500_transect = U500_day[:, i]
            U500max = np.where(U500_transect==U500_transect.max())[0][0]   
            # Find latitude of jet
            lat_jet[day, i] = lat_fr[U500max]
    del i
    # If field of interest has dimensions [time, lat, lng], output will be 
    # field north/south of jet for every timestep and will thus have dimensions 
    # of[time, 2*jetdistance+1, lng]. 
    if np.ndim(field) == 2:
        field_jet = np.empty(shape=(2*jetdistance+1, len(lng_fr)))
        field_jet[:] = np.nan
        # Mean jet location 
        lat_jet_mean = np.nanmean(lat_jet, axis=0)
        # Loop through longitude
        for i in np.arange(0, len(lng_fr), 1):    
            # Find mean latitude index of jet
            U500max = (np.abs(lat_fr-lat_jet_mean[i])).argmin()
            # Find value of field at a particular longitude at the jet's center
            field_jetcenter = field[U500max, i]       
            # Find values above/below jet's center
            field_transect_above = field[U500max+1:U500max+jetdistance+1, i]
            field_transect_below = field[U500max-jetdistance:U500max, i]        
            # Fill output grid
            field_jet[jetdistance, i] = field_jetcenter
            field_jet[jetdistance+1:jetdistance+1+len(field_transect_above), 
                      i] = field_transect_above            
            field_jet[jetdistance-len(field_transect_below):jetdistance, 
                      i] = field_transect_below
    # If field of interest is already time-averaged and has dimensions [lat, 
    # lng], the output will only have dimensions [2*jetdistance+1, lng]    
    if np.ndim(field) == 3:
        field_jet = np.empty(shape=(len(field), 2*jetdistance+1, len(lng_fr)))
        field_jet[:] = np.nan
        # Loop through time
        for day in np.arange(0, len(field), 1):    
            U500_day = U500_fr[day]
            # Loop through longitude
            for i in np.arange(0, len(lng_fr), 1):    
                U500max = (np.abs(lat_fr-lat_jet[day, i])).argmin()
                field_jetcenter = field[day, U500max, i]       
                # Find values above/below jet's center
                field_transect_above = field[day, U500max+1:U500max+
                                             jetdistance+1, i]
                field_transect_below = field[day, U500max-jetdistance:U500max, 
                                             i]        
                # Fill output grid
                field_jet[day, jetdistance, i] = field_jetcenter
                field_jet[day, jetdistance+1:jetdistance+1+len(
                          field_transect_above), i] = field_transect_above            
                field_jet[day, jetdistance-len(field_transect_below):
                          jetdistance, i] = field_transect_below
    # If anom==True, then value of the field at the jet's center will be 
    # subtracted off from the field, so that field[jetdistance, :] = 0.0 by 
    # definition            
    if anom==True: 
        if np.ndim(field) == 2:                    
            field_alongcenter = field_jet[jetdistance, :]
            field_alongcenter = np.tile(field_alongcenter, (len(field_jet), 1))
            field_jet = field_jet-field_alongcenter
        if np.ndim(field) == 3:
            field_alongcenter = field_jet[:, jetdistance, :]
            # Expand along axis by copying 2D array into 3D array 
            field_alongcenter = np.repeat(field_alongcenter[:,:,np.newaxis], 
                                          field_jet.shape[1], axis=2)
            field_alongcenter = np.swapaxes(field_alongcenter, 2, 1)
            field_jet = field_jet - field_alongcenter
    return lat_jet, field_jet

def calculate_o3jet_relationship(o3_fr, lat_fr, lng_fr, jetpos, lng_jet):
    """Given the O3 concentrations and the jet's latitude in a focus region, 
    function determines the relationship between O3 and distance from the 
    jet at each grid cell. This is quantified in terms of the slope (i.e., 
    change in O3 per degree shift in the jet stream) and the correlation (i.e.,
    the Pearson product-moment correlation coefficient calculated between O3 
    and a grid cell's distance from the jet.)
    
    Parameters
    ----------
    o3_fr : numpy.ndarray
        O3 concentrations in the focus region (i.e., Northern Hemisphere), 
        units of ppbv, [time, lat, lng]
    lat_fr : numpy.ndarray
        Latitude coordinates of the focus region (i.e., Northern Hemisphere), 
        units of degrees north, [lat,]
    lng_fr : numpy.ndarray
        Longitude coordinates of the focus region (i.e., Northern Hemisphere), 
        units of degrees east, [lng,]        
    jetpos : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in region, units of degrees north, [time, lngjet]
    lng_jet = numpy.ndarray
        The longitude grid over which the eddy-driven jet's position was 
        analyzed (i.e., since the EDJ is found by examining 500 hPa zonal winds
        in the mid-latitudes, the latitudes and longitudes defining the mid-
        latitudes could differ from the focus region's coordinates), 
        units of degrees east, [lngjet]

    Returns
    -------
    slope : numpy.ndarray
        The slope of the linear regression of O3 versus a grid cell's distance
        from the jet, units of ppbv deg-1, [lat, lng]
    correl : numpy.ndarray
        The Pearson product-moment correlation coefficient calculated between 
        O3 and a grid cell's distance from the jet, [lat, lng]
    """
    import numpy as np
    # Array slope is the slope of the linear regression of surface-level O3 and 
    # a grid cell's distance from the eddy-driven jet (positive distance is a 
    # jet north of grid cell); array correl is the correlation between the two 
    # times series
    slope = np.empty(shape=(lat_fr.shape[0],lng_fr.shape[0]))
    slope[:] = np.nan
    correl = np.empty(shape=(lat_fr.shape[0],lng_fr.shape[0]))
    correl[:] = np.nan
    # Loop through latitudes
    for loci in np.arange(0,len(lat_fr),1):
        # Loop through longitudes
        for locj in np.arange(0,len(lng_fr),1):
            o3ij = o3_fr[:, loci, locj]
            # Latitude of eddy-driven jet at a given grid cell, found by finding 
            # a timeseries of the jet location at the longitude of interest
            jetj = jetpos[:, np.abs(lng_jet-lng_fr[locj]).argmin()]
            # Difference in grid cell's latitude and the latitude of the jet
            diff = jetj-lat_fr[loci]
            slope[loci,locj] = np.polyfit(diff, o3ij, 1)[0]
            correl[loci,locj] = np.corrcoef(diff, o3ij)[0,1]
    return slope, correl

def convert_uv_tocardinal(U, V): 
    """calculate the wind direction using the typical meteorological convention
    where, for example, 0 deg denotes wind from the north (southerly flow) and 
    90 deg denotes wind from the east (westerly flow). 
    
    Parameters
    ----------
    U : numpy.ndarray
        The zonal wind component. Array size and shape must match V, units of 
        m s-1, [time, lat, lng]
    V : numpy.ndarray
        The meridional wind component. Array size and shape must match V, units 
        of m s-1, [time, lat, lng]
        
    Returns
    -------
    wind_dir_met : numpy.ndarray     
        Cardinal wind direction, units of degrees, [time, lat, lng]    
    """
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
          'Calculating wind direction...')    
    import time
    start_time = time.time()    
    import numpy as np
    # Find wind magntiude
    wind_abs = np.hypot(U, V)
    wind_dir_trig = np.arctan2(V/wind_abs, U/wind_abs) 
    wind_dir_degrees = wind_dir_trig*(180/np.pi)
    # np.arctan2 has a discontinunity at 180 deg, switching to (-180 to 0 deg)
    # going clockwise. Use modulo operator to convert these negative values. 
    # n.b., for positive (1 to 180 deg): if you mod any positive number between
    # 1 and 180 by 360, you will get the exact same number you put in. for 
    # negative (-180 to -1): using mod here will return values in the range 
    # of 180 and 359 degrees.
    wind_dir_degrees = (wind_dir_degrees+360) % 360
    # Convert from cartesian directions to cardinal conventions used in 
    # meteorology 
    wind_dir_met = (270-wind_dir_degrees) % 360
    print('Wind direction calculated in %.2f seconds!'
          %(time.time()-start_time)) 
    return wind_dir_met

def identify_SLPcenter(lat, lng, SLP, dx, dy, kind, pr_crit, years, 
    checkplot='no', fstr=''):
    """function identifies the centers of cyclones or anticyclones from sea 
    level pressure fields (SLP) over the specified domain. The center of a 
    cyclone (anticyclone) is defined as a local pressure minima (maxima). In 
    order for this to be true for a cyclone (anticyclone), the pressure at a 
    particular grid cell must be less (more) than the dx/dy grid cells on all
    sides of it and the pressure at the grid cell must be less (more) than the 
    criteria pressure specified in variable 'pr_crit.' If checkplot='yes', a
    map showing the locations of all (anti)cyclones and a binned 2D histogram
    indicating their frequency is plotted.)

    Parameters
    ----------
    lat : numpy.ndarray
        Latitude coordinates for SLP, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates for SLP, units of degrees east, [lng,]
    SLP : numpy.ndarray
        MERRA-2 Sea level pressure (SLP), units of Pa, [time, lat, lng]
    dx : int
        The search area in the x (longitudinal) direction over which pressure 
        values are examined to identify minima/maxima
    dy : int
        The search area in the y (latitudinal) direction over which pressure 
        values are examined to identify minima/maxima    
    kind : str
        cyclone or anticyclone
    pr_crit : float
        The critical pressure above (below) which a SLP minima (maxima) cannot
        be considered a cyclone (anticyclone); n.b. I have been using 1000 hPa
        (100000 Pa) and 1012 hPa (101200 Pa) as the critical pressures for 
        cyclones and anticylones, respectively, units of Pa
    years : list
        Years in measuring period
    checkplot : str
        If 'yes' a map indicating the frequency/location of cyclones/
        anticyclones is plotted
    fstr : str
        Output filename suffix (should indicate whether plot shows cyclone
        or anticylone)
        
    Returns
    -------
    center : numpy.ndarray 
        A value of 1 indicates the presence of a cyclone/anticylone for a 
        particular day and locaditon, [time, lat, lng]   
    where_center_daycoord : list
        The dates (day indices) of cyclones/anticyclones        
    where_center_ycoord : list 
        The y-coordinates (latitude indices) of cyclones/anticyclones
    where_center_xcoord : list 
        The x-coordinates (longitude indices) of cyclones/anticyclones        
    """
    import time 
    start_time = time.time()    
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
          'Calculating %s locations...'%kind)   
    import numpy as np
    center = np.empty(shape=SLP.shape)
    center[:] = np.nan    
    for day in np.arange(0, len(SLP), 1):
        SLP_day = SLP[day]
        # For given day, index pressure grid. The outer rows/columns of the 
        # grid since the neighboring eight grid points cannot be searched for 
        # a pressure minima; here i is the latitude index and j is the 
        # longitude index
        for i in np.arange(dy, len(lat)-dy, 1):
            for j in np.arange(dx, len(lng)-dx, 1):
                # SLP at point of interest and SLP at the eight grid cells
                # surrounding that grid cell
                SLP_day_atpoint = SLP_day[i, j]
                # Similar to Lang & Waugh (2011), each of the surrounding grid
                # points within a search radius defined by dx and dy (not 
                # counting the center grid point) must be greater than (for 
                # cyclones) or less than (for anticyclones) the pressure in the 
                # center)
                SLP_day_surround = SLP_day[i-dy:i+dy+1, j-dx:j+dx+1]
                # For cyclones
                if kind=='cyclone':
                    if (np.min(SLP_day_surround) >= SLP_day_atpoint) & \
                    (SLP_day_atpoint < pr_crit):
                        center[day,i,j] = 1.
                # For anticyclones
                else: 
                    if (np.max(SLP_day_surround) <= SLP_day_atpoint) & \
                    (SLP_day_atpoint > pr_crit):
                        center[day,i,j] = 1.       
    where_center_ycoord, where_center_xcoord = [], []
    where_center_daycoord = []
    where_center_lat, where_center_lng = [], []
    # Loop through days in measuring period
    for day in np.arange(0, len(center), 1):
        # Find position(s) of (anti)cyclone on day of interest
        where_center_day = center[day]
        where_center_day = np.where(where_center_day==1.)
        # Append latitude/longitude indices/coordinates to list containing
        # values for all days
        where_center_ycoord.append(where_center_day[0])
        where_center_xcoord.append(where_center_day[1])      
        where_center_daycoord.append(np.repeat(day, len(where_center_day[0])))
        where_center_lat.append(lat[where_center_day[0]])
        where_center_lng.append(lng[where_center_day[1]])
    where_center_daycoord = np.hstack(where_center_daycoord)
    where_center_ycoord = np.hstack(where_center_ycoord)
    where_center_xcoord = np.hstack(where_center_xcoord)
    where_center_lat = np.hstack(where_center_lat)
    where_center_lng = np.hstack(where_center_lng)
    if checkplot=='yes':
        import matplotlib.pyplot as plt
        import cartopy.util
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature  
        fig = plt.figure(figsize=(9,4))
        ax = plt.subplot2grid((1,1), (0,0), projection=ccrs.PlateCarree())
        # Do coordinate conversion of (x,y)
        xynps = ax.projection.transform_points(ccrs.PlateCarree(), 
            where_center_lng, where_center_lat)
        # Make 2D histogram; n.b., h:(counts, xedges, yedges, image)
        h = ax.hist2d(xynps[:,0], xynps[:,1], bins=15, zorder=10, alpha=1., 
            vmin=0, vmax=7, cmap=plt.get_cmap('Blues'), 
            transform=ccrs.PlateCarree())
        colorbar_axes = plt.gcf().add_axes([0.83,0.25,0.02,0.5])
        colorbar = plt.colorbar(h[3], colorbar_axes, orientation='vertical',
                                extend='max')
        colorbar.ax.tick_params(labelsize=12)
        colorbar.set_label('Frequency', fontsize=14)
        # Add scatterpoints for all 
        ax.plot(where_center_lng, where_center_lat, 'ko', markersize=2, zorder=12,
                transform=ccrs.PlateCarree())
        ax.coastlines(lw=0.25, color='k', zorder=10) 
#        ax.set_extent([lng[0]-360, lng[-1]-360, lat[0], lat[-1]])
#        if (lng[0]==0.) and (lng[-1]==360.): 
        ax.set_extent([-180., 180., 0., 85.]) 
#        else: 
#            ax.set_extent([lng[0]-360, lng[-1]-360, lat[0], lat[-1]])
        plt.gcf().subplots_adjust(right=0.8)
        plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
            'identify_SLPcenter_%s.eps' %fstr, dpi = 300)
    print('%s locations for %d-%d calculated in %.2f seconds!'%(
        kind.capitalize(), years[0], years[-1], time.time() - start_time)) 
    # OPTIONAL: Plot maps showing the pressure centers on given days (needs
    # mtime as an argument) 
    #import matplotlib.pyplot as plt
    #import cartopy.crs as ccrs
    #import cartopy.feature as cfeature        
    #day = 0
    #from datetime import datetime 
    #plt.figure(figsize=(9,4))
    #ax = plt.subplot2grid((1,1), (0,0), projection=ccrs.PlateCarree())
    #ax.set_title('%s' %(datetime.strftime(mtime[day], '%m/%d/%Y')))
    #ax.set_extent([lng[0]-360, lng[-1]-360, lat[0], lat[-1]])
    ## Contours of SLP 
    #mb = ax.contourf(lng,lat, SLP[day]/100., cmap=plt.get_cmap('rainbow'), 
    #                 transform=ccrs.PlateCarree())
    #colorbar = plt.colorbar(mb)
    #colorbar.ax.tick_params(labelsize=12)
    #colorbar.set_label('Pressure [hPa]', fontsize=14)
    #ax.coastlines(lw=0.25, color='k')    
    ## Locations of identified cyclones/anticyclones
    #lng_grid, lat_grid = np.meshgrid(lng, lat)
    #ax.scatter(lng_grid, lat_grid, s=cyclones[day]*50, facecolor='k', 
    #            lw=0, marker='.', transform=ccrs.PlateCarree())
    #ax.scatter(lng_grid, lat_grid, s=anticyclones[day]*50, facecolor='k', 
    #            lw=0, marker='.', transform=ccrs.PlateCarree())
    #plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
    #            'identify_SLPcenter_%s.eps'
    #            %datetime.strftime(mtime[day], '%m-%d-%Y'))    
    return (center, where_center_daycoord, where_center_ycoord, 
        where_center_xcoord)
    
def filter_center_bylat(center, lat_thresh, lat_gmi):
    """filter (anti)cyclones by latitude, separating them into systems at 
    low latitudes (high latitudes) such that any systems above (below)
    the specified latitude threshold are set to NaN. 
    
    Parameters
    ----------
    center : numpy.ndarray 
        A value of 1 indicates the presence of a cyclone/anticylone for a 
        particular day and locaditon, [time, lat, lng]       
    lat_gmi : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
        
    Returns
    -------
    center_lowlat : numpy.ndarray
        A value of 1 indicates the presence of a cyclone/anticylone for a 
        particular day and location with all cyclones/anticyclones poleward of 
        the specified latitude threshold set to NaN, [time, lat, lng]       
    center_highlat : numpy.ndarray
        A value of 1 indicates the presence of a cyclone/anticylone for a 
        particular day and location with all cyclones/anticyclones equatorward
        of the specified latitude threshold set to NaN, [time, lat, lng]           
    """
    import numpy as np
    import copy
    center_lowlat = copy.deepcopy(center)
    center_highlat = copy.deepcopy(center)
    # Find index of latitude threshold
    lat_thresh_idx = np.abs(lat_gmi-lat_thresh).argmin()
    # Set all pressure centers (which are denoted with 1s above or below the 
    # latitude threshold to NaNs).
    center_lowlat[:, lat_thresh_idx:] = np.nan
    center_highlat[:, :lat_thresh_idx] = np.nan
    return center_lowlat, center_highlat    

def filter_center_byjet(center, jet, lat_systemcoords, lng_systemcoords, 
    lat_jetcoords, lng_jetcoords): 
    """filter (anti)cyclones by their latitude with respect to the eddy-driven 
    jet.
    
    Parameters
    ----------
    center : numpy.ndarray 
        A value of 1 indicates the presence of a cyclone/anticylone for a 
        particular day and locaditon, [time, lat_systemcoords, 
        lng_systemcoords]      
    jet : numpy.ndarray 
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in region, units of degrees north[time, lng_jetcoords]        
    lat_systemcoords : numpy.ndarray 
        Latitude coordinates corresponding to the (anti)cyclone array, units of 
        degrees north, [lat_systemcoords,]        
    lng_systemcoords : numpy.ndarray 
        Longitude coordinates corresponding to the (anti)cyclone array, units of 
        degrees east, [lng_systemcoords,]
    lat_jetcoords : numpy.ndarray
        Latitude coordinates corresponding to the jet array, units of degrees
        north, [lat_jetcoords,]
    lng_jetcoords : numpy.ndarray
        Longitude coordinates corresponding to the jet array, units of degrees 
        east, [lng_jetcoords,]
    lat_jet : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in region, units of degrees north[time, lng_jetcoords]        
        
    Returns
    -------
    abovejet : numpy.ndarray
        A value of 1 indicates the presence of an (anti)cyclone above the eddy-
        driven jet on the day which the (anti)cyclone occurs, [time, 
        lat_systemcoords, lng_systemcoords]       
    belowjet : numpy.ndarray
        A value of 1 indicates the presence of an (anti)cyclone below the eddy-
        driven jet on the day which the (anti)cyclone occurs, [time, 
        lat_systemcoords, lng_systemcoords]
    """
    import numpy as np
    # Find indicies corresponding to positions (with respect to space and time)
    # of (anti)cyclones    
    where_system = np.where(center == 1)
    # Positions of systems above/below jet
    abovejet = np.empty(shape=center.shape)
    abovejet[:] = np.nan
    belowjet = np.empty(shape=center.shape)
    belowjet[:] = np.nan
    # Loop through all systems
    for systemday, systemlat_idx, systemlng_idx in zip(where_system[0],
                                                       where_system[1],
                                                       where_system[2]): 
        # Find longitude index within the jet array coordinates closest to 
        # the system 
        system_lng_jetspace = np.abs(lng_jetcoords-
                                     lng_systemcoords[systemlng_idx]).argmin()    
        # Find latitude of jet at longitude/on day of interest
        jet_lat = jet[systemday, system_lng_jetspace]
        # Add 1 to array denoting the position of systems above (below) jet if 
        # the latitude of the system is greater than (less than) 0
        if (lat_systemcoords[systemlat_idx]-jet_lat) > 0:
            abovejet[systemday, systemlat_idx, systemlng_idx] = 1.
        else: 
            belowjet[systemday, systemlat_idx, systemlng_idx] = 1.
    # OPTIONAL: to check to see if function works
    #import matplotlib.pyplot as plt
    #for d in np.arange(0, 276, 10):    
    #    plt.plot(lng_jetcoords, jet[d], 'ko')   
    #    plt.plot(lng_systemcoords[np.where(abovejet[d]==1.)[1]],
    #             lat_systemcoords[np.where(abovejet[d]==1.)[0]], 'ro')
    #    plt.show()            
    return abovejet, belowjet