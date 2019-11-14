#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module calculates statistical analysis of GMI output and MERRA-2 and emissions
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
    11062019 -- function 'find_grid_in_bb' edited to handle output from 
                GEOS-C1SD
    26062019 -- remove function 'filter_center_bylat' - this simple approach is 
                not needed if examining (anti)cyclones with respect to the jet
    03072019 -- edit function 'find_grid_overland' such that basemap is larger
                than domain to pick out land at domain edges
    12072019 -- edited 'calculate_fieldjet_relationship' description to make it
                clear that the difference was calculated as the difference 
                between the jet latitude and the latitude of the grid cell of 
                interest and to handle any field (i.e., tracer, gas, 
                meteorological field)
    17072019 -- changed maximum 500 hPa wind locator in 'find_field_atjet' to 
                reflect possible NaNs in 500 hPa (on 18/1/2008 there is a NaN
                in the U500 wind - this created issues with DJF analysis)
    22072019 -- function 'calculate_schnell_do3dt_rto3' added
    23072019 -- function 'calculate_obs_do3dt_rto3' added
    25072019 -- function 'ctm_obs_bias' added
    01082019 -- function 'calculate_obs_do3dt_rto3' changed to 
                'calculate_obs_o3_temp_jet' and code added to calculate O3-jet
                relationship from observed O3 and jet (from reanalysis)
    12082019 -- function 'calculate_r_significance' added
    18082019 -- function 'calculate_fieldjet_relationship' edited to return the
                daily difference in the latitude of the jet and the latitude of
                each grid cell for use to determine significance of O3-T-jet 
                distance calculations
    20082019 -- function 'field_binner' added
    25082019 -- function 'segregate_cyclones_bylat' added
    29082019 -- function 'field_binner' edited to bin not just observations 
                over the entire measuring period but on daily timescales
    07102019 -- functions 'sortfield_byjetlat' and 'sortfield_byjetlat_column' 
                added
    15102019 -- function 'bin_observations_bylat' added
    21102019 -- function 'calculate_aqpi' added
    27102019 -- function 'convolve_jet' added
    13112019 -- functions 'calculate_initial_compass_bearing' and 
                'o3anom_cyclone' added
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
    np.seterr(invalid='ignore')
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
    np.seterr(invalid='ignore')    
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
    elif (lat_dim == 1) and (lng_dim == 0):
        outgrid = ingrid[left:right+1, down:up+1]        
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
    bm = Basemap(llcrnrlon=lng.min()-5,llcrnrlat=lat.min()-5,
                 urcrnrlon=lng.max()+5,urcrnrlat=lat.max()+5)
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
    pressure using the ratio of vapor pressure to saturation vapor pressure.
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
    remaining values are filled with NaNs. If anom=True, then the anomaly field
    is found (here anomaly means that grid cells about the jet have their 
    grid cell climatology mean subtracted from them).
    
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
    import time
    start_time = time.time()    
    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
          'Finding field anomaly about the eddy-driven jet...')   
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
            U500max = np.where(U500_transect==np.nanmax(U500_transect))[0][0]   
            # Find latitude of jet
            lat_jet[day, i] = lat_fr[U500max]
    del i
    # If field of interest has dimensions [time, lat, lng] (i.e., ndims=3) the
    # output will be the field north/south of jet for every timestep and will 
    # thus have dimensions of[time, 2*jetdistance+1, lng]. However, if the 
    # field of interest is already time-averaged and has dimensions [lat, 
    # lng] (i.e., ndims=2), the output will only have dimensions 
    # [2*jetdistance+1, lng]. For anom==True, then values about the jet will 
    # their climatological means (at each grid cell) subtracted from them
    # Two dimensions
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
            # If the jet is in an extreme equatorward position, the operation
            # U500max-jetdistance will yield a negative index value and 
            # the following operation will yield no values
            if (U500max-jetdistance) < 0:
                field_transect_below = field[0:U500max, i]                    
            else:  
                field_transect_below = field[U500max-jetdistance:U500max, i]        
            # Fill output grid
            field_jet[jetdistance, i] = field_jetcenter
            field_jet[jetdistance+1:jetdistance+1+len(field_transect_above), 
                      i] = field_transect_above            
            field_jet[jetdistance-len(field_transect_below):jetdistance, 
                      i] = field_transect_below
    # Three dimensions
    if np.ndim(field) == 3:
        if anom is False:
            field_jet = np.empty(shape=(len(field), 2*jetdistance+1, 
                len(lng_fr)))
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
                    field_transect_below = field[day, 
                        U500max-jetdistance:U500max, i]        
                    # Fill output grid
                    field_jet[day, jetdistance, i] = field_jetcenter
                    field_jet[day, jetdistance+1:jetdistance+1+len(
                              field_transect_above), i] = field_transect_above            
                    field_jet[day, jetdistance-len(field_transect_below):
                              jetdistance, i] = field_transect_below
        if anom is True: 
            field_jet = np.empty(shape=(len(field), 2*jetdistance+1, 
                                        len(lng_fr)))
            field_jet[:] = np.nan
            # Loop through time
            for day in np.arange(0, len(field), 1):  
                U500_day = U500_fr[day]
                # Loop through longitude
                for i in np.arange(0, len(lng_fr), 1):    
                    U500max = (np.abs(lat_fr-lat_jet[day, i])).argmin()
                    # Find value of field at the jet's center on the day of 
                    # interest and the climatology of grid cell over all days 
                    # in measuring period
                    field_jetcenter = field[day, U500max, i]       
                    field_jetcenter_clim = np.nanmean(field[:, U500max, i])
                    # Find field above/below jet's center and climatology
                    field_transect_above = field[day, U500max+1:U500max+
                        jetdistance+1, i]
                    field_transect_above_clim = np.nanmean(field[:, 
                        U500max+1:U500max+jetdistance+1, i], axis=0)                
                    field_transect_below = field[day, U500max-jetdistance:
                        U500max, i]
                    field_transect_below_clim = np.nanmean(field[:, 
                        U500max-jetdistance:U500max, i], axis=0)    
                    # Fill output grid
                    field_jet[day, jetdistance, i] = (field_jetcenter-
                        field_jetcenter_clim)
                    field_jet[day, jetdistance+1:jetdistance+1+len(
                        field_transect_above), i] = (field_transect_above-
                        field_transect_above_clim)
                    field_jet[day, jetdistance-len(field_transect_below):
                        jetdistance, i] = (field_transect_below-
                        field_transect_below_clim)
    print('Values of field about jet calculated in %.2f seconds!'%(
            time.time() - start_time))      
    return lat_jet, field_jet

def convolve_jet(lat_jet, w):
    """Perform moving average on daily jet latitude ("boxcar smoothing"). This 
    process replaces each data value with the average of the W neighboring 
    values. In equation form, the moving average is calculated by: 
        \hat{x}[i] = 1/(2*M+1) * sum from -M to +M of x[i+j]
    This is done by convolving the input data with a box-shaped pulse of 2*M+1 
    values all equal to 1/(2*M+1). See https://waterprogramming.wordpress.com/
    # 2018/09/04/implementation-of-the-moving-average-filter-using-convolution/ 
    
    Parameters
    ----------
    lat_jet : numpy.ndarray
        The latitude of the eddy-driven jet, units of degrees north, [time, 
        lng]
    w : int
        Window size for smoothing; n.b., choosing the width of the moving 
        window is not trivial, and if the window is too big, the trends in the 
        data may be masked, while too small of a window may not show the larger 
        underlying trend
        
    Returns
    -------
    lat_jet_convolve : numpy.ndarray
        The convolved (smoothed) latitude of the eddy-driven jet, units of 
        degrees north, [time, lng]
    """
    import numpy as np
    from astropy.convolution import convolve    
    lat_jet_convolve = np.empty(lat_jet.shape)
    lat_jet_convolve[:] = np.nan
    # Loop through days
    for day in np.arange(0, len(lat_jet), 1): 
        series = lat_jet[day]
        # Define mask and store as an array
        mask = np.ones((1,w))/w
        mask = mask[0,:]
        # Convolve the mask with the raw data
        convolved_data = convolve(series, mask, boundary='extend') 
        # Save convolved data
        lat_jet_convolve[day] = convolved_data
    return lat_jet_convolve

def calculate_fieldjet_relationship(field_fr, lat_fr, lng_fr, jetpos, lng_jet):
    """Given a field of interest (i.e., O3, tracer, temperature) and the jet's 
    latitude in a focus region, function determines the relationship between 
    the field and distance from the jet at each grid cell. This is quantified 
    in terms of the slope (i.e., changea in the field per degree shift in the 
    jet stream) and the correlation (i.e., the Pearson product-moment 
    correlation coefficient calculated between the field and a grid cell's 
    distance from the jet). Positive values for the correlation and the slope 
    imply that poleward movement of the jet relative to the grid cell of 
    interest increases value of field (due to the convention that the 
    distance from the jet is defined as the latitude of the jet minus the 
    latitud of the grid cell of interest).
    
    Parameters
    ----------
    field_fr : numpy.ndarray
        Field of interest in the focus region (i.e., Northern Hemisphere), 
        [time, lat, lng]
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
        The slope of the linear regression of the field versus a grid cell's 
        distance from the jet, [lat, lng]
    correl : numpy.ndarray
        The Pearson product-moment correlation coefficient calculated between 
        the field and a grid cell's distance from the jet, [lat, lng]
    diff : numpy.ndarray
        The difference between the latitude of the jet and the latitude of 
        each grid cell, units of degrees, [time, lat, lng]
    """
    import numpy as np
    # Array slope is the slope of the linear regression of the field and the
    # grid cell's distance from the eddy-driven jet (positive distance is a 
    # jet north of grid cell); array correl is the correlation between the two 
    # times series; array diff is the difference in each grid cell's distance 
    # from the jet
    slope = np.empty(shape=(lat_fr.shape[0],lng_fr.shape[0]))
    slope[:] = np.nan
    correl = np.empty(shape=(lat_fr.shape[0],lng_fr.shape[0]))
    correl[:] = np.nan
    diff = np.empty(shape=(len(field_fr), lat_fr.shape[0],
                           lng_fr.shape[0]))
    diff[:] = np.nan
    # Loop through latitudes
    for loci in np.arange(0,len(lat_fr),1):
        # Loop through longitudes
        for locj in np.arange(0,len(lng_fr),1):
            fieldij = field_fr[:, loci, locj]
            # Latitude of eddy-driven jet at a given grid cell, found by finding 
            # a timeseries of the jet location at the longitude of interest
            jetj = jetpos[:, np.abs(lng_jet-lng_fr[locj]).argmin()]
            # Difference in grid cell's latitude and the latitude of the jet
            diffij = jetj-lat_fr[loci]
            # Add timeseries difference in grid cell
            diff[:,loci,locj] = diffij
            # Mask out NaNs to avoid ValueError; i.e., for additional 
            # information see https://stackoverflow.com/questions/13675912/
            # python-programming-numpy-polyfit-saying-nan
            notnan = np.isfinite(diffij) & np.isfinite(fieldij)
            if True in notnan:
                slope[loci,locj] = np.polyfit(diffij[notnan], 
                     fieldij[notnan], 1)[0]
                correl[loci,locj] = np.corrcoef(diffij[notnan], 
                      fieldij[notnan])[0,1]
    return slope, correl, diff

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

#def identify_SLPcenter(lat, lng, SLP, dx, dy, kind, pr_crit, years, 
#    checkplot='no', fstr=''):
#    """function identifies the centers of cyclones or anticyclones from sea 
#    level pressure fields (SLP) over the specified domain. The center of a 
#    cyclone (anticyclone) is defined as a local pressure minima (maxima). In 
#    order for this to be true for a cyclone (anticyclone), the pressure at a 
#    particular grid cell must be less (more) than the dx/dy grid cells on all
#    sides of it and the pressure at the grid cell must be less (more) than the 
#    criteria pressure specified in variable 'pr_crit.' If checkplot='yes', a
#    map showing the locations of all (anti)cyclones and a binned 2D histogram
#    indicating their frequency is plotted.)
#
#    Parameters
#    ----------
#    lat : numpy.ndarray
#        Latitude coordinates for SLP, units of degrees north, [lat,]
#    lng : numpy.ndarray
#        Longitude coordinates for SLP, units of degrees east, [lng,]
#    SLP : numpy.ndarray
#        MERRA-2 Sea level pressure (SLP), units of Pa, [time, lat, lng]
#    dx : int
#        The search area in the x (longitudinal) direction over which pressure 
#        values are examined to identify minima/maxima
#    dy : int
#        The search area in the y (latitudinal) direction over which pressure 
#        values are examined to identify minima/maxima    
#    kind : str
#        cyclone or anticyclone
#    pr_crit : float
#        The critical pressure above (below) which a SLP minima (maxima) cannot
#        be considered a cyclone (anticyclone); n.b. I have been using 1000 hPa
#        (100000 Pa) and 1012 hPa (101200 Pa) as the critical pressures for 
#        cyclones and anticylones, respectively, units of Pa
#    years : list
#        Years in measuring period
#    checkplot : str
#        If 'yes' a map indicating the frequency/location of cyclones/
#        anticyclones is plotted
#    fstr : str
#        Output filename suffix (should indicate whether plot shows cyclone
#        or anticylone)
#        
#    Returns
#    -------
#    center : numpy.ndarray 
#        A value of 1 indicates the presence of a cyclone/anticylone for a 
#        particular day and locaditon, [time, lat, lng]   
#    where_center_daycoord : list
#        The dates (day indices) of cyclones/anticyclones        
#    where_center_ycoord : list 
#        The y-coordinates (latitude indices) of cyclones/anticyclones
#    where_center_xcoord : list 
#        The x-coordinates (longitude indices) of cyclones/anticyclones        
#    """
#    import time 
#    start_time = time.time()    
#    print('# # # # # # # # # # # # # # # # # # # # # # # # # # \n'+
#          'Calculating %s locations...'%kind)   
#    import numpy as np
#    center = np.empty(shape=SLP.shape)
#    center[:] = np.nan    
#    for day in np.arange(0, len(SLP), 1):
#        SLP_day = SLP[day]
#        # For given day, index pressure grid. The outer rows/columns of the 
#        # grid since the neighboring eight grid points cannot be searched for 
#        # a pressure minima; here i is the latitude index and j is the 
#        # longitude index
#        for i in np.arange(dy, len(lat)-dy, 1):
#            for j in np.arange(dx, len(lng)-dx, 1):
#                # SLP at point of interest and SLP at the eight grid cells
#                # surrounding that grid cell
#                SLP_day_atpoint = SLP_day[i, j]
#                # Similar to Lang & Waugh (2011), each of the surrounding grid
#                # points within a search radius defined by dx and dy (not 
#                # counting the center grid point) must be greater than (for 
#                # cyclones) or less than (for anticyclones) the pressure in the 
#                # center)
#                SLP_day_surround = SLP_day[i-dy:i+dy+1, j-dx:j+dx+1]
#                # For cyclones
#                if kind=='cyclone':
#                    if (np.min(SLP_day_surround) >= SLP_day_atpoint) & \
#                    (SLP_day_atpoint < pr_crit):
#                        center[day,i,j] = 1.
#                # For anticyclones
#                else: 
#                    if (np.max(SLP_day_surround) <= SLP_day_atpoint) & \
#                    (SLP_day_atpoint > pr_crit):
#                        center[day,i,j] = 1.       
#    where_center_ycoord, where_center_xcoord = [], []
#    where_center_daycoord = []
#    where_center_lat, where_center_lng = [], []
#    # Loop through days in measuring period
#    for day in np.arange(0, len(center), 1):
#        # Find position(s) of (anti)cyclone on day of interest
#        where_center_day = center[day]
#        where_center_day = np.where(where_center_day==1.)
#        # Append latitude/longitude indices/coordinates to list containing
#        # values for all days
#        where_center_ycoord.append(where_center_day[0])
#        where_center_xcoord.append(where_center_day[1])      
#        where_center_daycoord.append(np.repeat(day, len(where_center_day[0])))
#        where_center_lat.append(lat[where_center_day[0]])
#        where_center_lng.append(lng[where_center_day[1]])
#    where_center_daycoord = np.hstack(where_center_daycoord)
#    where_center_ycoord = np.hstack(where_center_ycoord)
#    where_center_xcoord = np.hstack(where_center_xcoord)
#    where_center_lat = np.hstack(where_center_lat)
#    where_center_lng = np.hstack(where_center_lng)
#    if checkplot=='yes':
#        import matplotlib.pyplot as plt
#        import cartopy.util
#        import cartopy.crs as ccrs
#        import cartopy.feature as cfeature  
#        fig = plt.figure(figsize=(9,4))
#        ax = plt.subplot2grid((1,1), (0,0), projection=ccrs.PlateCarree())
#        # Do coordinate conversion of (x,y)
#        xynps = ax.projection.transform_points(ccrs.PlateCarree(), 
#            where_center_lng, where_center_lat)
#        # Make 2D histogram; n.b., h:(counts, xedges, yedges, image)
#        h = ax.hist2d(xynps[:,0], xynps[:,1], bins=15, zorder=10, alpha=1., 
#            vmin=0, vmax=7, cmap=plt.get_cmap('Blues'), 
#            transform=ccrs.PlateCarree())
#        colorbar_axes = plt.gcf().add_axes([0.83,0.25,0.02,0.5])
#        colorbar = plt.colorbar(h[3], colorbar_axes, orientation='vertical',
#                                extend='max')
#        colorbar.ax.tick_params(labelsize=12)
#        colorbar.set_label('Frequency', fontsize=14)
#        # Add scatterpoints for all 
#        ax.plot(where_center_lng, where_center_lat, 'ko', markersize=2, zorder=12,
#                transform=ccrs.PlateCarree())
#        ax.coastlines(lw=0.25, color='k', zorder=10) 
##        ax.set_extent([lng[0]-360, lng[-1]-360, lat[0], lat[-1]])
##        if (lng[0]==0.) and (lng[-1]==360.): 
#        ax.set_extent([-180., 180., 0., 85.]) 
##        else: 
##            ax.set_extent([lng[0]-360, lng[-1]-360, lat[0], lat[-1]])
#        plt.gcf().subplots_adjust(right=0.8)
#        plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#            'identify_SLPcenter_%s.eps' %fstr, dpi = 300)
#    print('%s locations for %d-%d calculated in %.2f seconds!'%(
#        kind.capitalize(), years[0], years[-1], time.time() - start_time)) 
#    # OPTIONAL: Plot maps showing the pressure centers on given days (needs
#    # mtime as an argument) 
#    #import matplotlib.pyplot as plt
#    #import cartopy.crs as ccrs
#    #import cartopy.feature as cfeature        
#    #day = 0
#    #from datetime import datetime 
#    #plt.figure(figsize=(9,4))
#    #ax = plt.subplot2grid((1,1), (0,0), projection=ccrs.PlateCarree())
#    #ax.set_title('%s' %(datetime.strftime(mtime[day], '%m/%d/%Y')))
#    #ax.set_extent([lng[0]-360, lng[-1]-360, lat[0], lat[-1]])
#    ## Contours of SLP 
#    #mb = ax.contourf(lng,lat, SLP[day]/100., cmap=plt.get_cmap('rainbow'), 
#    #                 transform=ccrs.PlateCarree())
#    #colorbar = plt.colorbar(mb)
#    #colorbar.ax.tick_params(labelsize=12)
#    #colorbar.set_label('Pressure [hPa]', fontsize=14)
#    #ax.coastlines(lw=0.25, color='k')    
#    ## Locations of identified cyclones/anticyclones
#    #lng_grid, lat_grid = np.meshgrid(lng, lat)
#    #ax.scatter(lng_grid, lat_grid, s=cyclones[day]*50, facecolor='k', 
#    #            lw=0, marker='.', transform=ccrs.PlateCarree())
#    #ax.scatter(lng_grid, lat_grid, s=anticyclones[day]*50, facecolor='k', 
#    #            lw=0, marker='.', transform=ccrs.PlateCarree())
#    #plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#    #            'identify_SLPcenter_%s.eps'
#    #            %datetime.strftime(mtime[day], '%m-%d-%Y'))    
#    return (center, where_center_daycoord, where_center_ycoord, 
#        where_center_xcoord)
    
#def filter_center_byjet(center, jet, lat_systemcoords, lng_systemcoords, 
#    lat_jetcoords, lng_jetcoords): 
#    """filter (anti)cyclones by their latitude with respect to the eddy-driven 
#    jet.
#    
#    Parameters
#    ----------
#    center : numpy.ndarray 
#        A value of 1 indicates the presence of a cyclone/anticylone for a 
#        particular day and locaditon, [time, lat_systemcoords, 
#        lng_systemcoords]      
#    jet : numpy.ndarray 
#        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
#        in region, units of degrees north[time, lng_jetcoords]        
#    lat_systemcoords : numpy.ndarray 
#        Latitude coordinates corresponding to the (anti)cyclone array, units of 
#        degrees north, [lat_systemcoords,]        
#    lng_systemcoords : numpy.ndarray 
#        Longitude coordinates corresponding to the (anti)cyclone array, units of 
#        degrees east, [lng_systemcoords,]
#    lat_jetcoords : numpy.ndarray
#        Latitude coordinates corresponding to the jet array, units of degrees
#        north, [lat_jetcoords,]
#    lng_jetcoords : numpy.ndarray
#        Longitude coordinates corresponding to the jet array, units of degrees 
#        east, [lng_jetcoords,]
#    lat_jet : numpy.ndarray
#        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
#        in region, units of degrees north[time, lng_jetcoords]        
#        
#    Returns
#    -------
#    abovejet : numpy.ndarray
#        A value of 1 indicates the presence of an (anti)cyclone above the eddy-
#        driven jet on the day which the (anti)cyclone occurs, [time, 
#        lat_systemcoords, lng_systemcoords]       
#    belowjet : numpy.ndarray
#        A value of 1 indicates the presence of an (anti)cyclone below the eddy-
#        driven jet on the day which the (anti)cyclone occurs, [time, 
#        lat_systemcoords, lng_systemcoords]
#    """
#    import numpy as np
#    # Find indicies corresponding to positions (with respect to space and time)
#    # of (anti)cyclones    
#    where_system = np.where(center == 1)
#    # Positions of systems above/below jet
#    abovejet = np.empty(shape=center.shape)
#    abovejet[:] = np.nan
#    belowjet = np.empty(shape=center.shape)
#    belowjet[:] = np.nan
#    # Loop through all systems
#    for systemday, systemlat_idx, systemlng_idx in zip(where_system[0],
#                                                       where_system[1],
#                                                       where_system[2]): 
#        # Find longitude index within the jet array coordinates closest to 
#        # the system 
#        system_lng_jetspace = np.abs(lng_jetcoords-
#                                     lng_systemcoords[systemlng_idx]).argmin()    
#        # Find latitude of jet at longitude/on day of interest
#        jet_lat = jet[systemday, system_lng_jetspace]
#        # Add 1 to array denoting the position of systems above (below) jet if 
#        # the latitude of the system is greater than (less than) 0
#        if (lat_systemcoords[systemlat_idx]-jet_lat) > 0:
#            abovejet[systemday, systemlat_idx, systemlng_idx] = 1.
#        else: 
#            belowjet[systemday, systemlat_idx, systemlng_idx] = 1.
#    # OPTIONAL: to check to see if function works
#    #import matplotlib.pyplot as plt
#    #for d in np.arange(0, 276, 10):    
#    #    plt.plot(lng_jetcoords, jet[d], 'ko')   
#    #    plt.plot(lng_systemcoords[np.where(abovejet[d]==1.)[1]],
#    #             lat_systemcoords[np.where(abovejet[d]==1.)[0]], 'ro')
#    #    plt.show()            
#    return abovejet, belowjet

def calculate_schnell_do3dt_rto3(years, months, domain): 
    """function opens O3 dataset from Schnell et al. (2014) for specified 
    domain and time period, commensurate 2-meter tempratures from MERRA-2, 
    interpolates to GMI CTM resolution, and calculates dO3/dT and r(O3, T). 
    Note that since the CTM and MERRA data are loaded for the Northern 
    Hemisphere, a portion of the interpolated fields need to be lobbed off. 

    Parameters
    ----------
    years : list
        Year or range of years in measuring period, [years,]
    months : list
        Three letter abbreviations (lowercase) for months in measuring period
    domain : str
        Either 'US' or 'EU'

    Returns\
    -------    
    do3dt2m : numpy.ndarray
        dO3/dT, units of ppbv K-1, [lat, lng]    
    r_t2mo3 : numpy.ndarray
        Pearson correlation coefficient calculated between O3 and 2-meter 
        temperature, [lat, lng]
    lat_gmi : numpy.ndarray
        GMI latitude coordinates corresponding to the Schnell et al. domain, 
        units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI longitude coordinates corresponding to the Schnell et al. domain, 
        units of degrees east, [lng,]

    References
    ----------
    [1] J. L. Schnell, C. D. Holmes, A. Jangam, and M. J. Prather, "Skill in 
    forecasting extreme ozone pollution episodes with a global atmospheric 
    chemistry model," Atmos. Chem. Phys., 14, 7721-7739, 2014.         
    """
    import numpy as np
    import sys
    sys.path.append('/Users/ghkerr/phd/globalo3/')
    import globalo3_open
    latmin_n, lngmin_n, latmax_n, lngmax_n = -1., 0., 90., 360.
    # Open Schnell et al. O3 data for given months, years
    lat_schnell, lng_schnell, o3_schnell = globalo3_open.open_schnello3(years, 
        months, domain)
    # Load 2-meter temperatures 
    lat_merra, lng_merra, t2m_merra = \
        globalo3_open.open_merra2t2m_specifieddomain(years, months, latmin_n, 
        latmax_n, lngmin_n, lngmax_n)
#        globalo3_open.open_merra2t2m_specifieddomain(years, months, 
#        lat_schnell.min(), lat_schnell.max(), lng_schnell.min(), 
#        lng_schnell.max())        
    # Load commensurate GMI O3 data (this is a "dummy" calculation; only the 
    # latitude and longitude coordinates are needed for interpolation)
    lat_gmi, lng_gmi, times, o3 = globalo3_open.open_overpass2_specifieddomain(
        years, months, latmin_n, latmax_n, lngmin_n, lngmax_n, 'O3', 
        'HindcastMR2')  
    if domain == 'eu':
        t2m_merra = np.roll(t2m_merra, 
            len(lng_merra)-np.abs(lng_merra-lng_schnell[0]).argmin())
        lng_merra = np.roll(lng_merra, 
            len(lng_merra)-np.abs(lng_merra-lng_schnell[0]).argmin())
        lng_gmi = np.roll(lng_gmi, 
            len(lng_gmi)-np.abs(lng_gmi-lng_schnell[0]).argmin())
    # Degrade O3 and temperature datasets to resolution of CTM 
    o3_schnell = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
        lng_gmi, lat_schnell, lng_schnell, o3_schnell)
    t2m_merra = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
        lng_gmi, lat_merra, lng_merra, t2m_merra)
    # Calculate r(O3, T) and dO3/dT
    do3dt2m = globalo3_calculate.calculate_do3dt(t2m_merra, o3_schnell, 
        lat_gmi, lng_gmi)
    r_t2mo3 = globalo3_calculate.calculate_r(t2m_merra, o3_schnell, 
        lat_gmi, lng_gmi)
    return do3dt2m, r_t2mo3, lat_gmi, lng_gmi

def calculate_obs_o3_temp_jet(obs, t2m, times_t2m, lat_t2m, lng_t2m, lat_jet, 
    lng_ml): 
    """function calculates the O3-temperature relationship from observational 
    O3 datasets and MERRA-2 2-meter temperatures. This function assumes that 
    the time coordinates for MERRA-2 are continuous. 
    
    Parameters
    ----------
    obs : pandas.core.frame.DataFrame
        Daily averaged O3 concentrations at individual observational sites for 
        the specified measuring period. DataFrame contains station ID and 
        station latitude and longitude
    t2m : numpy.ndarray
        Daily averaged MERRA-2 2-meter temperatures, units of K, [days, lat, 
        lng]
    times_t2m : pandas.core.indexes.datetimes.DatetimeIndex
        Timestamps of dates for which MERRA-2 data is desired, [days,]
    lat_t2m : numpy.ndarray
        MERRA-2 latitude coordinates, units of degrees north, [lat,]
    lng_t2m : numpy.ndarray
        MERRA-2 longitude coordinates, units of degrees east, [lng,]
    lat_jet : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in region, units of degrees north[days, lng]    
    lng_ml : numpy.ndarray
        Longitude coordinates corresponding to the jet dataset, units of 
        degrees east, [lng,]
        
    Returns
    -------
    r_all : list
        r(T, O3) at each observational station, [stations,]
    do3dt_all : list 
        dO3/dT at each observational station, units of ppbv K-1, [stations,]
    ro3jet : list
        r(jet lat - lat, O3) at each observational station, [stations,]
    djetdo3 : list
        dO3/d(jet lat - lat) at each observational station, units of ppbv 
        degree-1, [stations,]
    lat_all : list 
        Latitude coorindate at each observational stations, units of degrees 
        north, [stations,]
    lng_all : list 
        Longitude coorindate at observational stations, units of degrees east, 
        [stations,]
    """
    import time
    start_time = time.time()    
    import numpy as np
    import sys
    sys.path.append('/Users/ghkerr/phd/utils/')
    from geo_idx import geo_idx
    # Convert
    times_t2m = [x.strftime('%Y-%m-%d') for x in times_t2m]
    r_all, do3dt_all, lat_all, lng_all = [], [], [], [] 
    ro3jet, djetdo3 = [], []    
    # Loop through all stations represented in observaitons
    for station_id in np.unique(obs['Station ID'].values):
        obs_station = obs.loc[obs['Station ID'] == station_id]
        # Only calculate r(T, O3) and dO3/dT for stations with data 
        if (obs_station.shape[0] != 0) \
            & (True in np.isfinite(obs_station['Latitude'].values)) \
            & (True in np.isfinite(obs_station['Longitude'].values)) \
            & (True in np.isfinite(obs_station['O3 (ppbv)'].values)):
            station_lat = np.nanmean(obs_station['Latitude'].values)
            station_lng = np.nanmean(obs_station['Longitude'].values) % 360
            # Find closest MERRA-2 grid cell
            lat_closest = geo_idx(station_lat, lat_t2m)
            lng_closest = geo_idx(station_lng, lng_t2m)   
            # Find closest grid cell in jet longitude dataset 
            lng_jet_closest = geo_idx(station_lng, lng_ml)   
            # Timeseries of jet latitude for a given longitude
            lat_jet_closest = lat_jet[:, lng_jet_closest]
            # Find indices of dates in O3 observations included in MERRA data
            times_obs = np.in1d(times_t2m, obs_station.index).nonzero()[0]
            # Select temperature on days with observations
            t2m_inperiod = t2m[times_obs, lat_closest, lng_closest]
            # Calculate r(T, O3) and dO3/dT 
            idx = np.isfinite(t2m_inperiod) \
                & np.isfinite(obs_station['O3 (ppbv)'].values)
            do3dt_all.append(np.polyfit(t2m_inperiod[idx], 
                             obs_station['O3 (ppbv)'].values[idx], deg=1)[0])
            r_all.append(np.corrcoef(t2m_inperiod[idx], 
                             obs_station['O3 (ppbv)'].values[idx])[0,1])
            # Calculate r(jet lat-lat, O3) and dO3/d(jet lat-lat)
            lat_jet_closest = lat_jet_closest[times_obs]
            idx = np.isfinite(lat_jet_closest) \
                & np.isfinite(obs_station['O3 (ppbv)'].values)     
            djetdo3.append(np.polyfit(lat_jet_closest[idx] - station_lat, 
                             obs_station['O3 (ppbv)'].values[idx], deg=1)[0])
            ro3jet.append(np.corrcoef(lat_jet_closest[idx] - station_lat, 
                             obs_station['O3 (ppbv)'].values[idx])[0,1])
            lat_all.append(station_lat)
            lng_all.append(station_lng)
    print('Observational O3, temperature, jet metrics in calculated ' + 
          'in %.2f seconds' %((time.time()-start_time)))            
    return r_all, do3dt_all, ro3jet, djetdo3, lat_all, lng_all

def ctm_obs_bias(obs_lat, obs_lng, obs_field, ctm_lat, ctm_lng, ctm_field): 
    """function calculates mean bias between CTM and observational datasets
    by selecting CTM grid cell containing observational site and determining
    the simple difference (i.e., CTM - observations, so positive values would
    indicate that the CTM is biased high). Note that the input array for the 
    CTM must be 2D, so CTM output with a time dimensions should be averaged 
    before finding the mean bias. 
    
    Parameters
    ----------
    obs_lat : list
        Latitude coordinates of observational stations, units of degrees north, 
        [stations,]
    obs_lng : list
        Longitude coordinates of observational stations, units of degrees east, 
        [stations,]    
    obs_field : list
        Field of interest at each observational station, [stations,]
    ctm_lat : numpy.ndarray
        Gridded latitude coordinates of CTM, units of degrees north, [lat,]
    ctm_lng : numpy.ndarray
        Gridded longitude coordinates of CTM, units of degrees east, [lng,]        
    ctm_field : numpy.ndarray
        Gridded field of interest from CTM, [lat, lng]
        
    Returns
    -------
    bias_all : list
        CTM-observation bias; order of biases corresponds to the order of the 
        input latitude and longitude coordinates for the observational 
        stations, [stations,]
    """
    import sys
    sys.path.append('/Users/ghkerr/phd/utils/')
    from geo_idx import geo_idx
    # Determine and plot mean bias 
    i = 0
    bias_all = []
    for ilat, ilng in zip(obs_lat, obs_lng): 
        # loop through CASTNet sites and find nearest GMI grid cell
        ctm_lat_near = geo_idx(ilat, ctm_lat)
        ctm_lng_near = geo_idx(ilng, ctm_lng)
        # Field at each CTM grid cell 
        ctm_near_obs = ctm_field[ctm_lat_near, ctm_lng_near]
        obs_atsite = obs_field[i]
        bias_atsite = (ctm_near_obs-obs_atsite)
        bias_all.append(bias_atsite)
        i = i + 1
    return bias_all

def calculate_r_significance(x, y, r, lat, lng):
    """function calculates the lag-1 autocorrelation (rho1) of the dependent 
    variable (y), and thereafter calculates the effective sample size, n', 
    where n' = n * ((1-rho1)/(1+rho1)). 
    The critical value for the one-sided t-test is found for the given 
    effective sample size and compared with the test statistic given the 
    observed correlation.

    Parameters
    ----------
    x : numpy.ndarray
        Independent variable [time, lat, lng]
    y : numpy.ndarray
        Dependent variable [time, lat, lng]
    r : numpy.ndarray     
        Pearson correlation coefficient between x and y, [lat, lng]        
    lat : numpy.ndarray
        Gridded latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Gridded longitude coordinates, units of degrees east, [lng,]        
        
    Returns
    -------
    significance : numpy.ndarray
        NaN for significant correlation at the alpha = 0.05 level and 1 for 
        insigificant correlation at the alpha = 0.05 level, [lat, lng]
    """    
    import numpy as np
    from scipy import stats
    import pandas as pd
    np.seterr(invalid='ignore')    
    import time
    start_time = time.time()
    alpha = 0.05
    rho1 = np.empty(y.shape[1:])
    # Lag-1 autocorrelation coefficient
    for i,ilat in enumerate(lat):
        for j,ilng in enumerate(lng):
            rho1[i,j] = pd.Series.autocorr(pd.Series(y[:, i, j]), lag=1)
    # Effective sample size (from Wilks)
    neff = len(y)*((1-rho1)/(1+rho1))
    # Determine the appropriate t value 
    # (from http://janda.org/c10/Lectures/topic06/L24-significanceR.htm)
    t = r*np.sqrt((neff-2.)/(1.-r**2))
    # Get the critical t test statistic (adapted from 
    # https://stackoverflow.com/questions/19339305/python-function-to-get-the-t-statistic)
    cv = np.empty(y.shape[1:])
    # Second argument is the degrees of freedom for entering 
    # the t-distribution is N - 2
    for i,ilat in enumerate(lat):
        for j,ilng in enumerate(lng):
            cv[i,j] = stats.t.ppf(1-alpha, neff[i,j]-2.)
    # For the given degrees of freedom for a one-tailed test, see which 
    # grid nodes lie below the critical value. If this is the case, the null
    # hypothesis of no relationship (r = 0) cannot be rejected
    where_insignificant = np.where(np.abs(t) < cv)
    # Create grid for significance where np.nan refers to a grid cell that IS 
    # significant at the alpha = 0.05 level and 1 corresponds to a grid cell 
    # that is NOT significant at that level
    significance = np.empty(y.shape[1:])
    significance[:] = np.nan
    significance[where_insignificant[0], where_insignificant[1]] = 1.
    print('significance of correlation coefficient in calculated '+
          'in %.2f seconds'%((time.time()-start_time)))    
    return significance

def field_binner(lat_grid, lng_grid, times_grid, lat_obs, lng_obs, val_obs,
    times_obs, operation): 
    """for a given field (i.e., O3 observations, cyclone centers, etc.) 
    function bins values into grid boxes of the specified size and calculates
    the mean values, variability, or frequency of the field in each bin. 
    
    Parameters
    ----------    
    lat_grid : numpy.ndarray
        Latitude coordinates to which field/observations will be binned (n.b., 
        for a given set of coordinates, function will search to 1/2 * 
        resolution on either side of each grid node, units of degrees north, 
        [lat,]
    lng_grid : numpy.ndarray
        Longitude coordinates to which field/observations will be binned,
        units of degrees east, [lng]
    times_grid : numpy.ndarray
        datetime.date timestamps corresponding to gridded field, [time,]
    lat_obs : numpy.ndarray
        Latitude coordinates corresponding to observations, units of degrees
        north, [no. obs.,]
    lng_obs : numpy.ndarray
        Longitude coordinates corresponding to observations, units of degrees
        east, [no. obs.,]  
    val_obs : numpy.ndarray
        Observational values, [no. obs.,] (n.b., if the frequency of the field 
        is desired, 'val_field' should be an array of 1s and operation = 'sum')
    times_obs : numpy.ndarray
        Timestamps of field, %Y-%m-%d format
    operation : str
        The operation performed on 'val_field' in every bin ('mean', 'std', 
        or 'sum')
        
    Returns
    -------        
    val_grid : numpy.ndarray
        The specified operation performed over the binned field, [lat_bin, 
        lng_bin,]
    lat_grid : numpy.ndarray
        The latitude array corresponding to the binned field (n.b., same as 
        input field), units of degrees north, [lat_bin,]
    lng_grd : numpy.ndarray
        The longitude array corresponding to the binned field, units of degrees
        east, [lng_bin,]    
    """
    import time
    start_time = time.time()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Binning observations and calculating %s...' %operation)
    import numpy as np
    # Convert timestamps from gridded dataset to YYYY-MM-DD format
    times_grid = [x.strftime('%Y-%m-%d') for x in times_grid]
    times_grid = np.array(times_grid)
    lat_res = np.mean(np.diff(lat_grid)[1:-1])
    lng_res = np.mean(np.diff(lng_grid)[1:-1])
    # "Regridded" (binned coarser) arrays
    val_grid = np.empty(shape=(lat_grid.shape[0], lng_grid.shape[0]))
    val_grid[:] = np.nan
    if operation == 'daily':
        val_grid = np.empty(shape=(len(times_grid), lat_grid.shape[0], 
            lng_grid.shape[0]))
        val_grid[:] = np.nan
    # Loop through coarse latitude and longitude 
    for iidx in np.arange(0, len(lat_grid), 1):
        for jidx in np.arange(0, len(lng_grid), 1):
            # Bounding box 
            left = lng_grid[jidx]-(lng_res/2.)
            right = lng_grid[jidx]+(lng_res/2.) 
            up = lat_grid[iidx]+(lat_res/2.)
            down = lat_grid[iidx]-(lat_res/2.)
            in_bb = np.where((lat_obs > down) & (lat_obs <= up) & 
                             (lng_obs > left) & (lng_obs <= right))[0]
            if np.shape(in_bb)[0] is not 0: 
                if operation == 'mean':
                    val_grid[iidx, jidx] = np.nanmean(val_obs[in_bb])            
                elif operation == 'std':        
                    val_grid[iidx, jidx] = np.nanstd(val_obs[in_bb])
                elif operation == 'sum':        
                    val_grid[iidx, jidx] = np.nansum(val_obs[in_bb])
                elif operation == 'daily':
                    # If cyclone frequency on daily timescales is desired, 
                    # function will loop through all unique times during 
                    # which cyclones occurred in the selected grid cell and 
                    # find the days and occurrences of cyclones
                    for day in np.unique(times_obs[in_bb]):
                        time_idx = np.where(times_grid == day)[0]
                        # Find number of cyclones on day at grid cell
                        nocyclones = np.where(times_obs[in_bb] == day)[0]
                        val_grid[time_idx, iidx, jidx] = len(nocyclones)
    print('Observations binned in %.2f seconds!' %(time.time() - start_time))                
    return val_grid, lat_grid, lng_grid

def segregate_cyclones_bylat(cyclones, field, lng_jet, lat_jet, times): 
    """function identifies cyclone locations on days when the eddy-driven jet
    is "extreme" poleward (> 70th percentile) and equatorward (< 30th 
    percentile). These operations are done locally, meaning that the poleward 
    and equatorward positions of the jet are indentified using the timeseries 
    for each longitudinal band and thereafter finding cyclones within these 
    bands. The position and variability of the jet on these extreme days is 
    also saved. 

    Parameters
    ----------
    cyclones : pandas.core.frame.DataFrame
        DataFrame containing the date, hour, fraction of land cover, latitude, 
        longitude, SLP, and storm ID
    field : numpy.ndarray
        3D (tyx) field that will be parsed on days with "extreme" poleward
        and equatorward jet. 
    lng_jet : numpy.ndarray
        The longitude dataset corresponding to the latitude of the jet, [lng]
    lat_jet : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in region, units of degrees north[time, lng]
    times : numpy.ndarray
        datetime.date objects corresponding to every day in measuring period, 
        [time,]        

    Returns
    -------
    lowthresh_lat_cyclone : numpy.ndarray
        Latitudes of cylones on days where the jet is "extreme" equatorward 
        (n.b., "extreme" refers to < 30th percentile), units of degrees north
    lowthresh_lng_cyclone : numpy.ndarray
        Longitude of cylones on days where the jet is "extreme" equatorward, 
        units of degrees east
    lowthresh_time_cyclone : numpy.ndarray   
        Timestamps (%Y-%m-%d format) of cyclones on days where the jet is 
        "extreme" equatorward
    highthresh_lat_cyclone : numpy.ndarray   
        Latitudes of cylones on days where the jet is "extreme" poleward, 
        (n.b., "extreme" refers to > 70th percentile), units of degrees north    
    highthresh_lng_cyclone : numpy.ndarray   
        Longitude of cylones on days where the jet is "extreme" poleward, 
        units of degrees east    
    highthresh_time_cyclone : numpy.ndarray   
        Timestamps (%Y-%m-%d format) of cyclones on days where the jet is 
        "extreme" poleward
    lowthresh_lat_jet : list
        The mean latitude of the eddy-driven jet on days when it is 
        "extreme" equatorward, [lng]
    lowthresh_lat_jet_var : list
        The variability of the latitude of the eddy-driven jet on days when 
        it is "extreme" equatorward, [lng]    
    highthresh_lat_jet : list
        The mean latitude of the eddy-driven jet on days when it is 
        "extreme" poleward, [lng]    
    highthresh_lat_jet_var : list
        The variability of the latitude of the eddy-driven jet on days when 
        it is "extreme" poleward, [lng]     
    pwjet_field_anom : numpy.ndarray
        The anomaly of the field of interest on days when the jet is "extreme" 
        poleward, [lat, lng]    
    eqjet_field_anom : numpy.ndarray   
        The anomaly of the field of interest on days when the jet is "extreme" 
        equatorward, [lat, lng]    
    """
    import time
    start_time = time.time()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Separating cyclones on days with equator/poleward jet...')
    import numpy as np
    lowthresh_lng_cyclone, lowthresh_lat_cyclone = [], []
    highthresh_lng_cyclone, highthresh_lat_cyclone = [], []
    highthresh_time_cyclone, lowthresh_time_cyclone = [], []
    lowthresh_lat_jet, lowthresh_lat_jet_var = [], []
    highthresh_lat_jet, highthresh_lat_jet_var = [], []
    # Return empty arrays to fill with value of 3D (tyx) field on days with 
    # poleward/equatorward jet 
    pwjet_field_anom = np.empty(shape=field.shape[1:])
    pwjet_field_anom[:] = np.nan
    eqjet_field_anom = np.empty(shape=field.shape[1:])
    eqjet_field_anom[:] = np.nan    
    # Resolution of longitudes corresponding to jet dataset 
    res_lng = np.diff(lng_jet[1:-1]).mean()
    for ilng in np.arange(0, len(lng_jet), 1):
        # Find threshold (latitude) for what qualifies as an "extreme" poleward
        # or equatoward jet at a particular longitude
        highthresh_ilng = np.percentile(lat_jet[:,ilng], 70)
        lowthresh_ilng = np.percentile(lat_jet[:,ilng], 30)
        highthresh_days = np.where(lat_jet[:,ilng] > highthresh_ilng)[0]
        lowthresh_days = np.where(lat_jet[:,ilng] < lowthresh_ilng)[0]
        # Mean values of the field on days when the jet extreme poleward or 
        # equatorward at longitude 
        pwjet_field_anom[:, ilng] = (
            np.nanmean(field[highthresh_days, :, ilng], axis=0))#-
#            np.nanmean(field[:, :, ilng], axis=0))
        eqjet_field_anom[:, ilng] = (
            np.nanmean(field[lowthresh_days, :, ilng], axis=0))#-
#            np.nanmean(field[:, :, ilng], axis=0))    
        # Find mean jet latitude and its variability on days when jet is 
        # extreme poleward or equatorward at longitude of interest
        lowthresh_lat_jet.append(np.nanmean(lat_jet[lowthresh_days,ilng]))
        lowthresh_lat_jet_var.append(np.nanstd(lat_jet[lowthresh_days,ilng]))
        highthresh_lat_jet.append(np.nanmean(lat_jet[highthresh_days,ilng]))
        highthresh_lat_jet_var.append(np.nanstd(lat_jet[highthresh_days,ilng]))
        # Dates ('YYYY-mm-dd' format) when jet is extreme poleward or 
        # equatorward at longitude of interest
        highthresh_days = [x.strftime('%Y-%m-%d') for x in 
            times[highthresh_days]]
        lowthresh_days = [x.strftime('%Y-%m-%d') for x in 
            times[lowthresh_days]]
        # Find cyclones at (near) longitude (here near is defined as the 
        # area within +/- the longitude of interest and half its resolution) on 
        # day of interest
        highthresh_cyclones = cyclones.loc[
            (cyclones['Longitude'] > lng_jet[ilng]-(res_lng/2.)) &
            (cyclones['Longitude'] <= lng_jet[ilng]+(res_lng/2.)) & 
            (cyclones['Date'].isin(highthresh_days))]
        highthresh_lng_cyclone.append(
            highthresh_cyclones['Longitude'].values)
        highthresh_lat_cyclone.append(
            highthresh_cyclones['Latitude'].values)  
        highthresh_time_cyclone.append(highthresh_cyclones['Date'].values)                  
        lowthresh_cyclones = cyclones.loc[
            (cyclones['Longitude'] > lng_jet[ilng]-(res_lng/2.)) &
            (cyclones['Longitude'] <= lng_jet[ilng]+(res_lng/2.)) & 
            (cyclones['Date'].isin(lowthresh_days))]    
        lowthresh_lng_cyclone.append(
            lowthresh_cyclones['Longitude'].values)
        lowthresh_lat_cyclone.append(
            lowthresh_cyclones['Latitude'].values)
        lowthresh_time_cyclone.append(lowthresh_cyclones['Date'].values)
    lowthresh_lat_cyclone = np.hstack(lowthresh_lat_cyclone) 
    lowthresh_lng_cyclone = np.hstack(lowthresh_lng_cyclone) 
    highthresh_lat_cyclone = np.hstack(highthresh_lat_cyclone) 
    highthresh_lng_cyclone = np.hstack(highthresh_lng_cyclone) 
    highthresh_time_cyclone = np.hstack(highthresh_time_cyclone)
    lowthresh_time_cyclone = np.hstack(lowthresh_time_cyclone)
    print('Cyclones separated in %.2f seconds!' %(time.time() - start_time))                    
    return (lowthresh_lat_cyclone, lowthresh_lng_cyclone, 
        lowthresh_time_cyclone, highthresh_lat_cyclone, highthresh_lng_cyclone, 
        highthresh_time_cyclone, lowthresh_lat_jet, 
        lowthresh_lat_jet_var, highthresh_lat_jet, highthresh_lat_jet_var,
        pwjet_field_anom, eqjet_field_anom)
    
def reynolds_decomposition(x, dtime, lat, lng):
    """function computes eddy fluxes through Reynolds decomposition 
    into mean and eddy transports. 
    
    Parameters
    ----------
    x : numpy.ndarray 
        Field of interest, [time, lat, lng] 
    dtime : numpy.ndarray or pandas.core.indexes.datetimes.DatetimeIndex
        Time stamps corresponding to x, [time, lat, lng]
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]       
        
    Returns
    -------
    xbar : numpy.ndarray
        Time-averaged field, [lat, lng]
    xbar_star : numpy.ndarray
        The deviations of the time-averaged field from its zonal average, 
        [lat, lng]
    xbar_zm : numpy.ndarray
        The zonal average of the time-averaged field (zonal average 
        repeated for each longitude), [lat, lng]
    xprime : numpy.ndarray
        Deviations from the time average, [time, lat, lng]

    References
    ----------
    Hartmann, D. L., (2016). Global Physical Climatology (Second Ed). 
    Elsevier: San Diego. 
    """
    import numpy as np
    lngaxis = np.where(np.array(x.shape) == lng.shape[0])[0][0]
    lataxis = np.where(np.array(x.shape) == lat.shape[0])[0][0]
    timeaxis = np.where(np.array(x.shape) == dtime.shape[0])[0][0]
    # Time average
    xbar = np.nanmean(x, axis=timeaxis)
    # Transient component (deviation from time mean)
    xprime = x-xbar
    # Stationary component (deviation from zonal mean of the time average, 
    # given by \bar{x}* = \bar{x} - [\bar{x}], where \bar{x} is the 
    # instantaneous zonal mean field (at each timestep) and [\bar{x}] is the 
    # time-averaged zonal value field.
    # Begin by finding the zonal mean of the time average
    xbar_zm = np.nanmean(xbar, axis=np.where(np.array(xbar.shape)==
        lng.shape[0])[0][0])
    xbar_zm = np.repeat(xbar_zm[:, np.newaxis], len(lng), axis=1)
    xbar_star = xbar - xbar_zm
    return xbar, xbar_star, xbar_zm, xprime

def meridional_flux(v, f, dtime, lat, lng):
    """function calculates the meridional transport of a given field as the 
    sum of contributions from the mean circulation, the stationary eddies, and 
    the transient eddies. 
    
    Parameters
    ----------
    v : numpy.ndarray 
        Northward (meridional) wind field, units of m s-1, [time, lat, lng]
    f : numpy.ndarray 
        Field of interest, [time, lat, lng] 
    dtime : numpy.ndarray or pandas.core.indexes.datetimes.DatetimeIndex
        Time stamps corresponding to x, [time, lat, lng]
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]       
        
    Returns
    -------
    mean : numpy.ndarray
        Northward transport of field due to mean meridional circulation, [lat,]
    stationary : numpy.ndarray
        Northward transport of field due to stationary eddies, [lat,]    
    transient : numpy.ndarray
        Northward transport of field due to transient eddies, [lat,]        
    total : numpy.ndarray
        Total northward transport of field, [lat,]
    """
    import numpy as np
    vbar, vbar_star, vbar_zm, vprime = reynolds_decomposition(v, dtime, lat, 
        lng)
    fbar, fbar_star, fbar_zm, fprime = reynolds_decomposition(f, dtime, lat, 
        lng)
    # From Hartmann (2016), using the definition of time and zonal averages, 
    # the flux of field f by wind u can be written as the sum of contributions 
    # from the mean meridional circulation, stationary eddies, and the 
    # transient eddyes
    mean = np.nanmean(vbar, axis=1)*np.nanmean(fbar, axis=1)
    stationary = np.nanmean((vbar_star*fbar_star), axis=1)
    transient = np.nanmean((vprime*fprime), axis=tuple((0,2)))
    total = mean+stationary+transient
    return mean, stationary, transient, total

def verticallyintegrated_meridional_flux(v, f, dtime, lat, lng, column, levmax, 
    levmin, ratio):
    """integrate meridional flux of tracer f over the specified vertical 
    pressure levels by calculating the tracer mass flow across each latitude.

    Parameters
    ----------
    v : numpy.ndarray 
        Northward (meridional) wind field, units of m s-1, [time, lev, lat, 
        lng]
    f : numpy.ndarray 
        Field of interest, [time, lev, lat, lng] 
    dtime : numpy.ndarray or pandas.core.indexes.datetimes.DatetimeIndex
        Time stamps corresponding to x, [time, lat, lng]
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]    
    column : 
        Model pressure levels (must be evenly spaced!), units of hPa, [lev]       
    levmax : int/float
        Desired pressure level closest to the surface (must match exactly with 
        pressure levels!)
    levmin : int/float
        Desired pressure level aloft 
    ratio : float
        Ratio of the molecular mass weight between tracer and dry air (28.97 g
        mol-1)
        
    Returns
    -------
    total : numpy.ndarray
        Total vertically-integrated meridional flux of tracer f, units of kg 
        s-1, [lat,]
    mean : numpy.ndarray 
        Mean vertically-integrated meridional flux of tracer f, units of kg 
        s-1, [lat,]    
    eddy : numpy.ndarray
        Eddy vertically-integrated meridional flux of tracer f, units of kg 
        s-1, [lat,]    
    """
    import time
    start_time = time.time()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Calculating vertically-integrated tracer flux...')    
    import numpy as np
    mean_all, eddy_all, total_all = [], [], []
    # Find pressure difference between levels, convert from hPa to Pa
    dp = -np.diff(column).mean()*100.
    # Cosine of latitude 
    coslat = np.cos(np.deg2rad(lat))
    # Loop over specified pressure levels
    for lev in np.arange(np.where(column==levmax)[0][0], 
        np.where(column==levmin)[0][0], 1):
        # Calculate meridional transport
        mean, stationary, transient, total = meridional_flux(v[:,lev], 
            f[:,lev], dtime, lat, lng)
        # The northward vertically integrated flux averaged around a latitude
        # circle; see Yang et al. (2019)
        # <F> = int_p1^p2 dP (2 pi cos(lat) r_M)/(g) • F
        # where r_M is the ratio of the molecule mass weight of tracer and 
        # dry air
        eddy = stationary+transient
        total_all.append((dp*2*np.pi*6370000.*coslat*ratio*total)/9.81)
        mean_all.append((dp*2*np.pi*6370000.*coslat*ratio*mean)/9.81)
        eddy_all.append((dp*2*np.pi*6370000.*coslat*ratio*eddy)/9.81)
    # Integrate over pressure levels
    total = np.sum(np.vstack(total_all),axis=0)
    mean = np.sum(np.vstack(mean_all),axis=0)
    eddy = np.sum(np.vstack(eddy_all),axis=0)
    print('Flux calculated in %.2f seconds!' %(time.time() - start_time))                    
    return total, mean, eddy    

def sortfield_byjetlat(field1, field2, lat_jet, lng_jet, lat, psize=30):
    """for each longitude, function fetches the jet latitude time series at 
    that longitude and, using the jet latitude, sorts the specified field by 
    the jet latitude. The top and bottom ~30th percentiles of the sorted field 
    are returned. There were some potential issues since the jet latitude is 
    repeated (i.e., many days where the jet latitude is, say, 65 or 48˚N). 
    This caused sorting to be inconsistent across days where the jet is at the 
    same latitude (i.e., two fields could be sorted differently based on 
    jet latitude because of the repeats.) To ameliorate this, the function 
    was edited to sort two fields simultaneously by jet latitude. If two 
    sorted fields aren't desired, the same field could be passed in as input 
    parameters "field1" and "field2."
    
    Parameters
    ----------
    field1 : numpy.ndarray
        First field of interest, [time, lat, lng]  
    field2 : numpy.ndarray
        Second field of interest, [time, lat, lng]          
    lat_jet : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in region, units of degrees north[time, lng]
    lng_jet : numpy.ndarray
        Longitude coordinates corresponding to jet dataset, units of degrees 
        east, [lng,]
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    psize : int
        Percentiles that will define the equator- and poleward eddy-driven jet
        days, default 30th percentile

    Returns
    -------
    eqjet_field1 : numpy.ndarray
        First field of interest on days when the jet is in an equatorward 
        position, [~time*psize*0.01, lat, lng]
    pwjet_field1 : numpy.ndarray
        First field of interest on days when the jet is in an poleward 
        position, [~time*psize*0.01, lat, lng]    
    eqjet_field2 : numpy.ndarray
        Second field of interest on days when the jet is in an equatorward 
        position, [~time*psize*0.01, lat, lng]
    pwjet_field2 : numpy.ndarray
        Second field of interest on days when the jet is in an poleward 
        position, [~time*psize*0.01, lat, lng]            
    """
    import time
    start_time = time.time()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Determining field on days with equator/poleward jet...')
    import numpy as np
    # Given the total sample size find the nearest sample size 
    # corresponding to the top and bottom 30th %-ile
    psize = np.int(np.round(len(field1)*(psize*0.01)))
    # Empty arrays to fill with value of 3D (tyx) field on days with 
    # poleward/equatorward jet 
    pwjet_field1 = np.empty(shape=(psize, field1.shape[1], field1.shape[2]))
    pwjet_field1[:] = np.nan
    eqjet_field1 = np.empty(shape=(psize, field1.shape[1], field1.shape[2]))
    eqjet_field1[:] = np.nan
    pwjet_field2 = np.empty(shape=(psize, field2.shape[1], field2.shape[2]))
    pwjet_field2[:] = np.nan
    eqjet_field2 = np.empty(shape=(psize, field2.shape[1], field2.shape[2]))
    eqjet_field2[:] = np.nan    
    for ilng in np.arange(0,len(lng_jet),1):
        # Find indices of days at a given longitude that have a poleward-
        # shifted or equatorward-shift jet. Previously this was accomplished
        # using percentiles (i.e., np.percentile(lat_jet[:,ilng], 70) and
        # np.percentile(lat_jet[:,ilng], 30)); however, this doesn't 
        # work because it creates an uneven sample size depending on the 
        # longitude considered. 
        # For a given longitude, find field of interest at each latitude 
        for ilat in np.arange(0,len(lat),1):
            # Sort field by jet latitude
            sorted_fieldjet = sorted(zip(lat_jet[:,ilng],
                field1[:,ilat,ilng], field2[:,ilat,ilng]))
            # Find values of field and time on days with pole- and equatorward 
            # jet
            eqjet_field1[:,ilat,ilng] = np.array(sorted_fieldjet[:psize])[:,1]
            pwjet_field1[:,ilat,ilng] = np.array(sorted_fieldjet[-psize:])[:,1]
            eqjet_field2[:,ilat,ilng] = np.array(sorted_fieldjet[:psize])[:,2]
            pwjet_field2[:,ilat,ilng] = np.array(sorted_fieldjet[-psize:])[:,2]            
    print('Field determined in %.2f seconds!' %(time.time() - start_time))
    return eqjet_field1, pwjet_field1, eqjet_field2, pwjet_field2

def sortfield_byjetlat_column(field1, field2, lat_jet, lng_jet, lat, column, 
    psize=30):
    """same as function 'sortfield_byjetlat' but for columned fields (i.e., 
    tracer mixing ratio from 900-850 hPa). Function segregates fields on days
    where the jet is in a pole- and equatorward position. The exact number of 
    days considered to be pole- or equatorward days is determined by optional
    argument "psize," corresponding to the percentile of interest. 
    
    Parameters
    ----------
    field1 : numpy.ndarray
        First field of interest; shape MUST be [time, lev, lat, lng]  
    field2 : numpy.ndarray
        Second field of interest; shape MUST be [time, lev, lat, lng]  
    lat_jet : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in region, units of degrees north[time, lng]
    lng_jet : numpy.ndarray
        Longitude coordinates corresponding to jet dataset, units of degrees 
        east, [lng,]
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    column : numpy.ndarray
        Model pressure levels (must be evenly spaced!), units of hPa, [lev]       
    psize : int
        Percentiles that will define the equator- and poleward eddy-driven jet
        days, default 30th percentile

    Returns
    -------
    eqjet_field1 : numpy.ndarray
        First field of interest on days when the jet is in an equatorward 
        position, [~time*psize*0.01, lev, lat, lng]
    pwjet_field1 : numpy.ndarray
        First field of interest on days when the jet is in an poleward 
        position, [~time*psize*0.01, lev, lat, lng]  
    eqjet_field2 : numpy.ndarray
        Second field of interest on days when the jet is in an equatorward 
        position, [~time*psize*0.01, lev, lat, lng]
    pwjet_field2 : numpy.ndarray
        Second field of interest on days when the jet is in an poleward 
        position, [~time*psize*0.01, lev, lat, lng]          
    """
    import numpy as np
    # Given the total sample size find the nearest sample size 
    # corresponding to the top and bottom 30th %-ile
    # Given the total sample size find the nearest sample size 
    # corresponding to the top and bottom %-ile
    pwjet_field1 = np.empty(shape=(np.int(np.round(len(field1)*(psize*0.01))), 
        field1.shape[1], field1.shape[2], field1.shape[3]))
    pwjet_field1[:] = np.nan
    eqjet_field1 = np.empty(shape=(np.int(np.round(len(field1)*(psize*0.01))), 
        field1.shape[1], field1.shape[2], field1.shape[3]))
    eqjet_field1[:] = np.nan

    pwjet_field2 = np.empty(shape=(np.int(np.round(len(field2)*(psize*0.01))), 
        field2.shape[1], field2.shape[2], field2.shape[3]))
    pwjet_field2[:] = np.nan
    eqjet_field2 = np.empty(shape=(np.int(np.round(len(field2)*(psize*0.01))), 
        field2.shape[1], field2.shape[2], field2.shape[3]))
    eqjet_field2[:] = np.nan
    # Loop through atmospheric layers
    for layer in np.arange(0, field1.shape[1], 1):
        print('For %d hPa...'%column[layer]) 
        # Segregate on days when jet is pole- versus equatorward     
        eqjet_field_layer1, pwjet_field_layer1, eqjet_field_layer2, pwjet_field_layer2 = sortfield_byjetlat(
            field1[:,layer], field2[:,layer], lat_jet, lng_jet, lat, psize=psize)
        # Append to multi-layer arrays
        eqjet_field1[:,layer] = eqjet_field_layer1
        pwjet_field1[:,layer] = pwjet_field_layer1
        eqjet_field2[:,layer] = eqjet_field_layer2
        pwjet_field2[:,layer] = pwjet_field_layer2       
    return eqjet_field1, pwjet_field1, eqjet_field2, pwjet_field2

def bin_observations_bylat(lat_inregion, obs_inregion, lat_obs_inregion): 
    """function iterates through latitudes in variable 'lat_inregion' and 
    finds observations bounded by latitude(i) and latitude(i+1) and thereafter
    calculates the mean and standard deviation of observations within 
    these bounded. If no observations are within latitude(i, i+1), the mean 
    and standard deviation are set equal to NaN. 

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates in region, units of degrees north, [lat,]
    obs_inregion : list 
        Field of interest at each observational station in region, [stations,]
    lat_obs_inregion : list
        Latitude coorindates at each observational station in region, units of 
        degrees north, [stations,]

    Returns
    -------
    obs_binned_lat : numpy.ndarray
        The average latitude of the bin; i.e., items in this list are 
        calculated by averaging the elements from the input variable: 
        lat_gmi(i) and lat_gmi(i+1) 
    obs_binned_mean : numpy.ndarray
        The value of the field averaged over all stations in latitude bins,
        [lat_binned,]
    obs_binned_err : numpy.ndarray
        The standard deviation of the field calculated over all stations in 
        latitude bins, [lat_binned,]
    """
    import numpy as np
    # Lists will be filled with the mean value and 2 sigma of observations of 
    # interest within the latitude bins as well as the mean latitude 
    obs_binned_mean = []
    obs_binned_err = []
    obs_binned_lat = []
    # Loop through model latitudes and bin by every X degrees
    for i in np.arange(0, len(lat_inregion)-1, 1):
        latlower = lat_inregion[i]
        latupper = lat_inregion[i+1]
        obs_binned_lat.append(np.mean([latlower,latupper]))
        # Find observations within band of observations
        obs_inlatband = np.where((lat_obs_inregion>=latlower) &  
                               (lat_obs_inregion<latupper))[0]
        if obs_inlatband.shape[0] > 0:
            obs_inlatband = np.array(obs_inregion)[obs_inlatband]
            obs_binned_mean.append(np.nanmean(obs_inlatband))
            obs_binned_err.append(2*np.nanstd(obs_inlatband))
        else: 
            obs_binned_mean.append(np.nan)
            obs_binned_err.append(np.nan)    
    return (np.array(obs_binned_lat), np.array(obs_binned_mean), 
        np.array(obs_binned_err))
    
def calculate_aqpi(obs, model, lat_model, lng_model, time_model, method):
    """determine the Air Quality Performance Index (AQPI) using several 
    statistics (e.g., FAC2, NMB, FMB, R) as performance metrics of O3. The AQPI
    is defined as 
        100 * AVG[FAC2 + R + (1 - ABS(MFB/2))]/3
    The AQPI provides values ranging from -33 to 0 (no skill) to 100 (perfect 
    model). Note from ECCC: Statistics are calculated using maximum daily 
    concentrations (observed and forecasted). Additional details can be found 
    at atmosphere.copernicus.eu/sites/default/files/2018-11/2_3rd_ECCC_NOAA_ECMWF_v06.pdf        
    The AQPI is either determined at each observational station and its 
    co-located grid cell or by aggregating all stations' O3 time series in 
    model grid cells and thereafter determining the AQPI.

    Parameters
    ----------
    obs : pandas.core.frame.DataFrame
        Daily averaged O3 concentrations at individual observational sites for 
        the specified measuring period. DataFrame contains station ID and 
        station latitude and longitude
    model : numpy.ndarray
        Modeled O3 concentrations, units of ppbv, [time, lat, lng]
    time_model : numpy.ndarray
        Contains datetime.date objects corresponding to the timestamps of the
        model output, [time,]
    lat_model : numpy.ndarray
        Model latitude coordinates, units of degrees north, [lat,]
    lng_model : numpy.ndarray
        Model longitude coordinates, units of degrees east, [lng,]
    method : str
        Either 'pointwise' (at each observational station and its co-located 
        grid cell) or 'gridwise' (at each model grid cell with all O3 from 
        observational stations within the cell averaged to a single time 
        series).
        
    Returns
    -------
    aqpi_all : list
        The AQPI calculated between each observational station and the 
        co-located model grid cell [stations,]
    lat_all : list 
        Latitude coorindate at each observational stations, units of degrees 
        north, [stations,]
    lng_all : list 
        Longitude coorindate at observational stations, units of degrees east, 
        [stations,]
    aqpi : numpy.ndarray
        The AQPI calculated with O3 averaged over observational stations within
        each model grid cell, [lat, lng]
    fac2 : numpy.ndarray
        Factor-of-2 fraction (FAC2) calculated between O3 averaged over 
        observational stations within each model grid cell, [lat, lng]
    r : numpy.ndarray
        Pearson correlation coefficient (r) calculated between O3 averaged over 
        observational stations within each model grid cell, [lat, lng]    
    mfb : numpy.ndarray
        Mean fractional bias (MFB) calculated between O3 averaged over 
        observational stations within each model grid cell, [lat, lng]    
    """
    import time
    start_time = time.time()    
    import numpy as np
    import numpy.ma as ma
    import sys
    sys.path.append('/Users/ghkerr/phd/utils/')
    from geo_idx import geo_idx
    time_model = [x.strftime('%Y-%m-%d') for x in time_model]
    # for station-by-station analysis
    if method=='pointwise':
        aqpi_all, lat_all, lng_all = [], [], []
        # Loop through all stations represented in observaitons
        for station_id in np.unique(obs['Station ID'].values):
            obs_station = obs.loc[obs['Station ID'] == station_id]
            # Only calculate r(T, O3) and dO3/dT for stations with data 
            if (obs_station.shape[0] != 0) \
                & (True in np.isfinite(obs_station['Latitude'].values)) \
                & (True in np.isfinite(obs_station['Longitude'].values)) \
                & (True in np.isfinite(obs_station['O3 (ppbv)'].values)):
                station_lat = np.nanmean(obs_station['Latitude'].values)
                station_lng = np.nanmean(obs_station['Longitude'].values)
                # Find closest MERRA-2 grid cell
                lat_closest = geo_idx(station_lat, lat_model)
                lng_closest = geo_idx(station_lng, lng_model)   
                station_times = np.in1d(time_model, 
                    obs_station.index).nonzero()[0]
                # Select O3 on days with observations
                o3_inperiod = model[station_times, lat_closest, lng_closest]
                idx = np.isfinite(o3_inperiod) \
                        & np.isfinite(obs_station['O3 (ppbv)'].values)
                # Predicted and observed O3
                p = o3_inperiod[idx]
                o = obs_station['O3 (ppbv)'].values[idx]
                # Mask values of 0 in observations as these cause errors in the 
                # calculation of the FAC2
                o = ma.masked_where(o==0., o)
                # Factor-of-2 fraction (FAC2) (measure of error or scatter)
                # Provides fraction (0-1) of modelled and observed pairs 
                # meeting this criterion (dimensionless statistic, not 
                # sensitive to outliers)
                fac2 = (p/o)
                fac2 = np.where((fac2>0.5) & (fac2<2.))[0]
                fac2 = len(fac2)/len(p)
                # Correlation coefficient (r)
                # Measure of linearity of relationship (dimensionless, values 
                # between -1 and 1)
                r = np.corrcoef(o,p)[0, 1]
                # Mean fractional bias (MFB) (measure of bias or offset)
                # Where MFB = 2x[(model-obs)/(model+obs)] and 1-ABS(MFB/2) 
                # provides values in range 0-1 (dimensionless, symmetric and 
                # bounded statistic)
                mfb = 2*((p-o)/(p+o))
                # Calculate Air Quality Performance Index (AQPI)
                aqpi = 100*np.mean(fac2+r+(1-np.abs(mfb/2.)))/3.
                aqpi_all.append(aqpi)
                lat_all.append(station_lat)
                lng_all.append(station_lng)
        print('AQPI calculated in %.2f seconds' %((time.time()-start_time))) 
        return aqpi_all, lat_all, lng_all        
    # For grid cell-aggregate analysis
    if method=='gridwise':    
        aqpi = np.empty(shape=(lat_model.shape[0], lng_model.shape[0]))
        aqpi[:] = np.nan
        fac2 = np.empty(shape=(lat_model.shape[0], lng_model.shape[0]))
        fac2[:] = np.nan
        r = np.empty(shape=(lat_model.shape[0], lng_model.shape[0]))
        r[:] = np.nan
        mfb = np.empty(shape=(lat_model.shape[0], lng_model.shape[0]))
        mfb[:] = np.nan        
        # Resolution of model 
        latres = np.diff(lat_model[:-1]).mean()
        lngres = np.diff(lng_model[:-1]).mean()
        # Loop through model grid cells
        for i, ilat in enumerate(lat_model):
            for j, ilng in enumerate(lng_model):
                # Define bounding box
                bbleft = ilng-(lngres/2.)
                bbup = ilat+(latres/2.)
                bbright = ilng+(lngres/2.)
                bbdown = ilat-(latres/2.)
                # Find observations in bounding box
                obs_inbb = obs.loc[(obs['Latitude'] >= bbdown) & 
                        (obs['Latitude'] < bbup) & 
                        (obs['Longitude'] >= bbleft) & 
                        (obs['Longitude'] < bbright)]
                if (obs_inbb.shape[0] != 0) \
                    & (True in np.isfinite(obs_inbb['Latitude'].values)) \
                    & (True in np.isfinite(obs_inbb['Longitude'].values)) \
                    & (True in np.isfinite(obs_inbb['O3 (ppbv)'].values)):
                    # Average over all station in grid cell 
                    obs_inbb = obs_inbb.groupby(obs_inbb.index).mean()
                    station_times = np.in1d(time_model, 
                        obs_inbb.index).nonzero()[0]
                    # Select O3 on days with observations
                    o3_inperiod = model[station_times, i, j]
                    idx = np.isfinite(o3_inperiod) \
                            & np.isfinite(obs_inbb['O3 (ppbv)'].values)
                    # Predicted and observed O3
                    p = o3_inperiod[idx]
                    o = obs_inbb['O3 (ppbv)'].values[idx]            
                    # Mask values of 0 in observations as these cause errors in the 
                    # calculation of the FAC2
                    o = ma.masked_where(o==0., o)
                    # Factor-of-2 fraction (FAC2) (measure of error or scatter)
                    # Provides fraction (0-1) of modelled and observed pairs meeting 
                    # this criterion (dimensionless statistic, not sensitive to 
                    # outliers)
                    fac2_ij = (p/o)
                    fac2_ij = np.where((fac2_ij>0.5) & (fac2_ij<2.))[0]
                    fac2_ij = len(fac2_ij)/len(p)
                    fac2[i,j]=fac2_ij                    
                    # Correlation coefficient (r)
                    # Measure of linearity of relationship (dimensionless, values 
                    # between -1 and 1)
                    r_ij = np.corrcoef(o,p)[0, 1]
                    r[i,j]=r_ij
                    # Mean fractional bias (MFB) (measure of bias or offset)
                    # Where MFB = 2x[(model-obs)/(model+obs)] and 1-ABS(MFB/2) provides 
                    # values in range 0-1 (dimensionless, symmetric and bounded 
                    # statistic)
                    mfb_ij = 2*((p-o)/(p+o))
                    mfb[i,j]=np.nanmean(mfb_ij)
                    # Calculate Air Quality Performance Index (AQPI)
                    aqpi_ij = 100*np.mean(fac2_ij+r_ij+(1-
                        np.abs(mfb_ij/2.)))/3.
                    aqpi[i,j]=aqpi_ij
        print('AQPI calculated in %.2f seconds' %((time.time()-start_time))) 
        return aqpi, fac2, r, mfb

def calculate_initial_compass_bearing(pointA, pointB):
    """Calculates the bearing between two points. The formula used is the 
    following:
        θ = atan2(sin(Δlong).cos(lat2),
                  cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(Δlong))
    Taken from https://gist.github.com/jeromer/2005586
    
    Parameters
    ----------
    pointA : tuple
        The tuple representing the latitude/longitude for the first point. 
        Latitude and longitude must be in decimal degrees
    pointB : tuple
        The tuple representing the latitude/longitude for the second point. 
        Latitude and longitude must be in decimal degrees

    Returns
    -------
    compass_bearing : float
      The bearing in degrees
    """
    import math
    if (type(pointA) != tuple) or (type(pointB) != tuple):
        raise TypeError("Only tuples are supported as arguments")
    lat1 = math.radians(pointA[0])
    lat2 = math.radians(pointB[0])
    diffLong = math.radians(pointB[1] - pointA[1])
    x = math.sin(diffLong) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1)
            * math.cos(lat2) * math.cos(diffLong))
    initial_bearing = math.atan2(x, y)
    # Now we have the initial bearing but math.atan2 return values
    # from -180° to + 180° which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = math.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360
    return compass_bearing

def o3anom_cyclone(cyclones, time_model, lat_model, lng_model, o3_model):
    """Screen MCMS cyclone dataset to consider only cyclones over land (land 
    fraction > 0.5) with > 1 daily mean time step (to determine direction of 
    cyclone propogation). With this subset of cyclones, the O3 and O3 anomaly
    within +/- 5 grid cells of the cyclones' centers ares found. 
    
    Parameters
    ----------
    cyclones : pandas.core.frame.DataFrame
        DataFrame containing the date, hour, fraction of land cover, latitude, 
        longitude, SLP, and storm ID
    time_model : numpy.ndarray
        Contains datetime.date objects corresponding to the timestamps of the
        model output, [time,]        
    lat_model : numpy.ndarray
        Latitude coordinates for the Northern Hemisphere, units of degrees 
        north, [lat,]                
    lng_model : numpy.ndarray
        Longitude coordinates for the Northern Hemisphere, units of degrees 
        east, [lng,]          
    o3_model : numpy.ndarray
        Modeled O3 for all days, units of ppbv, [time, lat, lng]

    Returns
    -------
    o3_anom_all : list
        The daily O3 anomaly in the vicinity of cyclone for all cyclones over 
        land (land fraction > 0.5) and with > 1 daily mean timestamp, 
        [cyclones,]
    o3_all : list
        The daily O3 in the vicinity of cyclone for all cyclones over 
        land (land fraction > 0.5) and with > 1 daily mean timestamp, 
        [cyclones,]    
    o3_anom_rotated : list 
        Same as o3_anom_all but individual cyclones have been rotated such that 
        the upward direction is the direction of cyclone motion, [cyclones,]    
    o3_rotated : list 
        Same as o3_all but individual cyclones have been rotated such that 
        the upward direction is the direction of cyclone motion, [cyclones,]            
    """
    import time
    start_time = time.time()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Determining O3 and O3 anomaly in vicinity of cyclones...')   
    import numpy as np    
    from scipy.ndimage import rotate
    import pandas as pd
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    # Drop all cyclone detections that are not in the latitudinal range of the 
    # focus region 
    cyclones = cyclones[cyclones['Latitude'] > lat_model.min()]
    # Combine date and hour columns to a single column and drop hour column 
    cyclones['Date'] = pd.to_datetime(cyclones['Date']+' '+cyclones['Hour'])
    del cyclones['Hour']
    # List creation 
    o3_anom, o3_anom_rotated = [], []
    o3_all, o3_rotated = [], []
    brngs = []
    totnumb = 0
    tothours = []
    # Loop through unique cyclone IDofTracks 
    for id in np.unique(cyclones['Storm ID'].values):
        cycloneid = cyclones.loc[cyclones['Storm ID']==id]
        # Sort by time (although this shouldn't be needed)
        cycloneid = cycloneid.sort_values(by='Date')
        # Only consider cyclones that are detected for at least two 6-hourly 
        # periods; if this is not done, we can't find the bearing (degrees) 
        # which corresponds to the direction that the cyclone is tracking
        if ((len(np.unique(cycloneid.Date.dt.date.values))>1)  & 
            (cycloneid['Land Cover'].values.mean()>0.5)): 
            tothours.append(len(cycloneid))
            # Calculate daily mean value
            cycloneid = cycloneid.resample('D', on='Date').mean()
            totnumb = totnumb+1
            # Loop through daily mean values and find cyclone bearing and O3 in
            # vicinity of cyclone
            for day in np.arange(0, len(cycloneid), 1):
                # Determine bearing of cyclone (note this only works until 
                # the day before the last cyclone timestamp; on the final day
                # that a particular cyclone is detected, use the bearing 
                # from the previous day)
                lat1 = cycloneid.iloc[day].Latitude
                lng1 = cycloneid.iloc[day].Longitude
                try: 
                    lat2 = cycloneid.iloc[day+1].Latitude
                    lng2 = cycloneid.iloc[day+1].Longitude                
                except IndexError:
                    lat1 = cycloneid.iloc[day-1].Latitude
                    lng1 = cycloneid.iloc[day-1].Longitude                
                    lat2 = cycloneid.iloc[day].Latitude
                    lng2 = cycloneid.iloc[day].Longitude      
                # n.b., if cyclone is tracking north, bearing = 0
                brng = calculate_initial_compass_bearing((lat1,lng1), 
                    (lat2,lng2))
                brngs.append(brng)
                # Find model indices at center of cyclone            
                latidx = cycloneid.iloc[day].Latitude
                lngidx = cycloneid.iloc[day].Longitude
                latidx = geo_idx(latidx, lat_model)
                lngidx = geo_idx(lngidx, lng_model)
                timeidx = np.where(time_model==
                    cycloneid.iloc[day].name.date())[0][0]
                # Determine O3 anomaly (difference from all values in measuring 
                # period) in the vicinity of cyclone (+/- 10 grid cells)
                o3_anom_day = (o3_model[timeidx,latidx-6:latidx+7,
                    lngidx-6:lngidx+7]-np.nanmean(o3_model[:,latidx-6:latidx+7,
                    lngidx-6:lngidx+7],axis=0))
                o3_day = (o3_model[timeidx,latidx-6:latidx+7,
                    lngidx-6:lngidx+7])    
                # Exclude cyclones on the periphery of domain 
                if o3_anom_day.shape == (13, 13):
                    # Rotate O3 anomaly to orient the right/east direction as 
                    # the mean bearing of the cyclone using spline 
                    # interpolation of the third order
                    o3_anom.append(o3_anom_day)
                    o3_anom_rotated.append(rotate(o3_anom_day, -brng, 
                        mode='constant', cval=np.nan, reshape=False))
                    o3_all.append(o3_day)
                    o3_rotated.append(rotate(o3_day, -brng, mode='constant', 
                        cval=np.nan, reshape=False))
    print('Mean bearing = %.3f (0 deg is north)'%np.nanmean(brngs))
    print('Total land-based cyclones with >1 timestamp = %.1f'%totnumb)
    print('Mean number of 6-hourly timesteps for cyclones = %.2f timesteps'
          %np.nanmean(tothours))
    print('O3 and O3 anomaly determined in %.2f seconds!'
        %(time.time()-start_time))
    return o3_anom, o3_all, o3_anom_rotated, o3_rotated