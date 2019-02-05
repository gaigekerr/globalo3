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
            do3dt[i,j] = np.polyfit(t2m[:,i,j], o3[:,i,j], deg=1)[0]
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