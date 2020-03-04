#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains statistical modules supporting the O3-temperature-jet 
relationship, specifically bootstrapping and significance 

Revision History
    15112019 -- initial version created
    03032020 -- edited to include function 'mmb_grid' to conduct moving block
                bootstrap on gridded fields
"""

def movingblocks_bootstrap_r(x, y, nb=10000, L=10, alpha=0.05):
    """Function performs a moving-blocks bootstrap by resampling given 1D
    timeseries (x and y) in contiguous blocks of length L and thereafter 
    calculating the correlation coefficient. The function constructs the 95% 
    (or any specified value) confidence intervals using the percentile method. 
    Given the null hypothesis that r_xy = 0, the null hypothesis can be 
    rejected if the confidence interval contains zero. The moving-block 
    bootstrap is used here as the data are not independent. Since the data are 
    mutually correlated (time correlation/persistence), the random sampling 
    used in ordinary bootstrapping would destroying the ordering that produces 
    that autocorrelation. The confidence interval is constructed using a two-
    sided test: in our case, values of r_xy that are both very large (strongly
    correlated) or very small (strongly anticorrelated) would be unfavorable
    to the null hypothesis of no correlation. Thus the rejection region for 
    this two-sided test consists of both extreme left and extreme right tails 
    of the null distribution. 

    Parameters
    ----------
    x : numpy.ndarray
        Serially-correlated x time series, (no. obs,)
    y : numpy.ndarray
        Serially-correlated y time series, (no. obs,)        
    nb : int
        Number of bootstrap samples; according to Wilks (2011), the 
        bootstrapping process should be repeated a large number of times 
        (~10,000) to yield nb samples
    L : int
        Blocklength; vectors of length L are resampled from the data in order
        to build up a synthetic sample. L can be determined using an implicit 
        equation 
    alpha : float
        Significance level (default alpha = 0.05), where the confidence level 
        is equivalent to (1 – alpha) level. n.b., alpha is the probability that
        a confidence interval will not include the population parameter, and 
        (1 - alpha) is the probability the population parameter will be 
        in the interval. The 100(1 - alpya)% confidence interval will include 
        the true value of the population parameter with probability 1 - α, 
        i.e., if alpha = .05, the probability is about .95 that the 95% 
        confidence interval will include the true population parameter.
        
    Returns
    -------
    ub : numpy.float64 
        Upper bound for (1 - alpha)% confidence interval for r_xy
    lb : numpy.float64
        Lower bound for (1 - alpha)% confidence interval for r_xy    
    
    References
    ----------
    Mudelsee, M. (2003). Estimating Pearson's correlation coefficient with 
    bootstrap confidence interval from serially dependent time series. Math. 
    Geosci. 35(6):651-665.
    
    Thiébaux, H. J., & Zwiers, F. W. (1984). The interpretation and estimation
    of effect sample size. J. Clim. App. Met. 23:800-811.
    
    Wilks, D. S. (2011). Statistical methods in the atmospheric sciences (3rd 
    ed.). Oxford ; Waltham, MA: Academic Press.
    
    Notes
    -----
    To visualize the null distribution of r for a particular grid cell, compute
    the following with time series x and y (e.g., 2-meter temperature and 
    surface-level O3) of equal length
    import matplotlib.pyplot as plt
    # Time series of T, O3
    plt.plot(y, '-k', label='O$_{\mathregular{3}}$'); plt.legend() 
    plt.ylabel('O$_{3}$ [ppbv]');plt.twinx(); plt.plot(x, '-r', label='T')
    plt.ylabel('T [K]');plt.legend(loc=1)
    plt.title('r = %.3f'%(np.corrcoef(x, y)[0, 1]))
    plt.savefig('/Users/ghkerr/Desktop/o3ttimeseries.png', dpi=300)
    # Histogram of r(T, O3) from resampled values
    plt.hist(mbbr, bins=20), 
    plt.vlines(np.corrcoef(x, y)[0, 1], ymin=plt.ylim()[0], ymax=160)
    plt.vlines(upper, linestyles='--', ymin=plt.ylim()[0], ymax=plt.ylim()[1])
    plt.vlines(lower, linestyles='--', ymin=plt.ylim()[0], ymax=plt.ylim()[1])
    plt.xlabel('r(T, O$_{\mathregular{3}}$)')
    plt.ylabel('Frequency')
    plt.ylim([0, 160])
    plt.show()
    """
    import numpy as np
    # Blocks for moving block bootstrapping have length L; from Wilks 3rd. ed 
    # 2011 (pg 178), the blocklength can be chosen according to the implicit 
    # equation L = (n - L + 1)^((2/3)((1-n')/n)) where n' is given by 
    # n' = n ((1 - rho1)/(1 + rho1)) 
    # where rho1 is the lag-1 autocorrelation coefficient 
    L = 10
    # Using this blocklength, we find how many blocks we can pick given our 
    # sample size. Ideally the length of all our blocks would be identically 
    # equal to the sample size, n, but this isn't always the case, so we
    # round to the nearest integer. The following line draws m = n/L blocks
    m = np.round(len(x)/L)
    # Unlike ordinary bootstrapping where we resample from a collection of n
    # individual, independent values, the objects to be resampled with 
    # replacement are all the n - L + 1 contiguous subseries of length L
    # (Wilks 2011) 
    totm = len(x)-L+1
    # Will be filled with estimations of the Pearson correlastion coefficient 
    # for each moving-block bootstrapped sample
    mbbr = []
    # Conduct nb samples; using np.seed function provides an input for the 
    # pseudo-random number generator so that the bootstrapping methods and 
    # subsequent distribution of correlation coefficients should be preserved
    # if nb is the same
    for seedn in np.arange(0, nb, 1):
        np.random.seed(seedn)
        # Given the n - L + 1 contiguous subseries that exist, the following 
        # line finds the first index of m blocks
        startm = np.random.randint(0, totm, size=int(m))
        mbbx, mbby = [], []
        # Draw m blocks of length L (with replacement) and form resultant time
        # series
        for startmi in startm:
            mbbx.append(x[startmi:startmi+L])
            mbby.append(y[startmi:startmi+L])
        # Align blocks back to back in the order they were picked
        mbbx = np.hstack(mbbx)
        mbby = np.hstack(mbby)
        # Calculate correlation coefficient 
        mbbr.append(np.corrcoef(mbbx, mbby)[0,1])
    # Construct confidence regions for r_xy using the percentile method 
    # (Wilks 2011). To do so, the bootstrap replications are taken to construct
    # an equitailed 100% *(1 - alpha)/2 confidence interval find the 
    # values of the parameter, in this case r_xy, defining the largest and 
    # smallest nb * alpha/2 of the nb bootstrap estimates
    bounds = np.int(nb*(alpha/2.))
    ub = np.sort(mbbr)[-bounds]
    lb = np.sort(mbbr)[bounds]
    # One could also explore using bias-corrected and accelerated (BCa) 
    # intervals, which are more accurate that bootstrap confidence intervals 
    # (Wilks 2011)
    return ub, lb

def mmb_grid(x, y, lat, lng):
    """Calculate significance using moving block bootstrapping for gridded 
    fields. 

    Parameters
    ----------
    x : numpy.ndarray
        Independent variable, [time, lat, lng]
    y : numpy.ndarray
        Dependent variable, [time, lat, lng]
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]
        
    Returns
    -------
    sig : numpy.ndarray
        Significance at the significance level set in function 
        'movingblocks_bootstrap_r' denoted with grid cells equal to 1. are 
        insignificant at specified significance level, [lat, lng]
    """
    import numpy as np
    # Lower and upper bounds of confidence interval
    lb = np.empty(shape=y.shape[1:])
    lb[:] = np.nan
    ub = np.empty(shape=y.shape[1:])
    ub[:] = np.nan
    for ilat in np.arange(0, lat.shape[0], 1):
        print('For latitude circle %.1f deg N...'%lat[ilat])
        for jlng in np.arange(0, lng.shape[0], 1):
            ubxy, lbxy = movingblocks_bootstrap_r(x[:, ilat, jlng], 
                y[:, ilat, jlng], nb=10000)
            lb[ilat, jlng] = lbxy
            ub[ilat, jlng] = ubxy
    sig = np.empty(shape=y.shape[1:])
    sig[:] = np.nan
    for ilat in np.arange(0, lat.shape[0], 1):
        for jlng in np.arange(0, lng.shape[0], 1):
            if lb[ilat, jlng] <= 0.0 <= ub[ilat, jlng]:
                sig[ilat, jlng] = 1.  
    return sig
    
        
        
      
         
# datapath = '/Users/ghkerr/phd/globalo3/data/parsed/'      
datapath = '/mnt/scratch3/gaige/kerr_surface_2020/data/parsed/'
import netCDF4 as nc
import xarray as xr
o3_gmi = nc.Dataset(datapath+'gmi_O3control_JJA2008-2010.nc')['O3_control'][:].data
o3_gmi_transport = nc.Dataset(datapath+'gmi_O3transportonly_JJA2008-2010.nc')['O3_transportonly'][:].data
lat_gmi = nc.Dataset(datapath+'gmi_O3control_JJA2008-2010.nc')['lat'][:].data
lng_gmi = nc.Dataset(datapath+'gmi_O3control_JJA2008-2010.nc')['lng'][:].data
t2m_merra = nc.Dataset(datapath+'merra2_T2M_JJA2008-2010.nc')['T2M'][:].data

x = t2m_merra
y = o3_gmi
lat = lat_gmi
lng = lng_gmi
sig = mmb_grid(x, y, lat, lng)
# 
ds = xr.Dataset({'sig_T2M_O3_control': (('lat', 'lng'), sig)},
    coords={'lat': lat, 'lng': lng})
ds.attrs['title'] ='T2M-O3 significance'
ds.attrs['history'] ='Significance of correlation between 2-meter '+\
    'temperature and surface-level O3 determined with moving block '+\
    'bootstrapping  with alpha = 0.05 and 10 000 realizations. Grid cells '+\
    'with insigificant correlation given by 1., significant given by NaN.'
ds.attrs['author'] ='Gaige Hunter Kerr, gaige.kerr@jhu.edu'
ds.to_netcdf(datapath+'sig_merra2_t2m_gmi_O3control_JJA2008-2010.nc')


     








        

# lb = np.empty(shape=o3_gmi.shape[1:])
# lb[:] = np.nan
# ub = np.empty(shape=o3_gmi.shape[1:])
# ub[:] = np.nan
# for x in np.arange(0, lat.shape[0], 1):
#     for y in np.arange(0, lng.shape[0], 1):
#        ubxy, lbxy = movingblocks_bootstrap_r(t2m_merra[:, x, y], 
#                                              o3_gmi[:, x, y], nb=1000)
#        lb[x, y] = lbxy
#        ub[x, y] = ubxy
#sig = np.empty(shape=o3_gmi.shape[1:])
#sig[:] = np.nan
#for x in np.arange(0, lat_gmi.shape[0], 1):
#    for y in np.arange(0, lng_gmi.shape[0], 1):
#        if lb[x, y] <= 0.0 <= ub[x,y]:
#            sig[x,y] = 1.     