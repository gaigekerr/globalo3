#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot quantities related to the distribution of O3, r(T, O3), and dO3/dT 
and their covariance with the eddy-driven jet in the Northern Hemisphere.

Revision History
    23012019 -- initial version created
    28012019 -- function 'map_global' edited to include stippling for grid 
                cells with statistically significant values
    03022019 -- function 'timeseries_seasonalo3' added
    19022019 -- map functions changed to just focus on Northern Hemisphere
    23042019 North American r(T, O3)-- function 'timeseries_o3atabovebelow_jet' added
    26052019 -- function 'map_jet_centerdist' added
    17062019 -- function 'edjetlocation_fieldatedjet' added
    18062019 -- contour/label capability added to 'map_nh' 
    09072019 -- edit function 'edjetlocation_fieldatedjet' to handle and 
                number of grid cells/degrees from jet center
    17072019 -- function 'edjetlocation_fieldatedjet' name changed to 
                'fieldatjet'
"""
# Change font
import sys
if 'mpl' not in sys.modules:
    import matplotlib.font_manager
    prop = matplotlib.font_manager.FontProperties(
            fname='/Users/ghkerr/Library/Fonts/cmunbmr.ttf')
    matplotlib.rcParams['font.family'] = prop.get_name()
    prop = matplotlib.font_manager.FontProperties(
        fname='/Users/ghkerr/Library/Fonts/cmunbbx.ttf')
    matplotlib.rcParams['mathtext.bf'] = prop.get_name()
    # for unicode minus/negative sign implementation
    matplotlib.rcParams['axes.unicode_minus'] = False
    # change width and thickness of ticks/spines
    matplotlib.rcParams['axes.linewidth'] = 1.5
    matplotlib.rcParams['xtick.major.width'] = 1.5
    matplotlib.rcParams['xtick.minor.width'] = 1.5
    matplotlib.rcParams['ytick.major.width'] = 1.5
    matplotlib.rcParams['ytick.minor.width'] = 1.5
        

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap; adapted 
    from stackoverflow.com/questions/18704353/correcting-matplotlib-colorbar-ticks

    Parameters
    ----------
    cmap: matplotlib.colors.LinearSegmentedColormap
        colormap instance, e.g., plt.get_cmap('jet')
    N: int
        Number of colors

    Returns
    -------
    cmap_d : matplotlib.colors.LinearSegmentedColormap
        Discretized colormap

    Example
    -------    
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])
                       for i in np.arange(N+1) ]
    # Return colormap object
    cmap_d = mcolors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)
    return cmap_d

def map_do3dt_endminusbeginning(lng_gmi_n, lat_gmi_n, lng_gmi_s, lat_gmi_s, 
    do3dt_byyr_n, do3dt_byyr_s, years):
    """plot simple difference between O3-climate penalty during the last and 
    first years (and months considered within each year) in measuring period 
    
    Parameters
    ----------
    lng_gmi_n : numpy.ndarray 
        GMI longitude coordinates for the Northern Hemisphere, units of degrees 
        east, [lng,]      
    lat_gmi_n : numpy.ndarray 
        GMI latitude coordinates for the Northern Hemisphere, units of degrees 
        north, [lat,]        
    lng_gmi_s : numpy.ndarray 
        GMI longitude coordinates for the Southern Hemisphere, units of degrees 
        east, [lng,]      
    lat_gmi_s : numpy.ndarray
        GMI latitude coordinates for the Southern Hemisphere, units of degrees 
        north, [lat,]    
    do3dt_byyr_n : numpy.ndarray     
        The O3-climate penalty for each year for the Northern Hemisphere, units 
        of ppbv K-1, [years, lat, lng] 
    do3dt_byyr_s : numpy.ndarray     
        The O3-climate penalty for each year for the Southern Hemisphere, units 
        of ppbv K-1, [years, lat, lng]        
    years : list
        Year or range of years in measuring period
        
    Returns
    -------
    None    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    plt.figure()
    axt=plt.subplot2grid((2,2), (0,0), colspan=2,
                         projection=ccrs.Robinson(central_longitude=0.))
    axb=plt.subplot2grid((2,2), (1,0), colspan=2,
                         projection=ccrs.Robinson(central_longitude=0.))
    clevs = np.linspace(-1.5, 1.5, 11)
    cmap = plt.get_cmap('bwr')
    # Northern Hemisphere
    axt.contourf(lng_gmi_n, lat_gmi_n, (do3dt_byyr_n[-1]-do3dt_byyr_n[0]), 
                 clevs, cmap=cmap, transform=ccrs.PlateCarree(), extend='both')
    axt.add_feature(cfeature.OCEAN, zorder=10, lw = 0.0, color='lightgrey')
    axt.coastlines(resolution='50m', lw=0.25)
    axt.set_title('June-August', fontsize=16)
    # Southern Hemisphere
    mb = axb.contourf(lng_gmi_s, lat_gmi_s, (do3dt_byyr_s[-1]-do3dt_byyr_s[0]), 
                      clevs, cmap=cmap, transform=ccrs.PlateCarree(), 
                      extend='both')
    axb.add_feature(cfeature.OCEAN, zorder=10, lw = 0.0, color='lightgrey')
    axb.coastlines(resolution='50m', lw=0.25)
    axb.set_title('July-September', y=-0.3, fontsize=16)
    # Add colorbar
    colorbar_axes = plt.gcf().add_axes([0.83,0.15,0.03,0.7])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical')
    colorbar.set_label('$\mathregular{\partial}$O$_{\mathregular{3}}$'+
                       ' $\mathregular{\partial}$T$^{\mathregular{-1}}$'+
                       '(%s - %s)'%(years[-1], years[0])+
                       ' [ppbv K$^{\mathregular{-1}}$]', fontsize=16)
    plt.gcf().subplots_adjust(right=0.8, hspace=-0.2)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'map_do3dt_endminusbeginning_%d-%d.eps'%(years[0],years[-1]))
    return 

def map_trend(lat_n, lat_s, lng_n, lng_s, trend_n, trend_s, cbar_label, 
              clevs, fstr): 
    """plot the linear trend in yearly/seasonal data

    Parameters
    ----------
    lat_n : numpy.ndarray
        Latitude coordinates for the Northern Hemisphere, units of degrees 
        north, [lat,]                
    lat_s : numpy.ndarray
        Latitude coordinates for the Southern Hemisphere, units of degrees 
        north, [lat,]    
    lng_n : numpy.ndarray
        Longitude coordinates for the Northern Hemisphere, units of degrees 
        east, [lng,]          
    lng_s : numpy.ndarmray
        Longitude coordinates for the Southern Hemisphere, units of degrees 
        east, [lng,]      
    trend_n : numpy.ndarray
        The trend in the yearly/seasonally-averaged data for the Northern 
        Hemisphere, units are data units yr-1, [lat, lng]
    trend_s : numpy.ndarray
        The trend in the yearly/seasonally-averaged data for the Southern 
        Hemisphere, units are data units yr-1, [lat, lng]
    cbar_label : str
        Label for the colorbar that will proceed "Trend"
    clevs : numpy.ndarray
        Contour levels values
    fstr : str
        Output filename suffix

    Returns
    -------
    None             
    """
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    plt.figure()
    axt=plt.subplot2grid((2,2), (0,0), colspan=2,
                         projection=ccrs.Robinson(central_longitude=0.))
    axb=plt.subplot2grid((2,2), (1,0), colspan=2,
                         projection=ccrs.Robinson(central_longitude=0.))
    cmap = plt.get_cmap('bwr')
    # Northern Hemisphere
    axt.contourf(lng_n, lat_n, trend_n, clevs, cmap=cmap,
                 transform=ccrs.PlateCarree(), extend='both')
    axt.add_feature(cfeature.OCEAN, zorder=10, lw = 0.0, color='lightgrey')
    axt.coastlines(resolution='50m', lw=0.25)
    axt.set_title('June-August', fontsize=16)
    # Southern Hemisphere
    mb = axb.contourf(lng_s, lat_s, trend_s, clevs, cmap=cmap,
                      transform=ccrs.PlateCarree(), extend='both')
    axb.add_feature(cfeature.OCEAN, zorder=10, lw = 0.0, color='lightgrey')
    axb.coastlines(resolution='50m', lw=0.25)
    axb.set_title('July-September', y=-0.3, fontsize=16)
    # Add colorbar
    colorbar_axes = plt.gcf().add_axes([0.78,0.15,0.03,0.7])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical')
    colorbar.set_label('Trend %s' %cbar_label, 
                       fontsize=16)
    plt.gcf().subplots_adjust(right=0.75, hspace=-0.2)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'map_trend_%s.eps'%fstr)
    return

def map_hemisphere(lat_n, lng_n, field_n, title, cbar_label, clevs, cmap, 
    hemisphere, fstr, contour=None, contour_levs=None, p_n=None, e_n=None, 
    eerr_n=None, alpha=0.05, extent=None, quiver=None, oceanon='yes', 
    extend='both'):
    """plot desired quantity over the Northern or Southern hemisphere (however, 
    if variable 'extent' is specified, map can be tailored to any region).
    
    Parameters
    ----------
    lat_n : numpy.ndarray
        Latitude coordinates for the Northern Hemisphere, units of degrees 
        north, [lat,]                
    lng_n : numpy.ndarray
        Longitude coordinates for the Northern Hemisphere, units of degrees 
        east, [lng,]          
    field_n : numpy.ndarray
        Desired field/quantity in the Northern Hemisphere, [lat, lng]
    title : str
        Title for plot        
    cbar_label : str
        Label for the colorbar (field and units)
    clevs : numpy.ndarray
        Filled contour levels values
    cmap : str
        Colormap name
    hemisphere : str
        Focus region for output filename (i.e., 'nh' for northern hemisphere, 
        'sh' for southern hemisphere, 'northamerica', etc.)
    fstr : str
        Output filename suffix
    contour : numpy.ndarray
        If type(f_c) is numpy.ndarray, black contours and labels are plotted 
        atop filled contours
    contour_levs : numpy.ndarray
        Contour level values
    p_n : numpy.ndarray or NoneType
        If type(p_n) is numpy.ndarray, gridded values less than the 
        significance level (alpha) are plotted as scatterpoints
    e_n : numpy.ndarray or NoneType
        If type(e_n) is numpy.ndarray, values are plotted as scatterpoints
    eerr_n : numpy.ndarray or NoneType
        If type(eerr_n) is numpy.ndarray, values are plotted as errorbars 
        on the scatterpoints of e_n
    alpha : int
        The significance level, the probability of rejecting the null 
        hypothesis p-values are less than alpha
    extent : list or NoneType
        The extent (x0, x1, y0, y1) of the map in the given coordinate system    
    quiver : tuple or NoneType
        First (second) object, an numpy.ndarray, in tuple corresponds to 
        U-wind (V-wind)
    oceanon : str
        If 'yes', map adds ocean polygons feature        
    extend : str
        Extend settings for matplotlib colormap/colorbar (i.e., 'neither', 
        'both', 'min', 'max')

    Returns
    -------
    None              
    """
    import copy
    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure(figsize=(8,3.5))
    ax=plt.subplot2grid((1,2), (0,0), colspan=2,
                        projection=ccrs.Miller(central_longitude=0.))
    ax.set_title(title, fontsize=14, x=0.02, ha='left')    
    if oceanon == 'yes':
        ax.add_feature(cfeature.OCEAN, zorder=10, lw = 0.0, color='lightgrey')
    ax.coastlines(lw=0.25, color='k')
    if extent is None:
        ax.set_extent([lng_n.min()-180., lng_n.max()-180., 
                       lat_n.min(), lat_n.max()])    
    else: 
        ax.set_extent(extent)
    # Plot filled contours
    cmap = plt.get_cmap(cmap)
    mb = ax.contourf(lng_n, lat_n, field_n, clevs, cmap=cmap, extend='both', 
                     transform=ccrs.PlateCarree())
    # If specified, add contours and labels. 
    if contour is None: pass
    else:
        cs = ax.contour(lng_n, lat_n, contour, contour_levs, colors = 'k', 
                        transform=ccrs.PlateCarree())
        plt.clabel(cs, fontsize=10, inline=1, fmt = '%1.0f')
    # If specified, add stippling (e.g., for statistical significance)
    if p_n is None: pass
    else:
        p_n = copy.deepcopy(p_n)
        p_n[p_n >= alpha] = np.nan
        p_n[np.isnan(p_n) == False] = 5.
        lng_n_gridded, lat_n_gridded = np.meshgrid(lng_n, lat_n)
        ax.scatter(lng_n_gridded[1::2,::3], lat_n_gridded[1::2,::3], 
                   s=p_n[1::2,::3], facecolor='k', lw=0, marker='.', 
                   transform=ccrs.PlateCarree())
    # If specified add scatterplots (e.g., for location of eddy-driven jet)
    if e_n is None: pass
    else:
        # Plot only every X number of longitude values
        skiplng = 6
        ax.errorbar(lng_n[::skiplng], e_n[::skiplng], 
                    yerr=eerr_n[::skiplng], zorder=12, color='k', markersize=2, 
                    elinewidth=0.5, ecolor='k', fmt='o', 
                    transform=ccrs.PlateCarree())
    # If specified, add quiver (e.g., for wind vectors)
    if quiver is None: pass
    else:    
        fu = quiver[0]
        fv = quiver[1]
        skip=(slice(None,None,3),slice(None,None,3))
        Q = ax.quiver(lng_n[skip[0]], lat_n[skip[1]], fu[skip], fv[skip], 
                      units='inches', pivot='middle', scale=150,
                      transform=ccrs.PlateCarree())
        plt.quiverkey(Q, 1.05, 0.05, 15, '15 m s$^{\mathregular{-1}}$', 
                      labelpos='E',coordinates='axes')
    # Add colorbar
    colorbar_axes = plt.gcf().add_axes([0.78,0.25,0.02,0.5])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
                            extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label(cbar_label, fontsize=14)
    plt.gcf().subplots_adjust(right=0.75, hspace=-0.2)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'map_%s_%s.eps'%(hemisphere, fstr), transparent=True)
    return 

def timeseries_seasonalo3():
    """this is a stand-alone function, meaning it loads all inputs it needs. 
    Function loads modeled O3 for the Northern and Southern hemispheres for 
    the entire year for the specified measuring period (i.e., 2004-2012) and  
    averages over the years to produce an average seasonal cycle. The seasonal
    cycles of O3 in the North America, Africa, the Southern Hemisphere, and 
    South America are plotted. 

    Parameters
    ----------
    None

    Returns
    -------
    None    
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    sys.path.append('/Users/ghkerr/phd/globalo3/')
    import globalo3_open
    import globalo3_calculate
    # Constants    
    years = [2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012]
    latmin_n, lngmin_n, latmax_n, lngmax_n = 0., 0., 90., 360.
    latmin_s, lngmin_s, latmax_s, lngmax_s = -90., 0., 0., 360.
    months_all = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 
                  'sep', 'oct', 'nov', 'dec']
    lat_gmi_nall, lng_gmi_nall, times_nall, o3_nall = \
        globalo3_open.open_overpass2_specifieddomain(years, months_all, 
        latmin_n, latmax_n, lngmin_n, lngmax_n, 'O3', 'HindcastMR2')
    lat_gmi_sall, lng_gmi_sall, times_sall, o3_sall = \
        globalo3_open.open_overpass2_specifieddomain(years, months_all, 
        latmin_s, latmax_s, lngmin_s, lngmax_s, 'O3', 'HindcastMR2')
    o3_nall = o3_nall*1e9    
    o3_sall = o3_sall*1e9
    # Separate continuous daily O3 into yearly chunks
    o3_nall_byyr = globalo3_calculate.separate_years(o3_nall, lat_gmi_nall, 
        lng_gmi_nall, years)
    o3_sall_byyr = globalo3_calculate.separate_years(o3_sall, lat_gmi_sall, 
        lng_gmi_sall, years)    
    o3_nall_byyr = np.mean(o3_nall_byyr, axis = 0)
    o3_sall_byyr = np.mean(o3_sall_byyr, axis = 0)
    # Find North America
    nal = geo_idx(240., lng_gmi_nall)
    nar = geo_idx(285., lng_gmi_nall)
    nau = geo_idx(49., lat_gmi_nall)
    nad = geo_idx(30., lat_gmi_nall)
    # Find South America
    sal = geo_idx(285., lng_gmi_sall)
    sar = geo_idx(315., lng_gmi_sall)
    sau = geo_idx(0., lat_gmi_sall)
    sad = geo_idx(-32., lat_gmi_sall)
    # Find Africa 
    al = geo_idx(11., lng_gmi_sall)
    ar = geo_idx(39., lng_gmi_sall)
    au = geo_idx(0., lat_gmi_sall)
    ad = geo_idx(-21., lat_gmi_sall)
    # Sample O3 and average over South American, African, and Northern/Southern 
    # Hemisphere domains
    nao3 = o3_nall_byyr[:, nad:nau+1, nal:nar+1]
    nao3 = np.mean(nao3, axis=tuple((1,2)))
    sao3 = o3_sall_byyr[:, sad:sau+1, sal:sar+1]
    sao3 = np.mean(sao3, axis=tuple((1,2)))
    ao3 = o3_sall_byyr[:, ad:au+1, al:ar+1]
    ao3 = np.mean(ao3, axis=tuple((1,2)))
    sho3 = np.mean(o3_sall_byyr, axis=tuple((1,2)))
    # Plotting
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))
    ax.plot(nao3, color='#1b9e77', lw=2., label='North America')
    ax.plot(sho3, lw=2., color='#d95f02', label='Southern Hemisphere')
    ax.plot(ao3, lw=2., color='#7570b3', label='Africa')
    ax.plot(sao3, lw=2., color='#e7298a', label='South America')
    ax.set_xlim([0, len(sho3)])
    ax.set_xticks([0,31,59,90,120,151,181,212,243,273,304,335])
    ax.set_xticklabels(['Jan', '', 'Mar', '', 'May', '', 'Jul', '', 'Sep', 
                        '', 'Nov', ''])
    ax.set_ylim([10, 55])
    ax.set_ylabel('O$_{\mathregular{3}}$ [ppbv]')
    plt.legend(ncol=2, fontsize=8, loc=8)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+'timeseries_seasonalo3.eps')    
    return 

def timeseries_rado3dt(t2m, stdo3, emfixo3, stdnox, emfixnox, years, lat, lng, 
    frl, fru, frr, frd, fstr):
    """function finds regionally-averaged 2-meter temperatures and O3 from 
    Strode et al. (2015) Std and EmFix simulations over the specified focus 
    region. dO3/dT values for each year are thereafter calculated and plotted. 
    
    Parameters
    ----------
    t2m : numpy.ndarray 
        Daily maximum 2-meter temperatures interpolated to the resolution of 
        the CTM, [time, lat, lng]
    stdo3 : numpy.ndarray 
        O3 concentrations at 1300-1400 hours local time from the 
        HindcastFFIgac2 simulation, units of ppbv, [time, lat, lng]            
    emfixo3 : numpy.ndarray 
        O3 concentrations at 1300-1400 hours local time from the 
        Hindcast3Igac2 simulation, units of ppbv, [time, lat, lng] 
    stdnox : numpy.ndarray 
        NOx concentrations at 1300-1400 hours local time from the 
        HindcastFFIgac2 simulation, units of ppbv, [time, lat, lng]        
    emfixnox : numpy.ndarray 
        NOx concentrations at 1300-1400 hours local time from the 
        Hindcast3Igac2 simulation, units of ppbv, [time, lat, lng]
    years : list
        Year or range of years in measuring period
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]
    frl : float
        Longitude coordinate of the left side (minimum) of the bounding box 
        containing the focus region, units of degrees east        
    fru : float 
        Latitude coordinate of the top side (maximum) of the bounding box 
        containing the focus region, units of degrees north    
    frr : float 
        Longitude coordinate of the right side (maximum) of the bounding box 
        containing the focus region, units of degrees east        
    frd : float
        Latitude coordinate of the bottom side (minimum) of the bounding box 
        containing the focus region, units of degrees north    
    fstr : str
        Output filename suffix (i.e. should specify which focus region is 
        shown in time series)

    Returns
    -------
    None        
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    sys.path.append('/Users/ghkerr/phd/globalo3/')    
    import globalo3_calculate    
    # Define focus region (i.e., Eastern U.S., Europe, etc.)
    frl = geo_idx(frl, lng)
    frr = geo_idx(frr, lng)
    fru = geo_idx(fru, lat)
    frd = geo_idx(frd, lat)
    # Sample O3, NOx, and temperature; average over focus region 
    # For Std O3
    frstdo3 = stdo3[:, frd:fru+1, frl:frr+1]
    frstdo3 = np.mean(frstdo3, axis=tuple((1,2)))
    frstdo3 = np.expand_dims(frstdo3, axis=1)
    frstdo3 = np.expand_dims(frstdo3, axis=1)
    # For EmFix O3
    fremfixo3 = emfixo3[:, frd:fru+1, frl:frr+1]
    fremfixo3 = np.mean(fremfixo3, axis=tuple((1,2)))
    fremfixo3 = np.expand_dims(fremfixo3, axis=1)
    fremfixo3 = np.expand_dims(fremfixo3, axis=1)
    # For temperature
    frt2m = t2m[:, frd:fru+1, frl:frr+1]
    frt2m = np.mean(frt2m, axis=tuple((1,2)))
    frt2m = np.expand_dims(frt2m, axis=1)
    frt2m = np.expand_dims(frt2m, axis=1)
    # For Std NOx
    frstdnox = stdnox[:, frd:fru+1, frl:frr+1]
    frstdnox = np.mean(frstdnox, axis=tuple((1,2)))  
    frstdnox = globalo3_calculate.separate_years(frstdnox, np.array([0.]), 
        np.array([0.]), years)
    frstdnox = np.mean(frstdnox[:,:,0,0], axis = 1)
    # Find dO3/dT with regionally-averaged data
    (frstd_byyr, frstd_ls, frstd_lsp, frstd_mkz, frstd_mkp) = \
        globalo3_calculate.calculate_trends_do3dt(frt2m, frstdo3, 
        np.array([0.]), np.array([0.]), years)
    (fremfix_byyr, fremfix_ls, fremfix_lsp, fremfix_mkz, fremfix_mkp) = \
        globalo3_calculate.calculate_trends_do3dt(frt2m, fremfixo3, 
        np.array([0.]), np.array([0.]), years)
    # Plotting
    fig = plt.figure()
    ax = plt.subplot2grid((1,1),(0,0)) 
    ax.plot(frstd_byyr[:,0,0], color='#e41a1c', lw=2, label='')
    ax.plot(fremfix_byyr[:,0,0], color='#377eb8', lw=2, label='')
    ax.set_xlim([0, len(years)-1])
    ax.set_xticks(np.arange(0, len(years), 1))
    ax.set_xticklabels(years)
    ax.set_ylabel('dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]')
    axt = plt.twinx()
    axt.plot(frstdnox, color='k', lw=2)
    axt.set_ylabel('<NO$_{x}$> [ppbv]')
    # Insert table containing slope and significance  
    table_vals=[['', 'm', 'p'],
                ['Std', '%.2f' %frstd_ls, '%.2f' %frstd_lsp],
                ['EmFix', '%.2f' %fremfix_ls, '%.2f' %fremfix_lsp]]
    table = plt.table(cellText=table_vals, colWidths = [0.13]*3, 
        loc='best')
    # Change color of table's cells to act as a legend for the timeseries
    table._cells[(1,0)]._text.set_color('#e41a1c')
    table._cells[(1,1)]._text.set_color('#e41a1c')
    table._cells[(1,2)]._text.set_color('#e41a1c')
    table._cells[(2,0)]._text.set_color('#377eb8')
    table._cells[(2,1)]._text.set_color('#377eb8')
    table._cells[(2,2)]._text.set_color('#377eb8')
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'timeseries_rado3dt_%s.eps'%fstr)
    return 

def map_do3dt_ro3t_nh(do3dt, r_t2mo3, lat, lng, fstr): 
    """plot map of mean JJA dO3/dT (filled contours) and the Pearson 
    correlation coefficient calculated between temperature and O3 (contours)
    for the Northern Hemisphere. 
    
    Parameters
    ----------
    do3dt : numpy.ndarray     
        Northern Hemisphere dO3/dT, units of ppbv K-1, [lat, lng]
    r_t2mo3 : numpy.ndarray     
        Pearson correlation coefficient between temperature and O3, [lat, lng]        
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]
    fstr : str
        Output filename suffix (i.e. should specify which focus region is 
        shown in time series)

    Returns
    -------
    None  
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    plt.figure()
    axt = plt.subplot2grid((1,1), (0,0), 
                           projection=ccrs.Robinson(central_longitude=0.))
    # Filled contour
    mb = axt.contourf(lng, lat, do3dt, np.linspace(0, 3., 7), 
                      cmap=plt.get_cmap('PuBu'), transform=ccrs.PlateCarree(), 
                      extend='both')
    # Contour
    levels = [0., 0.25, 0.5, 0.75]
    colors =['#fecc5c', '#fd8d3c', '#f03b20', '#bd0026']
    cs = axt.contour(lng, lat, r_t2mo3, levels, colors=colors, 
                     linewidths = 0.5, transform=ccrs.PlateCarree())
    # Legend for contours
    labels = ['$r\:$=$\:$0.0', '$r\:$=$\:$0.25','$r\:$=$\:$0.5',
              '$r\:$=$\:$0.75']
    for i in range(len(labels)):
        cs.collections[i].set_label(labels[i])
    plt.legend(loc=8, bbox_to_anchor=(0.5,-0.6), ncol=2, frameon=False)
    # Add colorbar for contourf
    colorbar_axes = plt.gcf().add_axes([0.82,0.25,0.03,0.5])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical')
    colorbar.set_label('dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]', 
                       fontsize=12)
    plt.gcf().subplots_adjust(right=0.8, hspace=-0.2)
    # Coastlines/oceans
    axt.add_feature(cfeature.OCEAN, zorder=10, lw=0.0, color='lightgrey')
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'map_do3dt_ro3t_nh_%s.eps'%fstr)    
    return 
    
def map_extent(extent, title, cbar_label, cmap, clevs, region, fstr, 
    lat_ctm=None, lng_ctm=None, field_ctm=None, lat_obs=None, lng_obs=None, 
    field_obs=None, extend='neither'):
    """plot desired quantity over the Northern or Southern hemisphere (however, 
    if variable 'extent' is specified, map can be tailored to any region).
    
    Parameters
    ----------
    extent : list
        The extent (x0, x1, y0, y1) of the map in the given coordinate system; 
        note that the longitude coordinates should be in (-180-180˚).
    title : str
        Title for plot        
    cbar_label : str
        Label for the colorbar (units)
    cmap : str
        Colormap name
    clevs : numpy.ndarray
        Filled contour levels values     
    region : str
        Focus region for output filename (i.e., 'nh' for northern hemisphere, 
        'sh' for southern hemisphere, 'northamerica', etc.)
    fstr : str
        Output filename suffix
    lat_ctm : numpy.ndarray
        Gridded latitude coordinates of CTM, units of degrees north, [lat,]
    lng_ctm : numpy.ndarray
        Gridded longitude coordinates of CTM, units of degrees east, [lng,]        
    field_ctm : numpy.ndarray
        Gridded field of interest from CTM, [lat, lng]        
    lat_obs : list
        Latitude coordinates of observational stations, units of degrees north, 
        [stations,]
    lng_obs : list
        Longitude coordinates of observational stations, units of degrees east, 
        [stations,]    
    field_obs : list
        Field of interest at each observational station, [stations,]
    extend : str
        Extend settings for matplotlib colormap/colorbar (i.e., 'neither', 
        'both', 'min', 'max')

    Returns
    -------
    None              
    """
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure(figsize=(8,3.5))
    ax=plt.subplot2grid((1,1), (0,0), colspan=2,
                        projection=ccrs.Miller(central_longitude=0.))
    ax.set_title(title, fontsize=14, x=0.02, ha='left')    
    ax.add_feature(cfeature.OCEAN, lw = 0.0, color='lightgrey', zorder = 10)
    ax.coastlines(lw=0.5, color='k', zorder = 15)
    ax.set_extent(extent)
    # Create colormap with discrete colorlevels (for scatterplot)
    cmap = plt.get_cmap(cmap)
    # Extract all colors from the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    ## Force the first color entry to be grey
    #cmaplist[0] = (.5, .5, .5, 1.0)
    # Create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)
    norm = mpl.colors.BoundaryNorm(clevs, cmap.N)
    # Contourf for CTM/reanalysis field
    if field_ctm is None: pass
    else: 
        mb = ax.contourf(lng_ctm, lat_ctm, field_ctm, clevs, cmap = cmap,
                         transform = ccrs.PlateCarree(), extend = extend, 
                         zorder = 5)
    # Scatterpoints for observational field 
    if field_obs is None: pass
    else: 
        mb = ax.scatter(lng_obs, lat_obs, c = field_obs, cmap = cmap, 
                        edgecolor = 'k', linewidth = 0.5, s = 20, 
                        vmin = clevs[0], vmax = clevs[-1], norm = norm, 
                        transform = ccrs.PlateCarree(), zorder = 20)
    # Add colorbar
    colorbar = plt.colorbar(mb, orientation='vertical', extend = extend)
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label(cbar_label, fontsize=14)
    plt.gcf().subplots_adjust(right=0.8)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'map_%s_%s.eps'%(region, fstr), transparent=True)
    return

def scatter_latitude_o3u500(lat, lng, U500, o3, c, cmap, cbar_label, clevs, 
    region):
    """function plots O3 and the U500 within a given region as a function of
    latitude. U500 is zonally-averaged and scatterpoints corresponding to O3 
    represent all grid cells in the region. The colorscale of O3 corresponds 
    to the O3-climate penalty.
    
    Parameters
    ----------  
    lat : numpy.ndarray
        Latitude coordinates in region, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates in region, units of degrees east, [lng,]
    U500 : numpy.ndarry
        U wind at 500 hPa in region, units of m s-1, [time, lat, lng]
    o3 : numpy.ndarray
        Surface-level O3 in region, units of ppbv, [time, lat, lng]
    c : numpy.ndarray     
        The O3-climate penalty (dO3/dT) in region, units of ppbv K-1, [lat, 
        lng]
    cmap : str
        Colormap name        
    cbar_label : str
        Label for the colorbar (field and units)        
    clevs : numpy.ndarray
        Contour levels values        
    region : str
        Focus region, used for naming of output file

    Returns
    -------
    None       
    """
    import numpy as np
    import matplotlib.pyplot as plt
    lng_m, lat_m = np.meshgrid(lng,lat)
    # Calculate zonally-averaged U wind at 500 hPa
    U500_za = np.nanmean(U500, axis=tuple((0,2)))
    # n.b. Python crashes if the whole flattened grid is plotted (26000 
    # points), skip every N points in flattened grid
    skippoints = 1
    # Flatten latitude grid for plotting
    lat_f = lat_m.flatten()[::skippoints]
    o3_f = np.mean((o3),axis=0).flatten()[::skippoints]
    c_f = c.flatten()[::skippoints]
    # Plotting O3 against latitude
    fig = plt.figure()
    ax = plt.subplot2grid((1,1),(0,0))
    mb = ax.scatter(lat_f, o3_f, c=c_f, s=8, 
                    cmap=cmap_discretize(plt.get_cmap(cmap), len(clevs)-1), 
                    vmin=clevs[0], vmax=clevs[-1])
    #plt.colorbar(mb)
    colorbar_axes = plt.gcf().add_axes([0.83,0.15,0.03,0.7])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
                            extend='both')
    colorbar.set_label(cbar_label)
    colorbar.set_ticks(clevs)
    colorbar.set_ticklabels(clevs)
    ax.set_xlim([20,70])
#    ax.set_xlim([0,30])    
    ax.set_xlabel('Latitude [$^{\mathregular{\circ}}$]')
    ax.set_ylabel('O$_{\mathregular{3}}$ [ppbv]')
    # Plotting zonally-averaged U500 against latitude
    ax2 = ax.twinx()
    ax2.plot(lat, U500_za, '-k')
    ax2.set_ylabel('U$_{\mathregular{500\:hPa}}$ [m s$^{\mathregular{-1}}$]')
    plt.gcf().subplots_adjust(right=0.7, hspace=-0.2)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'scatter_latitude_o3u500_%s.eps'%region, dpi=300)
    return 

def map_toar(data, lat, lng, frl, fru, frr, frd, res, vmin, vmax, nlevs, 
    title, cbar_label, fstr):
    """function plots a map of TOAR data over a specified region. 

    Parameters
    ----------
    data : numpy.ndarray
        TOAR data to be plotted, [lat, lng]    
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]
    frl : float
        Longitude coordinate of the left side (minimum) of the bounding box 
        containing the focus region that map's extent will be set to, units of 
        degrees east        
    fru : float 
        Latitude coordinate of the top side (maximum) of the bounding box 
        containing the focus region  that map's extent will be set to, units of 
        degrees north    
    frr : float 
        Longitude coordinate of the right side (maximum) of the bounding box 
        containing the focus region that map's extent will be set to, units of 
        degrees east        
    frd : float
        Latitude coordinate of the bottom side (minimum) of the bounding box 
        containing the focus region  that map's extent will be set to, units of 
        degrees north  
    res : int 
        Resolution of TOAR gridded set (resolutions of 10°×10°, 5°×5°, and 
        2°×2° are available). Recommended use of the 5° longitude ×5° latitude 
        products is encouraged as they provide a reasonable compromise between 
        global coverage and regional differentiation.    
    vmin : int
        Lower colorbar limit
    vmax : int
        Upper colorbar limit
    nlevs : int
        Number of levels on colorbar 
    title : str
        Title for plot
    cbar_label : str
        Label for the colorbar (field and units)
    fstr : str
        Output filename suffix (should specify which TOAR metric and region 
        are plotted)     

    Returns
    -------
    None        
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.colors
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    # Set-up colormap
    cmap = cmap_discretize('YlGnBu', nlevs)
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    # Initialize figure, axis
    fig = plt.figure()
    ax = plt.subplot2grid((1,1), (0,0), colspan=1, 
        projection=ccrs.Miller(central_longitude=(frl+frr)/2.))    
    ax.set_title(title, fontsize=18, x=0.02, ha='left')
    ax.add_feature(cfeature.OCEAN, zorder=1, lw=0.0, color='lightgrey')
    ax.add_feature(cfeature.BORDERS, zorder=2, edgecolor='k', linewidth=.25)
    ax.coastlines(lw=0.25, color='k')
    ax.set_extent([frl,frr,frd,fru])
    cmap = plt.get_cmap(cmap)
    # Add patches for non-NaN grid cells
    for i, lat_gb in enumerate(lat):
        for j, lon_gb in enumerate(lng):
            if np.isnan(data[i, j]) != True:
                ax.add_patch(mpatches.Rectangle(xy=[lon_gb, lat_gb], width=res, 
                    height=res, facecolor=cmap(norm(data[i,j])), 
                    edgecolor='lightgrey', lw = 0.25, alpha=1.,
                    transform=ccrs.PlateCarree(), zorder=10))
    # Add colorbar 
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    colorbar_axes = plt.gcf().add_axes([0.82,0.25,0.03,0.5])
    colorbar = plt.colorbar(sm, colorbar_axes, extend='both', 
                            orientation='vertical')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label(cbar_label, fontsize=16)
    plt.subplots_adjust(right=0.8)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'map_toar_%s.eps' %fstr, dpi=300, transparent=True)
    return 

def timeseries_o3atabovebelow_jet(o3, lat_jet, lat, lng, lngl, lngr):
    """for O3 and the latitude of the eddy-driven jet in a given region, plot 
    a timeseries of O3 at (with +/- 2 degrees) of the mean jet location, O3 
    above (8-12 degrees above) the mean jet location, and O3 below (8-12 
    degrees below) the mean jet location. A timeseries of the regionally-
    averaged jet latitude is plotted for each subplot. 
    
    Parameters
    ----------
    o3 : numpy.ndarray
        O3 in region, units of ppbv, [time, lat, lng]  
    lat_jet : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in region, units of degrees north[time, lng]
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]
    lngl : float
        Longitude coordinate corresponding to the west side of the focus 
        region; n.b., lngl < lngr
    lngr : float
        Longitude coordinate corresponding to the east side of the focus region
        
    Returns
    -------
    None
    """
    import numpy as np
    import matplotlib.pyplot as plt
    # Find mean jet position in Eastern North America
    lat_jet_fr = lat_jet[:, np.abs(lng-lngl).argmin(): 
        np.abs(lng-lngr).argmin()]
    lat_jet_fr = np.nanmean(lat_jet_fr, axis=1) 
    lat_jet_mean = np.nanmean(lat_jet_fr)
    # O3 at the mean position of the jet (+/- 3 degrees)
    o3at = o3[:, np.abs(lat-lat_jet_mean).argmin()-2:
        np.abs(lat-lat_jet_mean).argmin()+2, np.abs(lng-lngl).argmin(): 
        np.abs(lng-lngr).argmin()]
    # O3 above the mean position of the jet (7-13 degrees above)    
    o3above = o3[:, np.abs(lat-lat_jet_mean).argmin()+8:
        np.abs(lat-lat_jet_mean).argmin()+12, np.abs(lng-lngl).argmin(): 
        np.abs(lng-lngr).argmin()]    
    # O3 below the mean position of the jet (7-13 degrees below)    
    o3below = o3[:, np.abs(lat-lat_jet_mean).argmin()-12:
        np.abs(lat-lat_jet_mean).argmin()-8, np.abs(lng-lngl).argmin(): 
        np.abs(lng-lngr).argmin()]        
    # Regional averaging 
    o3at = np.nanmean(o3at, axis=tuple((1,2)))
    o3above = np.nanmean(o3above, axis=tuple((1,2)))
    o3below = np.nanmean(o3below, axis=tuple((1,2)))
    # Plotting
    fig = plt.figure()
    ax1 = plt.subplot2grid((3,3), (0,0), colspan=3)
    ax1b = ax1.twinx()
    ax2 = plt.subplot2grid((3,3), (1,0), colspan=3)
    ax2b = ax2.twinx()
    ax3 = plt.subplot2grid((3,3), (2,0), colspan=3)
    ax3b = ax3.twinx()
    # Above jet
    ax1b.axhspan(lat_jet_mean+8, lat_jet_mean+12, xmin=ax1b.get_xlim()[0], 
        xmax=ax1b.get_xlim()[1], alpha=0.2, color='dodgerblue', zorder=0)
    ax1.plot(o3above[:92], '-k')
    ax1b.plot(lat_jet_fr[:92], ls='-', color='dodgerblue')
    # At jet
    ax2b.axhline(y=lat_jet_mean, xmin=ax2b.get_xlim()[0], 
        xmax=ax2b.get_xlim()[1], color='dodgerblue', ls='--', alpha=0.5, 
        zorder=2)
    ax2b.axhspan(lat_jet_mean-2, lat_jet_mean+2, xmin=ax2b.get_xlim()[0], 
        xmax=ax2b.get_xlim()[1], alpha=0.2, color='dodgerblue', zorder=0)
    ax2b.plot(lat_jet_fr[:92], ls='-', color='dodgerblue', zorder=5)
    ax2.plot(o3at[:92], '-k', zorder=5)
    # Below jet
    ax3b.axhspan(lat_jet_mean-12, lat_jet_mean-8, xmin=ax3b.get_xlim()[0], 
        xmax=ax3b.get_xlim()[1], alpha=0.2, color='dodgerblue', zorder=0)
    ax3b.plot(lat_jet_fr[:92], ls='-', color='dodgerblue')
    ax3.plot(o3below[:92], '-k')
    # Set limits, labels
    for ax in [ax1, ax1b, ax2, ax2b, ax3, ax3b]:
        ax.set_xlim([0, 91])
        ax.set_xticks([0, 14, 30, 44, 61, 75])  
        ax.set_xticklabels([''])
    for ax in [ax1, ax2, ax3]:
        ax.set_ylim([15, 60])
        ax.set_yticks([15, 30, 45, 60])
        ax.set_ylabel('O$_{\mathregular{3}}$ [ppbv]')
    for ax in [ax1b, ax2b, ax3b]:
        ax.set_ylim([30, 60])
        ax.set_yticks([30, 40, 50, 60])
        ax.set_ylabel('Jet position [$^{\mathregular{\circ}}$]')  
    ax3.set_xticklabels(['1 June 2008', '', '1 July', '', '1 Aug', ''], 
        ha = 'center', fontsize = 12)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'timeseries_o3atabovebelow_jet.png', dpi=300)
    return 

def contourf_var_atcenter(center, var, lat, lng, varname, exam_rad, 
    kind, cbar_label, clevs, cmap, fstr, SLP=None, H=None): 
    """function retrieves variable of interest in the vincinity of 
    (anti)cyclones and plots the mean O3 concentrations averaged over all 
    systems in region/measuring period. Only (anti)cyclones over land are 
    considered as part of the average.

    Parameters
    ----------
    center : numpy.ndarray 
        A value of 1 indicates the presence of a cyclone/anticylone for a 
        particular day and location, [time, lat, lng]       
    var : numpy.ndarray
        Variable of interest in region, units of ppbv, [time, lat, lng]
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]
    lng : numpy.ndarray
        Longitude coordinates, units of degrees east, [lng,]
    varname : str
        Variable name for output file
    exam_rad : int
        Search radius; the number of grid cells on each side of the (anti)
        cyclone's center over which O3 concentrations will be retrieved
    kind : str
        'cyclone' or 'anticyclone' (for axis labels)
    cbar_label : str
        Label for the colorbar (field and units)
    clevs : numpy.ndarray
        Contour levels values
    cmap : str
        Colormap name        
    fstr : str
        Output filename suffix (should specify the type of system in variable
        'kind' and the latitude) 
    SLP: numpy.ndarray
        Sea level pressure in region, units of Pa, [time, lat, lng]
    H : numpy.ndarray
        Geopotential height at 500 hPa in region, units of m, [time, lat, lng]
        
    Returns
    -------
    None
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    bm = Basemap()    
    # Identify locations of (anti)cyclones
    where_center = np.where(center==1.)
    # Find O3, SLP, and  surrounding (within +/- exam_rad) systems for each 
    # (anti)cyclone
    var_atcenter = np.empty(shape=(where_center[0].shape[0], exam_rad, 
                                   exam_rad))
    var_atcenter[:] = np.nan
    SLP_atcenter = np.empty(shape=(where_center[0].shape[0], exam_rad, 
                                   exam_rad))
    SLP_atcenter[:] = np.nan
    H_atcenter = np.empty(shape=(where_center[0].shape[0], exam_rad, 
                                 exam_rad))
    H_atcenter[:] = np.nan            
    for system in np.arange(0, len(where_center[0]), 1):
        system_coords = (where_center[0][system], 
                         where_center[1][system], 
                         where_center[2][system])
        # (xpt, ypt)
        if bm.is_land(lng[where_center[2][system]]-360., 
                      lat[where_center[1][system]]) == True:   
            var_atcenter[system]= var[system_coords[0], 
                     system_coords[1]-(int(np.round(exam_rad/2))-1):
                     system_coords[1]+int(np.round(exam_rad/2)),
                     system_coords[2]-(int(np.round(exam_rad/2))-1):
                     system_coords[2]+int(np.round(exam_rad/2))]
            if SLP is not None: 
                # Convert from Pa to hPa
                SLP_atcenter[system]= SLP[system_coords[0], 
                     system_coords[1]-(int(np.round(exam_rad/2))-1):
                     system_coords[1]+int(np.round(exam_rad/2)),
                     system_coords[2]-(int(np.round(exam_rad/2))-1):
                     system_coords[2]+int(np.round(exam_rad/2))]/100.
            if H is not None: 
                # Convert from m to dm 
                H_atcenter[system]= H[system_coords[0], 
                     system_coords[1]-(int(np.round(exam_rad/2))-1):
                     system_coords[1]+int(np.round(exam_rad/2)),
                     system_coords[2]-(int(np.round(exam_rad/2))-1):
                     system_coords[2]+int(np.round(exam_rad/2))]/10.            
    # Plotting 
    fig = plt.figure()
    ax = plt.subplot2grid((1,1),(0,0))
    mb = ax.contourf(np.nanmean(var_atcenter, axis=0), clevs, extend='both',
                     cmap=plt.get_cmap(cmap))
    
    # Add contours for SLP and H, if needed
    if SLP is None: pass
    else:
        cs = ax.contour(np.nanmean(SLP_atcenter, axis=0), colors='k')
        plt.clabel(cs, fontsize=12, inline=1, fmt='%1.0f')    
    if H is None: pass
    else:
        cs = ax.contour(np.nanmean(H_atcenter, axis=0), colors='w')
        plt.clabel(cs, fontsize=12, inline=1, fmt='%1.0f')    
    # Add colorbar
    colorbar_axes = plt.gcf().add_axes([0.83,0.25,0.02,0.5])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical',
                            extend='max')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label(cbar_label, fontsize=14)
    # Aesthetics (crosshairs for system's center, axes labels)
#    ax.axhline(y=(int(np.round(exam_rad/2))-1), xmin=ax.get_xlim()[0], 
#               xmax=ax.get_xlim()[1], ls='--', lw=2., color='k')
#    ax.axvline(x=(int(np.round(exam_rad/2))-1), ymin=ax.get_ylim()[0], 
#               ymax=ax.get_ylim()[1], ls='--', lw=2., color='k') 
    ax.set_xticks(np.arange(0, exam_rad, 1))
    ax.set_xticklabels(np.arange(-(int(np.round(exam_rad/2))-1), 
                       int(np.round(exam_rad/2)), 1), fontsize=12)
    ax.set_xlabel('E-W grid cells from %s center'%kind, fontsize=14)
    ax.set_yticks(np.arange(0, exam_rad, 1))
    ax.set_yticklabels(np.arange(-(int(np.round(exam_rad/2))-1), 
                       int(np.round(exam_rad/2)), 1), fontsize=12)
    ax.set_ylabel('N-S grid cells from %s center'%kind, fontsize=14)
    plt.gcf().subplots_adjust(right=0.8)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'contourf_%s_at%s.eps' %(varname,fstr), dpi=300)
    return

def map_jet_centerdist(centers, lat_jet, lng_jetcoords, lat_centercoords, 
    lng_centercoords, kind): 
    """for a given day/longitude, function finds all (anti)cyclones at that 
    longitude and calculates the difference between the latitude of the (anti)-
    cyclones and the jet latitude (n.b. positive values imply that the (anti)-
    cyclone is above the jet on a particular day). Plotted values are (a) the 
    location of the centers in map coordinates and the mean jet location and 
    variability and (b) the difference between centers' latitudes and the 
    jet versus longitude. 
    
    Parameters
    ----------
    center : numpy.ndarray 
        A value of 1 indicates the presence of a cyclone/anticylone for a 
        particular day and locaditon, [time, lat_centercoords, 
        lng_centercoords]  
    lat_jet : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in region, units of degrees north[time, lng_jetcoords]
    lng_jetcoords : numpy.ndarray
        Longitude coordinates corresponding to the jet array, units of degrees 
        east, [lng_jetcoords,]
    lat_centercoords : numpy.ndarray 
        Latitude coordinates corresponding to the (anti)cyclone array, units of 
        degrees north, [lat_centercoords,]        
    lng_centercoords : numpy.ndarray 
        Longitude coordinates corresponding to the (anti)cyclone array, units 
        of degrees east, [lng_centercoords,]
    kind : str
        'cyclone' or 'anticyclone' (for output filename suffix) 

    Returns
    -------
    None    
    """
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure()
    axt = plt.subplot2grid((2, 2), (0, 0), colspan = 2, 
                          projection=ccrs.Miller(central_longitude=0.))
    axb = plt.subplot2grid((2, 2), (1, 0), colspan = 2)
    # Create 2D array of lists such that daily values for fields can be 
    # appended to the appropriate date/longitude index
    diff = np.empty(shape=lat_jet.shape, dtype=object)
    for i in np.arange(0, diff.shape[0], 1):
        for j in np.arange(0, diff.shape[1], 1):
            diff[i,j]=[]
    # Loop through longitudes 
    for i, lng in enumerate(lng_jetcoords):
        for day in np.arange(0, len(lat_jet), 1):
            # Find the latitude of the jet at longitude/on day of interest
            jetlat_atlng = lat_jet[day,i]
            # Find the index cooresponding to the closest longitude in the 
            # (anti)cyclone center dataset (n.b., as long as the resolution/
            # longitudinal span is the same for both the jet/center dataset, 
            # this will be the same value)
            lngwhere_centercoords = np.abs(lng_centercoords-lng).argmin()
            # (Anticyclones) at longitude/on day of interest 
            centers_atlng = centers[day, :, lngwhere_centercoords]
            centers_atlng = np.where(centers_atlng==1.)[0]
            centerlats_atlng = lat_centercoords[centers_atlng]
            # Find the distance (in degrees) between the jet and center with 
            # positive (negative) distances defined as the center above (below)
            # the jet
            distance = centerlats_atlng - jetlat_atlng
            # Save off values
            diff[day,i].append(list(distance))
            # Plotting 
            if centerlats_atlng.shape[0] > 0:
                lng_forplotting = np.empty(shape=centerlats_atlng.shape)
                lng_forplotting[:] = np.mod((lng+180),360)-180
                axt.plot(lng_forplotting, centerlats_atlng, 'ko', markersize=1,
                         transform=ccrs.PlateCarree())
                axb.plot(lng_forplotting, distance, 'ko', markersize=1)
    # Add mean jet position and variability 
    skiplng = 6
    axt.errorbar(lng_jetcoords[::skiplng], np.nanmean(lat_jet, axis=0)[::skiplng], 
                yerr=np.std(lat_jet,axis=0)[::skiplng], zorder=12, color='r', 
                markersize=2, elinewidth=0.5, ecolor='r', fmt='o', 
                transform=ccrs.PlateCarree())
    axt.coastlines(lw=0.25, color='k')
    # Aesthetics
    axb.set_xlim([-180, 180])
    axb.set_xlabel('Longitude [$^{\circ}$E]')
    axb.set_ylabel('$\mathregular{\phi}_{\mathregular{%s}}$ - '%kind+
                   '$\mathregular{\phi}_{\mathregular{jet}}$ [$^{\circ}$]')
    plt.subplots_adjust(hspace=0.1)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'map_jet_centerdist_%s.eps'%kind, dpi=300)
    return 

def fieldatjet(lat_fr, lng_fr, lat_jet, field_jet, cmap, cbar_label, clevs, 
    fieldstr, fstr, skiplng=6):
    """function plots the mean position and variability (standard deviation) 
    of the eddy-driven jet on a map of the mid-latitudes of the Northern 
    Hemisphere (upper subplot) and a given field depicting a variable equator-
    and poleward of the latitude of the jet averaged over all days. 
    
    Parameters
    ----------
    lat_fr : numpy.ndarray
        Latitude coordinates for the focus region, units of degrees north, 
        [lat,]                
    lng_fr : numpy.ndarray
        Longitude coordinates for the focus region, units of degrees east, 
        [lng,]          
    lat_jet : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in region, units of degrees north[time, lng]
    field_jet : numpy.ndarray
        The value of the field at the jet and within +/- jetdistance of jet, 
        [lat, lng] or [time, lat, lng]
    cmap : str
        Colormap name    
    cbar_label : str
        Label for the colorbar (field and units)
    clevs : numpy.ndarray
        Contour levels values
    fieldstr : str
        For naming the output file; this variable should specify which 
        field is shown (i.e., o3anom, do3dt, etc.)
    fstr : str
        For naming the output file; this variable should specify the model 
        shown (i.e., GEOS-C1SD, etc. )
    skiplng : int
        This number of longitude coordinates will be skipped when the mean 
        position and variability of the eddy-driven jet is plotted (default is
        6, which looks visually nice for 1˚ resolution output)

    Returns
    -------
    None    
    """
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure()
    axt = plt.subplot2grid((2,2), (0,0), colspan=2, projection=
                           ccrs.Miller(central_longitude=0.))
    axt.set_title('(a)', fontsize=14, x=0.02, ha='left')    
    axb = plt.subplot2grid((2,2), (1,0), colspan=2)
    axb.set_title('(b)', fontsize=14, x=0.02, ha='left')    
    axt.add_feature(cfeature.OCEAN, zorder=1, lw=0.0, color='lightgrey')
    axt.coastlines(lw=0.25, color='k')
    axt.set_extent([lng_fr.min()-180, lng_fr.max()-180, lat_fr.min(), 
                    lat_fr.max()])    
    axt.set_xticks([-180,-135,-90,-45,0,45,90,135,180])
    axt.set_xticklabels([''])
    # Plot only every X number of longitude values
    axt.errorbar(lng_fr[::skiplng], np.nanmean(lat_jet, axis=0)[::skiplng], 
                 yerr=np.nanstd(lat_jet, axis=0)[::skiplng], zorder=10, 
                 color='k', markersize=2, elinewidth=0.5, ecolor='k', 
                 fmt='o', transform=ccrs.PlateCarree())
    # Roll O3 differences arond the jet so that they align with the top (map)
    # subplot
    field_jet_roll = np.roll(field_jet, int(lng_fr.shape[0]/2.), axis=1)
    mb = axb.contourf(lng_fr, np.arange(-int((field_jet_roll.shape[0]-1)/2.), 
                      int((field_jet_roll.shape[0]-1)/2.)+1, 1),
                      field_jet_roll, clevs, cmap=plt.get_cmap(cmap), 
                      extend='both')
    axb.set_xlabel('Longitude [$^{\circ}$]', fontsize=14)
    axb.set_xticks(np.linspace(0, 360, 9))
    axb.set_xticklabels([-180,-135,-90,-45,0,45,90,135,180], fontsize=12)
    axb.set_ylabel('Distance from jet center [$^{\mathregular{\circ}}$]', 
                   fontsize=14)
    axb.tick_params(labelsize=12)
    # Add colorbar
    colorbar_axes = plt.gcf().add_axes([0.78,0.135,0.02,0.375])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label(cbar_label, fontsize=14)
    ## Optional: this is clunky and should be improved to more clearly indicate
    ## the location of shorelines in (b)
    ## Find locations of shorelines by finding fraction of NaNs (ocean) grid 
    ## cells in a given longitudinal transect 
    #ocean_frac = np.empty(shape=(lng_fr.shape[0]))
    #ocean_frac[:] = np.nan
    #land_fr_roll = np.roll(land_fr, shift=int(len(lng_fr)/2), axis=1)
    #for i in np.arange(0, len(lng_fr), 1):
    #    land_fr_transect = land_fr_roll[:, i]
    #    frac = np.where(np.isnan(land_fr_transect)==
    #                    True)[0].shape[0]/len(land_fr_roll[:, i])
    #    ocean_frac[i] = frac
    #del i
    ## Find where at least 50% of transect are land-based grid cells
    #whereshore = np.where((ocean_frac>0.47) & (ocean_frac<0.5))[0]
    ## Manually remove duplicates and close doubles (this is clunky!)
    #whereshore = np.delete(whereshore, [1,3,5,6])
    #for ws in whereshore: 
    #    axb.axvline(x=lng_fr[ws], c='k', lw=0.75)
    plt.gcf().subplots_adjust(right=0.75)#, hspace=-0.1)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'fieldatjet_%s_%s.eps'%(fieldstr, fstr), dpi=300)
    return 

def zonalavg_byregion_rt2mo3_do3dt2m(r_t2mo3, do3dt2m, lat_gmi, lng_gmi, r_aqs, 
    do3dt2m_aqs, lat_aqs, lng_aqs, r_naps, do3dt2m_naps, lat_naps, lng_naps, 
    r_emep, do3dt2m_emep, lat_emep, lng_emep, lng_ml, lat_jet_ml, fstr): 
    """find zonally-averaged mean values and variability of dO3/dT and r(T, O3)
    in Western North America, Eastern North America, Europe, and East Asia
    from both CTM simulations and observations and plot as a functino of 
    latitude alongside the mean position and variability of the eddy-driven 
    jet. 
    
    Parameters
    ----------
    r_t2mo3 : numpy.ndarray     
        Pearson correlation coefficient between temperature and O3, [lat, lng]        
    do3dt2m : numpy.ndarray     
        Northern Hemisphere dO3/dT, units of ppbv K-1, [lat, lng]
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    r_aqs : list
        r(T, O3) at each AQS station, [stations,]
    do3dt_aqs : list 
        dO3/dT at each AQS station, units of ppbv K-1, [stations,]
    lat_aqs : list 
        Latitude coorindate at each AQS station, units of degrees north, 
        [stations,]
    lng_aqs : list 
        Longitude coorindate at each AQS station, units of degrees east, 
        [stations,]
    r_naps : list
        r(T, O3) at each NAPS station, [stations,]
    do3dt_naps : list 
        dO3/dT at each NAPS station, units of ppbv K-1, [stations,]
    lat_naps : list 
        Latitude coorindate at each NAPS station, units of degrees north, 
        [stations,]
    lng_naps : list 
        Longitude coorindate at each NAPS station, units of degrees east, 
        [stations,]        
    r_emep : list
        r(T, O3) at each EMEP station, [stations,]
    do3dt_emep : list 
        dO3/dT at each EMEP station, units of ppbv K-1, [stations,]
    lat_emep : list 
        Latitude coorindate at each EMEP station, units of degrees north, 
        [stations,]
    lng_emep : list 
        Longitude coorindate at each EMEP station, units of degrees east, 
        [stations,]                    
    lng_ml : numpy.ndarray
        The longitude coordinates corresponding to the jet latitude, 
        units of degrees east, [lng]
    lat_jet_ml : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in the Northern Hemisphere mid-latitudes, units of degrees north, 
        [time, lng]
    fstr : str
        Output filename suffix (should specify season, years)  

    Returns
    -------
    None
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    # Find quantities (model, model error, and observations) in regions and 
    # calculate land-based zonally-averaged quantities
    # Western North America
    do3dt2m_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(
        do3dt2m, lat_gmi, lng_gmi, 235., 260., 25., 70.) 
    r_t2mo3_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(
        r_t2mo3, lat_gmi, lng_gmi, 235., 260., 25., 70.) 
    wna_left = geo_idx(235., lng_ml)
    wna_right = geo_idx(260., lng_ml)
    lat_jet_wna = lat_jet_ml[:, wna_left:wna_right+1]
    lat_jet_err_wna = np.nanstd(np.nanmean(lat_jet_wna, axis = 1))
    lat_jet_mean_wna = np.nanmean(lat_jet_wna)
    land_wna = globalo3_calculate.find_grid_overland(lat_wna, lng_wna)
    do3dt2m_err_wna = np.nanstd(do3dt2m_wna*land_wna, axis=1)
    do3dt2m_wna = np.nanmean(do3dt2m_wna*land_wna, axis=1)
    r_t2mo3_err_wna = np.nanstd(r_t2mo3_wna*land_wna, axis=1)
    r_t2mo3_wna = np.nanmean(r_t2mo3_wna*land_wna, axis=1)
    aqs_wna = np.where((np.array(lat_aqs) < 70.) &
                       (np.array(lat_aqs) > 25.) &
                       (np.array(lng_aqs) > 235.) &
                       (np.array(lng_aqs) < 260.))[0]
    do3dt2m_aqs_wna = np.array(do3dt2m_aqs)[aqs_wna]
    r_aqs_wna = np.array(r_aqs)[aqs_wna]
    lat_aqs_wna = np.array(lat_aqs)[aqs_wna]
    naps_wna = np.where((np.array(lat_naps) < 70.) &
                        (np.array(lat_naps) > 25.) &
                        (np.array(lng_naps) > 235.) &
                        (np.array(lng_naps) < 260.))[0]
    do3dt2m_naps_wna = np.array(do3dt2m_naps)[naps_wna]
    r_naps_wna = np.array(r_naps)[naps_wna]
    lat_naps_wna = np.array(lat_naps)[naps_wna]
    # Eastern North America
    do3dt2m_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
        do3dt2m, lat_gmi, lng_gmi, 260., 295., 25., 70.)
    r_t2mo3_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
        r_t2mo3, lat_gmi, lng_gmi, 260., 295., 25., 70.) 
    ena_left = geo_idx(260., lng_ml)
    ena_right = geo_idx(295., lng_ml)
    lat_jet_ena = lat_jet_ml[:, ena_left:ena_right+1]
    lat_jet_err_ena = np.nanstd(np.nanmean(lat_jet_ena, axis = 1))
    lat_jet_mean_ena = np.nanmean(lat_jet_ena)
    land_ena = globalo3_calculate.find_grid_overland(lat_ena, lng_ena)
    do3dt2m_err_ena = np.nanstd(do3dt2m_ena*land_ena, axis=1)
    do3dt2m_ena = np.nanmean(do3dt2m_ena*land_ena, axis=1)
    r_t2mo3_err_ena = np.nanstd(r_t2mo3_ena*land_ena, axis=1)
    r_t2mo3_ena = np.nanmean(r_t2mo3_ena*land_ena, axis=1)
    aqs_ena = np.where((np.array(lat_aqs) < 70.) &
                       (np.array(lat_aqs) > 25.) &
                       (np.array(lng_aqs) > 260.) &
                       (np.array(lng_aqs) < 295.))[0]
    do3dt2m_aqs_ena = np.array(do3dt2m_aqs)[aqs_ena]
    r_aqs_ena = np.array(r_aqs)[aqs_ena]
    lat_aqs_ena = np.array(lat_aqs)[aqs_ena]
    naps_ena = np.where((np.array(lat_naps) < 70.) &
                        (np.array(lat_naps) > 25.) &
                        (np.array(lng_naps) > 260.) &
                        (np.array(lng_naps) < 295.))[0]
    do3dt2m_naps_ena = np.array(do3dt2m_naps)[naps_ena]
    r_naps_ena = np.array(r_naps)[naps_ena]
    lat_naps_ena = np.array(lat_naps)[naps_ena]
    # Europe (since Europe crosses the prime meridian, finding grid cells over
    # 0 deg longitude will not work, so find European domain in two steps)
    do3dt2m_eu1, lat_eu1, lng_eu1 = globalo3_calculate.find_grid_in_bb(
        do3dt2m, lat_gmi, lng_gmi, 350., 360., 25., 70.) 
    r_t2mo3_eu1, lat_eu1, lng_eu1 = globalo3_calculate.find_grid_in_bb(
        r_t2mo3, lat_gmi, lng_gmi, 350., 360., 25., 70.) 
    eu1_left = geo_idx(350., lng_ml)
    eu1_right = geo_idx(360., lng_ml)
    lat_jet_eu1 = lat_jet_ml[:, eu1_left:eu1_right+1]
    land_eu1 = globalo3_calculate.find_grid_overland(lat_eu1, lng_eu1)
    do3dt2m_eu2, lat_eu2, lng_eu2 = globalo3_calculate.find_grid_in_bb(
        do3dt2m, lat_gmi, lng_gmi, 0., 30., 25., 70.) 
    r_t2mo3_eu2, lat_eu2, lng_eu2 = globalo3_calculate.find_grid_in_bb(
        r_t2mo3, lat_gmi, lng_gmi, 0., 30., 25., 70.) 
    eu2_left = geo_idx(0., lng_ml)
    eu2_right = geo_idx(30., lng_ml)
    lat_jet_eu2 = lat_jet_ml[:, eu2_left:eu2_right+1]
    land_eu2 = globalo3_calculate.find_grid_overland(lat_eu2, lng_eu2)
    # Concatenate
    do3dt2m_eu = np.hstack([do3dt2m_eu1, do3dt2m_eu2])
    r_t2mo3_eu = np.hstack([r_t2mo3_eu1, r_t2mo3_eu2])
    lat_jet_eu = np.hstack([lat_jet_eu1, lat_jet_eu2])
    lat_jet_err_eu = np.nanstd(np.nanmean(lat_jet_eu, axis = 1))
    lat_jet_mean_eu = np.nanmean(lat_jet_eu)
    land_eu = np.hstack([land_eu1, land_eu2])
    lat_eu = lat_eu1
    do3dt2m_err_eu = np.nanstd(do3dt2m_eu*land_eu, axis=1)
    do3dt2m_eu = np.nanmean(do3dt2m_eu*land_eu, axis=1)
    r_t2mo3_err_eu = np.nanstd(r_t2mo3_eu*land_eu, axis=1)
    r_t2mo3_eu = np.nanmean(r_t2mo3_eu*land_eu, axis=1)
    emep_eu1 = np.where((np.array(lat_emep) < 70.) &
                        (np.array(lat_emep) > 25.) &
                        (np.array(lng_emep) > 350.) &
                        (np.array(lng_emep) < 360.))[0]
    emep_eu2 = np.where((np.array(lat_emep) < 70.) &
                        (np.array(lat_emep) > 25.) &
                        (np.array(lng_emep) > 0.) &
                        (np.array(lng_emep) < 30.))[0]
    emep_eu = np.hstack([emep_eu1, emep_eu2])
    do3dt2m_emep_eu = np.array(do3dt2m_emep)[emep_eu]
    r_emep_eu = np.array(r_emep)[emep_eu]
    lat_emep_eu = np.array(lat_emep)[emep_eu]
    # East Asia 
    do3dt2m_asia, lat_asia, lng_asia = globalo3_calculate.find_grid_in_bb(
        do3dt2m, lat_gmi, lng_gmi, 90., 125., 25., 70.) 
    r_t2mo3_asia, lat_asia, lng_asia = globalo3_calculate.find_grid_in_bb(
        r_t2mo3, lat_gmi, lng_gmi, 90., 125., 25., 70.) 
    asia_left = geo_idx(90., lng_ml)
    asia_right = geo_idx(125., lng_ml)
    lat_jet_asia = lat_jet_ml[:, asia_left:asia_right+1]
    lat_jet_err_asia = np.nanstd(np.nanmean(lat_jet_asia, axis = 1))
    lat_jet_mean_asia = np.nanmean(lat_jet_asia)
    land_asia = globalo3_calculate.find_grid_overland(lat_asia, lng_asia)
    do3dt2m_err_asia = np.nanstd(do3dt2m_asia*land_asia, axis=1)
    do3dt2m_asia = np.nanmean(do3dt2m_asia*land_asia, axis=1)
    r_t2mo3_err_asia = np.nanstd(r_t2mo3_asia*land_asia, axis=1)
    r_t2mo3_asia = np.nanmean(r_t2mo3_asia*land_asia, axis=1)
    # Plotting
    fig = plt.figure(figsize=(6.5, 6.))
    ax1 = plt.subplot2grid((4,2), (0,0), colspan = 2)
    ax2 = plt.subplot2grid((4,2), (1,0), colspan = 2)
    ax3 = plt.subplot2grid((4,2), (2,0), colspan = 2)
    ax4 = plt.subplot2grid((4,2), (3,0), colspan = 2)
    color_r = '#EF9802'
    color_do3dt2m = '#A2C3E0'
    color_jet = 'k' #'#3F79B7'
    # Western North America r(T, O3)
    ax1.text(0.98, 0.06, 'Western North America', ha = 'right', 
             transform=ax1.transAxes, fontsize=12)
    p1, = ax1.plot(lat_wna, r_t2mo3_wna, ls = '-', lw = 2, color = color_r, 
             zorder = 5)
    p2 = ax1.fill_between(lat_wna, r_t2mo3_wna-r_t2mo3_err_wna,
                     r_t2mo3_wna+r_t2mo3_err_wna, color = color_r, alpha = 0.2,
                     zorder = 5)
    ax1.plot(lat_aqs_wna, r_aqs_wna, 'o', color = color_r, markersize=1, 
             zorder = 2)
    p3, = ax1.plot(lat_naps_wna, r_naps_wna, 'o', color = color_r, markersize=1, 
             zorder = 2)
    # Western North America dO3/dT
    ax1b = ax1.twinx()
    p4, = ax1b.plot(lat_wna, do3dt2m_wna, ls = '-', lw = 2, color = color_do3dt2m, zorder = 5)
    p5 = ax1b.fill_between(lat_wna, do3dt2m_wna-do3dt2m_err_wna,
                     do3dt2m_wna+do3dt2m_err_wna, color = color_do3dt2m, alpha = 0.2, zorder = 5)
    p6, = ax1b.plot(lat_aqs_wna, do3dt2m_aqs_wna, 'o', color = color_do3dt2m, 
              markersize=1, zorder = 2)
    ax1b.plot(lat_naps_wna, do3dt2m_naps_wna, 'o', color = color_do3dt2m, 
              markersize=1, zorder = 2)
    # Western North America eddy-driven jet
    p7 = ax1.errorbar(lat_jet_mean_wna, 0., xerr=[lat_jet_err_wna], fmt='o', 
                 color=color_jet, ecolor=color_jet, elinewidth=2, capsize=3, 
                 zorder = 20)#, label = 'Eddy-driven jet')
    ax1.set_zorder(ax1b.get_zorder()+1)
    ax1.patch.set_visible(False)
    # Add legend
    ax1.legend([p3, (p1, p2), p6, (p4, p5), p7], 
        ['Observed $r\:$(T, O$_{\mathregular{3}}$)', 
         'Modeled $r\:$(T, O$_{\mathregular{3}}$) $\mathregular{\pm}$ 1$\mathregular{\sigma}$', 
         'Observed dO$_{\mathregular{3}}$/dT', 
         'Modeled dO$_{\mathregular{3}}$/dT $\mathregular{\pm}$ 1$\mathregular{\sigma}$',
         'Eddy-driven jet'], loc = 2, ncol = 2, bbox_to_anchor=(0.05, 1.75), 
         frameon=False)
    # Eastern North America r(T, O3)
    ax2.text(0.98, 0.06, 'Eastern North America', ha = 'right', 
             transform=ax2.transAxes, fontsize=12)
    ax2.plot(lat_ena, r_t2mo3_ena, ls = '-', lw = 2, color = color_r, zorder = 5)
    ax2.fill_between(lat_ena, r_t2mo3_ena-r_t2mo3_err_ena,
                     r_t2mo3_ena+r_t2mo3_err_ena, color = color_r, alpha = 0.2, zorder = 5)
    ax2.plot(lat_aqs_ena, r_aqs_ena, 'o', color = color_r, markersize=1, zorder = 2)
    ax2.plot(lat_naps_ena, r_naps_ena, 'o', color = color_r, markersize=1, zorder = 2)
    # Eastern North America dO3/dT
    ax2b = ax2.twinx()
    ax2b.plot(lat_ena, do3dt2m_ena, ls = '-', lw = 2, color = color_do3dt2m, 
              zorder = 15)
    ax2b.fill_between(lat_ena, do3dt2m_ena-do3dt2m_err_ena,
                     do3dt2m_ena+do3dt2m_err_ena, color = color_do3dt2m, 
                     alpha = 0.2, zorder = 5)
    ax2b.plot(lat_aqs_ena, do3dt2m_aqs_ena, 'o', color = color_do3dt2m, 
              markersize=1, zorder = 2)
    ax2b.plot(lat_naps_ena, do3dt2m_naps_ena, 'o', color = color_do3dt2m, 
              markersize=1, zorder = 2)
    # Eastern North America eddy-driven jet
    ax2.errorbar(lat_jet_mean_ena, 0., xerr=[lat_jet_err_ena], fmt='o', 
                 color=color_jet, ecolor=color_jet, elinewidth=2, capsize=3, 
                 zorder = 20)
    ax2.set_zorder(ax2b.get_zorder()+1)
    ax2.patch.set_visible(False)
    # Europe r(T, O3)
    ax3.text(0.98, 0.06, 'Europe', ha = 'right', 
             transform=ax3.transAxes, fontsize=12)
    ax3.plot(lat_eu, r_t2mo3_eu, ls = '-', lw = 2, color = color_r, zorder = 5)
    ax3.fill_between(lat_eu, r_t2mo3_eu-r_t2mo3_err_eu,
                     r_t2mo3_eu+r_t2mo3_err_eu, color = color_r, alpha = 0.2, 
                     zorder = 5)
    ax3.plot(lat_emep_eu, r_emep_eu, 'o', color = color_r, markersize=1, zorder = 2)
    # Europe dO3/dT
    ax3b = ax3.twinx()
    ax3b.plot(lat_eu, do3dt2m_eu, ls = '-', lw = 2, color = color_do3dt2m, 
              zorder = 15)
    ax3b.fill_between(lat_eu, do3dt2m_eu-do3dt2m_err_eu,
                     do3dt2m_eu+do3dt2m_err_eu, color = color_do3dt2m, 
                     alpha = 0.2, zorder = 5)
    ax3b.plot(lat_emep_eu, do3dt2m_emep_eu, 'o', color = color_do3dt2m, 
              markersize=1, zorder = 2)
    # Europe eddy-driven jet
    ax3.errorbar(lat_jet_mean_eu, 0., xerr=[lat_jet_err_eu], fmt='o', 
                 color=color_jet, ecolor=color_jet, elinewidth=2, capsize=3, 
                 zorder = 20)
    ax3.set_zorder(ax3b.get_zorder()+1)
    ax3.patch.set_visible(False)
    # East Asia r(T, O3)
    ax4.text(0.98, 0.06, 'East Asia', ha = 'right', 
             transform=ax4.transAxes, fontsize=12)
    ax4.plot(lat_asia, r_t2mo3_asia, ls = '-', lw = 2, color = color_r, zorder = 5)
    ax4.fill_between(lat_asia, r_t2mo3_asia-r_t2mo3_err_asia,
                     r_t2mo3_asia+r_t2mo3_err_asia, color = color_r, alpha = 0.2, 
                     zorder = 5)
    # East Asia dO3/dT
    ax4b = ax4.twinx()
    ax4b.plot(lat_asia, do3dt2m_asia, ls = '-', lw = 2, color = color_do3dt2m, 
              zorder = 15)
    ax4b.fill_between(lat_asia, do3dt2m_asia-do3dt2m_err_asia,
                     do3dt2m_asia+do3dt2m_err_asia, color = color_do3dt2m, 
                     alpha = 0.2, zorder = 5)
    # East Asia eddy-driven jet
    ax4.errorbar(lat_jet_mean_asia, 0., xerr=[lat_jet_err_asia], fmt='o', 
                 color=color_jet, ecolor=color_jet, elinewidth=2, capsize=2, 
                 zorder = 10)
    for ax in [ax1, ax2, ax3, ax4]: 
        ax.set_xlim([25, 70])
        ax.set_xticks(np.linspace(25, 70, 10))
        ax.set_xticklabels([])
        ax.set_ylim([-1, 1])
    ax4.set_xticklabels([ '25', '', '35', '', '45', '', '55', '', '65', ''])    
    ax4.set_xlabel('Latitude [$^{\mathregular{\circ}}$N]', fontsize=14)
    ax2.set_ylabel('$r\:$(T, O$_{\mathregular{3}}$) [$\cdot$]', fontsize=14, 
                   y = -0.15)
    ax2b.set_ylabel('dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]',
                    fontsize=14, labelpad = 20, y = -0.2, rotation=270)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'zonalavg_byregion_rt2mo3_do3dt2m_%s.png' %fstr, dpi = 300)
    return

import numpy as np
import sys
sys.path.append('/Users/ghkerr/phd/globalo3/')
import globalo3_open, globalo3_calculate, observations_open
# Load data the first iteration only 
try:lat_gmi
except NameError:
    years = [2008, 2009, 2010]
    hours = [14]
    season = 'SON'
    months = [9, 10, 11]
    months_str = ['sep', 'oct', 'nov']
    # Hemisphere definition
    latmin, lngmin, latmax, lngmax = -1., 0., 90., 360.
    # Load Northern Hemisphere GMI CTM O3
    lat_gmi, lng_gmi, times_gmi, o3_gmi = \
        globalo3_open.open_overpass2_specifieddomain(years, months_str, latmin, 
        latmax, lngmin, lngmax, 'O3', 'HindcastMR2')
    o3_gmi = o3_gmi*1e9        
    # Load Northern Hemisphere 2-meter temperatures
    lat_merra, lng_merra, t2m_merra = \
        globalo3_open.open_merra2t2m_specifieddomain(years, months_str, latmin, 
        latmax, lngmin, lngmax)        
    # Interpolate 2-meter temperature
    t2m_merra = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, lng_gmi, 
        lat_merra, lng_merra, t2m_merra)
    # Load Northern Hemisphere 500 hPa MERRA-2 data
    U500, mtime, lat_merra, lng_merra = \
        globalo3_open.open_merra2_specifieddomain(years, months_str,
        [0,3,6,9,12,15,18,21], 'U', 'inst3_3d_asm_Np_500hPa', lngmin, latmax, 
        lngmax, latmin, dailyavg='yes')
    # Interpolate 500 hPa u-wind
    U500 = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
        lng_gmi, lat_merra, lng_merra, U500, checkplot='yes')
    # Load observational datasets
    naps = observations_open.open_napso3(years, months, hours)
    aqs = observations_open.open_aqso3(years, months, hours)
    emep = observations_open.open_emepo3(years, months, hours)
    # Calculate dO3/dT and r(T, O3) from model and observational datasets
    do3dt2m = globalo3_calculate.calculate_do3dt(t2m_merra, o3_gmi, lat_gmi, 
        lng_gmi)
    r_t2mo3 = globalo3_calculate.calculate_r(t2m_merra, o3_gmi, lat_gmi, lng_gmi)
    r_naps, do3dt2m_naps, lat_naps, lng_naps = \
        globalo3_calculate.calculate_obs_do3dt_rto3(naps, t2m_merra, times_gmi, 
        lat_gmi, lng_gmi)
    r_aqs, do3dt2m_aqs, lat_aqs, lng_aqs = \
        globalo3_calculate.calculate_obs_do3dt_rto3(aqs, t2m_merra, times_gmi, 
        lat_gmi, lng_gmi)    
    r_emep, do3dt2m_emep, lat_emep, lng_emep = \
        globalo3_calculate.calculate_obs_do3dt_rto3(emep, t2m_merra, times_gmi, 
        lat_gmi, lng_gmi)    
    # Find bias 
    r_naps_bias = globalo3_calculate.ctm_obs_bias(lat_naps, lng_naps, r_naps, 
        lat_gmi, lng_gmi, r_t2mo3)
    do3dt2m_naps_bias = globalo3_calculate.ctm_obs_bias(lat_naps, lng_naps, 
        do3dt2m_naps, lat_gmi, lng_gmi, do3dt2m)
    r_aqs_bias = globalo3_calculate.ctm_obs_bias(lat_aqs, lng_aqs, r_aqs, 
        lat_gmi, lng_gmi, r_t2mo3)
    do3dt2m_aqs_bias = globalo3_calculate.ctm_obs_bias(lat_aqs, lng_aqs, 
        do3dt2m_aqs, lat_gmi, lng_gmi, do3dt2m)
    r_emep_bias = globalo3_calculate.ctm_obs_bias(lat_emep, lng_emep, r_emep, 
        lat_gmi, lng_gmi, r_t2mo3)
    do3dt2m_emep_bias = globalo3_calculate.ctm_obs_bias(lat_emep, lng_emep, 
        do3dt2m_emep, lat_gmi, lng_gmi, do3dt2m)
    # Subset fields in mid-latitudes
    U500_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(U500, lat_gmi, 
        lng_gmi, 0., 360., 20., 70.)
    o3_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(o3_gmi, lat_gmi, 
        lng_gmi, 0., 360., 20., 70.)
    t2m_ml, lat_ml, lng_nhml = globalo3_calculate.find_grid_in_bb(t2m_merra, 
        lat_gmi, lng_gmi, 0., 360., 20., 70.)
    do3dt2m_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(do3dt2m, 
        lat_gmi, lng_gmi, 0., 360., 20., 70.)
    r_t2mo3_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(r_t2mo3, 
        lat_gmi, lng_gmi, 0., 360., 20., 70.)
    # Find latitude of eddy-driven jet and fields at jet
    lat_jet_ml, o3_jet_ml = globalo3_calculate.find_field_atjet(o3_ml, U500_ml, 
        lat_ml, lng_ml, 20, anom = True)    
    lat_jet_ml, t2m_jet_ml = globalo3_calculate.find_field_atjet(t2m_ml, 
        U500_ml, lat_ml, lng_ml, 20, anom=True)
    lat_jet_ml, do3dt2m_jet_ml = globalo3_calculate.find_field_atjet(
        do3dt2m_ml, U500_ml, lat_ml, lng_ml, 20)
    lat_jet_ml, r_t2mo3_jet_ml = globalo3_calculate.find_field_atjet(
        r_t2mo3_ml, U500_ml, lat_ml, lng_ml, 20) 
    # Slope and correlation of O3/jet distance and 2-meter temperature and 
    # jet distance
    m_o3jetdist, r_o3jetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(o3_gmi, lat_gmi, 
        lng_gmi, lat_jet_ml, lng_ml)
    m_t2mjetdist, r_t2mjetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(t2m_merra, lat_gmi, 
        lng_gmi, lat_jet_ml, lng_ml)    

maparea = 'nh'
# Mean O3, mean eddy-driven jet position and variability
map_hemisphere(lat_gmi, 
    lng_gmi, 
    np.mean(o3_gmi, axis=0), 
    '%s O$_{\mathregular{3}}$' %season, 
    '[ppbv]', 
    np.linspace(25, 65, 11), 
    'PuBu', 
    maparea,
    'meano3_jet_%s_%d-%d'%(season, years[0],years[-1]), 
    e_n=np.nanmean(lat_jet_ml,axis=0), 
    eerr_n=np.std(lat_jet_ml,axis=0),
    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
            lat_gmi.min()+1, lat_gmi.max()-5])
 dO3/dT
map_hemisphere(lat_gmi, 
    lng_gmi,
    do3dt2m, 
    '%s dO$_{\mathregular{3}}$/dT' %season,
    '[ppbv K$^{\mathregular{-1}}$]', 
    np.linspace(-2, 2., 9), 
    'bwr', 
    maparea,
    'do3dt2m_jet_%s_%d-%d'%(season, years[0],years[-1]), 
    e_n=np.nanmean(lat_jet_ml,axis=0), 
    eerr_n=np.std(lat_jet_ml,axis=0),
    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
            lat_gmi.min()+1, lat_gmi.max()-5])
# r(T, O3)
map_hemisphere(lat_gmi, 
    lng_gmi, 
    r_t2mo3, 
    r'%s $\it{r}\:$(T, O$_\mathregular{3}$)' %season,  
    '[$\cdot$]', 
    np.linspace(-1., 1., 11), 
    'bwr', 
    maparea,
    'rt2mo3_jet_%s_%d-%d'%(season, years[0],years[-1]), 
    e_n=np.nanmean(lat_jet_ml,axis=0), 
    eerr_n=np.std(lat_jet_ml,axis=0),
    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
            lat_gmi.min()+1, lat_gmi.max()-5])
# dO3/d(jet lat - lat)
map_hemisphere(lat_gmi, 
    lng_gmi,
    m_o3jetdist, 
    '%s dO$_{\mathregular{3}}$/d($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$)' %season, 
    '[ppbv degree$^{\mathregular{-1}}$]', 
    np.linspace(-0.3, 0.3, 7), 
    'bwr', 
    maparea,
    'do3djetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
    e_n=np.nanmean(lat_jet_ml,axis=0), 
    eerr_n=np.std(lat_jet_ml,axis=0),
    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
            lat_gmi.min()+1, lat_gmi.max()-5])
# r(jet lat - lat, O3)
map_hemisphere(lat_gmi, 
    lng_gmi,
    r_o3jetdist, 
    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, O$_{\mathregular{3}}$)' %season, 
    '[$\cdot$]', 
    np.linspace(-0.7, 0.7, 8), 
    'bwr', 
    maparea,
    'ro3jetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
    e_n=np.nanmean(lat_jet_ml,axis=0), 
    eerr_n=np.std(lat_jet_ml,axis=0),
    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
            lat_gmi.min()+1, lat_gmi.max()-5])
# dT/d(jet lat - lat)
map_hemisphere(lat_gmi, 
    lng_gmi,
    m_t2mjetdist, 
    '%s dT/d($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$)' %season, 
    '[K degree$^{\mathregular{-1}}$]', 
    np.linspace(-0.1, 0.4, 11), 
    'gist_earth',
    maparea,
    'dt2mdjetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
    e_n=np.nanmean(lat_jet_ml,axis=0), 
    eerr_n=np.std(lat_jet_ml,axis=0),
    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
            lat_gmi.min()+1, lat_gmi.max()-5])
# r(jet lat - lat, O3)
map_hemisphere(lat_gmi, 
    lng_gmi,
    r_t2mjetdist, 
    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, T)' %season, 
    '[$\cdot$]', 
    np.linspace(-0.7, 0.7, 8), 
    'bwr',
    maparea,
    'rt2mjetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
    e_n=np.nanmean(lat_jet_ml,axis=0), 
    eerr_n=np.std(lat_jet_ml,axis=0),
    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
            lat_gmi.min()+1, lat_gmi.max()-5])
# O3 anomaly in vicinity of jet 
fieldatjet(lat_ml, 
    lng_ml, 
    lat_jet_ml, 
    np.nanmean(o3_jet_ml, axis=0),
    'bwr', 
    r'$\mathregular{\delta}$O$_{\mathregular{3}}$ [ppbv]',  
    np.linspace(-5., 5., 6), 
    '%s_o3anom' %maparea, 
    '%s_%d-%d'%(season, years[0],years[-1]), 
    skiplng=6)
# 2-meter temperature anomaly in vicinity of jet 
fieldatjet(lat_ml, 
    lng_ml,
    lat_jet_ml, 
    np.nanmean(t2m_jet_ml, axis=0),
    'bwr', 
    r'$\mathregular{\delta}$T [K]',  
    np.linspace(-3., 3., 7), 
    '%s_t2manom' %maparea, 
    '%s_%d-%d'%(season, years[0],years[-1]), 
    skiplng=6)
# dO3/dT in vicinity of jet
fieldatjet(lat_ml, 
    lng_ml,
    lat_jet_ml, 
    do3dt2m_jet_ml,
    'Reds', 
    'dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]',  
    np.linspace(0, 2, 6),
    '%s_do3dt' %maparea, 
    '%s_%d-%d'%(season, years[0],years[-1]), 
    skiplng=6)
# r(T, O3) in vicinity of jet
fieldatjet(lat_ml, 
    lng_ml,
    lat_jet_ml, 
    r_t2mo3_jet_ml,
    'bwr', 
    'dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]',  
    np.linspace(-1, 1, 11),
    '%s_rt2mo3' %maparea, 
    '%s_%d-%d'%(season, years[0],years[-1]), 
    skiplng=6)

"""MODEL-OBSERVATION COMPARISON"""
# North American r(T, O3)
map_extent([-170, -50, 24, 60], 
    '%s $r\:$(T, O$_{\mathregular{3}}$)' %season, 
    '[$\mathregular{\cdot}$]', 
    'bwr', 
    np.linspace(-1., 1, 9),
    'northamerica',
    'rt2mo3_gmi_aqsnaps_%s_%d-%d' %(season, years[0], years[-1]),
    lat_ctm=lat_gmi, 
    lng_ctm=lng_gmi, 
    field_ctm=r_t2mo3, 
    lat_obs=np.hstack((lat_naps, lat_aqs)), 
    lng_obs=np.hstack((lng_naps, lng_aqs)), 
    field_obs=np.hstack((r_naps, r_aqs)))
# North American r(T, O3) bias
map_extent([-170, -50, 24, 60], 
    '%s $r\:$(T, O$_{\mathregular{3}}$)' %season, 
    '[$\mathregular{\cdot}$]', 
    'bwr', 
    np.linspace(-.5, .5, 6),
    'northamerica',
    'rt2mo3bias_gmi_aqsnaps_%s_%d-%d' %(season, years[0], years[-1]),
    lat_obs=np.hstack((lat_naps, lat_aqs)), 
    lng_obs=np.hstack((lng_naps, lng_aqs)), 
    field_obs=np.hstack((r_naps_bias, r_aqs_bias)), 
    extend='both')
# North American dO3/dT
map_extent([-170, -50, 24, 60], 
    '%s dO$_{\mathregular{3}}$/dT' %season,
    '[ppbv K$^{\mathregular{-1}}$]',
    'Reds', 
    np.linspace(0., 2.5, 6),
    'northamerica',
    'do3dt2m_gmi_aqsnaps_%s_%d-%d' %(season, years[0], years[-1]),
    lat_ctm=lat_gmi, 
    lng_ctm=lng_gmi, 
    field_ctm=do3dt2m, 
    lat_obs=np.hstack((lat_naps, lat_aqs)), 
    lng_obs=np.hstack((lng_naps, lng_aqs)), 
    field_obs=np.hstack((do3dt2m_naps, do3dt2m_aqs)), 
    extend='both')
# North American dO3/dT bias
map_extent([-170, -50, 24, 60], 
    '%s dO$_{\mathregular{3}}$/dT' %season,
    '[ppbv K$^{\mathregular{-1}}$]',
    'bwr', 
    np.linspace(-1., 1., 6),
    'northamerica',
    'do3dt2mbias_gmi_aqsnaps_%s_%d-%d' %(season, years[0], years[-1]),
    lat_obs=np.hstack((lat_naps, lat_aqs)), 
    lng_obs=np.hstack((lng_naps, lng_aqs)), 
    field_obs=np.hstack((do3dt2m_naps_bias, do3dt2m_aqs_bias)), 
    extend='both')
# European r(T, O3)
map_extent([-14, 37, 28, 70], 
    '%s $r\:$(T, O$_{\mathregular{3}}$)' %season, 
    '[$\mathregular{\cdot}$]', 
    'bwr', 
    np.linspace(-1., 1., 9),
    'europe',
    'rt2mo3_gmi_emep_%s_%d-%d' %(season, years[0], years[-1]),
    lat_ctm=lat_gmi, 
    lng_ctm=lng_gmi, 
    field_ctm=r_t2mo3,
    lat_obs=lat_emep,
    lng_obs=lng_emep, 
    field_obs=r_emep)
# European r(T, O3) bias
map_extent([-14, 37, 28, 70], 
    '%s $r\:$(T, O$_{\mathregular{3}}$)' %season, 
    '[$\mathregular{\cdot}$]', 
    'bwr', 
    np.linspace(-.5, .5, 6),
    'europe',
    'rt2mo3bias_gmi_emep_%s_%d-%d' %(season, years[0], years[-1]),
    lat_obs=lat_emep,
    lng_obs=lng_emep, 
    field_obs=r_emep_bias,
    extend='both')
# European dO3/dT
map_extent([-14, 37, 28, 70], 
    '%s dO$_{\mathregular{3}}$/dT' %season,
    '[ppbv K$^{\mathregular{-1}}$]',
    'Reds', 
    np.linspace(0., 2.5, 6),
    'europe',
    'do3dt2m_gmi_emep_%s_%d-%d' %(season, years[0], years[-1]),
    lat_ctm=lat_gmi, 
    lng_ctm=lng_gmi, 
    field_ctm=do3dt2m, 
    lat_obs=lat_emep,
    lng_obs=lng_emep, 
    field_obs=do3dt2m_emep,
    extend='both')
# European dO3/dT bias
map_extent([-14, 37, 28, 70], 
    '%s dO$_{\mathregular{3}}$/dT' %season,
    '[ppbv K$^{\mathregular{-1}}$]',
    'bwr', 
    np.linspace(-1., 1., 6),
    'europe',
    'do3dt2mbias_gmi_emep_%s_%d-%d' %(season, years[0], years[-1]),
    lat_obs=lat_emep,
    lng_obs=lng_emep, 
    field_obs=do3dt2m_emep_bias,
    extend='both')
# Zonally-averaged model/observations comparison 
zonalavg_byregion_rt2mo3_do3dt2m(r_t2mo3, 
    do3dt2m, 
    lat_gmi, 
    lng_gmi, 
    r_aqs, 
    do3dt2m_aqs, 
    lat_aqs, 
    lng_aqs, 
    r_naps, 
    do3dt2m_naps, 
    lat_naps, 
    lng_naps, 
    r_emep, 
    do3dt2m_emep, 
    lat_emep, 
    lng_emep, 
    lng_ml, 
    lat_jet_ml, 
    '%s_%d-%d' %(season, years[0], years[-1]))



























"""O3-CYCLONE RELATIONSHIP FROM GISS MERRA-2 CYCLONE DATABASE""" 
#cyclones = globalo3_open.open_merra2_cyclones('06/01/2008', '08/31/2008')
#lat = lat_gmi_n
#lng = lng_gmi_n
#times = times_n_jja
#o3 = o3_n_jja
#import pandas as pd
#from datetime import date
#from datetime import datetime
#o3anomalies = []
## Select cyclones over North America with at least 25% land cover    
#cyclones_na = cyclones.loc[(cyclones['Latitude'] > 20.) &
#                     (cyclones['Latitude'] < 70.) & 
#                     (cyclones['Longitude'] > 205.) &
#                     (cyclones['Longitude'] < 304.) &
#                     (cyclones['Land Cover'] > 0.25)]
## Loop through cyclones, extract dates, latitude, and longitude in 
## order to find corresponding indices in the O3 dataset. As there are 
## multiple entries corresponding to each cyclone ID, values are averaged 
## over the lifespan of a cyclone
#for id in np.unique(cyclones_na['Storm ID']):
#    id = cyclones_na.loc[cyclones_na['Storm ID'] == id]
#    bye = []
#    for idx in np.arange(0, id.shape[0], 1):
#        entry = id.iloc[idx]
#        clat = entry['Latitude']
#        clng = entry['Longitude']
#        # Find mean time and round to nearest day 
#        ctime = entry['Date']
#        ctime = (np.array(ctime, dtype='datetime64[s]').view('i8').
#                 mean().astype('datetime64[s]'))
#        ctime = pd.to_datetime(ctime).round('1D')
#        ctime = datetime.date(ctime.to_pydatetime())
#        # Time, latitude, longitude indices corresponding to cyclone in 
#        # CTM dataset
#        t = np.where(times==ctime)[0][0]
#        x = np.abs(lng-clng).argmin()
#        y = np.abs(lat-clat).argmin()
#        lat_jet_today = lat_jet_ml[t, x]
#        if lat_jet_today < y:
#            # O3 anomaly within +/- 10 degrees of cyclone center
#            o3anom_atcyclone = (o3[t,y-10:y+11,x-10:x+11]-
#                                np.nanmean(o3[:,y-10:y+11,x-10:x+11]))
#            o3anomalies.append(o3anom_atcyclone)
#fig = plt.figure()
#ax = plt.subplot2grid((1,1),(0,0))
#ax.set_title('O$_\mathregular{3}$ Anomaly for all North American Cyclones, N=%d' %len(hi))
#mp = ax.contourf(np.nanmean(np.array(o3anomalies),axis=0), np.linspace(-6, 6, 11), 
#                 cmap=plt.get_cmap('bwr'), extend='both')
#plt.colorbar(mp, label = '[ppbv]')
#plt.savefig('/Users/ghkerr/Desktop/o3_allcyclones_below.eps', dpi=300)
#fig = plt.figure(figsize=(8, 10))
#ax = plt.subplot2grid((1,1),(0,0))
## Select cyclones over North America with at least 25% land cover    
#cyclones_na = cyclones.loc[(cyclones['Latitude'] > 20.) &
#                     (cyclones['Latitude'] < 70.) & 
#                     (cyclones['Longitude'] > 205.) &
#                     (cyclones['Longitude'] < 304.) &
#                     (cyclones['Land Cover'] > 0.25)]
## Loop through cyclones, extract dates, latitude, and longitude in 
## order to find corresponding indices in the O3 dataset. As there are 
## multiple entries corresponding to each cyclone ID, values are averaged 
## over the lifespan of a cyclone
#for id in np.unique(cyclones_na['Storm ID']):
#    id = cyclones_na.loc[cyclones_na['Storm ID'] == id]
#    displacement = []
#    o3_atdisplacement = []
#    for idx in np.arange(0, id.shape[0], 1):
#        entry = id.iloc[idx]
#        clat = entry['Latitude']
#        clng = entry['Longitude']
#        # Find mean time and round to nearest day 
#        ctime = entry['Date']
#        ctime = (np.array(ctime, dtype='datetime64[s]').view('i8').
#                 mean().astype('datetime64[s]'))
#        ctime = pd.to_datetime(ctime).round('1D')
#        ctime = datetime.date(ctime.to_pydatetime())
#        # Time, latitude, longitude indices corresponding to cyclone in 
#        # CTM dataset
#        t = np.where(times==ctime)[0][0]
#        x = np.abs(lng-clng).argmin()
#        y = np.abs(lat-clat).argmin()
#        lat_jet_now = lat_jet_ml[t, x]
#        displacement.append(clat-lat_jet_now)
#        o3_atdisplacement.append(o3[t,y,x]-np.nanmean(o3[:,y,x]))
#        ax.plot(np.arange(0, len(displacement), 1), displacement, '-ko', 
#                lw=0.5, markersize=5, zorder=1) 
#        p = ax.scatter(np.arange(0, len(displacement), 1), displacement, 
#                       c=o3_atdisplacement, s=5, cmap='bwr', vmin=-6, vmax=6, 
#                       zorder=2)        
#ax.set_xlim([0, 25])
#ax.set_xticks([0, 5, 10, 15, 20, 25])
#ax.set_xticklabels([0, 30, 60, 90, 120, 150])
#ax.set_xlabel('Hours since genesis', fontsize=16)
#ax.set_ylabel('$\mathregular{\phi_{\:cyclone}}$ - $\mathregular{\phi_{\:jet}}$', fontsize=16)
#cb = fig.colorbar(p, extend='both')
#cb.set_label('$\mathregular{\delta}$O$_{\mathregular{3}}$ [ppbv]', fontsize=16)
#plt.savefig('/Users/ghkerr/Desktop/cyclone_jet_diff.eps', dpi=300)


"""PLOT MEAN LOCATION OF JET STREAM (MAX U WINDS AT 500 HPA) AND 
   O3, dO3/dT, and r(O3, T) WITHIN +/- 10 DEGREES OF JET STREAM""" 
#season = 'JJA' 
#o3 = o3_n_jja
#t2m = t2m_n_jja
#r_t2mo3 = r_t2mo3_n_jja
#U500 = U500_n_jja
#do3dt2m = do3dt2m_n_jja
#lat = lat_gmi_n
#lng = lng_gmi_n
#latmin = 20.
#latmax = 70.
#latmin = -70.
#latmax = -20.
#lngmin = 0.
#lngmax = 360.
#maparea = 'nh'
#"""PLOT MEAN LOCATION OF JET STREAM (MAX U WINDS AT 500 HPA) AND 
#   O3, dO3/dT, and r(O3, T) WITHIN +/- 10 DEGREES OF JET STREAM""" 
# Subset fields in mid-latitudes
#U500_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(U500, lat, lng, 
#    lngmin, lngmax, latmin, latmax)
#o3_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(o3, lat, lng, 
#    lngmin, lngmax, latmin, latmax)
#t2m_ml, lat_ml, lng_nhml = globalo3_calculate.find_grid_in_bb(t2m, lat, lng, 
#    lngmin, lngmax, latmin, latmax)
#do3dt2m_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(do3dt2m, lat, 
#    lng, lngmin, lngmax, latmin, latmax)
#r_t2mo3_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(r_t2mo3, lat, 
#    lng, lngmin, lngmax, latmin, latmax)
#land_nhml = globalo3_calculate.find_grid_overland(lat_ml, lng_ml)   
#lat_jet_ml, o3_jet_ml = globalo3_calculate.find_field_atjet(o3_ml, U500_ml, 
#    lat_ml, lng_ml, 20, anom = True)
#lat_jet_ml, t2m_jet_ml = globalo3_calculate.find_field_atjet(t2m_ml, 
#    U500_ml, lat_ml, lng_ml, 20, anom=True)
#lat_jet_ml, do3dt2m_jet_ml = globalo3_calculate.find_field_atjet(do3dt2m_ml, 
#    U500_ml, lat_ml, lng_ml, 20)
#lat_jet_ml, r_t2mo3_jet_ml = globalo3_calculate.find_field_atjet(
#    r_t2mo3_ml, U500_ml, lat_ml, lng_ml, 20) 
## Slope and correlation of O3/jet distance and 2-meter temperature and 
## jet distance
#m_o3jetdist, r_o3jetdist = \
#    globalo3_calculate.calculate_fieldjet_relationship(o3, lat, lng, 
#    lat_jet_ml, lng_ml)
#m_t2mjetdist, r_t2mjetdist = \
#    globalo3_calculate.calculate_fieldjet_relationship(t2m, lat, lng, 
#    lat_jet_ml, lng_ml)
## Mean O3, mean eddy-driven jet position and variability
#map_hemisphere(lat, 
#    lng, 
#    np.mean(o3, axis=0), 
#    '%s O$_{\mathregular{3}}$' %season, 
#    '[ppbv]', 
#    np.linspace(25, 65, 11), 
#    'PuBu', 
#    maparea,
#    'meano3_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    extent=[lng.min()-180., lng.max()-180., lat.min()+1, lat.max()-5])
## dO3/dT
#map_hemisphere(lat, 
#    lng, 
#    do3dt2m, 
#    '%s dO$_{\mathregular{3}}$/dT' %season,
#    '[ppbv K$^{\mathregular{-1}}$]', 
#    np.linspace(0, 3., 7), 
#    'Reds', 
#    maparea,
#    'do3dt2m_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    extent=[lng.min()-180., lng.max()-180., lat.min()+1, lat.max()-5])
## r(T, O3)
#map_hemisphere(lat,
#    lng, 
#    r_t2mo3, 
#    r'%s $\it{r}\:$(T, O$_\mathregular{3}$)' %season,  
#    '[$\cdot$]', 
#    np.linspace(-1., 1., 11), 
#    'bwr', 
#    maparea,
#    'rt2mo3_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    extent=[lng.min()-180., lng.max()-180., lat.min()+1, lat.max()-5])
## dO3/d(jet lat - lat)
#map_hemisphere(lat, 
#    lng, 
#    m_o3jetdist, 
#    '%s dO$_{\mathregular{3}}$/d($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$)' %season, 
#    '[ppbv degree$^{\mathregular{-1}}$]', 
#    np.linspace(0.0, 0.3, 7), 
#    'gist_earth_r', 
#    maparea,
#    'do3djetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    extent=[lng.min()-180., lng.max()-180., lat.min()+1, lat.max()-5])
## r(jet lat - lat, O3)
#map_hemisphere(lat, 
#    lng, 
#    r_o3jetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, O$_{\mathregular{3}}$)' %season, 
#    '[$\cdot$]', 
#    np.linspace(-0.6, 0.6, 13), 
#    'bwr', 
#    maparea,
#    'ro3jetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    extent=[lng.min()-180., lng.max()-180., lat.min()+1, lat.max()-5])
## dT/d(jet lat - lat)
#map_hemisphere(lat, 
#    lng, 
#    m_t2mjetdist, 
#    '%s dT/d($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$)' %season, 
#    '[K degree$^{\mathregular{-1}}$]', 
#    np.linspace(0., 0.3, 7), 
#    'gist_earth_r',
#    maparea,
#    'dt2mdjetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    extent=[lng.min()-180., lng.max()-180., lat.min()+1, lat.max()-5])
## r(jet lat - lat, O3)
#map_hemisphere(lat, 
#    lng, 
#    r_t2mjetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, T)' %season, 
#    '[$\cdot$]', 
#    np.linspace(-0.6, 0.6, 13), 
#    'bwr',
#    maparea,
#    'rt2mjetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    extent=[lng.min()-180., lng.max()-180., lat.min()+1, lat.max()-5])
## O3 anomaly in vicinity of jet 
#fieldatjet(lat_ml, 
#    lng_ml, 
#    lat_jet_ml, 
#    np.nanmean(o3_jet_ml, axis=0),
#    'bwr', 
#    r'$\mathregular{\delta}$O$_{\mathregular{3}}$ [ppbv]',  
#    np.linspace(-5., 5., 6), 
#    '%s_o3anom' %maparea, 
#    '%s_%d-%d'%(season, years[0],years[-1]), 
#    skiplng=6)
## 2-meter temperature anomaly in vicinity of jet 
#fieldatjet(lat_ml, 
#    lng_ml,
#    lat_jet_ml, 
#    np.nanmean(t2m_jet_ml, axis=0),
#    'bwr', 
#    r'$\mathregular{\delta}$T [K]',  
#    np.linspace(-3., 3., 7), 
#    '%s_t2manom' %maparea, 
#    '%s_%d-%d'%(season, years[0],years[-1]), 
#    skiplng=6)
## dO3/dT in vicinity of jet
#fieldatjet(lat_ml, 
#    lng_ml,
#    lat_jet_ml, 
#    do3dt2m_jet_ml,
#    'Reds', 
#    'dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]',  
#    np.linspace(0, 2, 6),
#    '%s_do3dt' %maparea, 
#    '%s_%d-%d'%(season, years[0],years[-1]), 
#    skiplng=6)
## r(T, O3) in vicinity of jet
#fieldatjet(lat_ml, 
#    lng_ml,
#    lat_jet_ml, 
#    r_t2mo3_jet_ml,
#    'bwr', 
#    'dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]',  
#    np.linspace(-1, 1, 11),
#    '%s_rt2mo3' %maparea, 
#    '%s_%d-%d'%(season, years[0],years[-1]), 
#    skiplng=6)

"""PLOT MEAN FIELDS AND CORRELATIONS"""
## Mean EDGAR NOx emissions
#map_nh(lat_gmi_n, lng_gmi_n, np.mean(nox_edgar_n, axis=0), '', 
#       'NO$_{x\:\mathregular{, EDGAR}}$ [kg m$^{\mathregular{-2}}$ '+
#       's$^{\mathregular{-1}}$]', np.linspace(0e-10, 1e-10, 11), 'PuBu', 
#       'meanedgarnox_%d-%d_jet'%(years[0],years[-1]), 
#       e_n=np.nanmean(lat_jet_nhml, axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))

""" O3-JET RELATIONSHIP DIFFERENCES BETWEEN STD/EMFIX SIMULATIONS """
## Load GMI CTM O3 from EmFix simulation
#lat_emfix_n, lng_emfix_n, times_emfix_n, o3_emfix_n = \
#globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_n, 
#latmax_n, lngmin_n, lngmax_n, 'O3', 'Hindcast3Igac2')
#o3_emfix_n = o3_emfix_n*1e9
## Interpolate EmFix to the resolution of HindcastMR2 (not the best practice, 
## change in the future)
#o3_emfix_n = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi_n, 
#    lng_gmi_n, lat_emfix_n, lng_emfix_n, o3_emfix_n)
## Load GMI CTM O3 over North America from Std simulation
#lat_std_n, lng_std_n, times_std_n, o3_std_n = \
#globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_n, 
#latmax_n, lngmin_n, lngmax_n, 'O3', 'HindcastFFIgac2')
#o3_std_n = o3_std_n*1e9
#o3_std_n = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi_n, 
#    lng_gmi_n, lat_std_n, lng_std_n, o3_std_n)
## Subset fields in the Northern Hemisphere mid-latitudes
#U500_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(U500, 
#    lat_gmi_n, lng_gmi_n, 0., 360., 20., 70.)
#o3_emfix_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(
#    o3_emfix_n, lat_gmi_n, lng_gmi_n, 0., 360., 20., 70.)
#o3_std_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(o3_std_n, 
#    lat_gmi_n, lng_gmi_n, 0., 360., 20., 70.)
## Identify eddy-driven jet and EmFix/Std O3 in the vicinity of the jet 
#lat_jet_nhml, o3_emfix_jet_nhml = globalo3_calculate.find_field_atjet(
#    o3_emfix_nhml, U500_nhml, lat_nhml, lng_nhml, 10, anom=True)
#lat_jet_nhml, o3_std_jet_nhml = globalo3_calculate.find_field_atjet(
#    o3_std_nhml, U500_nhml, lat_nhml, lng_nhml, 10, anom=True)
## O3 anomaly about jet from EmFix simulation
#edjetlocation_fieldatedjet(lat_nhml, 
#    lng_nhml,
#    lat_jet_nhml, 
#    np.nanmean(o3_emfix_jet_nhml, axis=0),
#    'bwr', 
#    '$\mathregular{\Delta}$O$_{\mathregular{3,\:EmFix}}$ [ppbv]',  
#    np.linspace(-14, 14, 8),
#    'o3anom_jet',
#    'emfix_%d-%d'%(years[0], years[-1]),
#    skiplng=6)
## O3 anomaly about jet from Std simulation
#edjetlocation_fieldatedjet(lat_nhml, 
#    lng_nhml,
#    lat_jet_nhml, 
#    np.nanmean(o3_std_jet_nhml, axis=0),
#    'bwr', 
#    '$\mathregular{\Delta}$O$_{\mathregular{3,\:Std}}$ [ppbv]',  
#    np.linspace(-14, 14, 8),
#    'o3anom_jet',
#    'std_%d-%d'%(years[0], years[-1]),
#    skiplng=6)
## Difference in O3 anomalies about jet from EmFix and Std simulations
#edjetlocation_fieldatedjet(lat_nhml, 
#    lng_nhml,
#    lat_jet_nhml, 
#    np.nanmean(o3_emfix_jet_nhml, axis=0)-np.nanmean(o3_std_jet_nhml, axis=0),
#    'bwr', 
#    '($\mathregular{\Delta}$O$_{\mathregular{3,\:EmFix}}$)-'+
#    '($\mathregular{\Delta}$O$_{\mathregular{3,\:Std}}$) [ppbv]',
#    np.linspace(-5, 5, 6),
#    'o3anom_jet',
#    'emfix-std_%d-%d'%(years[0], years[-1]),
#    skiplng=6)

"""O3 ANOMALIES AND JET POSITION ON DAYS WITH EQUATOR/POLEWARD JET"""
## Find fields over North American mid-latitudes (20-70˚N, 225-300˚E)
#U500_naml, lat_naml, lng_naml = globalo3_calculate.find_grid_in_bb(U500, 
#    lat_gmi_n, lng_gmi_n, 225., 300., 20., 70.)
#o3_naml, lat_naml, lng_naml = globalo3_calculate.find_grid_in_bb(o3_n, 
#    lat_gmi_n, lng_gmi_n, 225., 300., 20., 70.)
#t2m_naml, lat_naml, lng_naml = globalo3_calculate.find_grid_in_bb(t2m_n, 
#    lat_gmi_n, lng_gmi_n, 225., 300., 20., 70.)
#land_naml = globalo3_calculate.find_grid_overland(lat_naml, lng_naml)
#lat_jet_naml, o3_jet_naml = globalo3_calculate.find_field_atjet(o3_naml, 
#    U500_naml, lat_naml, lng_naml, 20, anom=True)
#lat_jet_naml, t2m_jet_naml = globalo3_calculate.find_field_atjet(t2m_naml, 
#    U500_naml, lat_naml, lng_naml, 20, anom=True)
#o3_anom_za = np.nanmean(o3_jet_naml, axis=2) 
#t2m_anom_za = np.nanmean(t2m_jet_naml, axis=2) 
#lat_jet_naml_za = np.nanmean(lat_jet_naml, axis=1)
## Group days into percentile categories based on the jet latitude
#jet_0_30 = np.where(lat_jet_naml_za < np.percentile(lat_jet_naml_za, 30))[0]
#jet_30_70 = np.where((lat_jet_naml_za >= np.percentile(lat_jet_naml_za, 30)) &
#                     (lat_jet_naml_za < np.percentile(lat_jet_naml_za, 70)))[0]
#jet_70_100 = np.where(lat_jet_naml_za >= np.percentile(lat_jet_naml_za, 70))[0]
## Since the O3 anomaly (o3_anom_za) is calculated about the jet axis, 
## create latitude "spines" where the center value of the spine is the mean 
## jet latitude for each percentile category and values away from the center
## reflect the resolution of the data
#res = np.diff(lat_naml).mean()
#js_0_30c = lat_jet_naml_za[jet_0_30].mean()
#js_0_30 = np.linspace(js_0_30c-20, js_0_30c+20, o3_anom_za.shape[1])
#js_30_70c = lat_jet_naml_za[jet_30_70].mean()
#js_30_70 = np.linspace(js_30_70c-20, js_30_70c+20, o3_anom_za.shape[1])
#js_70_100c = lat_jet_naml_za[jet_70_100].mean()
#js_70_100 = np.linspace(js_70_100c-20, js_70_100c+20, o3_anom_za.shape[1])
## Plot zonally-averaged O3 anomaly about the jet and indicate mean jet 
## latitude for percentile categories with vertical line
#fig = plt.figure()
#ax = plt.subplot2grid((1,1), (0,0))
#ax.plot(js_0_30, np.nanmean(o3_anom_za[jet_0_30], axis=0), '-', lw=2, 
#        color='#1b9e77', zorder=10,
#        label='$\mathregular{\phi_{jet}}$ 0-30th percentile')
#ax.axvline(x=js_0_30c, color='#1b9e77', linestyle='--', lw=1, 
#           zorder=1)
#ax.plot(js_30_70, np.nanmean(o3_anom_za[jet_30_70], axis=0), '-', lw=2, 
#        color='#d95f02', zorder=11,
#        label='$\mathregular{\phi_{jet}}$ 30-70th percentile')
#ax.axvline(x=js_30_70c, color='#d95f02', linestyle='--', lw=1, 
#           zorder=1)
#ax.plot(js_70_100, np.nanmean(o3_anom_za[jet_70_100], axis=0), '-', lw=2,
#        color='#7570b3', zorder=12,
#        label='$\mathregular{\phi_{jet}}$ 70-100th percentile')
#ax.axvline(x=js_70_100c, color='#7570b3', linestyle='--', lw=1, 
#           zorder=1)
#plt.legend(loc=9,  ncol=3, bbox_to_anchor=(0.5, 1.15), frameon=False)
#ax.set_xlabel('Latitude [$\mathregular{^{\circ}}$]')
#ax.set_ylabel('$\mathregular{\delta}$O$_\mathregular{3}$ [ppbv]')
#plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+'anomo3_latitude_naml.eps',
#            dpi=300)
#plt.show()
## Plot zonally-averaged 2-meter temperature anomaly about the jet and indicate 
## mean jet latitude for percentile categories with vertical line
#fig = plt.figure()
#ax = plt.subplot2grid((1,1), (0,0))
#ax.plot(js_0_30, np.nanmean(t2m_anom_za[jet_0_30], axis=0), '-', lw=2, 
#        color='#1b9e77', zorder=10,
#        label='$\mathregular{\phi_{jet}}$ 0-30th percentile')
#ax.axvline(x=js_0_30c, color='#1b9e77', linestyle='--', lw=1, 
#           zorder=1)
#ax.plot(js_30_70, np.nanmean(t2m_anom_za[jet_30_70], axis=0), '-', lw=2, 
#        color='#d95f02', zorder=11,
#        label='$\mathregular{\phi_{jet}}$ 30-70th percentile')
#ax.axvline(x=js_30_70c, color='#d95f02', linestyle='--', lw=1, 
#           zorder=1)
#ax.plot(js_70_100, np.nanmean(t2m_anom_za[jet_70_100], axis=0), '-', lw=2,
#        color='#7570b3', zorder=12,
#        label='$\mathregular{\phi_{jet}}$ 70-100th percentile')
#ax.axvline(x=js_70_100c, color='#7570b3', linestyle='--', lw=1, 
#           zorder=1)
#plt.legend(loc=9,  ncol=3, bbox_to_anchor=(0.5, 1.15), frameon=False)
#ax.set_xlabel('Latitude [$\mathregular{^{\circ}}$]')
#ax.set_ylabel('$\mathregular{\delta}$T [K]')
#plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+'anomt2m_latitude_naml.eps',
#            dpi=300)
#plt.show()
## Map of O3 anomalies and jet latitude/variability for days where the jet 
## latitude is within the 0-30th percentile
#map_nh(lat_naml, 
#    lng_naml, 
#    np.nanmean(o3_naml[jet_0_30], axis=0)-np.nanmean(o3_naml, axis=0), 
#    '$\mathregular{\delta}$ O$_\mathregular{3}$ for $\mathregular{\phi}_'+
#    '{\mathregular{jet}}$ 0-30th percentile',    
#    '[ppbv]', 
#    np.linspace(-5, 5, 11), 
#    'bwr', 
#    'northamerica_anomo3_jet_jetp0p30',
#    e_n=np.nanmean(lat_jet_naml[jet_0_30],axis=0), 
#    eerr_n=np.nanstd(lat_jet_naml[jet_0_30],axis=0), 
#    extent=[lng_naml.min(), lng_naml.max(), lat_naml.min(), lat_naml.max()-5], 
#    oceanon='yes', 
#    extend='both')
## 30-70th percentile
#map_nh(lat_naml, 
#    lng_naml, 
#    np.nanmean(o3_naml[jet_30_70], axis=0)-np.nanmean(o3_naml, axis=0), 
#    '$\mathregular{\delta}$ O$_\mathregular{3}$ for $\mathregular{\phi}_'+
#    '{\mathregular{jet}}$ 30-70th percentile',    
#    '[ppbv]', 
#    np.linspace(-5, 5, 11), 
#    'bwr', 
#    'northamerica_anomo3_jet_jetp30p70',
#    e_n=np.nanmean(lat_jet_naml[jet_30_70],axis=0), 
#    eerr_n=np.nanstd(lat_jet_naml[jet_30_70],axis=0), 
#    extent=[lng_naml.min(), lng_naml.max(), lat_naml.min(), lat_naml.max()-5], 
#    oceanon='yes', 
#    extend='both')
## 70-100th percentile
#map_nh(lat_naml, 
#    lng_naml, 
#    np.nanmean(o3_naml[jet_70_100], axis=0)-np.nanmean(o3_naml, axis=0), 
#    '$\mathregular{\delta}$ O$_\mathregular{3}$ for $\mathregular{\phi}_'+
#    '{\mathregular{jet}}$ 70-100th percentile',    
#    '[ppbv]', 
#    np.linspace(-5, 5, 11), 
#    'bwr', 
#    'northamerica_anomo3_jet_jetp70p100',
#    e_n=np.nanmean(lat_jet_naml[jet_70_100],axis=0), 
#    eerr_n=np.nanstd(lat_jet_naml[jet_70_100],axis=0), 
#    extent=[lng_naml.min(), lng_naml.max(), lat_naml.min(), lat_naml.max()-5], 
#    oceanon='yes', 
#    extend='both')
    
"""LOAD AND PLOT OUTPUT FROM GEOS-C1SD"""
#import numpy as np
#import sys
#sys.path.append('/Users/ghkerr/phd/GMI/')
#from geo_idx import geo_idx
#sys.path.append('/Users/ghkerr/phd/globalo3/')
#import globalo3_open, globalo3_calculate
#sys.path.append('/Users/ghkerr/phd/transporto3/')
#years_gsd = [1990, 1991, 1992, 1993, 1994]
#latmin, latmax, lngmin, lngmax = -2., 90., 0., 360.
## Open GEOS-C1SD CO emissions with 25 day lifetime output
#co25_gsd, lat_gsd, lng_gsd, pressure_co25_gsd = globalo3_open.open_geos_c1sd(
#    years_gsd, 'co25', 1000., 800., lngmin, latmax, lngmax, latmin, 
#    columnmean=True)
## Open GEOS-C1SD CO emissions with 50 day lifetime output
#co50_gsd, lat_gsd, lng_gsd, pressure_co50_gsd = globalo3_open.open_geos_c1sd(
#    years_gsd, 'co50', 1000., 800., lngmin, latmax, lngmax, latmin, 
#    columnmean=True)
## Open GEOS-C1SD output with fixed mixing ratio at the surface between 
## 30 and 50˚N and 50 day lifetime
#nh50_gsd, lat_gsd, lng_gsd, pressure_nh50_gsd = globalo3_open.open_geos_c1sd(
#    years_gsd, 'nh50', 1000., 800., lngmin, latmax, lngmax, latmin, 
#    columnmean=True)
## Open GEOS-C1SD output with fixed mixing ratio in the stratosphere (> 80 hPa)
## and a 25 day lifetime in the troposphere
#st80_25_gsd, lat_gsd, lng_gsd, pressure_st80_25_gsd = \
#    globalo3_open.open_geos_c1sd(years_gsd, 'st80_25', 1000., 800., lngmin, 
#    latmax, lngmax, latmin, columnmean=True)
## Open GEOS-C1SD 500 hPa U wind
#U500_gsd, lat_gsd, lng_gsd, pressure_U_gsd = globalo3_open.open_geos_c1sd(
#    years_gsd, 'U', 500., 500., lngmin, latmax, lngmax, latmin, 
#    columnmean=True)
## Open GEOS-C1SD surface pressure
#PS_gsd, lat_gsd, lng_gsd = globalo3_open.open_geos_c1sd(years_gsd, 'PS', 0., 
#    0., lngmin, latmax, lngmax, latmin) # note that levmax, levmin parameters 
#    # are dummies 
## Open GEOS-C1SD temperature
#T_gsd, lat_gsd, lng_gsd, pressure_T_gsd = globalo3_open.open_geos_c1sd(
#    years_gsd, 'T', 1000., 800., lngmin, latmax, lngmax, latmin, 
#    columnmean=True)
## Add cyclic point to longitude coordinates so that it wraps around the Prime
## Meridian when plotting
#lng_gsd[-1]=360.
## Plot mean fields for GEOS-C1SD simulations
#map_nh(lat_gsd, lng_gsd, np.nanmean(co25_gsd, axis=0)*1e9, 
#   'Mean 1000-800 hPa column', 'CO25 [ppbv]', 
#   np.linspace(0., 100., 6), 'PuBu', 'co25_1000-800hPa_%d-%d'%(years_gsd[0],
#   years_gsd[-1]), oceanon='no', extend='max')
#map_nh(lat_gsd, lng_gsd, np.nanmean(co50_gsd, axis=0)*1e9, 
#   'Mean 1000-800 hPa column', 'CO50 [ppbv]', 
#   np.linspace(0., 100., 6), 'PuBu', 'co50_1000-800hPa_%d-%d'%(years_gsd[0],
#   years_gsd[-1]), oceanon='no', extend='max')
#map_nh(lat_gsd, lng_gsd, np.nanmean(nh50_gsd, axis=0)*1e9, 
#   'Mean 1000-800 hPa column', 'NH50 [ppbv]', 
#   np.linspace(0., 100000., 6), 'PuBu', 'nh50_1000-800hPa_%d-%d'%(years_gsd[0],
#   years_gsd[-1]), oceanon='no', extend='max')
#map_nh(lat_gsd, lng_gsd, np.nanmean(st80_25_gsd, axis=0)*1e9, 
#   'Mean 1000-800 hPa column', 'ST80_25 [ppbv]', 
#   np.linspace(0., 1., 6), 'PuBu', 'st80_25_1000-800hPa_%d-%d'%(years_gsd[0],
#   years_gsd[-1]), oceanon='no', extend='max')
#map_nh(lat_gsd, lng_gsd, np.nanmean(U500_gsd, axis=0), 
#   '500 hPa', 'U [m s$^{\mathregular{-1}}$]', np.linspace(-16, 16, 9), 
#   'bwr', 'u_500hPa_%d-%d'%(years_gsd[0], years_gsd[-1]), oceanon='no', extend='both')
## Calculate tracer-temperature correlation and plot 
#r_tco25 = globalo3_calculate.calculate_r(T_gsd, co25_gsd, lat_gsd, lng_gsd)
#map_nh(lat_gsd, lng_gsd, r_tco25, r'$\it{r}\:$(T, CO25)', '', 
#    np.linspace(-1., 1., 11), 'bwr', 'r_tco25gsd_%d-%d'%(years_gsd[0], 
#    years_gsd[-1]), oceanon='no', extend='both')
#r_tco50 = globalo3_calculate.calculate_r(T_gsd, co50_gsd, lat_gsd, lng_gsd)
#map_nh(lat_gsd, lng_gsd, r_tco50, r'$\it{r}\:$(T, CO50)', '', 
#    np.linspace(-1., 1., 11), 'bwr', 'r_tco50gsd_%d-%d'%(years_gsd[0], 
#    years_gsd[-1]), oceanon='no', extend='both')
#r_tnh50 = globalo3_calculate.calculate_r(T_gsd, nh50_gsd, lat_gsd, lng_gsd)
#map_nh(lat_gsd, lng_gsd, r_tnh50, r'$\it{r}\:$(T, NH50)', '', 
#    np.linspace(-1., 1., 11), 'bwr', 'r_tnh50gsd_%d-%d'%(years_gsd[0], 
#    years_gsd[-1]), oceanon='no', extend='both')
#r_tst80_25 = globalo3_calculate.calculate_r(T_gsd, st80_25_gsd, lat_gsd, 
#    lng_gsd)
#map_nh(lat_gsd, lng_gsd, r_tst80_25, r'$\it{r}\:$(T, ST80_25)', '', 
#    np.linspace(-1., 1., 11), 'bwr', 'r_tst80_25gsd_%d-%d'%(years_gsd[0], 
#    years_gsd[-1]), oceanon='no', extend='both')
## Subset fields of interest in the Northern Hemisphere mid-latitudes
#U500_gsd_nhml, lat_gsd_nhml, lng_gsd_nhml = \
#    globalo3_calculate.find_grid_in_bb(U500_gsd, lat_gsd, lng_gsd, 0., 360., 
#    23., 60.)
#co25_gsd_nhml, lat_gsd_nhml, lng_gsd_nhml = globalo3_calculate.find_grid_in_bb(
#    co25_gsd, lat_gsd, lng_gsd, 0., 360., 23., 60.)
#co50_gsd_nhml, lat_gsd_nhml, lng_gsd_nhml = globalo3_calculate.find_grid_in_bb(
#    co50_gsd, lat_gsd, lng_gsd, 0., 360., 23., 60.)
#nh50_gsd_nhml, lat_gsd_nhml, lng_gsd_nhml = globalo3_calculate.find_grid_in_bb(
#    nh50_gsd, lat_gsd, lng_gsd, 0., 360., 23., 60.)
#st80_25_gsd_nhml, lat_gsd_nhml, lng_gsd_nhml = \
#    globalo3_calculate.find_grid_in_bb(st80_25_gsd, lat_gsd, lng_gsd, 0., 
#    360., 23., 60.)
#land_nhml = globalo3_calculate.find_grid_overland(lat_gsd, lng_gsd)
## Identify eddy-driven jet and tracers in the vicinity of the jet 
#lat_jet_gsd_nhml, co25_gsd_nhml = globalo3_calculate.find_field_atjet(
#    co25_gsd_nhml, U500_gsd_nhml, lat_gsd_nhml, lng_gsd_nhml, 10, anom=True)
#lat_jet_gsd_nhml, co50_gsd_nhml = globalo3_calculate.find_field_atjet(
#    co50_gsd_nhml, U500_gsd_nhml, lat_gsd_nhml, lng_gsd_nhml, 10, anom=True)
#lat_jet_gsd_nhml, nh50_gsd_nhml = globalo3_calculate.find_field_atjet(
#    nh50_gsd_nhml, U500_gsd_nhml, lat_gsd_nhml, lng_gsd_nhml, 10, anom=True)
#lat_jet_gsd_nhml, st80_25_gsd_nhml = globalo3_calculate.find_field_atjet(
#    st80_25_gsd_nhml, U500_gsd_nhml, lat_gsd_nhml, lng_gsd_nhml, 10, anom=True)
## Calculate CO25-jet relationship
#m_gsd_co25jetdist, r_gsd_co25jetdist = \
#    globalo3_calculate.calculate_o3jet_relationship(co25_gsd, lat_gsd, lng_gsd, 
#    lat_jet_gsd_nhml, lng_gsd_nhml)
## Plot tracer anomalies about the eddy-driven jet
#edjetlocation_fieldatedjet(lat_gsd_nhml, lng_gsd_nhml, lat_jet_gsd_nhml, 
#    co25_gsd_nhml*1e9, 'bwr', r'$\mathregular{\Delta}$ CO25 [ppbv]', 
#    np.linspace(-50., 50., 6), 'co25anom', 
#    '1000-800hPa_%d-%d'%(years_gsd[0], years_gsd[-1]), skiplng=4)
#edjetlocation_fieldatedjet(lat_gsd_nhml, lng_gsd_nhml, lat_jet_gsd_nhml, 
#    co50_gsd_nhml*1e9, 'bwr', r'$\mathregular{\Delta}$ CO50 [ppbv]', 
#    np.linspace(-50., 50., 6), 'co50anom', 
#    '1000-800hPa_%d-%d'%(years_gsd[0], years_gsd[-1]), skiplng=4)
#edjetlocation_fieldatedjet(lat_gsd_nhml, lng_gsd_nhml, lat_jet_gsd_nhml, 
#    nh50_gsd_nhml*1e9, 'bwr', r'$\mathregular{\Delta}$ NH50 [ppbv]',  
#    np.linspace(-50000., 50000., 6), 'nh50anom', 
#    '1000-800hPa_%d-%d'%(years_gsd[0], years_gsd[-1]), skiplng=4)
#edjetlocation_fieldatedjet(lat_gsd_nhml, lng_gsd_nhml, lat_jet_gsd_nhml, 
#    st80_25_gsd_nhml*1e9, 'bwr', r'$\mathregular{\Delta}$ ST80_25 [ppbv]',  
#    np.linspace(-0.5, 0.5, 6), 'st80_25anom', 
#    '1000-800hPa_%d-%d'%(years_gsd[0], years_gsd[-1]), skiplng=4)
## Find indices corresponding to the focus region (Eastern North America)
#left = geo_idx(269., lng_gsd_nhml)
#right = geo_idx(294., lng_gsd_nhml)
## For CO25 tracer
#co25_gsd_ena = co25_gsd_nhml[:, :, left:right+1]
#lng_gsd_ena = lng_gsd[left:right+1]
#field_vj = co25_gsd_ena*1e9
## Half of latitudinal span of field
#splitint = int(len(lat_gsd_nhml)/2.)
## Split field in the vicinity of jet to the field above/below the jet
#field_vj_below = field_vj[:, :splitint]
#field_vj_above = field_vj[:, -splitint:]
## Find regional average of field above/below jet
#field_vj_below_ra = np.nanmean(field_vj_below, axis=tuple((1,2)))
## Determine days with the largest increases (tracer anomaly > 0) equatorward 
## of the jet and largest decreases (tracer anomaly < 0) equatorward 
## of the jet
#below_increase = np.where(field_vj_below_ra > 
#                     np.nanpercentile(field_vj_below_ra, 80))[0]
#below_decrease = np.where(field_vj_below_ra < 
#                     np.nanpercentile(field_vj_below_ra, 20))[0]
## Jet location, CO25 anomaly in the vicinity of jet and the 500 hPa u-wind 
## anomaly on days where the CO25 anomaly below the jet has largest increase
## (> 80th percentile)
#edjetlocation_fieldatedjet(lat_gsd_nhml, lng_gsd_nhml, 
#    lat_jet_gsd_nhml[below_increase], co25_gsd_nhml[below_increase]*1e9, 
#    'bwr', r'$\mathregular{\Delta}$ CO25 [ppbv]',  np.linspace(-50., 50., 6), 
#    'co25anom_positivebelowjet_', '1000-800hPa_%d-%d'%(years_gsd[0], years_gsd[-1]), 
#    skiplng=4)
#map_nh(lat_gsd, lng_gsd, (np.nanmean(U500_gsd[below_increase], axis=0)-
#    np.nanmean(U500_gsd, axis=0)), '500 hPa | Sfc. Pressure', 
#    '$\mathregular{\Delta}$ U [m s$^{\mathregular{-1}}$]', 
#    np.linspace(-10, 10, 6), 'bwr', 
#    'uanomalyforco25_positivebelowjet_500hPa_%d-%d'%(years_gsd[0], years_gsd[-1]), 
#    contour=(np.nanmean(PS_gsd[below_increase], axis=0)-
#    np.nanmean(PS_gsd, axis=0))/100., contour_levs=[-3,-2,-1,1,2,3], 
#    oceanon='no', extend='both')
## Jet location, CO25 anomaly in the vicinity of jet and the 500 hPa u-wind 
## anomaly on days where the CO25 anomaly below the jet has largest decrease
## (< 20th percentile)
#edjetlocation_fieldatedjet(lat_gsd_nhml, lng_gsd_nhml, 
#    lat_jet_gsd_nhml[below_decrease], 
#    co25_gsd_nhml[below_decrease]*1e9, 
#    'bwr', r'$\mathregular{\Delta}$ CO25 [ppbv]',  np.linspace(-50., 50., 6), 
#    'co25anom_negativebelowjet_', '1000-800hPa_%d-%d'%(years_gsd[0], years_gsd[-1]), 
#    skiplng=4)
#map_nh(lat_gsd, lng_gsd, (np.nanmean(U500_gsd[below_decrease], axis=0)-
#    np.nanmean(U500_gsd, axis=0)), '500 hPa | Sfc. Pressure', 
#    '$\mathregular{\Delta}$ U [m s$^{\mathregular{-1}}$]', 
#    np.linspace(-10, 10, 6), 'bwr', 
#    'uanomalyforco25_negativebelowjet_500hPa_%d-%d'%(years_gsd[0], years_gsd[-1]), 
#    contour=(np.nanmean(PS_gsd[below_decrease], axis=0)-
#    np.nanmean(PS_gsd, axis=0))/100., contour_levs=[-3,-2,-1,1,2,3],
#    oceanon='no', extend='both')
## Correlation between CO25 tracer and distance from eddy-driven jet from 
## GEOS-C1SD    
#map_nh(lat_gsd, lng_gsd, r_gsd_co25jetdist, '', 
#    r'$\it{r}\:$(CO25, $\mathregular{\phi}_{\mathregular{jet}}$ $-$ '+
#    '$\mathregular{\phi}$)', np.linspace(-0.7, 0.7, 8), 'bwr', 
#    'r_co25jetdist_%d-%d_jet'%(years_gsd[0], years_gsd[-1]),
#    e_n=np.nanmean(lat_jet_gsd_nhml, axis=0), 
#    eerr_n=np.zeros(lat_jet_gsd_nhml.shape[1]), oceanon='yes')
## Change in CO25 with change in eddy-driven jet position from GMI CTM 
#map_nh(lat_gsd, lng_gsd, m_gsd_co25jetdist*1e9, '', 
#   r'$\mathregular{\Delta}$ CO25 [ppbv $^{\circ \mathregular{-1}}$]',
#   np.linspace(-0.2, 0.4, 7), 'gist_earth_r', 'dco25djet_%d-%d_jet'
#   %(years_gsd[0],years_gsd[-1]), e_n=np.nanmean(lat_jet_gsd_nhml, axis=0), 
#   eerr_n=np.zeros(lat_jet_gsd_nhml.shape[1]), oceanon='yes')

"""RELATIONSHIP OF (ANTI)CYCLONES AND EDDY-DRIVEN JET"""
#latmin_n = 25.
#latmax_n = 90.
#lngmin_n = 0.
#lngmax_n = 360.
#lat_gmi_n, lng_gmi_n, times_n, o3_n = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_n, 
#    latmax_n, lngmin_n, lngmax_n, 'O3', 'HindcastMR2')
## Load daily-averaged sea level pressure
#hours = [0,3,6,9,12,15,18,21]
#SLP, mtime, lat_merra_n, lng_merra_n = transporto3_open.open_merra2(years, 
#    hours, 'SLP', 'inst3_3d_asm_Np', 'JJA_500mb.nc', lngmin_n, latmax_n, 
#    lngmax_n, latmin_n, dailyavg='yes')
## Interpolate MERRA-2 to the resolution of the CTM
#SLP = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi_n, lng_gmi_n, 
#    lat_merra_n, lng_merra_n, SLP)
## Identify anticyclones/cyclones over the hemispheric domain and their 
## coordinates
#cyclones, cyclones_daycoord, cyclones_ycoord, cyclones_xcoord = \
#    globalo3_calculate.identify_SLPcenter(lat_gmi_n, lng_gmi_n, SLP, 10, 10, 
#    'cyclone', 100000., years, checkplot='yes', fstr='cyclones_%d-%d'%(
#    years[0], years[-1]))
#anticyclones, anticyclones_daycoord, anticyclones_ycoord, anticyclones_xcoord = \
#    globalo3_calculate.identify_SLPcenter(lat_gmi_n, lng_gmi_n, SLP, 10, 10, 
#    'anticyclone', 101200., years, checkplot='yes', fstr='anticyclones_%d-%d'%(
#    years[0], years[-1]))
#map_jet_centerdist(cyclones, lat_jet_nhml, lng_nhml, lat_gmi_n, lng_gmi_n, 
#    'cyclone')
#map_jet_centerdist(anticyclones, lat_jet_nhml, lng_nhml, lat_gmi_n, lng_gmi_n, 
#    'anticyclone')

"""IDENTIFY (ANTI)CYCLONES AND COINCIDENT O3"""
#latmin_na = 30.
#latmax_na = 90.
#lngmin_na = 200.
#lngmax_na = 350.
## Load GMI CTM O3 over North America
#lat_gmi_na, lng_gmi_na, times_na, o3_na = \
#globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_na, 
#latmax_na, lngmin_na, lngmax_na, 'O3', 'HindcastMR2')
#o3_na = o3_na*1e9
## Load GMI CTM O3 over North America from EmFix simulation
#lat_gmi_na_emfix, lng_gmi_na_emfix, times_na_emfix, o3_na_emfix = \
#globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_na, 
#latmax_na, lngmin_na, lngmax_na, 'O3', 'Hindcast3Igac2')
#o3_na_emfix = o3_na_emfix*1e9
## Interpolate EmFix to the resolution of HindcastMR2 (not the best practice, 
## change in the future)
#o3_na_emfix = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi_na, 
#    lng_gmi_na, lat_gmi_na_emfix, lng_gmi_na_emfix, o3_na_emfix)
## Load GMI CTM O3 over North America from Std simulation
#lat_gmi_na_std, lng_gmi_na_std, times_na_std, o3_na_std = \
#globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_na, 
#latmax_na, lngmin_na, lngmax_na, 'O3', 'HindcastFFIgac2')
#o3_na_std = o3_na_std*1e9
#o3_na_std = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi_na, 
#    lng_gmi_na, lat_gmi_na_std, lng_gmi_na_std, o3_na_std)
## Load daily-averaged sea level pressure
#hours = [0,3,6,9,12,15,18,21]
#SLP_na, mtime, lat_merra_na, lng_merra_na = transporto3_open.open_merra2(years, 
#    hours, 'SLP', 'inst3_3d_asm_Np', 'JJA_500mb.nc', lngmin_na, latmax_na, 
#    lngmax_na, latmin_na, dailyavg='yes')
## Interpolate MERRA-2 to the resolution of the CTM
#SLP_na = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi_na, lng_gmi_na, 
#    lat_merra_na, lng_merra_na, SLP_na)
#H_na, mtime, lat_merra_na, lng_merra_na = transporto3_open.open_merra2(years, 
#    hours, 'H', 'inst3_3d_asm_Np', 'JJA_500mb.nc', lngmin_na, latmax_na, 
#    lngmax_na, latmin_na, dailyavg='yes')
## Interpolate MERRA-2 to the resolution of the CTM
#H_na = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi_na, lng_gmi_na, 
#    lat_merra_na, lng_merra_na, H_na)        
## Identify anticyclones/cyclones in North America and their coordinates
#cyclones, cyclones_daycoord, cyclones_ycoord, cyclones_xcoord = \
#    globalo3_calculate.identify_SLPcenter(lat_gmi_na, lng_gmi_na, SLP_na, 10, 
#    10, 'cyclone', 100000., years, checkplot='yes', fstr='cyclones_%d-%d'%(
#    years[0], years[-1]))
#anticyclones, anticyclones_daycoord, anticyclones_ycoord, anticyclones_xcoord = \
#    globalo3_calculate.identify_SLPcenter(lat_gmi_na, lng_gmi_na, SLP_na, 10, 
#    10, 'anticyclone', 101200., years, checkplot='yes', 
#    fstr='anticyclones_%d-%d'%(years[0], years[-1]))
## Segregate (anti)cyclones based their position above or below the jet
#cyclones_abovejet, cyclones_belowjet = globalo3_calculate.filter_center_byjet(
#    cyclones, lat_jet_nhml, lat_gmi_na, lng_gmi_na, lat_nhml, lng_nhml)
#anticyclones_abovejet, anticyclones_belowjet = \
#    globalo3_calculate.filter_center_byjet(anticyclones, lat_jet_nhml, 
#    lat_gmi_na, lng_gmi_na, lat_nhml, lng_nhml)
## O3, SLP and Z500 around all cyclones in North America
#contourf_var_atcenter(cyclones, o3_na-np.mean(o3_na,axis=0), lat_gmi_na, 
#    lng_gmi_na, 'o3', 15, 'cyclone', 'O$_{\mathregular{3}}$ [ppbv]', 
#    np.linspace(-4, 4, 9), 'bwr', 'cyclones', 
#    SLP=SLP_na-np.nanmean(SLP_na, axis=0), H=H_na-np.nanmean(H_na,axis=0))
## For cyclones below the jet
#contourf_var_atcenter(cyclones_belowjet, o3_na-np.mean(o3_na,axis=0), 
#    lat_gmi_na, lng_gmi_na, 'o3', 15, 'cyclone', 
#    'O$_{\mathregular{3}}$ [ppbv]', np.linspace(-4, 4, 9), 'bwr', 
#    'cyclones_belowjet', SLP=SLP_na-np.nanmean(SLP_na, axis=0), 
#    H=H_na-np.nanmean(H_na,axis=0))
## For cyclones below the jet with and without emissions controls (EmFix/Std)
#contourf_var_atcenter(cyclones_belowjet, o3_na_std-np.mean(o3_na_std, axis=0),
#    lat_gmi_na, lng_gmi_na, 'o3', 15, 'cyclone', 
#    'O$_{\mathregular{3}}$ [ppbv]', np.linspace(-4, 4, 9), 'bwr', 
#    'cyclones_belowjet_std', SLP=SLP_na-np.nanmean(SLP_na, axis=0), 
#    H=H_na-np.nanmean(H_na,axis=0))
#contourf_var_atcenter(cyclones_belowjet, o3_na_emfix-
#    np.mean(o3_na_emfix,axis=0), lat_gmi_na, lng_gmi_na, 'o3', 15, 'cyclone', 
#    'O$_{\mathregular{3}}$ [ppbv]', np.linspace(-4, 4, 9), 'bwr', 
#    'cyclones_belowjet_emfix', SLP=SLP_na-np.nanmean(SLP_na, axis=0), 
#    H=H_na-np.nanmean(H_na,axis=0))
## For cyclones above the jet
#contourf_var_atcenter(cyclones_abovejet, o3_na-np.mean(o3_na,axis=0), 
#    lat_gmi_na, lng_gmi_na, 'o3', 15, 'cyclone', 
#    'O$_{\mathregular{3}}$ [ppbv]', np.linspace(-4, 4, 9), 'bwr', 
#    'cyclones_abovejet', SLP=SLP_na-np.nanmean(SLP_na, axis=0), 
#    H=H_na-np.nanmean(H_na,axis=0))
## O3, SLP and Z500 around all anticyclones in North America
#contourf_var_atcenter(anticyclones, o3_na-np.mean(o3_na,axis=0), lat_gmi_na, 
#    lng_gmi_na, 'o3', 15, 'anticyclone', 'O$_{\mathregular{3}}$ [ppbv]', 
#    np.linspace(-4, 4, 9), 'bwr', 'cyclones', 
#    SLP=SLP_na-np.nanmean(SLP_na, axis=0), H=H_na-np.nanmean(H_na,axis=0))
## For anticyclones below the jet
#contourf_var_atcenter(anticyclones_belowjet, o3_na-np.mean(o3_na,axis=0), 
#    lat_gmi_na, lng_gmi_na, 'o3', 15, 'anticyclone', 
#    'O$_{\mathregular{3}}$ [ppbv]', np.linspace(-4, 4, 9), 'bwr', 
#    'anticyclones_belowjet', SLP=SLP_na-np.nanmean(SLP_na, axis=0), 
#    H=H_na-np.nanmean(H_na,axis=0))
## For anticyclones below the jet with and without emissions controls 
## (EmFix/Std)
#contourf_var_atcenter(anticyclones_belowjet, o3_na_std-
#    np.mean(o3_na_std,axis=0), lat_gmi_na, lng_gmi_na, 'o3', 15, 'anticyclone', 
#    'O$_{\mathregular{3}}$ [ppbv]', np.linspace(-4, 4, 9), 'bwr', 
#    'anticyclones_belowjet_std', SLP=SLP_na-np.nanmean(SLP_na, axis=0), 
#    H=H_na-np.nanmean(H_na,axis=0))
#contourf_var_atcenter(anticyclones_belowjet, o3_na_emfix-
#    np.mean(o3_na_emfix,axis=0), lat_gmi_na, lng_gmi_na, 'o3', 15, 
#    'anticyclone', 'O$_{\mathregular{3}}$ [ppbv]', np.linspace(-4, 4, 9), 
#    'bwr', 'anticyclones_belowjet_emfix', 
#    SLP=SLP_na-np.nanmean(SLP_na, axis=0), H=H_na-np.nanmean(H_na,axis=0))
## For anticyclones above the jet
#contourf_var_atcenter(anticyclones_abovejet, o3_na-np.mean(o3_na,axis=0), 
#    lat_gmi_na, lng_gmi_na, 'o3', 15, 'anticyclone', 
#    'O$_{\mathregular{3}}$ [ppbv]', np.linspace(-4, 4, 9), 'bwr', 
#    'anticyclones_abovejet', SLP=SLP_na-np.nanmean(SLP_na, axis=0), 
#    H=H_na-np.nanmean(H_na,axis=0))

"""SLP AND O3 OVER NORTH ATLANTIC ON DAYS WITH CYCLONES"""
#from datetime import datetime
#mtimes = [13, 14, 272, 273, 274]
#for day in mtimes: 
#    fig = plt.figure()#figsize=(8,3.5))
#    axl = plt.subplot2grid((1,2), (0,0), projection=ccrs.PlateCarree())
#    axl.set_title('%s' %datetime.strftime(mtime[day], '%m/%d/%Y'), 
#                  x=0.05, ha='left')
#    axr = plt.subplot2grid((1,2), (0,1), projection=ccrs.PlateCarree())    
#    # SLP
#    clevs = np.arange(992, 1024, 4)
#    mbl = axl.contourf(lng_gmi_na, lat_gmi_na, SLP[day]/100., clevs, 
#                       cmap=plt.get_cmap('rainbow'), extend='both', 
#                       transform=ccrs.PlateCarree())  
#    plt.colorbar(mbl, fraction=0.03, pad=0.04, ax=axl, label='Pressure [hPa]')
#    axl.coastlines(lw=0.25, color='k')
#    axl.set_extent([-75, -25, 30, 55])
#    # O3
#    clevs = np.linspace(20, 40, 21) 
#    mbr = axr.contourf(lng_gmi_na, lat_gmi_na, o3_na[day], clevs,
#                     cmap=plt.get_cmap('gist_earth_r'), extend='both',
#                     transform=ccrs.PlateCarree())
#    axr.coastlines(lw=0.25, color='k')
#    plt.colorbar(mbr, fraction=0.03, pad=0.04, ax=axr, 
#                 label='O$_{\mathregular{3}}$ [ppbv]')
#    axr.coastlines(lw=0.25, color='k')
#    axr.set_extent([-75, -25, 30, 55])
#    plt.subplots_adjust(wspace=0.4)
#    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#                'map_SLP_o3_northatlantic_%s.eps'
#                %datetime.strftime(mtime[day], '%m-%d-%Y'), dpi=300)  
#    plt.show()

"""ANIMATION OF 500 HPA WINDS AND O3/O3 ANOMALY FOR JJA 2008"""
#for day in np.arange(92+92, 92+92+92, 1):
#    map_nh(lat_gmi_n, lng_gmi_n, o3_n[day], 
#           datetime.strftime(mtime[day], '%m/%d/%Y'), 
#           'O$_{\mathregular{3}}$ [ppbv]', np.linspace(35, 60, 6), 'PuBu', 
#           'o3_wind500hPa_%.2d' %day, extent=[-135., -55, 20., 50.], 
#           quiver=(U500[day],V500[day]), oceanon='no')
#    plt.show()
#    map_nh(lat_gmi_n, lng_gmi_n, o3_n[day]-np.mean(o3_n, axis=0),
#           datetime.strftime(mtime[day], '%m/%d/%Y'),
#           'O$_{\mathregular{3}}$ [ppbv]', np.linspace(-15, 15, 11), 'bwr', 
#           'o3anom_wind500hPa_%.2d' %day, extent=[-135., -55, 20., 50.], 
#           quiver=(U500[day],V500[day]), oceanon='no')
    
"""PLOT MEAN LOCATION OF JET STREAM (MAX U WINDS AT 500 HPA) AND 
   2-METER TEMPERATURE WITHIN +/- 10 DEGREES OF JET STREAM"""
#t2m_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(t2m_n, 
#    lat_gmi_n, lng_gmi_n, 0., 360., 23., 60.)   
#t2m_fr = t2m_nhml
#t2m_jet = np.empty(shape=(len(t2m_fr), 2*jetdistance+1, len(lng_fr)))
#t2m_jet[:] = np.nan
#for day in np.arange(0, len(U500_fr), 1):
#    U500_day = U500_fr[day]
#    # Loop through longitude
#    for i in np.arange(0, len(lng_fr), 1):    
#        # U wind at 500 hPa for longitude/day of interest 
#        U500_transect = U500_day[:, i]
#        U500max = np.where(U500_transect==U500_transect.max())[0][0]   
#        # RH at jet 
#        t2m_jetdayi = t2m_fr[day, U500max, i]       
#        # Longitudinal cross-section of O3 within the number of grid cells 
#        # specified by variable jetdistance on each side ("above" and "below")
#        # of the jet
#        t2m_transect_above = t2m_fr[day, U500max+1:
#            U500max+jetdistance+1, i]-t2m_fr[day, U500max, i]
#        t2m_transect_below = t2m_fr[day, U500max-jetdistance:
#            U500max, i]-t2m_fr[day, U500max, i]     
#        # O3 at jet maximum 
#        t2m_jet[day, jetdistance, i] = 0.0
#        # O3 above jet maximum
#        t2m_jet[day, jetdistance+1:jetdistance+1+len(t2m_transect_above), 
#                    i] = t2m_transect_above            
#        # O3 below jet maximum 
#        t2m_jet[day, jetdistance-len(t2m_transect_below):jetdistance, 
#                    i] = t2m_transect_below
#for day in np.arange(0, 92, 1):
#    fig = plt.figure()
#    axt = plt.subplot2grid((2,2), (0,0), colspan=2, projection=
#                           ccrs.Miller(central_longitude=0.))
#    axt.set_title('%s'%datetime.datetime.strftime(mtime[day], '%d/%m/%y'))
#    axb = plt.subplot2grid((2,2), (1,0), colspan=2)
#    axt.add_feature(cfeature.OCEAN, zorder=1, lw=0.0, color='lightgrey')
#    axt.coastlines(lw=0.25, color='k')
#    axt.set_extent([lng_fr.min()-180, lng_fr.max()-180, lat_fr.min(), 
#                    lat_fr.max()])    
#    axt.set_xticks([-180,-135,-90,-45,0,45,90,135,180])
#    axt.set_xticklabels([''])
#    # Plot only every X number of longitude values
#    skiplng = 6
#    axt.scatter(lng_fr[::skiplng], jet_lat[day][::skiplng], zorder=10, color='k',
#                s=2, transform=ccrs.PlateCarree())
#    # Roll O3 differences arond the jet so that they align with the top (map)
#    # subplot
#    t2m_jet_roll = np.roll(t2m_jet[day], int(lng_fr.shape[0]/2.), axis=1)
#    mb = axb.contourf(lng_fr, np.arange(-10, 11, 1), 
#                      t2m_jet_roll, np.linspace(-10, 10, 11), 
#                      cmap=plt.get_cmap('bwr'), extend='both')
#    axb.set_xlabel('Longitude [$^{\circ}$]', fontsize=16)
#    axb.set_xticks(np.linspace(0, 360, 9))
#    axb.set_xticklabels([-180,-135,-90,-45,0,45,90,135,180])
#    axb.set_ylabel('Grid cells from jet', fontsize=16)
#    colorbar_axes = plt.gcf().add_axes([0.78,0.18,0.03,0.5])
#    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical')
#    colorbar.ax.tick_params(labelsize=12)
#    plt.gcf().subplots_adjust(right=0.75, hspace=-0.1)
#    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#                'edjetlocation_t2medjet_%.2d.jpg'%day, dpi=200)
#    plt.show()

"""PLOT DO3/DT AND MEAN U500"""
#clevs = np.linspace(0, 3, 7)
#cmap = 'Reds'
#lng_n = lng_gmi_n
#lat_n = lat_gmi_n
#import numpy as np
#import matplotlib.pyplot as plt
#import cartopy.util
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
#plt.figure(figsize=(9,4))
#ax=plt.subplot2grid((1,1), (0,0), colspan=1,
#                    projection=ccrs.Miller(central_longitude=0.))
#ax.add_feature(cfeature.OCEAN, zorder=1, lw=0.0, color='lightgrey')
#ax.coastlines(lw=0.25, color='k')
#ax.set_extent([-180., 180., 0., 85.])    
#cmap = plt.get_cmap(cmap)
#mb = ax.contourf(lng_n, lat_n, do3dt_n, clevs, cmap=cmap,
#                 transform=ccrs.PlateCarree(), extend='both')
## Add cyclic point to longitude/data so that it wraps around the Prime 
## meridian 
#lng_n_cp = cartopy.util.add_cyclic_point(lng_n, coord=None, axis=-1).data
#data = cartopy.util.add_cyclic_point(np.abs(np.mean(U500,axis=0)), 
#                                            coord=None, axis=-1).data
#cs = ax.contour(lng_n_cp, lat_n, data+0.01, [6, 9, 12, 15], colors='k', 
#                linewidths=0.5, transform=ccrs.PlateCarree(), 
#                extend ='both')
#labels = plt.clabel(cs, fontsize=5, inline=True, inline_spacing=2, fmt='%d')
#for l in labels:
#    l.set_rotation(0)
##X, Y = np.meshgrid(lng_n, lat_n)
##lw = np.abs(np.mean(U500,axis=0))/10.
##ax.streamplot(X, Y, np.mean(U500,axis=0), np.mean(V500,axis=0), 
##              color='k', linewidth=lw, transform=ccrs.PlateCarree())
## Add colorbar
#colorbar_axes = plt.gcf().add_axes([0.82,0.25,0.03,0.5])
#colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical')
#colorbar.ax.tick_params(labelsize=12)
#colorbar.set_label('dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]', 
#                   fontsize=16)
#plt.gcf().subplots_adjust(right=0.75, hspace=-0.2)
#plt.savefig('/Users/ghkerr/Desktop/trial.eps', dpi=300)


""" PLOT O3, DO3/DT, AND U500 """
## Find fields over Eastern North American 
#latmin_ena, latmax_ena, lngmin_ena, lngmax_ena = 20., 70., 260., 294.
#o3_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(o3_n, lat_gmi_n, 
#    lng_gmi_n, lngmin_ena, lngmax_ena, latmin_ena, latmax_ena)
#do3dt_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(do3dt_n, 
#    lat_gmi_n, lng_gmi_n, lngmin_ena, lngmax_ena, latmin_ena, latmax_ena)
#r_t2mo3_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(r_t2mo3_n, 
#    lat_gmi_n, lng_gmi_n, lngmin_ena, lngmax_ena, latmin_ena, latmax_ena)
#U500_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(U500, 
#    lat_gmi_n, lng_gmi_n, lngmin_ena, lngmax_ena, latmin_ena, latmax_ena)
#nox_edgar_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
#    nox_edgar_n, lat_gmi_n, lng_gmi_n, lngmin_ena, lngmax_ena, latmin_ena, 
#    latmax_ena)
#land_ena = globalo3_calculate.find_grid_overland(lat_ena, lng_ena)
## Find zonally-averaged fields
#o3_ena_za = np.nanmean(o3_ena*land_ena, axis=tuple((0,2)))
#do3dt_ena_za = np.nanmean(do3dt_ena*land_ena, axis=1)
#r_t2mo3_ena_za = np.nanmean(r_t2mo3_ena*land_ena, axis=1) 
#nox_edgar_ena_za = np.nanmean(nox_edgar_ena*land_ena, axis=tuple((0,2)))
#U500_ena_za = np.nanmean(U500_ena*land_ena, axis=tuple((0,2)))
## Plotting
#fig = plt.figure(figsize=(8,8))
#ax1 = plt.subplot2grid((5, 5), (0, 0), rowspan=5)
#ax2 = plt.subplot2grid((5, 5), (0, 1), rowspan=5)
#ax3 = plt.subplot2grid((5, 5), (0, 2), rowspan=5)
#ax4 = plt.subplot2grid((5, 5), (0, 3), rowspan=5)
#ax5 = plt.subplot2grid((5, 5), (0, 4), rowspan=5)
## U-wind at 500 hPa
#ax1.plot(U500_ena_za, lat_ena, '-k', lw=2.)
#ax1.set_xticks(np.linspace(-4, 14, 5))
#ax1.set_xticklabels(['-4', '', '5', '', '14'], fontsize=12)
#ax1.set_title('U$_{\mathregular{500\:hPa}}$', fontsize=16, ha='center', y=1.02)
#ax1.set_xlabel('[m s$^{\mathregular{-1}}$]', fontsize=14)
## Surface-level O3
#ax2.plot(o3_ena_za , lat_ena, '-k', lw=2.)
#ax2.set_xlim([20, 50])
#ax2.set_xticks(np.linspace(20, 50, 5))
#ax2.set_xticklabels(['20', '', '35', '', '50'], fontsize=12)
#ax2.set_title('O$_{\mathregular{3}}$', fontsize=16, ha='center', y=1.02)
#ax2.set_xlabel('[ppbv]', fontsize=14, labelpad=7)
## EDGAR NOx emissions (n.b., multiply by no. sec/day and 1000 m/km to 
## convert from /s/m2 to /day/km2)
#ax3.plot(nox_edgar_ena_za*86400.*1000*1000., lat_ena, '-k', lw=2., 
#         clip_on=False)
#ax3.set_xlim([0, 36])
#ax3.set_xticks(np.linspace(0, 36, 5))
#ax3.set_xticklabels(['0', '', '18', '', '36'], fontsize=12)
#ax3.set_title('NO$_{{x}}$', fontsize=16, ha='center', y=1.02)
#ax3.set_xlabel('[kg km$^{\mathregular{-2}}$ day$^{\mathregular{-1}}$]', fontsize=14)
## r(T, O3)
#ax4.plot(r_t2mo3_ena_za, lat_ena, '-k', lw=2.)
#ax4.set_xlim([-0.25, 0.75])
#ax4.set_xticks(np.linspace(-0.25, 0.75, 5))
#ax4.set_xticklabels(['-0.25', '', '0.25', '', '0.75'], fontsize=12)
#ax4.set_title('$r$(T, O$_{\mathregular{3}}$)', fontsize=16, ha='center', y=1.02)
#ax4.set_xlabel('[$\mathregular{\cdot}$]', fontsize=14)
## dO3/dT
#ax5.plot(do3dt_ena_za, lat_ena, '-k', lw=2., clip_on=False)
#ax5.set_xlim([-0.25, 1.75])
#ax5.set_xticks(np.linspace(-0.25, 1.75, 5))
#ax5.set_xticklabels(['-0.25', '', '0.75', '', '1.75'], fontsize=12)
#ax5.set_title('dO$_{\mathregular{3}}$/dT$^{\mathregular{-1}}$', fontsize=16, ha='center', y=1.02)
#ax5.set_xlabel('[ppbv K$^{\mathregular{-1}}$]', fontsize=14)
## Aesthetics
#for ax in [ax1, ax2, ax3, ax4, ax5]:
#    ax.spines['right'].set_color(None)   
#    ax.tick_params(top=True, labeltop=False)
#    ax.set_ylim([20, 70])    
#for ax in [ax2, ax3, ax4, ax5]:
#    ax.spines['left'].set_color(None)
#    ax.yaxis.set_ticks_position('none')
#    ax.tick_params(labelleft=False)    
#ax1.set_yticks(np.linspace(20, 70, 11))
#ax1.set_yticklabels(['20', '', '30', '', '40', '', '50', '', '60',
#    '', '70'], fontsize=12)
#ax1.set_ylabel('Latitude [$^{\mathregular{\circ}}$N]', fontsize=14)
#plt.subplots_adjust(wspace=0.4)
#plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#            'zonalavg_u500_o3_nox_ro3t_do3dt_ena.eps', dpi=300)
#plt.show()
## Find fields over Eastern Asia
#latmin_ea, latmax_ea, lngmin_ea, lngmax_ea = 20., 70., 100., 130.
#o3_ea, lat_ea, lng_ea = globalo3_calculate.find_grid_in_bb(o3_n, lat_gmi_n, 
#    lng_gmi_n, lngmin_ea, lngmax_ea, latmin_ea, latmax_ea)
#do3dt_ea, lat_ea, lng_ea = globalo3_calculate.find_grid_in_bb(do3dt_n, 
#    lat_gmi_n, lng_gmi_n, lngmin_ea, lngmax_ea, latmin_ea, latmax_ea)
#r_t2mo3_ea, lat_ea, lng_ea = globalo3_calculate.find_grid_in_bb(r_t2mo3_n, 
#    lat_gmi_n, lng_gmi_n, lngmin_ea, lngmax_ea, latmin_ea, latmax_ea)
#U500_ea, lat_ea, lng_ea = globalo3_calculate.find_grid_in_bb(U500, 
#    lat_gmi_n, lng_gmi_n, lngmin_ea, lngmax_ea, latmin_ea, latmax_ea)
#nox_edgar_ea, lat_ea, lng_ea = globalo3_calculate.find_grid_in_bb(
#    nox_edgar_n, lat_gmi_n, lng_gmi_n, lngmin_ea, lngmax_ea, latmin_ea, 
#    latmax_ea)
#land_ea = globalo3_calculate.find_grid_overland(lat_ea, lng_ea)
## Find zonally-averaged fields
#o3_ea_za = np.nanmean(o3_ea*land_ea, axis=tuple((0,2)))
#do3dt_ea_za = np.nanmean(do3dt_ea*land_ea, axis=1)
#r_t2mo3_ea_za = np.nanmean(r_t2mo3_ea*land_ea, axis=1) 
#nox_edgar_ea_za = np.nanmean(nox_edgar_ea*land_ea, axis=tuple((0,2)))
#U500_ea_za = np.nanmean(U500_ea*land_ea, axis=tuple((0,2)))
## Plotting
#fig = plt.figure(figsize=(8,8))
#ax1 = plt.subplot2grid((5, 5), (0, 0), rowspan=5)
#ax2 = plt.subplot2grid((5, 5), (0, 1), rowspan=5)
#ax3 = plt.subplot2grid((5, 5), (0, 2), rowspan=5)
#ax4 = plt.subplot2grid((5, 5), (0, 3), rowspan=5)
#ax5 = plt.subplot2grid((5, 5), (0, 4), rowspan=5)
## U-wind at 500 hPa
#ax1.plot(U500_ea_za , lat_ea, '-k', lw=2.)
#ax1.set_xticks(np.linspace(-1, 9, 5))
#ax1.set_xticklabels(['-1', '', '4', '', '9'], fontsize=12)
#ax1.set_title('U$_{\mathregular{500\:hPa}}$', fontsize=16, ha='center', 
#    y=1.02)
#ax1.set_xlabel('[m s$^{\mathregular{-1}}$]', fontsize=14)
## Surface-level O3
#ax2.plot(o3_ea_za , lat_ea, '-k', lw=2.)
#ax2.set_xlim([20, 62])
#ax2.set_xticks(np.linspace(20, 62, 5))
#ax2.set_xticklabels(['20', '', '41', '', '62'], fontsize=12)
#ax2.set_title('O$_{\mathregular{3}}$', fontsize=16, ha='center', y=1.02)
#ax2.set_xlabel('[ppbv]', fontsize=14, labelpad=7)
## EDGAR NOx emissions (n.b., multiply by no. sec/day and 1000 m/km to 
## convert from /s/m2 to /day/km2)
#ax3.plot(nox_edgar_ea_za*86400.*1000*1000., lat_ea, '-k', lw=2., 
#         clip_on=False)
#ax3.set_xlim([0, 132])
#ax3.set_xticks(np.linspace(0, 132, 5))
#ax3.set_xticklabels(['0', '', '66', '', '132'], fontsize=12)
#ax3.set_title('NO$_{{x}}$', fontsize=16, ha='center', y=1.02)
#ax3.set_xlabel('[kg km$^{\mathregular{-2}}$ day$^{\mathregular{-1}}$]', 
#    fontsize=14)
## r(T, O3)
#ax4.plot(r_t2mo3_ea_za, lat_ea, '-k', lw=2., clip_on=False)
#ax4.set_xlim([-0.1, 0.7])
#ax4.set_xticks(np.linspace(-0.1, 0.7, 5))
#ax4.set_xticklabels(['-0.1', '', '0.3', '', '0.7'], fontsize=12)
#ax4.set_title('$r$(T, O$_{\mathregular{3}}$)', fontsize=16, ha='center', y=1.02)
#ax4.set_xlabel('[$\mathregular{\cdot}$]', fontsize=14)
## dO3/dT
#ax5.plot(do3dt_ea_za, lat_ea, '-k', lw=2.)
#ax5.set_xlim([-0.5, 1.6])
#ax5.set_xticks(np.linspace(-0.5, 1.6, 5))
#ax5.set_xticklabels(['-0.5', '', '0.55', '', '1.6'], fontsize=12)
#ax5.set_title('dO$_{\mathregular{3}}$/dT$^{\mathregular{-1}}$', fontsize=16, 
#    ha='center', y=1.02)
#ax5.set_xlabel('[ppbv K$^{\mathregular{-1}}$]', fontsize=14)
## Aesthetics
#for ax in [ax1, ax2, ax3, ax4, ax5]:
#    ax.spines['right'].set_color(None)   
#    ax.tick_params(top=True, labeltop=False)
#    ax.set_ylim([20, 70])    
#for ax in [ax2, ax3, ax4, ax5]:
#    ax.spines['left'].set_color(None)
#    ax.yaxis.set_ticks_position('none')
#    ax.tick_params(labelleft=False)    
#ax1.set_yticks(np.linspace(20, 70, 11))
#ax1.set_yticklabels(['20', '', '30', '', '40', '', '50', '', '60',
#    '', '70'], fontsize=12)
#ax1.set_ylabel('Latitude [$^{\mathregular{\circ}}$N]', fontsize=14)
#plt.subplots_adjust(wspace=0.4)
#plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#            'zonalavg_u500_o3_nox_ro3t_do3dt_ea.eps', dpi=300)
#plt.show()


"""CALCULATE THE STANDARD DEVIATION OF 2-METER TEMPERATURE WITH LATITUDE
   OVER THE EASTERN U.S. AFTER BARNES AND FIORE [2013] """
## Find 2-meter temperatures over Eastern North America 
#t2m_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(t2m_n, 
#    lat_gmi_n, lng_gmi_n, 269., 294., 39., 58.)
#land_ena = globalo3_calculate.find_grid_overland(lat_ena, lng_ena)
#lng_ena_m, lat_ena_m = np.meshgrid(lng_ena, lat_ena)
#lat_ena_f = (lat_ena_m*land_ena).flatten()
#fig = plt.figure()
#ax = plt.subplot2grid((1,1),(0,0))
#ax.plot(lat_ena_f, (np.std(t2m_ena,axis=0)*land_ena).flatten(), 'ko')
#ax.set_xlabel('Latitude [$^{\circ}$]')
#ax.set_ylabel('$\mathregular{\sigma}_{\mathregular{T}}$ [K]')
#plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#            'scatter_latitude_t2m_easternnorthamerica.eps', dpi=300)

"""PLOT GLOBAL TRENDS"""
#map_nh(lat_gmi_n, lat_gmi_s, lng_gmi_n, lng_gmi_s, nox_ls_n, 
#    nox_ls_s, 'L-S Trend NO$_{x}$ [ppbv yr$^{\mathregular{-1}}$]', 
#    np.linspace(-.05, .05, 9), 'bwr', 'lstrendnox', 
#    p_n = nox_lsp_n, p_s = nox_lsp_s)
#map_nh(lat_gmi_n, lat_gmi_s, lng_gmi_n, lng_gmi_s, nox_mk_n, 
#    nox_mk_s, 'M-K Trend NO$_{x}$ [$\cdot$]', 
#    np.linspace(-3.5, 3.5, 9), 'bwr', 'mktrendnox', p_n = nox_mkp_n, 
#    p_s = nox_mkp_s)
#map_nh(lat_gmi_n, lat_gmi_s, lng_gmi_n, lng_gmi_s, do3dt_ls_n, 
#    do3dt_ls_s, 'dO$_{\mathregular{3}}$/dT [ppbv '+
#    'K$^{\mathregular{-1}}$ yr$^{\mathregular{-1}}$]',
#    np.linspace(-.2, .2, 9), 'bwr', 'lstrenddo3dt', p_n = do3dt_lsp_n,
#    p_s = do3dt_lsp_s)
#map_nh(lat_gmi_n, lat_gmi_s, lng_gmi_n, lng_gmi_s, do3dt_mk_n, 
#    do3dt_mk_s, 'M-K Trend dO$_{\mathregular{3}}$/dT [$\cdot$]',
#    np.linspace(-2., 2., 9), 'bwr', 'mktrenddo3dt', p_n = do3dt_mkp_n,
#    p_s = do3dt_mkp_s)

"""CHECK STD/EMFIX FIELDS FROM STRODE ET AL. (2015)"""
#stdlat_gmi_n, stdlng_gmi_n, times_n, stdo3_n = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_n, 
#    latmax_n, lngmin_n, lngmax_n, 'O3', 'HindcastFFIgac2')
#stdlat_gmi_s, stdlng_gmi_s, times_s, stdo3_s = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_s, latmin_s, 
#    latmax_s, lngmin_s, lngmax_s, 'O3', 'HindcastFFIgac2')
#stdlat_gmi_n, stdlng_gmi_n, times_n, stdno_n = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_n, 
#    latmax_n, lngmin_n, lngmax_n, 'NO', 'HindcastFFIgac2')
#stdlat_gmi_s, stdlng_gmi_s, times_s, stdno_s = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_s, latmin_s, 
#    latmax_s, lngmin_s, lngmax_s, 'NO', 'HindcastFFIgac2')
#stdlat_gmi_n, stdlng_gmi_n, times_n, stdno2_n = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_n, 
#    latmax_n, lngmin_n, lngmax_n, 'NO2', 'HindcastFFIgac2')
#stdlat_gmi_s, stdlng_gmi_s, times_s, stdno2_s = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_s, latmin_s, 
#    latmax_s, lngmin_s, lngmax_s, 'NO2', 'HindcastFFIgac2')
## Load EmFix O3 from Strode et al. (2015)
#emfixlat_gmi_n, emfixlng_gmi_n, times_n, emfixo3_n = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_n, 
#    latmax_n, lngmin_n, lngmax_n, 'O3', 'Hindcast3Igac2')
#emfixlat_gmi_s, emfixlng_gmi_s, times_s, emfixo3_s = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_s, latmin_s, 
#    latmax_s, lngmin_s, lngmax_s, 'O3', 'Hindcast3Igac2')
#emfixlat_gmi_n, emfixlng_gmi_n, times_n, emfixno_n = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_n, 
#    latmax_n, lngmin_n, lngmax_n, 'NO', 'Hindcast3Igac2')
#emfixlat_gmi_s, emfixlng_gmi_s, times_s, emfixno_s = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_s, latmin_s, 
#    latmax_s, lngmin_s, lngmax_s, 'NO', 'Hindcast3Igac2')
#emfixlat_gmi_n, emfixlng_gmi_n, times_n, emfixno2_n = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_n, 
#    latmax_n, lngmin_n, lngmax_n, 'NO2', 'Hindcast3Igac2')
#emfixlat_gmi_s, emfixlng_gmi_s, times_s, emfixno2_s = \
#    globalo3_open.open_overpass2_specifieddomain(years, months_s, latmin_s, 
#    latmax_s, lngmin_s, lngmax_s, 'NO2', 'Hindcast3Igac2')    
## Convert trace gases from volume mixing ratio to ppbv
#stdo3_n = stdo3_n*1e9
#stdo3_s = stdo3_s*1e9
#stdno_n = stdno_n*1e9
#stdno_s = stdno_s*1e9
#stdno2_n = stdno2_n*1e9
#stdno2_s = stdno2_s*1e9
#emfixo3_n = emfixo3_n*1e9
#emfixo3_s = emfixo3_s*1e9
#emfixno_n = emfixno_n*1e9
#emfixno_s = emfixno_s*1e9
#emfixno2_n = emfixno2_n*1e9
#emfixno2_s = emfixno2_s*1e9
## Load MERRA-2 T2m
#lat_merra_n, lng_merra_n, stdt2m_n = \
#    globalo3_open.open_merra2t2m_specifieddomain(years, months_n, latmin_n, 
#    latmax_n, lngmin_n, lngmax_n)
#lat_merra_s, lng_merra_s, stdt2m_s = \
#    globalo3_open.open_merra2t2m_specifieddomain(years, months_s, latmin_s, 
#    latmax_s, lngmin_s, lngmax_s)
## Interpolate T2m
#stdt2m_n = globalo3_open.interpolate_merra_to_ctmresolution(stdlat_gmi_n, 
#    stdlng_gmi_n, lat_merra_n, lng_merra_n, stdt2m_n)
#stdt2m_s = globalo3_open.interpolate_merra_to_ctmresolution(stdlat_gmi_s, 
#    stdlng_gmi_s, lat_merra_s, lng_merra_s, stdt2m_s)
## Calculate dO3/dT
#stddo3dt_n = globalo3_calculate.calculate_do3dt(stdt2m_n, stdo3_n, 
#    stdlat_gmi_n, stdlng_gmi_n)    
#stddo3dt_s = globalo3_calculate.calculate_do3dt(stdt2m_s, stdo3_s, 
#    stdlat_gmi_s, stdlng_gmi_s)
#emfixdo3dt_n = globalo3_calculate.calculate_do3dt(stdt2m_n, emfixo3_n, 
#    stdlat_gmi_n, stdlng_gmi_n)    
#emfixdo3dt_s = globalo3_calculate.calculate_do3dt(stdt2m_s, emfixo3_s, 
#    stdlat_gmi_s, stdlng_gmi_s)
## Calculate least-squares/Mann-Kendall trends
#(stddo3dt_byyr_n, stddo3dt_ls_n, stddo3dt_lsp_n, stddo3dt_mkz_n, 
#    stddo3dt_mkp_n) = globalo3_calculate.calculate_trends_do3dt(stdt2m_n, 
#    stdo3_n, stdlat_gmi_n, stdlng_gmi_n, years)
#(stddo3dt_byyr_s, stddo3dt_ls_s, stddo3dt_lsp_s, stddo3dt_mkz_s, 
#    stddo3dt_mkp_s) = globalo3_calculate.calculate_trends_do3dt(stdt2m_s, 
#    stdo3_s, stdlat_gmi_s, stdlng_gmi_s, years)
#(emfixdo3dt_byyr_n, emfixdo3dt_ls_n, emfixdo3dt_lsp_n, emfixdo3dt_mkz_n, 
#    emfixdo3dt_mkp_n) = globalo3_calculate.calculate_trends_do3dt(stdt2m_n, 
#    emfixo3_n, stdlat_gmi_n, stdlng_gmi_n, years)
#(emfixdo3dt_byyr_s, emfixdo3dt_ls_s, emfixdo3dt_lsp_s, emfixdo3dt_mkz_s, 
#    emfixdo3dt_mkp_s) = globalo3_calculate.calculate_trends_do3dt(stdt2m_s, 
#    emfixo3_s, stdlat_gmi_s, stdlng_gmi_s, years)
## Plot mean global fields from Std and Emfix simulation and differences
## for O3
#map_global(stdlat_gmi_n, stdlat_gmi_s, stdlng_gmi_n, stdlng_gmi_s, 
#    np.mean(stdo3_n, axis=0), np.mean(stdo3_s, axis=0), 
#    '<O$_{\mathregular{3, Std}}$> [ppbv]', np.linspace(25, 65, 11), 'PuBu', 
#    'meano3_std') 
#map_global(stdlat_gmi_n, stdlat_gmi_s, stdlng_gmi_n, stdlng_gmi_s, 
#    np.mean(emfixo3_n, axis=0), np.mean(emfixo3_s, axis=0), 
#    '<O$_{\mathregular{3, EmFix}}$> [ppbv]', np.linspace(25, 65, 11), 'PuBu', 
#    'meano3_emfix')
#map_global(stdlat_gmi_n, stdlat_gmi_s, stdlng_gmi_n, stdlng_gmi_s, 
#    np.mean(emfixo3_n-stdo3_n, axis=0), np.mean(emfixo3_s-stdo3_s, axis=0), 
#    '<O$_{\mathregular{3}}$>$_{\mathregular{EmFix - Std}}$ [ppbv]', 
#    np.linspace(-8, 8, 9), 'bwr', 'meano3_emfixstddiff')
## For NOx
#map_global(stdlat_gmi_n, stdlat_gmi_s, stdlng_gmi_n, stdlng_gmi_s, 
#    np.mean((stdno_n+stdno2_n), axis=0), np.mean((stdno_s+stdno2_s), axis=0), 
#    '<NO$_{x\mathregular{, Std}}$> [ppbv]', np.linspace(0.5, 1.5, 11), 
#    'PuBu', 'meannox_std')
#map_global(stdlat_gmi_n, stdlat_gmi_s, stdlng_gmi_n, stdlng_gmi_s, 
#    np.mean((emfixno_n+emfixno2_n), axis=0), np.mean((emfixno_s+emfixno2_s), 
#    axis=0), '<NO$_{x\mathregular{, EmFix}}$> [ppbv]', 
#    np.linspace(0.5, 1.5, 11), 'PuBu', 'meannox_emfix')
#map_global(stdlat_gmi_n, stdlat_gmi_s, stdlng_gmi_n, stdlng_gmi_s, 
#    np.mean((emfixno_n+emfixno2_n)-(stdno_n+stdno2_n), axis=0), 
#    np.mean((emfixno_s+emfixno2_s)-(stdno_s+stdno2_s), 
#    axis=0), '<NO$_{x}$>$_{\mathregular{EmFix - Std}}$ [ppbv]', 
#    np.linspace(-0.5, 0.5, 11), 'bwr', 'meannox_emfixstddiff')
## For dO3/dT
#map_global(stdlat_gmi_n, stdlat_gmi_s, stdlng_gmi_n, stdlng_gmi_s, stddo3dt_n, 
#    stddo3dt_s,'<dO$_{\mathregular{3, Std}}$/dT> [ppbv K$^{\mathregular{-1}}$]', 
#    np.linspace(0, 3., 7), 'PuBu', 'meando3dt_std')
#map_global(stdlat_gmi_n, stdlat_gmi_s, stdlng_gmi_n, stdlng_gmi_s, 
#    emfixdo3dt_n, stddo3dt_s,'<dO$_{\mathregular{3, EmFix}}$/dT> '+
#    '[ppbv K$^{\mathregular{-1}}$]', np.linspace(0, 3., 7), 'PuBu', 
#    'meando3dt_emfix')    
#map_global(stdlat_gmi_n, stdlat_gmi_s, stdlng_gmi_n, stdlng_gmi_s, 
#    (emfixdo3dt_n-stddo3dt_n), (emfixdo3dt_s-stddo3dt_s),
#    '<dO$_{\mathregular{3}}$/dT>$_{\mathregular{EmFix - Std}}$ '+
#    '[ppbv K$^{\mathregular{-1}}$]', np.linspace(-1, 1., 9), 'bwr', 
#    'meando3dt_emfixstddiff')  
## Plot regionally-averaged trends
## Eastern U.S. 
#timeseries_rado3dt(stdt2m_n, stdo3_n, emfixo3_n, (stdno_n+stdno2_n), 
#    (emfixno_n+emfixno2_n), years, stdlat_gmi_n, stdlng_gmi_n, 270., 49., 290., 
#    32., 'eastus')
## China
#timeseries_rado3dt(stdt2m_n, stdo3_n, emfixo3_n, (stdno_n+stdno2_n), 
#    (emfixno_n+emfixno2_n), years, stdlat_gmi_n, stdlng_gmi_n, 110., 28., 118., 
#    22., 'china')