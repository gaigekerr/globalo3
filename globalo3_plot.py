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
    23042019 -- function 'timeseries_o3atabovebelow_jet' added
    26052019 -- function 'map_jet_centerdist' added
    17062019 -- function 'edjetlocation_fieldatedjet' added
    18062019 -- contour/label capability added to 'map_nh' 
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

def map_nh(lat_n, lng_n, field_n, title, cbar_label, clevs, cmap, fstr, 
    contour=None, contour_levs=None, p_n=None, e_n=None, eerr_n=None, 
    alpha=0.05, extent=None, quiver=None, oceanon='yes', extend='both'):
    """plot desired quantity over the Northern Hemisphere (however, if variable
    'extent' is specified, map can be tailored to any region).
    
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
        ax.set_extent([-180., 180., 0., 85.])    
    else: 
        ax.set_extent(extent)
    # Plot filled contours
    cmap = plt.get_cmap(cmap)
    mb = ax.contourf(lng_n, lat_n, field_n, clevs, cmap=cmap,
                     transform=ccrs.PlateCarree(), extend=extend)
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
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label(cbar_label, fontsize=14)
    plt.gcf().subplots_adjust(right=0.75, hspace=-0.2)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'map_nh_%s.eps'%fstr, transparent=True)
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
    
def map_extent(f, lat, lng, frl, fru, frr, frd, cbar_label, clevs, cmap, fstr):
    """function plots a filled contour map (Mercator projection) of the 
    specified field over a specified domain. 
    
    Parameters
    ----------  
    f : numpy.ndarray        
        Field of interest that is plotted on map, [lat, lng]        
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
    cbar_label : str
        Label for the colorbar (field and units)
    clevs : numpy.ndarray
        Contour levels values
    cmap : str
        Colormap name
    fstr : str
        Output filename suffix

    Returns
    -------
    None          
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')    
    from geo_idx import geo_idx
    # Set the extent (x0, x1, y0, y1) of the map in the given coordinate system    
    extent = [frl, frr, frd, fru-1]
    # Define focus region (i.e., Eastern U.S., Europe, etc.)
    frl = geo_idx(frl, lng)
    frr = geo_idx(frr, lng)
    fru = geo_idx(fru, lat)
    frd = geo_idx(frd, lat)
    # Find field in region 
    frf = f[frd:fru+1, frl:frr+1]
    lng = lng[frl:frr+1]
    lat = lat[frd:fru+1]
    # Plotting
    ax = plt.axes(projection=ccrs.Mercator())
    ax.set_extent(extent)
    ax.coastlines(resolution='50m')
    ax.add_feature(cfeature.OCEAN, zorder=10, lw = 0.0, color='lightgrey')
    ax.add_feature(cfeature.COASTLINE, zorder=10, lw = 1.)
    ax.add_feature(cfeature.BORDERS, zorder=10, lw = 1.)
    mb = ax.contourf(lng, lat, frf, clevs, cmap=plt.get_cmap(cmap), 
        transform=ccrs.PlateCarree(), extend='both')
    # Add colorbar    
    colorbar_axes = plt.gcf().add_axes([0.78,0.15,0.03,0.7])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical')
    colorbar.set_label('%s' %cbar_label)
    plt.gcf().subplots_adjust(right=0.75)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+'map_extent_%s.eps'%(fstr))
    plt.show()
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
            # Find the latitud of the jet at longitude/on day of interest
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

def edjetlocation_fieldatedjet(lat_fr, lng_fr, lat_jet, field_jet, cmap, 
    cbar_label, clevs, fieldstr, fstr, skiplng=6):
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
    mb = axb.contourf(lng_fr, np.arange(-10, 11, 1), field_jet_roll, clevs, 
                      cmap=plt.get_cmap(cmap), extend='both')
    axb.set_xlabel('Longitude [$^{\circ}$]', fontsize=14)
    axb.set_xticks(np.linspace(0, 360, 9))
    axb.set_xticklabels([-180,-135,-90,-45,0,45,90,135,180], fontsize=12)
    axb.set_ylabel('Grid cells from jet', fontsize=14)
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
    plt.gcf().subplots_adjust(right=0.75, hspace=-0.1)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'edjetlocation_%satedjet_%s.eps'%(fieldstr, fstr), dpi=300)
    return 

import numpy as np
from datetime import datetime
import sys
sys.path.append('/Users/ghkerr/phd/globalo3/')
import globalo3_open, globalo3_calculate
sys.path.append('/Users/ghkerr/phd/transporto3/')
import transporto3_open
years = [2008, 2009, 2010]
# Define hemispheres and their O3 seasons
latmin_n, lngmin_n, latmax_n, lngmax_n = -1., 0., 90., 360.
months_n = ['jun', 'jul', 'aug']
# Load data the first iteration only 
try:o3_n
except NameError:
    # Load trace gases from + Chemistry and Transport simulations
    lat_gmi_n, lng_gmi_n, times_n, o3_n = \
        globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_n, 
        latmax_n, lngmin_n, lngmax_n, 'O3', 'HindcastMR2')
    lat_gmi_n, lng_gmi_n, times_n, o3_dat_n = \
        globalo3_open.open_overpass2_specifieddomain(years, months_n, latmin_n, 
        latmax_n, lngmin_n, lngmax_n, 'O3', 'HindcastMR2-DiurnalAvgT')
    # Convert trace gases from volume mixing ratio to ppbv
    o3_n = o3_n*1e9
    o3_dat_n = o3_dat_n*1e9  
    # Load TOAR MDA8 O3
    urban_mean, toartime, toarlat, toarlng = globalo3_open.open_toar(years, 
       months_n, 'urban_mean', 2)
    rural_mean, toartime, toarlat, toarlng = globalo3_open.open_toar(years, 
       months_n, 'rural_mean', 2)
    # Load MERRA-2 T2m
    lat_merra_n, lng_merra_n, t2m_n = globalo3_open.open_merra2t2m_specifieddomain(
        years, months_n, latmin_n, latmax_n, lngmin_n, lngmax_n)
    # Interpolate T2m
    t2m_n = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi_n, lng_gmi_n, 
        lat_merra_n, lng_merra_n, t2m_n)
#    # Calculate relative humidity
#    rh_n = globalo3_calculate.calculate_rh_from_q(years, 
#        list(np.arange(0,24,1)), lngmin_n, latmax_n, lngmax_n, latmin_n)
    # Load MERRA-2 data at 500 hPa
    years = [2008,2009,2010]
    hours = [0,3,6,9,12,15,18,21]
    U500, mtime, lat_merra, lng_merra = transporto3_open.open_merra2(years, 
        hours, 'U', 'inst3_3d_asm_Np', 'JJA_500mb.nc', lngmin_n, latmax_n, 
        lngmax_n, latmin_n, dailyavg='yes')
    V500, mtime, lat_merra, lng_merra = transporto3_open.open_merra2(years, 
        hours, 'V', 'inst3_3d_asm_Np', 'JJA_500mb.nc', lngmin_n, latmax_n, 
        lngmax_n, latmin_n, dailyavg='yes')
    QV500, mtime, lat_merra, lng_merra = transporto3_open.open_merra2(years, 
        hours, 'QV', 'inst3_3d_asm_Np', 'JJA_500mb.nc', lngmin_n, latmax_n, 
        lngmax_n, latmin_n, dailyavg='yes')
    # Interpolate MERRA-2 data to resolution of CTM
    U500 = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi_n, 
        lng_gmi_n, lat_merra, lng_merra, U500, checkplot='yes')
    V500 = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi_n, 
        lng_gmi_n, lat_merra, lng_merra, V500, checkplot='yes')
    QV500 = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi_n, 
        lng_gmi_n, lat_merra, lng_merra, QV500, checkplot='yes')    
#    # Calculate wind direction at 500 hPa
#    DIR500 = globalo3_calculate.convert_uv_tocardinal(U500, V500) 
    # Load EDGAR inventory     
    lat_edgar_n, lng_edgar_n, nox_edgar_n = \
        globalo3_open.open_edgar_specifieddomain(years, latmin_n, latmax_n, 
        lngmin_n, lngmax_n, 'NOx')
    # Interpolate EDGAR 
    nox_edgar_n = globalo3_open.interpolate_edgar_to_ctmresolution(lat_gmi_n, 
        lng_gmi_n, lat_edgar_n, lng_edgar_n, nox_edgar_n)
    ## Load Schnell et al. (2014) O3 dataset
    #lat_js, lng_js, o3_js = globalo3_open.open_schnello3(years, months_n, 'US')
    # Calculate dO3/dT
    do3dt_n = globalo3_calculate.calculate_do3dt(t2m_n, o3_n, lat_gmi_n, 
        lng_gmi_n)
    do3dt_dat_n = globalo3_calculate.calculate_do3dt(t2m_n, o3_dat_n, 
        lat_gmi_n, lng_gmi_n)
#    # Calculate dO3/dRH
#    do3drh_n = globalo3_calculate.calculate_do3dt(rh_n*100., o3_n, lat_gmi_n, 
#        lng_gmi_n)
    # Calculate correlation coefficients
    r_t2mo3_n = globalo3_calculate.calculate_r(t2m_n, o3_n, lat_gmi_n, 
        lng_gmi_n)
    r_t2mo3_n_dat = globalo3_calculate.calculate_r(t2m_n, o3_dat_n, lat_gmi_n, 
        lng_gmi_n)    
#    r_rho3_n = globalo3_calculate.calculate_r(rh_n*100., o3_n, lat_gmi_n, 
#        lng_gmi_n)
#    r_t2mrh_n = globalo3_calculate.calculate_r(rh_n, t2m_n, lat_gmi_n, 
#        lng_gmi_n)
#    r_o3u500_n = globalo3_calculate.calculate_r(U500, o3_dat_n, lat_gmi_n, 
#        lng_gmi_n)
#    r_dirrh_n = globalo3_calculate.calculate_r(DIR500, rh_n, lat_gmi_n, 
#        lng_gmi_n)
#    r_dirq500_n = globalo3_calculate.calculate_r(DIR500, QV500, lat_gmi_n, 
#        lng_gmi_n)
#    r_diro3_n = globalo3_calculate.calculate_r(DIR500, o3_n, lat_gmi_n, 
#        lng_gmi_n)    
#    # Calculate trends
#    do3dt_byyr_n, do3dt_ls_n, do3dt_lsp_n, do3dt_mk_n, do3dt_mkp_n = \
#    globalo3_calculate.calculate_trends_do3dt(t2m_n, o3_n, lat_gmi_n, 
#        lng_gmi_n, years)
#    edgar_byyr_n, edgar_ls_n, edgar_lsp_n, edgar_mk_n, edgar_mkp_n = \
#    globalo3_calculate.calculate_trends(nox_edgar_n, lat_gmi_n, lng_gmi_n, 
#        years)
    # Calculate eddy-driven jet position in the mid-latitudes
    U500_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(U500, 
        lat_gmi_n, lng_gmi_n, 0., 360., 23., 60.)
    o3_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(o3_n, 
        lat_gmi_n, lng_gmi_n, 0., 360., 23., 60.)
    #do3dt_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(do3dt_n, 
    #    lat_gmi_n, lng_gmi_n, 0., 360., 23., 60.)   
    #r_t2mo3_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(
    #    r_t2mo3_n, lat_gmi_n, lng_gmi_n, 0., 360., 23., 60.)      
    #land_nhml = globalo3_calculate.find_grid_overland(lat_nhml, lng_nhml)   
    lat_jet_nhml, o3_jet_nhml = globalo3_calculate.find_field_atjet(o3_nhml, 
        U500_nhml, lat_nhml, lng_nhml, 10, anom=True)
    m_o3jetdist, r_o3jetdist = \
        globalo3_calculate.calculate_o3jet_relationship(o3_n, lat_gmi_n,
        lng_gmi_n, lat_jet_nhml, lng_nhml)

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
#land_naml = globalo3_calculate.find_grid_overland(lat_naml, lng_naml)
#lat_jet_naml, o3_jet_naml = globalo3_calculate.find_field_atjet(o3_naml, 
#    U500_naml, lat_naml, lng_naml, 10, anom=True)
## Find the daily mean jet latitude averaged over focus region
#jet_lat_mean = np.nanmean(lat_jet_naml, axis=1)
## Days where the mean jet position is 2 sigma above/below average
#where_high2sigma = np.where(jet_lat_mean > 
#    jet_lat_mean.mean()+(1.25*np.std(jet_lat_mean)))[0]
#where_low2sigma = np.where(jet_lat_mean < 
#    jet_lat_mean.mean()-(1.25*np.std(jet_lat_mean)))[0]
## Jet and mean O3 anomaly for the jet when poleward than 2 sigma its mean
## latitude
#map_nh(lat_naml, 
#    lng_naml, 
#    np.nanmean(o3_naml[where_high2sigma], axis=0)-np.nanmean(o3_naml, axis=0), 
#    '$\mathregular{\Delta}$ O$_\mathregular{3}$ for $\mathregular{\phi}_'+
#    '{\mathregular{jet}}$>($\mathregular{\mu}_{\mathregular{\phi_{jet}}}$+'+
#    '1.25$\mathregular{\sigma}_{\mathregular{\phi_{jet}}}$)',    
#    '[ppbv]', 
#    np.linspace(-10, 10, 11), 
#    'bwr', 
#    'northamerica_anomo3_jetlat_high2sigma',
#    e_n=np.nanmean(lat_jet_naml[where_high2sigma],axis=0), 
#    eerr_n=np.nanstd(lat_jet_naml[where_high2sigma],axis=0), 
#    extent=[lng_naml.min(), lng_naml.max(), lat_naml.min(), lat_naml.max()-5], 
#    oceanon='yes', 
#    extend='both')
## Jet and mean O3 anomaly for the jet when equatorward than 2 sigma its mean
## latitude
#map_nh(lat_naml, 
#    lng_naml, 
#    np.nanmean(o3_naml[where_low2sigma], axis=0)-np.nanmean(o3_naml, axis=0), 
#    '$\mathregular{\Delta}$ O$_\mathregular{3}$ for $\mathregular{\phi}_'+
#    '{\mathregular{jet}}$<($\mathregular{\mu}_{\mathregular{\phi_{jet}}}$-'+
#    '1.25$\mathregular{\sigma}_{\mathregular{\phi_{jet}}}$)',
#    '[ppbv]', 
#    np.linspace(-10, 10, 11), 
#    'bwr', 
#    'northamerica_anomo3_jetlat_low2sigma',
#    e_n=np.nanmean(lat_jet_naml[where_low2sigma],axis=0), 
#    eerr_n=np.nanstd(lat_jet_naml[where_low2sigma],axis=0), 
#    extent=[lng_naml.min(), lng_naml.max(), lat_naml.min(), lat_naml.max()-5], 
#    oceanon='yes', 
#    extend='both')    
## Maps for individual days/case studies with extreme poleward and equatorward
## jet streams positions
#map_nh(lat_naml, 
#    lng_naml, 
#    o3_naml[where_high2sigma[1]]-np.nanmean(o3_naml, axis=0), 
#    '$\mathregular{\Delta}$ O$_\mathregular{3}$ for %s' 
#    %datetime.strftime(mtime[where_high2sigma[1]], '%m/%d/%Y'),
#    '[ppbv]', 
#    np.linspace(-10, 10, 11), 
#    'bwr',
#    'northamerica_anomo3_jetlat-%s'
#    %datetime.strftime(mtime[where_high2sigma[1]], '%m-%d-%Y'),
#    e_n=jet_lat[where_high2sigma[1]],
#    eerr_n=np.zeros(shape=jet_lat[where_high2sigma[1]].shape), 
#    extent=[lng_fr.min(), lng_fr.max(), lat_fr.min(), lat_fr.max()-5], 
#    oceanon='yes', 
#    extend='both') 
#map_nh(lat_naml, 
#    lng_naml, 
#    o3_naml[where_high2sigma[8]]-np.nanmean(o3_naml, axis=0), 
#    '$\mathregular{\Delta}$ O$_\mathregular{3}$ for %s' 
#    %datetime.strftime(mtime[where_high2sigma[8]], '%m/%d/%Y'),    
#    '[ppbv]', 
#    np.linspace(-10, 10, 11), 
#    'bwr', 
#    'northamerica_anomo3_jetlat-%s'
#    %datetime.strftime(mtime[where_high2sigma[8]], '%m-%d-%Y'),
#    e_n=jet_lat[where_high2sigma[8]],
#    eerr_n=np.zeros(shape=jet_lat[where_high2sigma[8]].shape), 
#    extent=[lng_fr.min(), lng_fr.max(), lat_fr.min(), lat_fr.max()-5], 
#    oceanon='yes', 
#    extend='both')     
#map_nh(lat_naml, 
#    lng_naml, 
#    o3_naml[where_low2sigma[0]]-np.nanmean(o3_naml, axis=0), 
#    '$\mathregular{\Delta}$ O$_\mathregular{3}$ for %s' 
#    %datetime.strftime(mtime[where_low2sigma[0]], '%m/%d/%Y'),
#    '[ppbv]', 
#    np.linspace(-10, 10, 11), 
#    'bwr', 
#    'northamerica_anomo3_jetlat-%s'
#    %datetime.strftime(mtime[where_low2sigma[0]], '%m-%d-%Y'),
#    e_n=jet_lat[where_low2sigma[0]],
#    eerr_n=np.zeros(shape=jet_lat[where_low2sigma[0]].shape), 
#    extent=[lng_fr.min(), lng_fr.max(), lat_fr.min(), lat_fr.max()-5], 
#    oceanon='yes', 
#    extend='both')     
#map_nh(lat_naml, 
#    lng_naml, 
#    o3_naml[where_low2sigma[3]]-np.nanmean(o3_naml, axis=0), 
#    '$\mathregular{\Delta}$ O$_\mathregular{3}$ for %s' 
#    %datetime.strftime(mtime[where_low2sigma[3]], '%m/%d/%Y'),
#    '[ppbv]',
#    np.linspace(-10, 10, 11), 
#    'bwr', 
#    'northamerica_anomo3_jetlat-%s'
#    %datetime.strftime(mtime[where_low2sigma[3]], '%m-%d-%Y'),
#    e_n=jet_lat[where_low2sigma[3]],
#    eerr_n=np.zeros(shape=jet_lat[where_low2sigma[3]].shape), 
#    extent=[lng_fr.min(), lng_fr.max(), lat_fr.min(), lat_fr.max()-5], 
#    oceanon='yes', 
#    extend='both')    
    
"""PLOT MEAN GLOBAL FIELDS AND CORRELATIONS"""
## Mean O3
#map_nh(lat_gmi_n, lng_gmi_n, np.mean(o3_n, axis=0), '', 
#       'O$_{\mathregular{3}}$ [ppbv]', np.linspace(25, 65, 11), 'PuBu', 
#       'meano3chem_%d-%d_jet'%(years[0],years[-1]), 
#       e_n=np.nanmean(lat_jet_nhml,axis=0), 
#       eerr_n=np.std(lat_jet_nhml,axis=0))
## Mean EDGAR NOx emissions
#map_nh(lat_gmi_n, lng_gmi_n, np.mean(nox_edgar_n, axis=0), '', 
#       'NO$_{x\:\mathregular{, EDGAR}}$ [kg m$^{\mathregular{-2}}$ '+
#       's$^{\mathregular{-1}}$]', np.linspace(0e-10, 1e-10, 11), 'PuBu', 
#       'meanedgarnox_%d-%d_jet'%(years[0],years[-1]), 
#       e_n=np.nanmean(lat_jet_nhml, axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
## dO3/dT from +Chemistry
#map_nh(lat_gmi_n, lng_gmi_n, do3dt_n, '', 
#       'dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]', 
#       np.linspace(0, 3., 7), 'Reds', 'do3dtchem_%d-%d_jet'%(years[0],years[-1]),
#       e_n=np.nanmean(lat_jet_nhml, axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
## r(T, O3) from +Chemistry
#map_nh(lat_gmi_n, lng_gmi_n, r_t2mo3_n, '', 
#    r'$\it{r}\:$(T, O$_\mathregular{3}$)', np.linspace(-1., 1., 11), 'bwr', 
#    'r_t2mo3chem_%d-%d_jet'%(years[0],years[-1]), 
#    e_n=np.nanmean(lat_jet_nhml, axis=0), 
#    eerr_n=np.zeros(lat_jet_nhml.shape[1]))
## dO3/dT from Transport
#map_nh(lat_gmi_n, lng_gmi_n, do3dt_dat_n, '', 
#       'dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]', 
#       np.linspace(0, 3., 7), 'Reds', 'do3dtdat_%d-%d_jet'%(years[0],
#       years[-1]), e_n=np.nanmean(lat_jet_nhml, axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
## r(T, O3) from Transport
#map_nh(lat_gmi_n, lng_gmi_n, r_t2mo3_n_dat, '', 
#    r'$\it{r}\:$(T, O$_\mathregular{3}$)', np.linspace(-1., 1., 11), 'bwr', 
#    'r_t2mo3dat_%d-%d_jet'%(years[0],years[-1]), 
#    e_n=np.nanmean(lat_jet_nhml, axis=0), 
#    eerr_n=np.zeros(lat_jet_nhml.shape[1]))
## O3 variability from +Chemistry 
#map_nh(lat_gmi_n, lng_gmi_n, np.nanstd(o3_n, axis=0), '', 
#       '$\mathregular{\sigma}_{\mathregular{O}_{\mathregular{3}}}$ [ppbv]', 
#       np.linspace(5, 11, 7), 'PuBu', 'sigmao3chem_%d-%d_jet'%(years[0],
#       years[-1]), e_n=np.nanmean(lat_jet_nhml, axis=0),
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
## Temperature variability 
#map_nh(lat_gmi_n, lng_gmi_n, np.nanstd(t2m_n, axis=0), '', 
#       '$\mathregular{\sigma}_{\mathregular{T}}$ [K]', 
#       np.linspace(2, 6, 9), 'PuBu', 'sigmat2m_%d-%d_jet'%(years[0],years[-1]), 
#       e_n=np.nanmean(lat_jet_nhml, axis=0),
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
## Correlation between O3 and distance from eddy-driven jet
#map_nh(lat_gmi_n, lng_gmi_n, r_o3jetdist, '', 
#       r'$\it{r}\:$(O$_\mathregular{3}$, jet distance)', 
#       np.linspace(-0.7, 0.7, 8), 'bwr', 'r_o3jet_%d-%d_jet'%(years[0],
#       years[-1]), e_n=np.nanmean(lat_jet_nhml, axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
## Change in O3 with change in eddy-driven jet position 
#map_nh(lat_gmi_n, lng_gmi_n, m_o3jetdist, '', 
#       '[ppbv degree$^{\mathregular{-1}}$]', np.linspace(-0.2, 0.4, 7), 
#       'gist_earth_r', 'do3djet_%d-%d_jet'%(years[0],years[-1]), 
#       e_n=np.nanmean(lat_jet_nhml, axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
## dO3/dRH from +Chemistry
#map_nh(lat_gmi_n, lng_gmi_n, do3drh_n, '', 
#       'dO$_{\mathregular{3}}$/dRH [ppbv %$^{\mathregular{-1}}$]', 
#       np.linspace(-0.5, 0.5, 6), 'bwr', 'do3drhchem_%d-%d_jet'%(years[0],
#       years[-1]), e_n=np.nanmean(lat_jet_nhml, axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
## Correlation between O3 and distance from eddy-driven jet from GMI CTM 
#map_nh(lat_gmi_n, lng_gmi_n, r_o3jetdist, '', 
#   r'$\it{r}\:$(O$_\mathregular{3}$, $\mathregular{\phi}_{\mathregular'+
#   '{jet}}$ $-$ $\mathregular{\phi}$)', np.linspace(-0.7, 0.7, 8), 
#   'bwr', 'r_o3jetdist_%d-%d_jet'%(years[0], years[-1]), 
#   e_n=np.nanmean(lat_jet_nhml, axis=0), 
#   eerr_n=np.zeros(lat_jet_nhml.shape[1]), oceanon='yes')
## Change in O3 with change in eddy-driven jet position from GMI CTM 
#map_nh(lat_gmi_n, lng_gmi_n, m_o3jetdist, '', 
#   r'$\mathregular{\Delta}$ O$_{\mathregular{3}}$ [ppbv $^{\circ \mathregular{-1}}$]',
#   np.linspace(-0.2, 0.4, 7), 'gist_earth_r', 'do3djet_%d-%d_jet'%(years[0],years[-1]), 
#       e_n=np.nanmean(lat_jet_nhml, axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]), oceanon='yes')
## r(RH, O3) from +Chemistry
#map_nh(lat_gmi_n, lng_gmi_n, r_rho3_n, '', 
#    r'$\it{r}\:$(RH, O$_\mathregular{3}$)', np.linspace(-1., 1., 11), 'bwr', 
#    'r_rho3chem_%d-%d_jet'%(years[0],years[-1]), 
#    e_n=np.nanmean(lat_jet_nhml, axis=0), 
#    eerr_n=np.zeros(lat_jet_nhml.shape[1]))
## Old 
#map_nh(lat_gmi_n, lng_gmi_n, np.mean(DIR500, axis=0), '', 
#       'Wind direction$_{\mathregular{500\:hPa}}$ [$^{\circ}$]', 
#       np.linspace(0, 360, 16), 'twilight', 
#       'meanwinddir_%d-%d_jet'%(years[0],years[-1])) 
#map_nh(lat_gmi_n, lng_gmi_n, do3drh_n,
#       'dO$_{\mathregular{3}}$/dRH [ppbv %$^{\mathregular{-1}}$]', 
#       np.linspace(-0.75, 0.75, 7), 'bwr', 'do3drh_%d-%d'%(years[0],years[-1]))
## O3, T, RH, wind direction correlations in + Chemistry simulation
#map_nh(lat_gmi_n, lng_gmi_n, r_rho3_n, '', 
#    r'$\it{r}\:$(RH, O$_{\mathregular{3}}$)', np.linspace(-0.9, 0.9, 10), 
#    'bwr', 'r_rho3_%d-%d_jet'%(years[0],years[-1]))
#map_nh(lat_gmi_n, lng_gmi_n, r_dirrh_n, '', 
#       r'$\it{r}\:$(Wind direction$_{\mathregular{500\:hPa}}$, RH)', 
#       np.linspace(-0.9, 0.9, 10), 'bwr', 
#       'r_dirrh_%d-%d_jet'%(years[0],years[-1])) 
#map_nh(lat_gmi_n, lng_gmi_n, r_dirq500_n, '', 
#       r'$\it{r}\:$(Wind direction$_{\mathregular{500\:hPa}}$, '+
#       'q$_{\mathregular{500\:hPa}}$)', np.linspace(-0.9, 0.9, 10), 'bwr', 
#       'r_dirq500_%d-%d_jet'%(years[0],years[-1])) 
#map_nh(lat_gmi_n, lng_gmi_n, r_diro3_n, '', 
#       r'$\it{r}\:$(Wind direction$_{\mathregular{500\:hPa}}$, '+
#       'O$_{\mathregular{3}}$)', np.linspace(-0.9, 0.9, 10), 'bwr', 
#       'r_diro3_%d-%d_jet'%(years[0],years[-1]))
## dO3/dT from + Chemistry overlaid with O3-T correlation contours
#map_do3dt_ro3t_nh(do3dt_n, r_t2mo3_n, lat_gmi_n, lng_gmi_n, '2008-2010')
## Plot O3-T correlation over India for + Chemistry, Transport, and difference
#map_extent(r_t2mo3_n, lat_gmi_n, lng_gmi_n, 65., 35., 95., 5., 
#    'r(T, O$_\mathregular{3,\:+\:Chemistry}$)', np.linspace(-1.,1.,11),
#    'bwr', 'india_r_t2mo3chem_%d-%d.eps'%(years[0],years[-1]))
#map_extent(r_t2mo3_n_dat, lat_gmi_n, lng_gmi_n, 65., 35., 95., 5., 
#    'r(T, O$_\mathregular{3,\:Transport}$)', np.linspace(-1.,1.,11),
#    'bwr', 'india_r_t2mo3transport_%d-%d.eps'%(years[0],years[-1]))
#map_extent((r_t2mo3_n-r_t2mo3_n_dat), lat_gmi_n, lng_gmi_n, 65., 35., 95., 5., 
#    'r(T, O$_\mathregular{3,\:+\:Chemistry}$) $-$ r(T, O$_\mathregular{3,'+
#    '\:Transport}$)', np.linspace(-0.2,0.2,11), 'bwr', 
#    'india_r_t2mo3_diffchemtransport_%d-%d.eps'%(years[0],years[-1]))
## Plot mean U500 (eddy-driven jet) winds 
#map_nh(lat_gmi_n, lng_gmi_n, np.mean(U500, axis=0),
#       '<U$_{\mathregular{500\:hPa}}$> [m s$^{\mathregular{-1}}$]', 
#       np.linspace(-10,15,11), 'PuBu', 
#       'meanU500_%d-%d'%(years[0],years[-1]))
## Plot 2-meter temperature variability 
#map_nh(lat_gmi_n, lng_gmi_n, np.std(t2m_n, axis=0),
#   '$\mathregular{\sigma}_{\mathregular{T}}$ [K]', 
#   np.linspace(0,6,11), 'PuBu', 
#   'stdt2m_%d-%d'%(years[0],years[-1]))
## Plot O3 at, above, and below the mean jet position in Eastern North America
#timeseries_o3atabovebelow_jet(o3_nhml, lat_jet_nhml, lat_nhml, lng_nhml, 
#    269, 294)        
        
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


"""PLOT MEAN LOCATION OF JET STREAM (MAX U WINDS AT 500 HPA) AND 
   O3, dO3/dT, and r(O3, T) WITHIN +/- 10 DEGREES OF JET STREAM""" 
#U500_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(U500, 
#    lat_gmi_n, lng_gmi_n, 0., 360., 20., 70.)
#o3_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(o3_n, 
#    lat_gmi_n, lng_gmi_n, 0., 360., 20., 70.)
#do3dt_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(do3dt_n, 
#    lat_gmi_n, lng_gmi_n, 0., 360., 20., 70.)   
#r_t2mo3_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(
#    r_t2mo3_n, lat_gmi_n, lng_gmi_n, 0., 360., 20., 70.)      
#land_nhml = globalo3_calculate.find_grid_overland(lat_nhml, lng_nhml)   
#lat_jet_nhml, o3_jet_nhml = globalo3_calculate.find_field_atjet(o3_nhml, 
#    U500_nhml, lat_nhml, lng_nhml, 10, anom=True)
#lat_jet_nhml, do3dt_jet_nhml = globalo3_calculate.find_field_atjet(do3dt_nhml, 
#    U500_nhml, lat_nhml, lng_nhml, 10)
#lat_jet_nhml, r_t2mo3_jet_nhml = globalo3_calculate.find_field_atjet(
#    r_t2mo3_nhml, U500_nhml, lat_nhml, lng_nhml, 10) 
## O3 anomaly in vicinity of jet 
#edjetlocation_fieldatedjet(lat_nhml, 
#    lng_nhml,
#    lat_jet_nhml, 
#    np.nanmean(o3_jet_nhml, axis=0),
#    'bwr', 
#    r'$\mathregular{\Delta}$ O$_{\mathregular{3}}$ [ppbv]',  
#    np.linspace(-14., 14., 8), 
#    'o3anom_jet', 
#    'gmictm_%d-%d'%(years[0], years[-1]), 
#    skiplng=6)
## dO3/dT in vicinity of jet
#edjetlocation_fieldatedjet(lat_nhml, 
#    lng_nhml,
#    lat_jet_nhml, 
#    do3dt_jet_nhml,
#    'Reds', 
#    'dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]',  
#    np.linspace(0, 2, 6),
#    'do3dt_jet', 
#    'gmictm_%d-%d'%(years[0], years[-1]), 
#    skiplng=6)
## r(T, O3) in vicinity of jet
#edjetlocation_fieldatedjet(lat_nhml, 
#    lng_nhml,
#    lat_jet_nhml, 
#    r_t2mo3_jet_nhml,
#    'bwr', 
#    'dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]',  
#    np.linspace(-1, 1, 11),
#    'rt2mo3_jet', 
#    'gmictm_%d-%d'%(years[0], years[-1]), 
#    skiplng=6)


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