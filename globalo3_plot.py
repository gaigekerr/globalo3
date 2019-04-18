#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
xxx

Revision History
    23012019 -- initial version created
    28012019 -- function 'map_global' edited to include stippling for grid 
                cells with statistically significant values
    03022019 -- function 'timeseries_seasonalo3' added
    19022019 -- map functions changed to just focus on Northern Hemisphere
"""
# Change font
import sys
if 'mpl' not in sys.modules:
    import matplotlib as mpl
    prop = mpl.font_manager.FontProperties(
        fname='/Users/ghkerr/Library/Fonts/cmunbmr.ttf')
    mpl.rcParams['font.family'] = prop.get_name()
    prop = mpl.font_manager.FontProperties(
        fname='/Users/ghkerr/Library/Fonts/cmunbbx.ttf')
    mpl.rcParams['mathtext.bf'] = prop.get_name()
    # for unicode minus/negative sign implementation
    mpl.rcParams['axes.unicode_minus'] = False


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
    p_n=None, e_n=None, eerr_n=None, alpha=0.05, oceanon='yes'):
    """plot desired quantity over North America
    
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
        Contour levels values
    cmap : str
        Colormap name
    fstr : str
        Output filename suffix
    p_n : numpy.ndarray or NoneType
        If type(p_n) is numpy.ndarray, gridded values less than the 
        significance level (alpha) are plotted as scatterpoints
    e_n : 
        If type(e_n) is numpy.ndarray, values are plotted as scatterpoints
    eerr_n : 
        If type(eerr_n) is numpy.ndarray, values are plotted as errorbars 
        on the scatterpoints of e_n
    alpha = int
        The significance level, the probability of rejecting the null 
        hypothesis p-values are less than alpha
    oceanon : str
        If 'yes', map adds ocean polygons feature        

    Returns
    -------
    None                    
    """
    import copy
    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure(figsize=(8,3))
    ax=plt.subplot2grid((1,2), (0,0), colspan=2,
                        projection=ccrs.Miller(central_longitude=0.))
    ax.set_title(title, fontsize=14, x=0.02, ha='left')    
    if oceanon == 'yes':
        ax.add_feature(cfeature.OCEAN, zorder=10, lw = 0.0, color='lightgrey')
    ax.coastlines(lw=0.25, color='k')
    ax.set_extent([-180., 180., 0., 85.])    
    cmap = plt.get_cmap(cmap)
    mb = ax.contourf(lng_n, lat_n, field_n, clevs, cmap=cmap,
                     transform=ccrs.PlateCarree(), extend='both')
    if p_n is None: pass
    else:
        p_n = copy.deepcopy(p_n)
        p_n[p_n >= alpha] = np.nan
        p_n[np.isnan(p_n) == False] = 5.
        lng_n_gridded, lat_n_gridded = np.meshgrid(lng_n, lat_n)
        ax.scatter(lng_n_gridded[1::2,::3], lat_n_gridded[1::2,::3], 
                   s=p_n[1::2,::3], facecolor='k', lw=0, marker='.', 
                   transform=ccrs.PlateCarree())
    if e_n is None: pass
    else:
        # Plot only every X number of longitude values
        skiplng = 6
        ax.errorbar(lng_n[::skiplng], e_n[::skiplng], 
                    yerr=eerr_n[::skiplng], zorder=12, color='k', markersize=2, 
                    elinewidth=0.5, ecolor='k', fmt='o', 
                    transform=ccrs.PlateCarree())
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

import numpy as np
import sys
sys.path.append('/Users/ghkerr/phd/globalo3/')
import globalo3_open, globalo3_calculate
sys.path.append('/Users/ghkerr/phd/transporto3/')
import transporto3_open
years = [2008, 2009, 2010]
# Define hemispheres and their O3 seasons
latmin_n, lngmin_n, latmax_n, lngmax_n = -1., 0., 90., 360.
months_n = ['jun', 'jul', 'aug']
# Load data the 
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
    # Calculate relative humidity
    rh_n = globalo3_calculate.calculate_rh_from_q(years, 
        list(np.arange(0,24,1)), lngmin_n, latmax_n, lngmax_n, latmin_n)
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
    # Calculate wind direction at 500 hPa
    DIR500 = globalo3_calculate.convert_uv_tocardinal(U500, V500) 
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
    # Calculate dO3/dRH
    do3drh_n = globalo3_calculate.calculate_do3dt(rh_n*100., o3_n, lat_gmi_n, 
        lng_gmi_n)
    # Calculate correlation coefficients
    r_t2mo3_n = globalo3_calculate.calculate_r(t2m_n, o3_n, lat_gmi_n, 
        lng_gmi_n)
    r_rho3_n = globalo3_calculate.calculate_r(rh_n, o3_n, lat_gmi_n, 
        lng_gmi_n)
    r_t2mrh_n = globalo3_calculate.calculate_r(rh_n, t2m_n, lat_gmi_n, 
        lng_gmi_n)
    r_t2mo3_n_dat = globalo3_calculate.calculate_r(t2m_n, o3_dat_n, lat_gmi_n, 
        lng_gmi_n)
    r_o3u500_n = globalo3_calculate.calculate_r(U500, o3_dat_n, lat_gmi_n, 
        lng_gmi_n)
    r_dirrh_n = globalo3_calculate.calculate_r(DIR500, rh_n, lat_gmi_n, 
        lng_gmi_n)
    r_dirq500_n = globalo3_calculate.calculate_r(DIR500, QV500, lat_gmi_n, 
        lng_gmi_n)
    r_diro3_n = globalo3_calculate.calculate_r(DIR500, o3_n, lat_gmi_n, 
        lng_gmi_n)    
    # Calculate trends
    do3dt_byyr_n, do3dt_ls_n, do3dt_lsp_n, do3dt_mk_n, do3dt_mkp_n = \
    globalo3_calculate.calculate_trends_do3dt(t2m_n, o3_n, lat_gmi_n, 
        lng_gmi_n, years)
    edgar_byyr_n, edgar_ls_n, edgar_lsp_n, edgar_mk_n, edgar_mkp_n = \
    globalo3_calculate.calculate_trends(nox_edgar_n, lat_gmi_n, lng_gmi_n, 
        years)


"""PLOT MEAN GLOBAL FIELDS AND CORRELATIONS"""
#map_nh(lat_gmi_n, lng_gmi_n, np.mean(o3_n, axis=0), '', 
#       'O$_{\mathregular{3}}$ [ppbv]', np.linspace(25, 65, 11), 'PuBu', 
#       'meano3_%d-%d_jet'%(years[0],years[-1]), 
#       e_n=np.nanmean(lat_jet_nhml,axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
#map_nh(lat_gmi_n, lng_gmi_n, np.nanstd(o3_n, axis=0), '(a)', 
#       '$\mathregular{\sigma}_{\mathregular{O}_{\mathregular{3}}}$ [ppbv]', 
#       np.linspace(5, 11, 7), 'PuBu', 'sigmao3_%d-%d_jet'%(years[0],years[-1]), 
#       e_n=np.nanmean(lat_jet_nhml, axis=0),
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
#map_nh(lat_gmi_n, lng_gmi_n, np.mean(nox_edgar_n, axis=0), '', 
#       'NO$_{x\:\mathregular{, EDGAR}}$ [kg m$^{\mathregular{-2}}$ '+
#       's$^{\mathregular{-1}}$]', np.linspace(0e-10, 1e-10, 11), 'PuBu', 
#       'meanedgarnox_%d-%d_jet'%(years[0],years[-1]), 
#       e_n=np.nanmean(lat_jet_nhml, axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
#map_nh(lat_gmi_n, lng_gmi_n, np.mean(DIR500, axis=0), '', 
#       'Wind direction$_{\mathregular{500\:hPa}}$ [$^{\circ}$]', 
#       np.linspace(0, 360, 16), 'twilight', 
#       'meanwinddir_%d-%d_jet'%(years[0],years[-1])) 
#map_nh(lat_gmi_n, lng_gmi_n, do3dt_n, '(a) Transport and Chemistry', 
#       'dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]', 
#       np.linspace(0, 3., 7), 'Reds', 'do3dt_%d-%d_jet'%(years[0],years[-1]),
#       e_n=np.nanmean(lat_jet_nhml, axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
#map_nh(lat_gmi_n, lng_gmi_n, do3dt_dat_n, '(b) Transport only', 
#       'dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]', 
#       np.linspace(0, 3., 7), 'Reds', 'do3dt_transport_%d-%d_jet' %(years[0],
#       years[-1]), e_n=np.nanmean(lat_jet_nhml, axis=0),
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
#map_nh(lat_gmi_n, lng_gmi_n, do3drh_n,
#       'dO$_{\mathregular{3}}$/dRH [ppbv %$^{\mathregular{-1}}$]', 
#       np.linspace(-0.75, 0.75, 7), 'bwr', 'do3drh_%d-%d'%(years[0],years[-1]))
## Change in O3 with change in eddy-driven jet position 
#map_nh(lat_nhml, lng_nhml, m_o3jetdist, '', 
#       '[ppbv degree$^{\mathregular{-1}}$]', np.linspace(-0.2, 0.4, 7), 
#       'gist_earth_r', 'do3djet_%d-%d_jet'%(years[0],years[-1]), 
#       e_n=np.nanmean(lat_jet_nhml, axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
## O3, T, RH, wind direction correlations in + Chemistry simulation
#map_nh(lat_gmi_n, lng_gmi_n, r_t2mo3_n, '(c)', 
#    r'$\it{r}\:$(T, O$_\mathregular{3}$)', np.linspace(-1., 1., 11), 'bwr', 
#    'r_t2mo3chem_%d-%d_jet'%(years[0],years[-1]), 
#    e_n=np.nanmean(lat_jet_nhml, axis=0), 
#    eerr_n=np.zeros(lat_jet_nhml.shape[1]))
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
## Difference in O3-T correlation in + Chemistry and Transport 
#map_nh(lat_gmi_n, lng_gmi_n, (r_t2mo3_n-r_t2mo3_n_dat), 
#    'r(T, O$_\mathregular{3,\:+\:Chemistry}$) $-$ r(T, O$_\mathregular{3,'+
#    '\:Tr;ansport}$)', np.linspace(-.4, .4, 9), 'bwr', 
#    'r_t2mo3_diffchemtransport_%d-%d'%(years[0],years[-1]))
## Correlation between O3 and distance from eddy-driven jet
#map_nh(lat_nhml, lng_nhml, r_o3jetdist, '', 
#       r'$\it{r}\:$(O$_\mathregular{3}$, jet distance)', 
#       np.linspace(-0.7, 0.7, 8), 'bwr', 'r_o3jet_%d-%d_jet'%(years[0],
#       years[-1]), e_n=np.nanmean(lat_jet_nhml, axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))
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
#    lat_gmi_n, lng_gmi_n, 0., 360., 23., 60.)
#o3_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(o3_n, 
#    lat_gmi_n, lng_gmi_n, 0., 360., 23., 60.)
#do3dt_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(do3dt_n, 
#    lat_gmi_n, lng_gmi_n, 0., 360., 23., 60.)   
#r_t2mo3_nhml, lat_nhml, lng_nhml = globalo3_calculate.find_grid_in_bb(
#    r_t2mo3_n, lat_gmi_n, lng_gmi_n, 0., 360., 23., 60.)      
#land_nhml = globalo3_calculate.find_grid_overland(lat_nhml, lng_nhml)   
#lat_jet_nhml, o3_jet_nhml = globalo3_calculate.find_field_atjet(o3_nhml, 
#    U500_nhml, lat_nhml, lng_nhml, 10, anom=True)
#lat_jet_nhml, do3dt_jet_nhml = globalo3_calculate.find_field_atjet(do3dt_nhml, 
#    U500_nhml, lat_nhml, lng_nhml, 10)
#lat_jet_nhml, r_t2mo3_jet_nhml = globalo3_calculate.find_field_atjet(
#    r_t2mo3_nhml, U500_nhml, lat_nhml, lng_nhml, 10) 
#m_o3jetdist, r_o3jetdist = globalo3_calculate.calculate_o3jet_relationship(
#    o3_nhml, lat_jet_nhml, lat_nhml, lng_nhml)
#lng_fr = lng_nhml
#lat_fr = lat_nhml
#land_fr = land_nhml
#jet_lat = lat_jet_nhml
#r_t2mo3_jet = r_t2mo3_jet_nhml
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
#import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
#fig = plt.figure()
#axt = plt.subplot2grid((2,2), (0,0), colspan=2, projection=
#                       ccrs.Miller(central_longitude=0.))
#axt.set_title('(a)', fontsize=14, x=0.02, ha='left')    
#axb = plt.subplot2grid((2,2), (1,0), colspan=2)
#axb.set_title('(b)', fontsize=14, x=0.02, ha='left')    
#axt.add_feature(cfeature.OCEAN, zorder=1, lw=0.0, color='lightgrey')
#axt.coastlines(lw=0.25, color='k')
#axt.set_extent([lng_fr.min()-180, lng_fr.max()-180, lat_fr.min(), 
#                lat_fr.max()])    
#axt.set_xticks([-180,-135,-90,-45,0,45,90,135,180])
#axt.set_xticklabels([''])
## Plot only every X number of longitude values
#skiplng = 6
#axt.errorbar(lng_fr[::skiplng], np.nanmean(jet_lat, axis=0)[::skiplng], 
#             yerr=np.nanstd(jet_lat, axis=0)[::skiplng], zorder=10, color='k',
#             markersize=2, elinewidth=0.5, ecolor='k', fmt='o', 
#             transform=ccrs.PlateCarree())
## Roll O3 differences arond the jet so that they align with the top (map)
## subplot
#o3_anom_jet_roll = np.roll(o3_jet_nhml, int(lng_fr.shape[0]/2.), axis=2)
#mb = axb.contourf(lng_fr, np.arange(-10, 11, 1), 
#                  np.nanmean(o3_anom_jet_roll, axis=0), 
#                  np.linspace(-14, 14, 8), cmap=plt.get_cmap('bwr'), 
#                  extend='both')
#axb.set_xlabel('Longitude [$^{\circ}$]', fontsize=14)
#axb.set_xticks(np.linspace(0, 360, 9))
#axb.set_xticklabels([-180,-135,-90,-45,0,45,90,135,180], fontsize=12)
#axb.set_ylabel('Grid cells from jet', fontsize=14)
#axb.tick_params(labelsize=12)
## Add colorbar
#colorbar_axes = plt.gcf().add_axes([0.78,0.135,0.02,0.375])
#colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical')
#colorbar.ax.tick_params(labelsize=12)
#colorbar.set_label(r'$\mathregular{\Delta}$ O$_{\mathregular{3}}$', 
#                   fontsize=14)
## Find where at least 50% of transect are land-based grid cells
#whereshore = np.where((ocean_frac>0.47) & (ocean_frac<0.5))[0]
## Manually remove duplicates and close doubles (this is clunky!)
#whereshore = np.delete(whereshore, [1,3,5,6])
#for ws in whereshore: 
#    axb.axvline(x=lng_fr[ws], c='k', lw=0.75)
#plt.gcf().subplots_adjust(right=0.75, hspace=-0.1)
#plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#            'edjetlocation_o3edjet.eps', dpi=300)
#fig = plt.figure()
#axt = plt.subplot2grid((2,2), (0,0), colspan=2, projection=
#                       ccrs.Miller(central_longitude=0.))
#axb = plt.subplot2grid((2,2), (1,0), colspan=2)
#axb.set_title('(d)', fontsize=14, x=0.02, ha='left')    
#axt.add_feature(cfeature.OCEAN, zorder=1, lw=0.0, color='lightgrey')
#axt.coastlines(lw=0.25, color='k')
#axt.set_extent([lng_fr.min()-180, lng_fr.max()-180, lat_fr.min(), 
#                lat_fr.max()])    
#axt.set_xticks([-180,-135,-90,-45,0,45,90,135,180])
#axt.set_xticklabels([''])
## Plot only every X number of longitude values
#skiplng = 6
#axt.errorbar(lng_fr[::skiplng], np.nanmean(jet_lat, axis=0)[::skiplng], 
#             yerr=np.nanstd(jet_lat, axis=0)[::skiplng], zorder=10, color='k',
#             markersize=2, elinewidth=0.5, ecolor='k', fmt='o', 
#             transform=ccrs.PlateCarree())
## Roll O3 differences arond the jet so that they align with the top (map)
## subplot
#do3dt_jet_roll = np.roll(do3dt_jet_nhml, int(lng_fr.shape[0]/2.), axis=1)
#mb = axb.contourf(lng_fr, np.arange(-10, 11, 1), do3dt_jet_roll, 
#                  np.linspace(0, 2, 6), cmap=plt.get_cmap('Reds'),
#                  extend='both')
#axb.set_xlabel('Longitude [$^{\circ}$]', fontsize=14)
#axb.set_xticks(np.linspace(0, 360, 9))
#axb.set_xticklabels([-180,-135,-90,-45,0,45,90,135,180], fontsize=12)
#axb.set_ylabel('Grid cells from jet', fontsize=14)
#axb.tick_params(labelsize=12)
## Add colorbar
#colorbar_axes = plt.gcf().add_axes([0.78,0.135,0.02,0.375])
#colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical')
#colorbar.ax.tick_params(labelsize=12)
#colorbar.set_label('dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]',
#                   fontsize=14)
## Find where at least 50% of transect are land-based grid cells
#whereshore = np.where((ocean_frac>0.47) & (ocean_frac<0.5))[0]
## Manually remove duplicates and close doubles (this is clunky!)
#whereshore = np.delete(whereshore, [1,3,5,6])
#for ws in whereshore: 
#    axb.axvline(x=lng_fr[ws], c='k', lw=0.75)
#plt.gcf().subplots_adjust(right=0.75, hspace=-0.1)
#plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#            'edjetlocation_do3dtedjet.eps', dpi=300)
#fig = plt.figure()
#axt = plt.subplot2grid((2,2), (0,0), colspan=2, projection=
#                       ccrs.Miller(central_longitude=0.))
#axb = plt.subplot2grid((2,2), (1,0), colspan=2)
#axb.set_title('(c)', fontsize=14, x=0.02, ha='left')    
#axt.add_feature(cfeature.OCEAN, zorder=1, lw=0.0, color='lightgrey')
#axt.coastlines(lw=0.25, color='k')
#axt.set_extent([lng_fr.min()-180, lng_fr.max()-180, lat_fr.min(), 
#                lat_fr.max()])    
#axt.set_xticks([-180,-135,-90,-45,0,45,90,135,180])
#axt.set_xticklabels([''])
## Plot only every X number of longitude values
#skiplng = 6
#axt.errorbar(lng_fr[::skiplng], np.nanmean(jet_lat, axis=0)[::skiplng], 
#             yerr=np.nanstd(jet_lat, axis=0)[::skiplng], zorder=10, color='k',
#             markersize=2, elinewidth=0.5, ecolor='k', fmt='o', 
#             transform=ccrs.PlateCarree())
## Roll O3 differences arond the jet so that they align with the top (map)
## subplot
#r_t2mo3_jet_roll = np.roll(r_t2mo3_jet, int(lng_fr.shape[0]/2.), axis=1)
#mb = axb.contourf(lng_fr, np.arange(-10, 11, 1), r_t2mo3_jet_roll, 
#                  np.linspace(-1, 1, 11), cmap=plt.get_cmap('bwr'))
#axb.set_xlabel('Longitude [$^{\circ}$]', fontsize=14)
#axb.set_xticks(np.linspace(0, 360, 9))
#axb.set_xticklabels([-180,-135,-90,-45,0,45,90,135,180])
#axb.set_ylabel('Grid cells from jet', fontsize=14)
#axb.tick_params(labelsize=12)
## Add colorbar
#colorbar_axes = plt.gcf().add_axes([0.78,0.135,0.02,0.375])
#colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical')
#colorbar.ax.tick_params(labelsize=12)
#colorbar.set_label(r'$\it{r}\:$(T, O$_\mathregular{3}$)', fontsize=14)
## Find where at least 50% of transect are land-based grid cells
#whereshore = np.where((ocean_frac>0.47) & (ocean_frac<0.5))[0]
## Manually remove duplicates and close doubles (this is clunky!)
#whereshore = np.delete(whereshore, [1,3,5,6])
#for ws in whereshore: 
#    axb.axvline(x=lng_fr[ws], c='k', lw=0.75)
#plt.gcf().subplots_adjust(right=0.75, hspace=-0.1)
#plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#            'edjetlocation_rt2mo3edjet.eps', dpi=300)


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
#o3_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(o3_n, lat_gmi_n, 
#    lng_gmi_n, 269., 294., 30., 60.)
#do3dt_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(do3dt_n, 
#    lat_gmi_n, lng_gmi_n, 269., 294., 30., 60.)
#r_t2mo3_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(r_t2mo3_n, 
#    lat_gmi_n, lng_gmi_n, 269., 294., 30., 60.)
#do3drh_na, lat_na, lng_na = globalo3_calculate.find_grid_in_bb(do3drh_n, 
#    lat_gmi_n, lng_gmi_n, 235., 290., 23., 60.)
#U500_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(U500, 
#    lat_gmi_n, lng_gmi_n, 269., 294., 30., 60.)
#land_ena = globalo3_calculate.find_grid_overland(lat_ena, lng_ena)
## Find fields over Western North America
#o3_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(o3_n, lat_gmi_n, 
#    lng_gmi_n, 230., 250., 30., 60.)
#do3dt_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(do3dt_n, 
#    lat_gmi_n, lng_gmi_n, 230., 250., 30., 60.)
#r_t2mo3_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(r_t2mo3_n, 
#    lat_gmi_n, lng_gmi_n, 230., 250., 30., 60.)
#do3drh_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(do3drh_n, 
#    lat_gmi_n, lng_gmi_n, 230., 250., 28., 60.)
#U500_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(U500, 
#    lat_gmi_n, lng_gmi_n, 230., 250., 30., 60.)
#land_wna = globalo3_calculate.find_grid_overland(lat_wna, lng_wna)
## Find fields over Northern Africa
#o3_af, lat_af, lng_af = globalo3_calculate.find_grid_in_bb(o3_n, lat_gmi_n, 
#    lng_gmi_n, 1., 35., 1., 30.)
#do3dt_af, lat_af, lng_af = globalo3_calculate.find_grid_in_bb(do3dt_n, 
#    lat_gmi_n, lng_gmi_n, 1., 35., 1., 30.)
#U500_af, lat_af, lng_af = globalo3_calculate.find_grid_in_bb(U500, 
#    lat_gmi_n, lng_gmi_n, 1., 35., 1., 30.)
#land_af = globalo3_calculate.find_grid_overland(lat_af, lng_af)
## Over Europe
#o3_eu, lat_eu, lng_eu = globalo3_calculate.find_grid_in_bb(o3_n, lat_gmi_n, 
#    lng_gmi_n, 0., 35., 35., 70.)
#do3dt_eu, lat_eu, lng_eu = globalo3_calculate.find_grid_in_bb(do3dt_n, 
#    lat_gmi_n, lng_gmi_n, 0., 35., 35., 70.)
#r_t2mo3_eu, lat_eu, lng_eu = globalo3_calculate.find_grid_in_bb(
#    r_t2mo3_n, lat_gmi_n, lng_gmi_n, 0., 35., 35., 70.)
#do3drh_eu, lat_eu, lng_eu = globalo3_calculate.find_grid_in_bb(do3drh_n, 
#    lat_gmi_n, lng_gmi_n, 0., 35., 35., 70.)
#U500_eu, lat_eu, lng_eu = globalo3_calculate.find_grid_in_bb(U500, 
#    lat_gmi_n, lng_gmi_n, 0., 35., 35., 70.)
#land_eu = globalo3_calculate.find_grid_overland(lat_eu, lng_eu)
## Over China
#o3_china, lat_china, lng_china = globalo3_calculate.find_grid_in_bb(o3_n, 
#    lat_gmi_n, lng_gmi_n, 100., 130., 20., 55.)
#do3dt_china, lat_china, lng_china = globalo3_calculate.find_grid_in_bb(do3dt_n, 
#    lat_gmi_n, lng_gmi_n, 100., 130., 20., 55.)
#r_t2mo3_china, lat_china, lng_china = globalo3_calculate.find_grid_in_bb(
#    r_t2mo3_n, lat_gmi_n, lng_gmi_n, 100., 130., 20., 55.)
#do3drh_china, lat_china, lng_china = globalo3_calculate.find_grid_in_bb(
#    do3drh_n, lat_gmi_n, lng_gmi_n, 85., 130., 20., 55.)
#U500_china, lat_china, lng_china = globalo3_calculate.find_grid_in_bb(U500, 
#    lat_gmi_n, lng_gmi_n, 100., 130., 20., 55.)
#land_china = globalo3_calculate.find_grid_overland(lat_china, lng_china)
## Plot O3 (colorcoded by dO3/dT) and U500 versus latitude
#scatter_latitude_o3u500(lat_na, lng_na, U500_na, o3_na*land_na, 
#    do3dt_na*land_na, 'Reds', '$\mathregular{\partial}$O$_{\mathregular{3}}$'+
#    ' $\mathregular{\partial}$T$^{\mathregular{-1}}$ [ppbv K$^{'
#    '\mathregular{-1}}$]', np.linspace(0, 2, 6), 'do3dt_northamerica')
#scatter_latitude_o3u500(lat_wna, lng_wna, U500_wna, o3_wna*land_wna, 
#    do3dt_wna*land_wna, 'Reds', '$\mathregular{\partial}$O$_{\mathregular{3'+
#    '}}$ $\mathregular{\partial}$T$^{\mathregular{-1}}$ [ppbv K$^{'
#    '\mathregular{-1}}$]', np.linspace(0, 2, 6), 'do3dt_westernnorthamerica')
#scatter_latitude_o3u500(lat_af, lng_af, np.abs(U500_af), o3_af*land_af, 
#    do3dt_af*land_af, '$\mathregular{\partial}$O$_{\mathregular{3}}$'+
#    ' $\mathregular{\partial}$T$^{\mathregular{-1}}$ [ppbv K$^{'
#    '\mathregular{-1}}$]', np.linspace(0, 3, 7), 'do3dt_africa')
#scatter_latitude_o3u500(lat_eu, lng_eu, U500_eu, o3_eu*land_eu, 
#    do3dt_eu*land_eu, 'Reds', '$\mathregular{\partial}$O$_{\mathregular{3}}$'+
#    ' $\mathregular{\partial}$T$^{\mathregular{-1}}$ [ppbv K$^{'+
#    '\mathregular{-1}}$]', np.linspace(0, 2, 6), 'do3dt_europe')
#scatter_latitude_o3u500(lat_china, lng_china, U500_china, o3_china*land_china, 
#    do3dt_china*land_china, 'Reds', '$\mathregular{\partial}$O$_{'+
#    '\mathregular{3}}$ $\mathregular{\partial}$T$^{\mathregular{-1}}$ [ppbv '+
#    'K$^{\mathregular{-1}}$]', np.linspace(0, 2, 6), 'do3dt_china')
## Plot O3 (colorcoded by dO3/dRH) and U500 versus latitude
#scatter_latitude_o3u500(lat_na, lng_na, U500_na, o3_na*land_na, 
#    do3drh_na*land_na, 'PuBu_r', '$\mathregular{\partial}$O$_{\mathregular'+
#    '{3}}$ $\mathregular{\partial}$RH$^{\mathregular{-1}}$ [ppbv %$^{'
#    '\mathregular{-1}}$]', np.linspace(-0.75, 0.25, 6), 'do3drh_northamerica')
#scatter_latitude_o3u500(lat_wna, lng_wna, U500_wna, o3_wna*land_wna, 
#    do3drh_wna*land_wna, 'PuBu_r', '$\mathregular{\partial}$O$_{\mathregular{'+
#    '3}}$ $\mathregular{\partial}$RH$^{\mathregular{-1}}$ [ppbv %$^{'
#    '\mathregular{-1}}$]', np.linspace(-0.75, 0.25, 6), 
#    'do3drh_westernnorthamerica')    
#scatter_latitude_o3u500(lat_eu, lng_eu, U500_eu, o3_eu*land_eu, 
#    do3drh_eu*land_eu, 'PuBu_r', '$\mathregular{\partial}$O$_{\mathregular{'+
#    '3}}$ $\mathregular{\partial}$RH$^{\mathregular{-1}}$ [ppbv %$^{'+
#    '\mathregular{-1}}$]', np.linspace(-.75, 0.25, 6), 'do3drh_europe')    
#scatter_latitude_o3u500(lat_china, lng_china, U500_china, o3_china*land_china, 
#    do3drh_china*land_china, 'PuBu_r', '$\mathregular{\partial}$O$_{'+
#    '\mathregular{3}}$ $\mathregular{\partial}$RH$^{\mathregular{-1}}$ [ppbv '+
#    '%$^{\mathregular{-1}}$]', np.linspace(-.75, 0.25, 6), 'do3drh_china')
#lng = lng_ena
#lat = lat_ena
#o3 = o3_ena*land_ena
#o3 = do3dt_ena*land_ena
#c = r_t2mo3_ena*land_ena
#U500_fr = U500_ena
#region = 'ena'
#ylabel = '$\mathregular{\sigma}_{\mathregular{O}_{\mathregular{3}}}$ [ppbv]'
#ylabel = '$\mathregular{O}_{\mathregular{3}}}$ [ppbv]'
#ylabel = 'dO$_{\mathregular{3}}$/dT [ppbv K$^{\mathregular{-1}}$]'
#cbar_label = r'$\it{r}\:$(T, O$_\mathregular{3}$)'
#clevs = np.linspace(-0.9, 0.9, 10)
#cmap = 'bwr'
#import numpy as np
#import matplotlib.pyplot as plt
#lng_m, lat_m = np.meshgrid(lng,lat)
## Calculate zonally-averaged U wind at 500 hPa
#U500_za = np.nanmean(U500_fr, axis=tuple((0,2)))
## n.b. Python crashes if the whole flattened grid is plotted (26000 
## points), skip every N points in flattened grid
#skippoints = 1
## Flatten latitude grid for plotting
#lat_f = lat_m.flatten()[::skippoints]
#o3_f = np.nanmean((o3),axis=0).flatten()[::skippoints]
#o3_f = (o3).flatten()[::skippoints]
#c_f = c.flatten()[::skippoints]
## Plotting O3 against latitude
#fig = plt.figure()
#ax = plt.subplot2grid((1,1),(0,0))
#ax.set_title('(a)', fontsize=14, x=0.02, ha='left')    
#mb = ax.scatter(lat_f, o3_f, s=8, c=c_f,
#                cmap=cmap_discretize(plt.get_cmap(cmap), len(clevs)-1), 
#                vmin=clevs[0], vmax=clevs[-1])
#colorbar_axes = plt.gcf().add_axes([0.83,0.15,0.03,0.7])
#colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
#                        extend='both')
#colorbar.set_label(cbar_label, fontsize=14)
#colorbar.set_ticks(clevs)
#colorbar.set_ticklabels(clevs)
#colorbar.ax.tick_params(labelsize=12)
#ax.set_xlim([20,70])
#ax.set_xlabel('Latitude [$^{\mathregular{\circ}}$]', fontsize=14)
#ax.tick_params(labelsize=12)
#ax.set_ylabel(ylabel, fontsize=14)
## Plotting zonally-averaged U500 against latitude
#ax2 = ax.twinx()
#ax2.plot(lat, U500_za, '-k')
#ax2.set_ylabel('U$_{\mathregular{500\:hPa}}$ [m s$^{\mathregular{-1}}$]',
#               fontsize=14)
#ax2.tick_params(labelsize=12)
#plt.gcf().subplots_adjust(right=0.7, hspace=-0.2)
#plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#            'scatter_latitude_meano3_%s.eps'%region, dpi=300)


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


"""FIND DRIVER OF dO3/dT RELATIONSHIP BY FIXING COMPONENTS"""
## Global mean O3-T correlation, and O3/T standard deviation
#gm_r = np.empty(shape=(lat_gmi_n.shape[0],lng_gmi_n.shape[0]))
#gm_r[:] = np.mean(r_t2mo3_n)
#gm_sy = np.empty(shape=(lat_gmi_n.shape[0],lng_gmi_n.shape[0]))
#gm_sy[:] = np.std(o3_n)
#gm_sx = np.empty(shape=(lat_gmi_n.shape[0],lng_gmi_n.shape[0]))
#gm_sx[:] = np.std(t2m_n)
## Calculate dO3/dT holding factors constant; i.e., dO3/dT is calculated by
## beta=r(O3,T)*[std(O3)/std(T)] (or, generically, as beta=r*(sy/sx))
#beta = r_t2mo3_n * (np.std(o3_n, axis=0)/np.std(t2m_n, axis=0)) 
#map_nh(lat_gmi_n, lng_gmi_n, beta,
#       r'$\mathregular{\beta}$ [ppbv K$^{\mathregular{-1}}$]', 
#       np.linspace(0, 3., 7), 'PuBu', 'beta')
## Fix r(O3,T)
#beta_rconstant = gm_r * (np.std(o3_n, axis=0)/np.std(t2m_n, axis=0))
#map_nh(lat_gmi_n, lng_gmi_n, beta_rconstant,
#       r'$\mathregular{\beta}_{\overline{\mathregular{r}}}$ '+
#       '[ppbv K$^{\mathregular{-1}}$]', np.linspace(-1.8, 0.0, 7), 'PuBu', 
#       'beta_rconstant')
## Fix r(O3,T), std(O3)
#beta_rconstant_syconstant = gm_r * (gm_sy/np.std(t2m_n, axis=0))
#map_nh(lat_gmi_n, lng_gmi_n, beta_rconstant_syconstant,
#       r'$\mathregular{\beta}_{\overline{\mathregular{r}},\:'+
#       '\overline{\mathregular{O_{3}}}}$ [ppbv K$^{\mathregular{-1}}$]', 
#       np.linspace(-2.8, 0.0, 8), 'PuBu', 'beta_rconstant_syconstant')
## Fix r(O3,T), std(T)
#beta_rconstant_sxconstant = gm_r * (np.std(o3_n, axis=0)/gm_sx)
#map_nh(lat_gmi_n, lng_gmi_n, beta_rconstant_sxconstant,
#       r'$\mathregular{\beta}_{\overline{\mathregular{r}},\:'+
#       '\overline{\mathregular{T}}}$ [ppbv K$^{\mathregular{-1}}$]', 
#       np.linspace(-0.105, 0.0, 8), 'PuBu', 'beta_rconstant_sxconstant')
## Fix std(O3)
#beta_syconstant = r_t2mo3_n * (gm_sy/np.std(t2m_n, axis=0))
#map_nh(lat_gmi_n, lng_gmi_n, beta_syconstant,
#       r'$\mathregular{\beta}_{\overline{\mathregular{O_3}}}$ '+
#       '[ppbv K$^{\mathregular{-1}}$]', np.linspace(-20, 16, 7), 'PuBu', 
#       'beta_syconstant')
## Fix std(T)
#beta_sxconstant = r_t2mo3_n * (np.std(o3_n, axis=0)/gm_sx)
#map_nh(lat_gmi_n, lng_gmi_n, beta_sxconstant,
#       r'$\mathregular{\beta}_{\overline{\mathregular{T}}}$ '+
#       '[ppbv K$^{\mathregular{-1}}$]', np.linspace(-0.75, 1.0, 7), 'PuBu', 
#       'beta_sxconstant')
## Fix std(O3), std(T)
#beta_sxconstant_syconstant = r_t2mo3_n * (gm_sy/gm_sx) # fix sx, sy
#map_nh(lat_gmi_n, lng_gmi_n, beta_sxconstant_syconstant,
#       r'$\mathregular{\beta}_{\overline{\mathregular{O_3}},\:'+
#       '\overline{\mathregular{T}}}$ [ppbv K$^{\mathregular{-1}}$]', 
#       np.linspace(-1, 1, 9), 'PuBu', 'beta_sxconstant_syconstant')


"""PRIMARY AIR CARE POSTER"""
#vmin=26
#vmax=56
#nlevs=10
#map_toar(np.nanmean([np.nanmean(urban_mean, axis=0), 
#    np.nanmean(rural_mean, axis=0)], axis=0), toarlat, toarlng, 225, 70, 304, 
#    20, 2, vmin, vmax, nlevs, '', 'O$_{\mathregular{3}}$ [ppbv]', 
#    'ruralurbanmean_northamerica_PCA')
#map_toar(np.nanmean([np.nanmean(urban_mean, axis=0), 
#    np.nanmean(rural_mean, axis=0)], axis=0), toarlat, toarlng, -10, 72, 36, 
#    28, 2, vmin, vmax, nlevs, '', 'O$_{\mathregular{3}}$ [ppbv]', 
#    'ruralurbanmean_europe_PCA')
#map_toar(np.nanmean([np.nanmean(urban_mean, axis=0), 
#    np.nanmean(rural_mean, axis=0)], axis=0), toarlat, toarlng, 114, 50, 145, 
#    18, 2, vmin, vmax, nlevs, '', 'O$_{\mathregular{3}}$ [ppbv]', 
#    'ruralurbanmean_asia_PCA')
#map_nh(lat_gmi_n, lng_gmi_n, np.mean(o3_n, axis=0), '', 
#       'O$_{\mathregular{3}}$ [ppbv]', np.linspace(25, 65, 11), 'YlGnBu', 
#       'meano3_%d-%d_jet'%(years[0],years[-1]), 
#       e_n=np.nanmean(lat_jet_nhml,axis=0), 
#       eerr_n=np.zeros(lat_jet_nhml.shape[1]))


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