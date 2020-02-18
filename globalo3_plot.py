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
    09072019 -- edit function 'edjetlocation_fieldatedjet' to handle and 
                number of grid cells/degrees from jet center
    17072019 -- function 'edjetlocation_fieldatedjet' name changed to 
                'fieldatjet'
    13082019 -- change scatterpoints/hatching in function 'map_hemisphere' to 
                reflect significance
    26082019 -- updated 'map_hemisphere' to overlap pcolor atop of contourf
    27082019 -- function 'rt2mo3_iav' added
    26092019 -- code to open 2016-2017 Chinese O3 observations added and 
                function 'zonalavg_byregion' edited to handle observations 
                from China
    30092019 -- function 'zonalavg_verticallyintegrated_meridional_flux' added
    09102019 -- code to process the eddy and mean meridional tracer fluxes on 
                days with a pole- versus equatorward jet added
    07112019 -- function 'map_o3anomcyclones_pweqjet' added
    13112019 -- function 'o3anom_atcyclones' added
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
    prop = matplotlib.font_manager.FontProperties(
        fname='/Users/ghkerr/Library/Fonts/cmunbmr.ttf')
    matplotlib.rcParams['mathtext.it'] = prop.get_name()
    # for unicode minus/negative sign implementation
    matplotlib.rcParams['axes.unicode_minus'] = False
    # change width and thickness of ticks/spines
    matplotlib.rcParams['axes.linewidth'] = 1.0
    matplotlib.rcParams['xtick.major.width'] = 1.0
    matplotlib.rcParams['xtick.minor.width'] = 1.0
    matplotlib.rcParams['ytick.major.width'] = 1.0
    matplotlib.rcParams['ytick.minor.width'] = 1.0

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
    hemisphere, fstr, contour=None, contour_levs=None, e_lng=None, e_n=None, 
    eerr_n=None, hatch=None, hatch_freq=None, extent=None, quiver=None, 
    lat_n_binned=None, lng_n_binned=None, oceanon='yes', extend='both', 
    pcolorflag=False, lbound=None, ubound=None):
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
    e_lng : numpy.ndarray or NoneType
        If type(e_lng) is numpy.ndarray, this optional argument corresponds 
        to the jet (or other scatterpoint) dataset
    e_n : numpy.ndarray or NoneType
        If type(e_n) is numpy.ndarray, values are plotted as scatterpoints
    eerr_n : numpy.ndarray or NoneType
        If type(eerr_n) is numpy.ndarray, values are plotted as errorbars 
        on the scatterpoints of e_n
    hatch : numpy.ndarray or NoneType
        If type(hatch) is numpy.ndarray, grid cells equal to 1 are hatched and 
        values equal to NaN are left empty, [lat, lng]
    hatch_freq : list
        Hatching type (i.e., +, /, o, etc.) and frequency
    extent : list or NoneType
        The extent (x0, x1, y0, y1) of the map in the given coordinate system    
    quiver : tuple or NoneType
        First (second) object, an numpy.ndarray, in tuple corresponds to 
        U-wind (V-wind)
    lat_n_binned : numpy.ndarray or NoneType
        Binned latitude coordinates corresponding to the pcolor array, units
        of degrees north, [lat_binned,]
    lng_n_binned : numpy.ndarray or NoneType
        Binned longitude coordinates corresponding to the pcolor array, units
        of degrees east, [lng_binned,]
    oceanon : str
        If 'yes', map adds ocean polygons feature        
    extend : str
        Extend settings for matplotlib colormap/colorbar (i.e., 'neither', 
        'both', 'min', 'max')
    pcolorflag : bool
        If False, field is plotted as filled contour; if True, field is plotted
        as pcolor. 
    lbound : float or NoneType
        Lower bound for pcolor masking; any pcolor values greater than or equal 
        to this value will be masked
    ubound : float or NoneType
        Upper bound for pcolor masking; any pcolor values less than or equal 
        to this value will be masked

    Returns
    -------
    None              
    """
    import numpy as np
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = 0.3     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure(figsize=(9,3.5))
    ax=plt.subplot2grid((1,2), (0,0), colspan=2,
                        projection=ccrs.Miller(central_longitude=0.))
    ax.set_title(title, fontsize=16, x=0.02, ha='left')    
    if oceanon == 'yes':
        ax.add_feature(cfeature.OCEAN, zorder=2, lw=0.0, color='lightgrey')
    ax.coastlines(lw=0.25, color='k', zorder=3)
    if extent is None:
        ax.set_extent([lng_n.min()-180., lng_n.max()-180., 
                       lat_n.min(), lat_n.max()])    
    else: 
        ax.set_extent(extent)
    # Plot filled contours; if pcolorflag == True, then plot field as pcolor
    cmap = plt.get_cmap(cmap)
    if pcolorflag == False: 
        mb = ax.contourf(lng_n, lat_n, field_n, clevs, cmap=cmap, 
            extend=extend, transform=ccrs.PlateCarree(), zorder=1)
    else: 
        if lbound is not None: 
            field_n = np.ma.masked_array(field_n, ((field_n >= lbound) & 
                (field_n <= ubound)))
        mb = ax.pcolor(lng_n, lat_n, field_n, 
            cmap=plt.get_cmap(cmap, len(clevs)-1), vmin=clevs[0], 
            vmax=clevs[-1], norm = mpl.colors.BoundaryNorm(clevs, 
            ncolors=cmap.N, clip=False), transform=ccrs.PlateCarree(), 
            zorder=6)
    # If specified, add contours and labels. 
    if contour is None: pass
    else:
        csthick = ax.contour(lng_n, lat_n, contour, [10.], colors='k', 
            linewidths=1.5, transform=ccrs.PlateCarree(), zorder=15)
        csmedium = ax.contour(lng_n, lat_n, contour, [8.], colors='k', 
            linestyles='--', linewidths=0.75, transform=ccrs.PlateCarree(), 
            zorder=15)
    # If specified, add hatching (e.g., for significance)
    if hatch is None: pass
    else:
        # Hatching for significance
        ax.contourf(lng_n, lat_n, hatch, hatches=hatch_freq, colors='none', 
                    transform=ccrs.PlateCarree())
    # If specified add scatterplots (e.g., for location of eddy-driven jet)
    if e_n is None: pass
    else:
        # Plot only every X number of longitude values
        skiplng = 6
        ax.errorbar(e_lng[::skiplng], e_n[::skiplng], yerr=eerr_n[::skiplng], 
            zorder=10, color='k', markersize=3, elinewidth=1.25,
            ecolor='k', fmt='o', transform=ccrs.PlateCarree())
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
    plt.gcf().subplots_adjust(left=0.02, right=0.86)#, hspace=-0.2)        
    colorbar_axes = plt.gcf().add_axes([0.89, ax.get_position().y0, 
        0.02, (ax.get_position().y1-ax.get_position().y0)]) 
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=clevs, extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label(cbar_label, fontsize=16)
    ax.outline_patch.set_zorder(20)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'map_%s_%s.png'%(hemisphere, fstr), #transparent = True, 
        dpi = 350)
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

#def timeseries_o3atabovebelow_jet(o3, lat_jet, lat, lng, lngl, lngr):
#    """for O3 and the latitude of the eddy-driven jet in a given region, plot 
#    a timeseries of O3 at (with +/- 2 degrees) of the mean jet location, O3 
#    above (8-12 degrees above) the mean jet location, and O3 below (8-12 
#    degrees below) the mean jet location. A timeseries of the regionally-
#    averaged jet latitude is plotted for each subplot. 
#    
#    Parameters
#    ----------
#    o3 : numpy.ndarray
#        O3 in region, units of ppbv, [time, lat, lng]  
#    lat_jet : numpy.ndarray
#        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
#        in region, units of degrees north[time, lng]
#    lat : numpy.ndarray
#        Latitude coordinates, units of degrees north, [lat,]
#    lng : numpy.ndarray
#        Longitude coordinates, units of degrees east, [lng,]
#    lngl : float
#        Longitude coordinate corresponding to the west side of the focus 
#        region; n.b., lngl < lngr
#    lngr : float
#        Longitude coordinate corresponding to the east side of the focus region
#        
#    Returns
#    -------
#    None
#    """
#    import numpy as np
#    import matplotlib.pyplot as plt
#    # Find mean jet position in Eastern North America
#    lat_jet_fr = lat_jet[:, np.abs(lng-lngl).argmin(): 
#        np.abs(lng-lngr).argmin()]
#    lat_jet_fr = np.nanmean(lat_jet_fr, axis=1) 
#    lat_jet_mean = np.nanmean(lat_jet_fr)
#    # O3 at the mean position of the jet (+/- 3 degrees)
#    o3at = o3[:, np.abs(lat-lat_jet_mean).argmin()-2:
#        np.abs(lat-lat_jet_mean).argmin()+2, np.abs(lng-lngl).argmin(): 
#        np.abs(lng-lngr).argmin()]
#    # O3 above the mean position of the jet (7-13 degrees above)    
#    o3above = o3[:, np.abs(lat-lat_jet_mean).argmin()+8:
#        np.abs(lat-lat_jet_mean).argmin()+12, np.abs(lng-lngl).argmin(): 
#        np.abs(lng-lngr).argmin()]    
#    # O3 below the mean position of the jet (7-13 degrees below)    
#    o3below = o3[:, np.abs(lat-lat_jet_mean).argmin()-12:
#        np.abs(lat-lat_jet_mean).argmin()-8, np.abs(lng-lngl).argmin(): 
#        np.abs(lng-lngr).argmin()]        
#    # Regional averaging 
#    o3at = np.nanmean(o3at, axis=tuple((1,2)))
#    o3above = np.nanmean(o3above, axis=tuple((1,2)))
#    o3below = np.nanmean(o3below, axis=tuple((1,2)))
#    # Plotting
#    fig = plt.figure()
#    ax1 = plt.subplot2grid((3,3), (0,0), colspan=3)
#    ax1b = ax1.twinx()
#    ax2 = plt.subplot2grid((3,3), (1,0), colspan=3)
#    ax2b = ax2.twinx()
#    ax3 = plt.subplot2grid((3,3), (2,0), colspan=3)
#    ax3b = ax3.twinx()
#    # Above jet
#    ax1b.axhspan(lat_jet_mean+8, lat_jet_mean+12, xmin=ax1b.get_xlim()[0], 
#        xmax=ax1b.get_xlim()[1], alpha=0.2, color='dodgerblue', zorder=0)
#    ax1.plot(o3above[:92], '-k')
#    ax1b.plot(lat_jet_fr[:92], ls='-', color='dodgerblue')
#    # At jet
#    ax2b.axhline(y=lat_jet_mean, xmin=ax2b.get_xlim()[0], 
#        xmax=ax2b.get_xlim()[1], color='dodgerblue', ls='--', alpha=0.5, 
#        zorder=2)
#    ax2b.axhspan(lat_jet_mean-2, lat_jet_mean+2, xmin=ax2b.get_xlim()[0], 
#        xmax=ax2b.get_xlim()[1], alpha=0.2, color='dodgerblue', zorder=0)
#    ax2b.plot(lat_jet_fr[:92], ls='-', color='dodgerblue', zorder=5)
#    ax2.plot(o3at[:92], '-k', zorder=5)
#    # Below jet
#    ax3b.axhspan(lat_jet_mean-12, lat_jet_mean-8, xmin=ax3b.get_xlim()[0], 
#        xmax=ax3b.get_xlim()[1], alpha=0.2, color='dodgerblue', zorder=0)
#    ax3b.plot(lat_jet_fr[:92], ls='-', color='dodgerblue')
#    ax3.plot(o3below[:92], '-k')
#    # Set limits, labels
#    for ax in [ax1, ax1b, ax2, ax2b, ax3, ax3b]:
#        ax.set_xlim([0, 91])
#        ax.set_xticks([0, 14, 30, 44, 61, 75])  
#        ax.set_xticklabels([''])
#    for ax in [ax1, ax2, ax3]:
#        ax.set_ylim([15, 60])
#        ax.set_yticks([15, 30, 45, 60])
#        ax.set_ylabel('O$_{\mathregular{3}}$ [ppbv]')
#    for ax in [ax1b, ax2b, ax3b]:
#        ax.set_ylim([30, 60])
#        ax.set_yticks([30, 40, 50, 60])
#        ax.set_ylabel('Jet position [$^{\mathregular{\circ}}$]')  
#    ax3.set_xticklabels(['1 June 2008', '', '1 July', '', '1 Aug', ''], 
#        ha = 'center', fontsize = 12)
#    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#                'timeseries_o3atabovebelow_jet.png', dpi=300)
#    return 

#def contourf_var_atcenter(center, var, lat, lng, varname, exam_rad, 
#    kind, cbar_label, clevs, cmap, fstr, SLP=None, H=None): 
#    """function retrieves variable of interest in the vincinity of 
#    (anti)cyclones and plots the mean O3 concentrations averaged over all 
#    systems in region/measuring period. Only (anti)cyclones over land are 
#    considered as part of the average.
#
#    Parameters
#    ----------
#    center : numpy.ndarray 
#        A value of 1 indicates the presence of a cyclone/anticylone for a 
#        particular day and location, [time, lat, lng]       
#    var : numpy.ndarray
#        Variable of interest in region, units of ppbv, [time, lat, lng]
#    lat : numpy.ndarray
#        Latitude coordinates, units of degrees north, [lat,]
#    lng : numpy.ndarray
#        Longitude coordinates, units of degrees east, [lng,]
#    varname : str
#        Variable name for output file
#    exam_rad : int
#        Search radius; the number of grid cells on each side of the (anti)
#        cyclone's center over which O3 concentrations will be retrieved
#    kind : str
#        'cyclone' or 'anticyclone' (for axis labels)
#    cbar_label : str
#        Label for the colorbar (field and units)
#    clevs : numpy.ndarray
#        Contour levels values
#    cmap : str
#        Colormap name        
#    fstr : str
#        Output filename suffix (should specify the type of system in variable
#        'kind' and the latitude) 
#    SLP: numpy.ndarray
#        Sea level pressure in region, units of Pa, [time, lat, lng]
#    H : numpy.ndarray
#        Geopotential height at 500 hPa in region, units of m, [time, lat, lng]
#        
#    Returns
#    -------
#    None
#    """
#    import numpy as np
#    import matplotlib.pyplot as plt
#    from mpl_toolkits.basemap import Basemap
#    bm = Basemap()    
#    # Identify locations of (anti)cyclones
#    where_center = np.where(center==1.)
#    # Find O3, SLP, and  surrounding (within +/- exam_rad) systems for each 
#    # (anti)cyclone
#    var_atcenter = np.empty(shape=(where_center[0].shape[0], exam_rad, 
#                                   exam_rad))
#    var_atcenter[:] = np.nan
#    SLP_atcenter = np.empty(shape=(where_center[0].shape[0], exam_rad, 
#                                   exam_rad))
#    SLP_atcenter[:] = np.nan
#    H_atcenter = np.empty(shape=(where_center[0].shape[0], exam_rad, 
#                                 exam_rad))
#    H_atcenter[:] = np.nan            
#    for system in np.arange(0, len(where_center[0]), 1):
#        system_coords = (where_center[0][system], 
#                         where_center[1][system], 
#                         where_center[2][system])
#        # (xpt, ypt)
#        if bm.is_land(lng[where_center[2][system]]-360., 
#                      lat[where_center[1][system]]) == True:   
#            var_atcenter[system]= var[system_coords[0], 
#                     system_coords[1]-(int(np.round(exam_rad/2))-1):
#                     system_coords[1]+int(np.round(exam_rad/2)),
#                     system_coords[2]-(int(np.round(exam_rad/2))-1):
#                     system_coords[2]+int(np.round(exam_rad/2))]
#            if SLP is not None: 
#                # Convert from Pa to hPa
#                SLP_atcenter[system]= SLP[system_coords[0], 
#                     system_coords[1]-(int(np.round(exam_rad/2))-1):
#                     system_coords[1]+int(np.round(exam_rad/2)),
#                     system_coords[2]-(int(np.round(exam_rad/2))-1):
#                     system_coords[2]+int(np.round(exam_rad/2))]/100.
#            if H is not None: 
#                # Convert from m to dm 
#                H_atcenter[system]= H[system_coords[0], 
#                     system_coords[1]-(int(np.round(exam_rad/2))-1):
#                     system_coords[1]+int(np.round(exam_rad/2)),
#                     system_coords[2]-(int(np.round(exam_rad/2))-1):
#                     system_coords[2]+int(np.round(exam_rad/2))]/10.            
#    # Plotting 
#    fig = plt.figure()
#    ax = plt.subplot2grid((1,1),(0,0))
#    mb = ax.contourf(np.nanmean(var_atcenter, axis=0), clevs, extend='both',
#                     cmap=plt.get_cmap(cmap))
#    
#    # Add contours for SLP and H, if needed
#    if SLP is None: pass
#    else:
#        cs = ax.contour(np.nanmean(SLP_atcenter, axis=0), colors='k')
#        plt.clabel(cs, fontsize=12, inline=1, fmt='%1.0f')    
#    if H is None: pass
#    else:
#        cs = ax.contour(np.nanmean(H_atcenter, axis=0), colors='w')
#        plt.clabel(cs, fontsize=12, inline=1, fmt='%1.0f')    
#    # Add colorbar
#    colorbar_axes = plt.gcf().add_axes([0.83,0.25,0.02,0.5])
#    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical',
#                            extend='max')
#    colorbar.ax.tick_params(labelsize=12)
#    colorbar.set_label(cbar_label, fontsize=14)
#    # Aesthetics (crosshairs for system's center, axes labels)
##    ax.axhline(y=(int(np.round(exam_rad/2))-1), xmin=ax.get_xlim()[0], 
##               xmax=ax.get_xlim()[1], ls='--', lw=2., color='k')
##    ax.axvline(x=(int(np.round(exam_rad/2))-1), ymin=ax.get_ylim()[0], 
##               ymax=ax.get_ylim()[1], ls='--', lw=2., color='k') 
#    ax.set_xticks(np.arange(0, exam_rad, 1))
#    ax.set_xticklabels(np.arange(-(int(np.round(exam_rad/2))-1), 
#                       int(np.round(exam_rad/2)), 1), fontsize=12)
#    ax.set_xlabel('E-W grid cells from %s center'%kind, fontsize=14)
#    ax.set_yticks(np.arange(0, exam_rad, 1))
#    ax.set_yticklabels(np.arange(-(int(np.round(exam_rad/2))-1), 
#                       int(np.round(exam_rad/2)), 1), fontsize=12)
#    ax.set_ylabel('N-S grid cells from %s center'%kind, fontsize=14)
#    plt.gcf().subplots_adjust(right=0.8)
#    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#                'contourf_%s_at%s.eps' %(varname,fstr), dpi=300)
#    return

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

def zonalavg_byregion(field, lat_gmi, lng_gmi, obs_aqs, lat_aqs, lng_aqs, 
    obs_naps, lat_naps, lng_naps, obs_emep, lat_emep, lng_emep, 
    obs_china, lat_china, lng_china, lng_ml, lat_jet_ml, label, units, 
    yticks, fstr): 
    """find zonally-averaged mean values of field of interest (i.e., r(T, O3) 
    or dO3/dT in Western North America, Eastern North America, Europe, and East 
    Asia from both CTM simulations and observations and plot as a function of 
    latitude alongside the mean position and variability of the eddy-driven 
    jet. 
    
    Parameters
    ----------
    field : numpy.ndarray     
        Field of interest from CTM simulation, [lat, lng]        
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    obs_aqs : list
        Field of interest at each AQS station, [stations,]
    lat_aqs : list 
        Latitude coorindate at each AQS station, units of degrees north, 
        [stations,]
    lng_aqs : list 
        Longitude coorindate at each AQS station, units of degrees east, 
        [stations,]
    obs_naps : list
        Field of interest at each NAPS station, [stations,]
    lat_naps : list 
        Latitude coorindate at each NAPS station, units of degrees north, 
        [stations,]
    lng_naps : list 
        Longitude coorindate at each NAPS station, units of degrees east, 
        [stations,]        
    obs_emep : list
        Field of interest at each EMEP station, [stations,]
    lat_emep : list 
        Latitude coorindate at each EMEP station, units of degrees north, 
        [stations,]
    lng_emep : list 
        Longitude coorindate at each EMEP station, units of degrees east, 
        [stations,]
    obs_china : list
        Field of interest at each Chinese station, [stations,]
    lat_china : list 
        Latitude coorindate at each Chinese station, units of degrees north, 
        [stations,]
    lng_china : list 
        Longitude coorindate at each Chinese station, units of degrees east, 
        [stations,]
    lng_ml : numpy.ndarray
        The longitude coordinates corresponding to the jet latitude, 
        units of degrees east, [lng]
    lat_jet_ml : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in the Northern Hemisphere mid-latitudes, units of degrees north, 
        [time, lng]
    label : str
        Name of field for y-axis label and legend
    units : str
        Units of quantity for y-axis label
    yticks : numpy.ndarray
        Y-ticks
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
    # Find quantity (model, model error, and observations) in regions and 
    # calculate land-based zonally-averaged quantities
    # Western North America
    # Model
    field_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(
        field, lat_gmi, lng_gmi, 235., 260., 25., 70.) 
    wna_left = geo_idx(235., lng_ml)
    wna_right = geo_idx(260., lng_ml)
    lat_jet_wna = lat_jet_ml[:, wna_left:wna_right+1]
    lat_jet_err_wna = np.nanstd(np.nanmean(lat_jet_wna, axis = 1))
    lat_jet_mean_wna = np.nanmean(lat_jet_wna)
    land_wna = globalo3_calculate.find_grid_overland(lat_wna, lng_wna)
    field_err_wna = 2*np.nanstd(field_wna*land_wna, axis=1)
    field_wna = np.nanmean(field_wna*land_wna, axis=1)
    # Observations
    obs_wna = np.where((np.array(lat_naps+lat_aqs) < 70.) &
                         (np.array(lat_naps+lat_aqs) > 25.) &
                         (np.array(lng_naps+lng_aqs) > 235.) &
                         (np.array(lng_naps+lng_aqs) < 260.))[0]
    # Combine NAPS and AQS observations 
    lat_obs_wna = np.array(lat_naps+lat_aqs)[obs_wna]
    obs_wna = np.array(obs_naps+obs_aqs)[obs_wna]
    # Bin observations by latitude 
    lat_obs_wna, obs_mean_wna, obs_std_wna = \
        globalo3_calculate.bin_observations_bylat(lat_wna[::2], obs_wna, 
        lat_obs_wna)
    # Eastern North America
    # Model 
    field_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
        field, lat_gmi, lng_gmi, 260., 295., 25., 70.)
    ena_left = geo_idx(260., lng_ml)
    ena_right = geo_idx(295., lng_ml)
    lat_jet_ena = lat_jet_ml[:, ena_left:ena_right+1]
    lat_jet_err_ena = np.nanstd(np.nanmean(lat_jet_ena, axis = 1))
    lat_jet_mean_ena = np.nanmean(lat_jet_ena)
    land_ena = globalo3_calculate.find_grid_overland(lat_ena, lng_ena)
    field_err_ena = 2*np.nanstd(field_ena*land_ena, axis=1)
    field_ena = np.nanmean(field_ena*land_ena, axis=1)
    # Observations
    obs_ena = np.where((np.array(lat_naps+lat_aqs) < 70.) &
                         (np.array(lat_naps+lat_aqs) > 25.) &
                         (np.array(lng_naps+lng_aqs) > 260.) &
                         (np.array(lng_naps+lng_aqs) < 295.))[0]
    lat_obs_ena = np.array(lat_naps+lat_aqs)[obs_ena]
    obs_ena = np.array(obs_naps+obs_aqs)[obs_ena]
    lat_obs_ena, obs_mean_ena, obs_std_ena = \
        globalo3_calculate.bin_observations_bylat(lat_ena[::2], obs_ena, 
        lat_obs_ena)
    # Europe (since Europe crosses the prime meridian, finding grid cells over
    # 0 deg longitude will not work, so find European domain in two steps)
    # Model 
    field_eu1, lat_eu1, lng_eu1 = globalo3_calculate.find_grid_in_bb(
        field, lat_gmi, lng_gmi, 350., 360., 25., 70.) 
    eu1_left = geo_idx(350., lng_ml)
    eu1_right = geo_idx(360., lng_ml)
    lat_jet_eu1 = lat_jet_ml[:, eu1_left:eu1_right+1]
    land_eu1 = globalo3_calculate.find_grid_overland(lat_eu1, lng_eu1)
    field_eu2, lat_eu2, lng_eu2 = globalo3_calculate.find_grid_in_bb(
        field, lat_gmi, lng_gmi, 0., 30., 25., 70.) 
    eu2_left = geo_idx(0., lng_ml)
    eu2_right = geo_idx(30., lng_ml)
    lat_jet_eu2 = lat_jet_ml[:, eu2_left:eu2_right+1]
    land_eu2 = globalo3_calculate.find_grid_overland(lat_eu2, lng_eu2)
    # Concatenate
    field_eu = np.hstack([field_eu1, field_eu2])
    lat_jet_eu = np.hstack([lat_jet_eu1, lat_jet_eu2])
    lat_jet_err_eu = np.nanstd(np.nanmean(lat_jet_eu, axis = 1))
    lat_jet_mean_eu = np.nanmean(lat_jet_eu)
    land_eu = np.hstack([land_eu1, land_eu2])
    lat_eu = lat_eu1
    field_err_eu = 2*np.nanstd(field_eu*land_eu, axis=1)
    field_eu = np.nanmean(field_eu*land_eu, axis=1)
    # Observations
    emep_eu1 = np.where((np.array(lat_emep) < 70.) &
                        (np.array(lat_emep) > 25.) &
                        (np.array(lng_emep) > 350.) &
                        (np.array(lng_emep) < 360.))[0]
    emep_eu2 = np.where((np.array(lat_emep) < 70.) &
                        (np.array(lat_emep) > 25.) &
                        (np.array(lng_emep) > 0.) &
                        (np.array(lng_emep) < 30.))[0]
    emep_eu = np.hstack([emep_eu1, emep_eu2])
    obs_eu = np.array(obs_emep)[emep_eu]
    lat_obs_eu = np.array(lat_emep)[emep_eu]
    lat_obs_eu, obs_mean_eu, obs_std_eu = \
        globalo3_calculate.bin_observations_bylat(lat_eu[::2], obs_eu, 
        lat_obs_eu)
    # China
    # Model 
    field_asia, lat_asia, lng_asia = globalo3_calculate.find_grid_in_bb(
        field, lat_gmi, lng_gmi, 90., 125., 25., 70.) 
    asia_left = geo_idx(90., lng_ml)
    asia_right = geo_idx(125., lng_ml)
    lat_jet_asia = lat_jet_ml[:, asia_left:asia_right+1]
    lat_jet_err_asia = np.nanstd(np.nanmean(lat_jet_asia, axis = 1))
    lat_jet_mean_asia = np.nanmean(lat_jet_asia)
    land_asia = globalo3_calculate.find_grid_overland(lat_asia, lng_asia)
    field_err_asia = 2*np.nanstd(field_asia*land_asia, axis=1)
    field_asia = np.nanmean(field_asia*land_asia, axis=1)
    # Observations
    obs_asia = np.where((np.array(lat_china) < 70.) &
        (np.array(lat_china) > 25.) & (np.array(lng_china) > 90.) &
        (np.array(lng_china) < 125.))[0]
    lat_obs_asia = np.array(lat_china)[obs_asia]
    obs_asia = np.array(obs_china)[obs_asia]
    lat_obs_china, obs_mean_china, obs_std_china = \
        globalo3_calculate.bin_observations_bylat(lat_asia[::2], obs_asia, 
        lat_obs_asia)
    # Plotting
    fig = plt.figure(figsize=(6.5, 6.))
    ax1 = plt.subplot2grid((4,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((4,2), (1,0), colspan=2)
    ax3 = plt.subplot2grid((4,2), (2,0), colspan=2)
    ax4 = plt.subplot2grid((4,2), (3,0), colspan=2)
    color_r = '#EF9802' #'#A2C3E0'
    color_jet = 'k'
    for ax in [ax1, ax2, ax3, ax4]: 
        ax.set_xlim([25, 70])
        ax.set_xticks(np.linspace(25, 70, 10))
        ax.set_xticklabels([])
        ax.set_ylim([yticks[0], yticks[-1]])
        ax.set_yticks(yticks)
        # Add horizontal line for field = 0
        ax.axhline(y = 0.0, color = 'k', lw = 1., linestyle = '--')        
    ax4.set_xticklabels([ '25', '', '35', '', '45', '', '55', '', '65', ''])    
    ax4.set_xlabel('Latitude [$^{\mathregular{\circ}}$N]', fontsize = 14)
    ax2.set_ylabel('%s %s'%(label, units), fontsize = 14, y = -0.15)    
    # Western North America
    ax1.text(0.98, 0.8, 'Western North America', ha = 'right', 
             transform=ax1.transAxes, fontsize=12, zorder=20)
    p1, = ax1.plot(lat_wna, field_wna, ls='-', lw=2, color=color_r, 
             zorder=5)
    p2 = ax1.fill_between(lat_wna, field_wna-field_err_wna,
        field_wna+field_err_wna, color = color_r, alpha=0.2, zorder=5)
    p3, = ax1.plot(lat_obs_wna, obs_mean_wna, '-', lw=2, color='grey',     
        zorder=4)
    p4 = ax1.fill_between(lat_obs_wna, (obs_mean_wna-obs_std_wna),
        (obs_mean_wna+obs_std_wna), facecolor='None', hatch='////', 
        edgecolor='grey', linewidth=0.0, zorder=3)
    # Western North America eddy-driven jet
    p7 = ax1.errorbar(lat_jet_mean_wna, 0., xerr=[lat_jet_err_wna], 
        fmt = 'o', color = color_jet, ecolor = color_jet, elinewidth = 2, 
        capsize = 3, zorder = 20)
    # Add legend
    leg = ax1.legend([(p3, p4), (p1, p2), p7], 
        ['Observed $\mathregular{\pm}$ 2$\mathregular{\sigma}$',
         'Modeled $\mathregular{\pm}$ 2$\mathregular{\sigma}$',
         'Eddy-driven jet'], loc=6, ncol=3, bbox_to_anchor=(0.05, 1.3), 
         frameon=False)
    # Change the marker size manually for observations
    leg.legendHandles[0]._legmarker.set_markersize(5)
    # Eastern North America
    ax2.text(0.98, 0.8, 'Eastern North America', ha='right', 
        transform=ax2.transAxes, fontsize=12, zorder=20)
    ax2.plot(lat_ena, field_ena, ls='-', lw=2, color=color_r, zorder=5)
    ax2.fill_between(lat_ena, field_ena-field_err_ena, field_ena+field_err_ena, 
        color=color_r, alpha=0.2, zorder=5)
    ax2.plot(lat_obs_ena, obs_mean_ena, '-', lw=2, color='grey', zorder=4)
    ax2.fill_between(lat_obs_ena, (obs_mean_ena-obs_std_ena),
        (obs_mean_ena+obs_std_ena), facecolor='None', hatch='////', 
        edgecolor='grey', linewidth=0.0, zorder=3)
    # Eastern North America eddy-driven jet
    ax2.errorbar(lat_jet_mean_ena, 0., xerr=[lat_jet_err_ena], fmt='o', 
        color=color_jet, ecolor=color_jet, elinewidth=2, capsize=3, zorder=20)
    # Europe
    ax3.text(0.98, 0.8, 'Europe', ha='right', transform=ax3.transAxes, 
        fontsize=12, zorder=20)
    ax3.plot(lat_eu, field_eu, ls='-', lw=2, color=color_r, zorder=5)
    ax3.fill_between(lat_eu, field_eu-field_err_eu, field_eu+field_err_eu, 
        color=color_r, alpha=0.2, zorder=5)
    ax3.plot(lat_obs_eu, obs_mean_eu, '-', lw=2, color='grey', zorder=4)
    ax3.fill_between(lat_obs_eu, (obs_mean_eu-obs_std_eu),
        (obs_mean_eu+obs_std_eu), facecolor='None', hatch='////', 
        edgecolor='grey', linewidth=0.0, zorder=3)
    # Europe eddy-driven jet
    ax3.errorbar(lat_jet_mean_eu, 0., xerr=[lat_jet_err_eu], fmt='o', 
        color=color_jet, ecolor=color_jet, elinewidth=2, capsize=3, zorder=20)
    # East Asia
    ax4.text(0.98, 0.8, 'China', ha='right', transform=ax4.transAxes, 
        fontsize=12, zorder=20)
    ax4.plot(lat_asia, field_asia, ls='-', lw=2, color=color_r, zorder=5)
    ax4.fill_between(lat_asia, field_asia-field_err_asia,
        field_asia+field_err_asia, color=color_r, alpha=0.2, zorder=5)
    ax4.plot(lat_obs_china, obs_mean_china, '-', lw=2, color='grey', zorder=4)
    ax4.fill_between(lat_obs_china, (obs_mean_china-obs_std_china),
        (obs_mean_china+obs_std_china), facecolor='None', hatch='////', 
        edgecolor='grey', linewidth=0.0, zorder=3)
    # East Asia eddy-driven jet
    ax4.errorbar(lat_jet_mean_asia, 0., xerr=[lat_jet_err_asia], fmt='o', 
        color=color_jet, ecolor=color_jet, elinewidth=2, capsize=2, zorder=10)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
                'zonalavg_byregion_smoothedobs_%s.png' %fstr, dpi = 300,
                transparent = True)
    return

def rt2mo3_iav(): 
    """function loads JJA GMI CTM O3 from the HindcastMR2 simulation and 
    MERRA-2 2-meter temperatures for 2000-2010 and finds the mean summertime
    values for the fields. The correlation coefficient is calculated with these
    mean fields and plotted. 
    
    Parameters
    ----------
    None  

    Returns
    -------
    None    
    """
    import numpy as np
    import pandas as pd
    years_iav = np.arange(2000, 2011, 1)
    months_str = ['jun', 'jul', 'aug']
    # Hemisphere definition
    latmin, lngmin, latmax, lngmax = -1., 0., 90., 360.
    # Load Northern Hemisphere HindcastMR2 GMI CTM O3
    lat_gmi, lng_gmi, times_iav_gmi, o3_iav_gmi = \
        globalo3_open.open_overpass2_specifieddomain(list(years_iav), months_str, 
        latmin, latmax, lngmin, lngmax, 'O3', 'HindcastMR2')
    o3_iav_gmi = o3_iav_gmi*1e9
    # Load Northern Hemisphere 2-meter temperatures
    lat_merra, lng_merra, t2m_iav_merra = \
        globalo3_open.open_merra2t2m_specifieddomain(list(years_iav), months_str,  
        latmin, latmax, lngmin, lngmax)        
    # Interpolate 2-meter temperature
    t2m_iav_merra = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, lng_gmi, 
        lat_merra, lng_merra, t2m_iav_merra)
    season = 'JJA'
    times_iav_gmi = pd.to_datetime(times_iav_gmi)
    # Lists will be filled with mean JJA values for each year
    t2m_jja_mean, o3_jja_mean = [], []
    for year in years_iav:
        idx_ty = np.where(times_iav_gmi.year == year)[0]
        # Find O3 and temperature during summer of interest
        t2m_ty = np.nanmean(t2m_iav_merra[idx_ty], axis = 0)
        o3_ty = np.nanmean(o3_iav_gmi[idx_ty], axis = 0)
        # Append to list
        t2m_jja_mean.append(t2m_ty)
        o3_jja_mean.append(o3_ty)
    t2m_jja_mean = np.stack(t2m_jja_mean)    
    o3_jja_mean = np.stack(o3_jja_mean)    
    # Calculate correlation coefficient with JJA mean fields
    r_t2mo3_iav = globalo3_calculate.calculate_r(t2m_jja_mean, o3_jja_mean, 
        lat_gmi, lng_gmi)    
    # Plot
    map_hemisphere(lat_gmi, lng_gmi, r_t2mo3_iav, 
        r'$\overline{\mathregular{%s}}$ $\it{r}\:$(T, O$_\mathregular{3}$)' 
        %season, '[$\cdot$]', np.linspace(-1., 1, 9), 'bwr', 'nh', 
        'rt2mo3_iav_%s_%d-%d'%(season, years_iav[0], years_iav[-1]), 
        extent=[lng_gmi.min()-180., lng_gmi.max()-180., lat_gmi.min()+1, 
        lat_gmi.max()-5], extend = 'neither')
    return

def zonalavg_verticallyintegrated_meridional_flux(total, mean, eddy, lat, 
    title, fstr):
    """plot the total, mean, and eddy components of the zonally-averaged 
    meridional flux of tracer f over the Northern Hemisphere

    Parameters
    ----------
    total : numpy.ndarray
        Total vertically-integrated meridional flux of tracer f, units of kg 
        s-1, [lat,]
    mean : numpy.ndarray 
        Mean vertically-integrated meridional flux of tracer f, units of kg 
        s-1, [lat,]    
    eddy : numpy.ndarray
        Eddy vertically-integrated meridional flux of tracer f, units of kg 
        s-1, [lat,]    
    lat : numpy.ndarray
        Latitude coordinates, units of degrees north, [lat,]        
    title : str
        Title for plot        
    fstr : str
        Output file suffix; should specified levels of vertical integration 
        and tracer name
    
    Returns
    -------        
    None    
    """    
    import matplotlib.pyplot as plt
    # Plotting
    fig = plt.figure()
    ax = plt.subplot2grid((1,1),(0,0))
    ax.plot(lat, total, ls='-', color='#1b9e77', lw=2, label='Total', zorder=2)
    ax.plot(lat, mean, ls='-', color='#d95f02', lw=2, label='Mean', zorder=3)
    ax.plot(lat, eddy, ls='-', color='#7570b3', lw=2, label='Eddy', zorder=4)        
    # Aesthetics
    ax.set_xlim([0,90])
    ax.set_xticks([0,30,60,90])
    ax.set_xlabel('Latitude [$^{\circ}$N]', fontsize=14)
#    ax.set_ylabel('800-950 hPa Flux [kg s$^{\mathregular{-1}}$]', fontsize=14)
    ax.set_ylabel('Surface Flux [ppbv m s$^{\mathregular{-1}}$]', fontsize=14)    
    ax.set_title(title, fontsize=14)
    ax.hlines(0, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], zorder=1, 
        linestyles='--', linewidths=0.75)
    plt.legend(ncol=3, frameon=False)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'zonalavg_verticallyintegrated_meridional_%sflux.png' %fstr, dpi=300)
    return

def map_obso3ctmo3_performancemetric(field_aqs, field_naps, field_emep, 
    field_china, lat_field, lng_field, cbar_label, clevs, cmap, extend, 
    fstr):
    """plot gridded modeled O3-observed O3 performance metrics over North 
    America (NAPS, AQS), Europe (EMEP), and China (Li et al., 2019).
    
    Parameters
    ----------
    field_aqs : numpy.ndarray
        Gridded performance metric from AQS, [lat, lng]
    field_naps : numpy.ndarray
        Gridded performance metric from NAPS, [lat, lng]
    field_emep : numpy.ndarray
        Gridded performance metric from EMEP, [lat, lng]
    field_china : numpy.ndarray
        Gridded performance metric from China, [lat, lng]        
    lat_field : numpy.ndarray
        Modeled latitude coordinates, units of degrees north, [lat,]        
    lng_field : numpy.ndarray
        Modeled longitude coordinates, units of degrees east, [lng,]
    cbar_label : str
        Label for the colorbar (field and units)
    clevs : numpy.ndarray
        Filled contour levels values
    cmap : str
        Colormap name
    extend : str
        Extend settings for matplotlib colormap/colorbar (i.e., 'neither', 
        'both', 'min', 'max')
    fstr : str
        Output filename suffix

    Returns
    -------
    None 
    """
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = 0.3     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure(figsize=(5.5,6))
    # Enforce equal size (height) of subplots 
    # stackoverflow.com/questions/37767026/matplotlib-enforce-equal-size-
    # height-of-subplots
    ax1 = plt.subplot2grid((2,2),(0,0), colspan=2,
        projection=ccrs.Miller(central_longitude=0.))
    ax2 = plt.subplot2grid((2,2),(1,0), aspect='auto', adjustable='box',
        projection=ccrs.Miller(central_longitude=0.))
    ax3 = plt.subplot2grid((2,2),(1,1), aspect='auto', adjustable='box',
        projection=ccrs.Miller(central_longitude=0.))
    # Set subplot labels and titles
    ax1.set_title('(a) NAPS, AQS', fontsize=15, x=0.02, ha='left')  
    ax2.set_title('(b) EMEP', fontsize=15, x=0.02, ha='left')
    ax3.set_title('(c) MEE', fontsize=15, x=0.02, ha='left')
    # Add oceans and coasts
    for ax in [ax1, ax2, ax3]:
        ax.add_feature(cfeature.OCEAN, zorder=2, lw=0.0, color='lightgrey')
        ax.coastlines(lw=0.25, color='k', zorder=3)
    # Set extent
    ax1.set_extent([-165, -50, 24, 55])    
    ax2.set_extent([-14, 37, 34, 70])
    ax3.set_extent([95, 139, 15, 53])
    # Define colormap
    cmap = plt.get_cmap(cmap)
    # Extract all colors from the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)
    norm = mpl.colors.BoundaryNorm(clevs, cmap.N)
    # North America
    mb = ax1.pcolor(lng_field, lat_field, field_aqs, cmap=cmap, norm=norm, 
        transform=ccrs.PlateCarree(), zorder=20)
    ax1.pcolor(lng_field, lat_field, field_naps, cmap=cmap, norm=norm, 
        transform=ccrs.PlateCarree(), zorder=20)
    # Europe
    ax2.pcolor(lng_field, lat_field, field_emep, cmap=cmap, norm=norm, 
        transform=ccrs.PlateCarree(), zorder=20)
    # China
    ax3.pcolor(lng_field, lat_field, field_china, cmap=cmap, norm=norm, 
        transform=ccrs.PlateCarree(), zorder=20)
    # Add colorbar
    colorbar_axes = plt.gcf().add_axes([0.83,0.15,0.03,0.7])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        extend=extend)
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label(cbar_label, fontsize=15)
    plt.gcf().subplots_adjust(left=0.04, right=0.78)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'map_obso3ctmo3_%s.png'%fstr, dpi=300)
    return

def map_nh_rpblho3_diurnalcycle(ilat, ilng, lat_model, lng_model, pblh_daily, 
    pblh_hrly):
    """Plot map of the correlation coefficient calculated between mean daily 
    PBL height across the Northern Hemisphere with hatching indicating in-
    significant correlations (top subplot) and the diurnal cycle of PBL height
    and O3 for the specified grid cell on days with high PBL heights (70th 
    percentile) and low PBL heights (30th percentile) (bottom subplot).

    Parameters
    ----------
    ilat : float 
        The latitude of interest for which the PBLH and O3 from AQS stations 
        will be plotted in the bottom subplot
    ilng : float 
        The longitude of interest for which the PBLH from MERRA and O3 from AQS 
        stations will be plotted in the bottom subplot        
    lat_model    
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_model : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    pblh_daily : numpy.ndarray
        Daily mean PBL height, units of m, [time, lat, lng]        
    pblh_hrly : numpy.ndarray
        Hourly PBL height, units of m, [time, lat, lng]    

    Returns
    -------
    None    
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from timezonefinder import TimezoneFinder
    import pytz
    import datetime
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = 0.3     
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx
    PATH_AQS = '/Users/ghkerr/phd/observations/o3/AQS/'
    dtype = {'State Code' : np.str, 'County Code' : np.str,
           'Site Num' : np.str, 'Latitude' : np.float64, 
           'Longitude' : np.float64, 'Date Local' : np.str,
           'Time Local' : np.str, 'Sample Measurement' : np.float}
    O3files = []
    for year in years:
        # Fetch file names for measuring period    
        O3files.append(PATH_AQS+'hourly_44201_%s_reduced.csv' %year)
    O3files.sort()
    # Read multiple CSV files (yearly) into Pandas dataframe 
    aqs_allyr = observations_open.get_merged_csv(O3files, dtype=dtype, index_col = None, 
                          usecols = list(dtype.keys())) 
    # Select months in measuring period     
    aqs_allyr = aqs_allyr.loc[pd.to_datetime(aqs_allyr['Date Local']).
                         dt.month.isin(months)]    
    aqs_allyr['Station ID'] = np.nan
    # Select hours of interest
    hours_str = ['%s:00'%str(x).zfill(2) for x in np.arange(0,24,1)]
    aqs_allyr = aqs_allyr.loc[aqs_allyr['Time Local'].isin(hours_str)]
    # Combine State Code, County Code, and Site Number into single site 
    # ID column 
    aqs_allyr['Station ID'] = (aqs_allyr['State Code'] + '-' +
        aqs_allyr['County Code'] + '-' + aqs_allyr['Site Num'])
    # Rename O3 column and convert from ppm to ppb
    aqs_allyr = aqs_allyr.rename(columns=
                                 {'Sample Measurement' : 'O3 (ppbv)'})
    aqs_allyr['O3 (ppbv)'] = aqs_allyr['O3 (ppbv)']*1000.
    # Convert longitude from (-180-180) to (0-360)
    aqs_allyr['Longitude'] = aqs_allyr['Longitude']%360.
    # Drop individual code columns
    del aqs_allyr['State Code'], aqs_allyr['County Code'], \
        aqs_allyr['Site Num']
    aqs_allyr.index = pd.to_datetime(aqs_allyr['Date Local'])
    # Find grid cell corresponding to city/point of interest and PBLH 
    # in grid cells 
    ilat = geo_idx(ilat, lat_model)
    ilng = geo_idx(ilng, lng_model)
    # Based on O3 observations from AQS (above) find O3 near to grid cell 
    # where "near" is defined by the searchradius; n.b., this assumes that 
    # the resolution of the model is ~1 deg x ~1 deg
    searchradius = 2
    ipblh_daily = pblh_daily[:,ilat-searchradius:ilat+searchradius,
        ilng-searchradius:ilng+searchradius]
    ipblh_daily = np.nanmean(ipblh_daily, axis=tuple((1,2)))
    ipblh_hrly = pblh_hrly[:,ilat-searchradius:ilat+searchradius,
        ilng-searchradius:ilng+searchradius]
    ipblh_hrly = np.nanmean(ipblh_hrly, axis=tuple((1,2)))
    iaqs_hrly = aqs_allyr.loc[
        (aqs_allyr['Latitude'] > (lat_model[ilat]-searchradius)) &
        (aqs_allyr['Latitude'] < (lat_model[ilat]+searchradius)) &
        (aqs_allyr['Longitude'] > (lng_model[ilng]-searchradius)) &
        (aqs_allyr['Longitude'] < (lng_model[ilng]+searchradius))]
    # Days that have PBLH in the top and bottom 30th percentiles of all days
    # (using daily mean )
    highpblhidx = np.where(ipblh_daily > np.percentile(ipblh_daily, 70))[0]
    lowpblhidx = np.where(ipblh_daily < np.percentile(ipblh_daily, 30))[0]
    highpblhtime = mtime[highpblhidx]
    lowpblhtime = mtime[lowpblhidx]
    # Find O3 observations on days with high/low PBLH
    highpblh_iaqs_hrly = iaqs_hrly.loc[highpblhtime]
    lowpblh_iaqs_hrly = iaqs_hrly.loc[lowpblhtime]
    # Group AQS O3 by hour (column 'Time Local') and average to produce diel 
    # cycle
    highpblh_iaqs_hrly = highpblh_iaqs_hrly.groupby(['Time Local']).mean()
    lowpblh_iaqs_hrly = lowpblh_iaqs_hrly.groupby(['Time Local']).mean()
    # Find GMI O3 on days with high/low PBLH
    highpblh_o3_model = o3_model[highpblhidx, ilat, ilng].mean()
    lowpblh_o3_model = o3_model[lowpblhidx, ilat, ilng].mean()
    # Hourly PBLH on days with high/low PBLH 
    highpblh_hrly, lowpblh_hrly = [], []
    for idx in highpblhidx:
        highpblh_hrly.append(pblh_hrly[idx*24:idx*24+24, ilat, ilng])
    for idx in lowpblhidx:
        lowpblh_hrly.append(pblh_hrly[idx*24:idx*24+24, ilat, ilng])
    lowpblh_hrly = np.stack(lowpblh_hrly)
    lowpblh_hrly = np.nanmean(lowpblh_hrly, axis=0)
    highpblh_hrly = np.stack(highpblh_hrly)
    highpblh_hrly = np.nanmean(highpblh_hrly, axis=0)
    # Shift diel cycle of PBLH based on location of grid cell (UTC offset)
    tf = TimezoneFinder()
    timezone_str = tf.certain_timezone_at(lat=lat_model[ilat], 
        lng=(lng_model[ilng]+180)%360-180)
    tz = pytz.timezone(timezone_str)
    # Dummy summer date (to pick up on DST)
    dt = datetime.datetime.strptime('2010-06-01', '%Y-%m-%d')
    offset = int(tz.utcoffset(dt).total_seconds() / 3600.0)
    # Roll PBLH diel cycle by UTC offset 
    highpblh_hrly = np.roll(highpblh_hrly, offset)
    lowpblh_hrly = np.roll(lowpblh_hrly, offset)
    # Plotting
    fig = plt.figure()
    ax=plt.subplot2grid((2,2), (0,0), colspan=2, aspect="auto",
        projection=ccrs.Miller(central_longitude=0.))
    ax.set_title(r'(a) $\it{r}\:$(PBLH, O$_{\mathregular{3}}$)', fontsize=16, 
        x=0.02, ha='left')    
    ax.add_feature(cfeature.OCEAN, zorder=2, lw=0.0, color='lightgrey')
    ax.coastlines(lw=0.25, color='k', zorder=3)
    ax.set_extent([lng_model.min()-180., lng_model.max()-180., 
        lat_model.min()+1, lat_model.max()-5])
    # Plot filled contours
    cmap = plt.get_cmap('bwr')
    mb = ax.contourf(lng_model, lat_model, r_pblho3, np.linspace(-1., 1, 9), 
        cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance
    ax.contourf(lng_model, lat_model, significance_r_pblho3, hatches=['//////'], 
        colors='none', transform=ccrs.PlateCarree())
    # Draw box to indicate region over which O3 observations/PBL height is
    # examed 
    x = [lng_model[ilng]-searchradius, lng_model[ilng]+searchradius, 
        lng_model[ilng]+searchradius, lng_model[ilng]-searchradius, 
        lng_model[ilng]-searchradius]
    y = [lat_model[ilat]-searchradius, lat_model[ilat]-searchradius, 
        lat_model[ilat]+searchradius, lat_model[ilat]+searchradius, 
        lat_model[ilat]-searchradius]
    ax.fill(x, y, transform=ccrs.PlateCarree(),facecolor="none",
            edgecolor="k")
    # PBLH diurnal cycle
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=2)
    ax2.set_title('(b)', fontsize=16, x=0.02, ha='left')    
    ln1 = ax2.plot(highpblh_iaqs_hrly['O3 (ppbv)'].values, ls='-', lw=2., 
        color='#e41a1c', 
        label='O$_{\mathregular{3}}$(PBLH$_{\mathregular{70}}$)')
    ln2 = ax2.plot(lowpblh_iaqs_hrly['O3 (ppbv)'].values, ls='-', lw=2., 
        color='#377eb8', 
        label='O$_{\mathregular{3}}$(PBLH$_{\mathregular{30}}$)')
    ax2.set_ylim([10, 60])
    ax2.set_yticks([10, 20, 30, 40, 50, 60])
    ax2.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize=16)
    # PBL height 
    ax2b = ax2.twinx()
    ln3 = ax2b.plot(highpblh_hrly, ls='--', lw=2., color='#e41a1c', 
        label='PBLH$_{\mathregular{70}}$')
    ln4 = ax2b.plot(lowpblh_hrly, ls='--', lw=2., color='#377eb8',
        label='PBLH$_{\mathregular{30}}$')          
    ax2b.set_xlim([0, 23])
    ax2b.set_xticks([0, 3, 6, 9, 12, 15, 18, 21])
    ax2b.set_xticklabels(['', '3', '6', '9', '12', '15', '18', '21'])
    ax2b.set_ylim([0, 2300])
    ax2b.set_yticks(np.linspace(0, 2300, 6))
    ax2b.set_ylabel('PBLH [m]', fontsize=16, rotation=270, position=(1.1, 0.5))
    ax2b.yaxis.set_label_coords(1.16, 0.5)
    ax2.set_xlabel('Hour [LT]', fontsize=16)
    # Add legend
    lns = ln1+ln2+ln3+ln4
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc=0, frameon=False, fontsize=8)
    plt.gcf().subplots_adjust(left=0.1, right=0.8, hspace=0.3)
    # Add colorbar
    colorbar_axes = plt.gcf().add_axes([ax.get_position().x1+0.04,
        ax.get_position().y0,
        0.02,
        ax.get_position().y1-ax.get_position().y0]) 
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        extend='neither')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[$\mathregular{\cdot}$]', fontsize=16)
    ax.outline_patch.set_zorder(20)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'map_nh_rpblho3_diurnalcycle.png', dpi = 450)
    return 
 
def map_o3anomcyclones_pweqjet(lat_model, lng_model, o3_model, o3_pwjet, 
    o3_eqjet, lat_binned, lng_binned, pwjet_cyclones_binned, 
    eqjet_cyclones_binned, pwjet_lat, pwjet_lat_var, eqjet_lat, eqjet_lat_var): 
    """plot maps of the O3 anomaly (filled contours) calculated as the 
    difference of O3 on days with a poleward/equatorward jet and all 
    days and the cyclone frequency (pcolor) on days with a poleward/equatorward
    jet. The mean position and variability of the poleward and equatorward 
    jets are also indicated. 

    Parameters
    ----------
    lat_model : numpy.ndarray
        Latitude coordinates for the Northern Hemisphere, units of degrees 
        north, [lat,]                
    lng_model : numpy.ndarray
        Longitude coordinates for the Northern Hemisphere, units of degrees 
        east, [lng,]          
    o3_model : numpy.ndarray
        Modeled O3 for all days, units of ppbv, [time, lat, lng]
    o3_pwjet : numpy.ndarray
        Modeled O3 for days where the jet is in a poleward position, units of
        ppbv, [lat, lng]
    o3_eqjet : numpy.ndarray
        Modeled O3 for days where the jet is in a poleward position, units of
        ppbv, [lat, lng]
    lat_binned : numpy.ndarray
        Latitude coordinates for binned cyclone dataset, units of degrees
        north, [lat_b,]        
    lng_binned : numpy.ndarray
        Longitude coordinates for binned cyclone dataset, units of degrees
        east, [lng_b,]
    pwjet_cyclones_binned : numpy.ndarray
        Binned cyclone frequency on days with a poleward jet, [lat_b, lng_b]
    eqjet_cyclones_binned : numpy.ndarray
        Binned cyclone frequency on days with a equatorward jet, [lat_b, 
        lng_b]    
    pwjet_lat : numpy.ndarray
        Mean latitude of the eddy-driven jet on days when the jet is in a 
        poleward position, units of degrees north, [lng,]
    pwjet_lat_var : numpy.ndarray
        Variability of the eddy-driven jet on days when the jet is in a 
        poleward position, units of degrees, [lng,]
    eqjet_lat : numpy.ndarray
        Mean latitude of the eddy-driven jet on days when the jet is in an
        equatorward position, units of degrees north, [lng,]
    eqjet_lat_var : numpy.ndarray
        Variability of the eddy-driven jet on days when the jet is in an 
        equatorward position, units of degrees, [lng,]
    """
    import numpy as np
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = 0.3     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure(figsize=(8,5))
    ax1=plt.subplot2grid((2,2), (0,0), colspan=2,
                        projection=ccrs.Miller(central_longitude=0.))
    ax2=plt.subplot2grid((2,2), (1,0), colspan=2,
                        projection=ccrs.Miller(central_longitude=0.))
    ax1.set_title('(a) Poleward jet', fontsize=16, x=0.02, ha='left')    
    ax2.set_title('(b) Equatorward jet', fontsize=16, x=0.02, ha='left')    
    for ax in [ax1, ax2]:
        ax.add_feature(cfeature.OCEAN, zorder=2, lw=0.0, color='lightgrey')
        ax.coastlines(lw=0.25, color='k', zorder=3)
        ax.set_extent([lng_model.min()-180., lng_model.max()-180., 
            lat_model.min()+1, lat_model.max()-5])
    # Plot difference in O3 between poleward and equatorward jet days
    cmap = plt.get_cmap('bwr')
    mb = ax1.contourf(lng_model, lat_model, o3_pwjet-np.nanmean(o3_model, axis=0),
        np.linspace(-4, 4, 9), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_model, lat_model, o3_eqjet-np.nanmean(o3_model, axis=0),
        np.linspace(-4, 4, 9), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    # Cyclone anomaly calculated as difference in cyclones on all days 
    # versus cyclones on days with a poleward jet
    pwjet_cyclones_binned_ma = np.ma.masked_array(pwjet_cyclones_binned, 
        pwjet_cyclones_binned < 4)
    pcolor_clevs = np.linspace(0, 20, 6)
    mb2 = ax1.pcolor(lng_binned, lat_binned, pwjet_cyclones_binned_ma, 
        cmap=plt.get_cmap('binary', len(pcolor_clevs)-1), 
        vmin=pcolor_clevs[0], vmax=pcolor_clevs[-1], alpha=0.5, 
        transform=ccrs.PlateCarree(), zorder=5)
    eqjet_cyclones_binned_ma = np.ma.masked_array(eqjet_cyclones_binned, 
        eqjet_cyclones_binned < 5)
    ax2.pcolor(lng_binned, lat_binned, eqjet_cyclones_binned_ma,
        cmap=plt.get_cmap('binary', len(pcolor_clevs)-1), 
        vmin=pcolor_clevs[0], vmax=pcolor_clevs[-1], alpha=0.5, 
        transform=ccrs.PlateCarree(), zorder=5)
    # Poleward and equatorward jet and variability
    skiplng = 6
    ax1.errorbar(lng_model[::skiplng], pwjet_lat[::skiplng], 
        yerr=pwjet_lat_var[::skiplng], zorder=12, color='k', markersize=2, 
        elinewidth=0.5, ecolor='k', fmt='o', transform=ccrs.PlateCarree())
    ax2.errorbar(lng_model[::skiplng], eqjet_lat[::skiplng], 
        yerr=eqjet_lat_var[::skiplng], zorder=12, color='k', markersize=2, 
        elinewidth=0.5, ecolor='k', fmt='o', transform=ccrs.PlateCarree())
    # Add colorbars; colorbar on left (right) corresponds to cyclone (O3) 
    # anomaly. Align colorbars with top and bottom of subplots 
    colorbar_axes = plt.gcf().add_axes([0.87, ax2.get_position().y0, 0.02,
        (ax1.get_position().y1-ax2.get_position().y0)])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('$\mathregular{\delta}$ O$_\mathregular{3}$ [ppbv]', 
        fontsize=14)
    colorbar_axes = plt.gcf().add_axes([0.1, ax2.get_position().y0, 0.02,
        (ax1.get_position().y1-ax2.get_position().y0)])
    colorbar = plt.colorbar(mb2, colorbar_axes, orientation='vertical', 
        extend='both', ticks=pcolor_clevs)
    # Shift ticks, label to left hand side
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('Cyclones', fontsize=14)
    colorbar_axes.yaxis.set_label_position('left')
    colorbar_axes.yaxis.set_ticks_position('left')
    plt.gcf().subplots_adjust(left=0.12, right=0.88)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'map_o3anomcyclones_pweqjet.png', dpi = 350)
    return

import numpy as np
import sys
sys.path.append('/Users/ghkerr/phd/globalo3/')
import globalo3_open, globalo3_calculate, observations_open
sys.path.append('/Users/ghkerr/phd/transporto3/')
import transporto3_open
# Load data the first iteration only 
try:lat_gmi
except NameError:
    years = [2008, 2009, 2010]
    hours = [14]
    season = 'JJA'
    months = [6, 7, 8]
    months_str = ['jun', 'jul', 'aug']
    # Hemisphere definition
    latmin, lngmin, latmax, lngmax = -1., 0., 90., 360.
    # Load Northern Hemisphere HindcastMR2 GMI CTM O3 
    lat_gmi, lng_gmi, times_gmi, o3_gmi = \
        globalo3_open.open_overpass2_specifieddomain(years, months_str, latmin, 
        latmax, lngmin, lngmax, 'O3', 'HindcastMR2')
    o3_gmi = o3_gmi*1e9
    # Load Northern Hemisphere HindcastMR2-DiurnalAvgT GMI CTM O3 
    lat_gmi, lng_gmi, times_gmi, o3_transport_gmi = \
        globalo3_open.open_overpass2_specifieddomain(years, months_str, latmin, 
        latmax, lngmin, lngmax, 'O3', 'HindcastMR2-DiurnalAvgTQ')
    o3_transport_gmi = o3_transport_gmi*1e9    
    # Load MERRA-2 Northern Hemisphere 2-meter temperatures, specific humidity
    lat_merra, lng_merra, t2m_merra = \
        globalo3_open.open_merra2t2m_specifieddomain(years, months_str, latmin, 
        latmax, lngmin, lngmax)        
    qv2m_merra, mtime, lat_merra, lng_merra = transporto3_open.open_merra2(
        years, [0, 3, 9, 12, 15, 18, 21], 'QV2M', 'tavg1_2d_slv_Nx', 
        'JJA_rh.nc', lngmin, latmax, lngmax, latmin, dailyavg='yes')
    # Interpolate MERRA-2 fields
    t2m_merra = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
        lng_gmi, lat_merra, lng_merra, t2m_merra)
    qv2m_merra = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
        lng_gmi, lat_merra, lng_merra, qv2m_merra)   
    qv2m_merra = qv2m_merra*1000. # Convert from kg kg-1 to g kg-1
    # Load Northern Hemisphere 500 hPa MERRA-2 data
    U500, mtime, lat_merra, lng_merra = \
        globalo3_open.open_merra2_specifieddomain(years, months_str,
        [0,3,6,9,12,15,18,21], 'U', 'inst3_3d_asm_Np_500hPa', lngmin, latmax, 
        lngmax, latmin, dailyavg='yes')
    # Interpolate 500 hPa u-wind
    U500 = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
        lng_gmi, lat_merra, lng_merra, U500, checkplot='yes')
    # Load 10-meter U and V wind and interpolate to native CTM resolution
    U10M, mtime, lat_merra, lng_merra = transporto3_open.open_merra2(years, 
        np.arange(0, 24, 1), 'U10M', 'tavg1_2d_slv_Nx', 'JJA_wind.nc', lngmin, 
        latmax, lngmax, latmin, dailyavg='yes')
    U10M = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, lng_gmi, 
        lat_merra, lng_merra, U10M)
    V10M, mtime, lat_merra, lng_merra = transporto3_open.open_merra2(years, 
        np.arange(0, 24, 1), 'V10M', 'tavg1_2d_slv_Nx', 'JJA_wind.nc', lngmin, 
        latmax, lngmax, latmin, dailyavg='yes')
    V10M = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, lng_gmi, 
        lat_merra, lng_merra, V10M)
    # Calculate dO3/dT, dO3/dq, r(T, O3), and r(q, O3) from model
    do3dt2m = globalo3_calculate.calculate_do3dt(t2m_merra, o3_gmi, lat_gmi, 
        lng_gmi)
    do3dq = globalo3_calculate.calculate_do3dt(qv2m_merra, o3_gmi, lat_gmi, 
        lng_gmi)
    r_t2mo3 = globalo3_calculate.calculate_r(t2m_merra, o3_gmi, lat_gmi, 
        lng_gmi)
    r_t2mo3_transport = globalo3_calculate.calculate_r(t2m_merra, 
        o3_transport_gmi, lat_gmi, lng_gmi)
    r_qv2mo3 = globalo3_calculate.calculate_r(qv2m_merra, o3_gmi, lat_gmi, 
        lng_gmi)
    r_qv2mo3_transport = globalo3_calculate.calculate_r(qv2m_merra, 
        o3_transport_gmi, lat_gmi, lng_gmi)
#    r_U10Mo3 = globalo3_calculate.calculate_r(U10M, o3_gmi, lat_gmi, lng_gmi)
#    r_V10Mo3 = globalo3_calculate.calculate_r(V10M, o3_gmi, lat_gmi, lng_gmi)
    # Subset fields in mid-latitudes
    U500_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(U500, lat_gmi, 
        lng_gmi, 0., 360., 20., 70.)
    o3_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(o3_gmi, lat_gmi, 
        lng_gmi, 0., 360., 20., 70.)
    t2m_ml, lat_ml, lng_nhml = globalo3_calculate.find_grid_in_bb(t2m_merra, 
        lat_gmi, lng_gmi, 0., 360., 20., 70.)
    qv2m_ml, lat_ml, lng_nhml = globalo3_calculate.find_grid_in_bb(qv2m_merra, 
        lat_gmi, lng_gmi, 0., 360., 20., 70.)    
    do3dt2m_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(do3dt2m, 
        lat_gmi, lng_gmi, 0., 360., 20., 70.)
    r_t2mo3_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(r_t2mo3, 
        lat_gmi, lng_gmi, 0., 360., 20., 70.)
    # Find latitude of eddy-driven jet and fields at jet
    lat_jet_ml, o3_jet_ml = globalo3_calculate.find_field_atjet(o3_ml, U500_ml, 
        lat_ml, lng_ml, 20, anom = True)
    scratch, t2m_jet_ml = globalo3_calculate.find_field_atjet(t2m_ml, 
        U500_ml, lat_ml, lng_ml, 20, anom=True)
    scratch, do3dt2m_jet_ml = globalo3_calculate.find_field_atjet(
        do3dt2m_ml, U500_ml, lat_ml, lng_ml, 20)
    scratch, r_t2mo3_jet_ml = globalo3_calculate.find_field_atjet(
        r_t2mo3_ml, U500_ml, lat_ml, lng_ml, 20) 
    # Smooth jet latitude with window size of ~10 deg
    lat_jet_ml = globalo3_calculate.convolve_jet(lat_jet_ml, 9)    
    # Slope and correlation of O3/jet distance and 2-meter temperature and 
    # jet distance
    m_o3jetdist, r_o3jetdist, diff_o3jetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(o3_gmi, lat_gmi, 
        lng_gmi, lat_jet_ml, lng_ml)
    m_t2mjetdist, r_t2mjetdist, diff_t2mjetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(t2m_merra, lat_gmi, 
        lng_gmi, lat_jet_ml, lng_ml)
    m_qv2mjetdist, r_qv2mjetdist, diff_qv2mjetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(qv2m_merra, lat_gmi, 
        lng_gmi, lat_jet_ml, lng_ml)
    m_U10Mjetdist, r_U10Mjetdist, diff_U10Mjetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(U10M, lat_gmi, 
        lng_gmi, lat_jet_ml, lng_ml)
    m_V10Mjetdist, r_V10Mjetdist, diff_V10Mjetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(V10M, lat_gmi, 
        lng_gmi, lat_jet_ml, lng_ml)  
    # Load observational datasets
    promptobs = input('Load observational datasets? [y/n]\n').lower()    
    if promptobs == 'y':    
        naps = observations_open.open_napso3(years, months, hours)
        aqs = observations_open.open_aqso3(years, months, hours)
        emep = observations_open.open_emepo3(years, months, hours)
        china = observations_open.open_chinao3([2016,2017], months, hours, 
            cityavg=True)
        # Since Chinese observations span a different measuring period than 
        # the other observations, load different GMI and MERRA-2 data
        lat_gmi_china, lng_gmi_china, times_gmi_china, o3_gmi_china = \
            globalo3_open.open_overpass2_specifieddomain([2016,2017], 
            months_str, latmin, latmax, lngmin, lngmax, 'O3', 'HindcastMR2')
        o3_gmi_china = o3_gmi_china*1e9
        lat_merra_china, lng_merra_china, t2m_merra_china = \
            globalo3_open.open_merra2t2m_specifieddomain([2016,2017], 
            months_str, latmin, latmax, lngmin, lngmax)        
        t2m_merra_china = globalo3_open.interpolate_merra_to_ctmresolution(
            lat_gmi, lng_gmi, lat_merra_china, lng_merra_china, 
            t2m_merra_china)
        qv2m_merra_china, mtime_china, lat_merra_china, lng_merra_china = \
            transporto3_open.open_merra2([2016,2017], 
            [0, 3, 9, 12, 15, 18, 21], 'QV2M', 'tavg1_2d_slv_Nx', 'JJA_rh.nc', 
            lngmin, latmax, lngmax, latmin, dailyavg='yes')
        qv2m_merra_china = globalo3_open.interpolate_merra_to_ctmresolution(
            lat_gmi, lng_gmi, lat_merra_china, lng_merra_china, 
            qv2m_merra_china)
        qv2m_merra_china = qv2m_merra_china*1000. # Convert from kg kg-1 to g kg-1        
        import pandas as pd
        times_gmi_china = pd.date_range(start='01/01/2016', end='12/31/2017')
        where_months = np.where(np.in1d(times_gmi_china.month, 
            np.array(months))==True)[0]
        times_gmi_china = times_gmi_china[where_months]
        times_gmi_china = [x.to_pydatetime().date() for x in times_gmi_china]
        # Calculate O3-meteorology relationships from observational datasets
        (r_t2mo3_naps, do3dt2m_naps, ro3jet_naps, do3djet_naps, lat_naps, 
            lng_naps) = globalo3_calculate.calculate_obs_o3_temp_jet(naps, 
            t2m_merra, times_gmi, lat_gmi, lng_gmi, lat_jet_ml, lng_ml)
        (r_t2mo3_aqs, do3dt2m_aqs, ro3jet_aqs, do3djet_aqs, lat_aqs, 
            lng_aqs) = globalo3_calculate.calculate_obs_o3_temp_jet(aqs, 
            t2m_merra, times_gmi, lat_gmi, lng_gmi, lat_jet_ml, lng_ml)    
        (r_t2mo3_emep, do3dt2m_emep, ro3jet_emep, do3djet_emep, lat_emep, 
            lng_emep) = globalo3_calculate.calculate_obs_o3_temp_jet(emep, 
            t2m_merra, times_gmi, lat_gmi, lng_gmi, lat_jet_ml, lng_ml)
        (r_qv2mo3_naps, do3dqv2m_naps, ro3jet_naps, do3djet_naps, lat_naps, 
            lng_naps) = globalo3_calculate.calculate_obs_o3_temp_jet(naps, 
            qv2m_merra, times_gmi, lat_gmi, lng_gmi, lat_jet_ml, lng_ml)
        (r_qv2mo3_aqs, do3dqv2m_aqs, ro3jet_aqs, do3djet_aqs, lat_aqs, 
            lng_aqs) = globalo3_calculate.calculate_obs_o3_temp_jet(aqs, 
            qv2m_merra, times_gmi, lat_gmi, lng_gmi, lat_jet_ml, lng_ml)    
        (r_qv2mo3_emep, do3dqv2m_emep, ro3jet_emep, do3djet_emep, lat_emep, 
            lng_emep) = globalo3_calculate.calculate_obs_o3_temp_jet(emep, 
            qv2m_merra, times_gmi, lat_gmi, lng_gmi, lat_jet_ml, lng_ml)        
        # n.b., I haven't downloaded 2016-2017 MERRA-2 data yet, so the 
        # jet-O3-temperature connections for China is made up
        (r_t2mo3_china, do3dt2m_china, ro3jet_china, do3djet_china, lat_china,
            lng_china) = globalo3_calculate.calculate_obs_o3_temp_jet(china, 
            t2m_merra_china, np.array(times_gmi_china), lat_gmi, lng_gmi, 
            lat_jet_ml[:len(t2m_merra_china)], lng_ml)
        (r_qv2mo3_china, do3dqv2m_china, ro3jet_china, do3djet_china, 
            lat_china, lng_china) = \
            globalo3_calculate.calculate_obs_o3_temp_jet(china, 
            qv2m_merra_china, np.array(times_gmi_china), lat_gmi, lng_gmi, 
            lat_jet_ml[:len(t2m_merra_china)], lng_ml)
#        # Find bias 
#        r_naps_bias = globalo3_calculate.ctm_obs_bias(lat_naps, lng_naps, r_naps, 
#            lat_gmi, lng_gmi, r_t2mo3)
#        do3dt2m_naps_bias = globalo3_calculate.ctm_obs_bias(lat_naps, lng_naps, 
#            do3dt2m_naps, lat_gmi, lng_gmi, do3dt2m)
#        r_aqs_bias = globalo3_calculate.ctm_obs_bias(lat_aqs, lng_aqs, r_aqs, 
#            lat_gmi, lng_gmi, r_t2mo3)
#        do3dt2m_aqs_bias = globalo3_calculate.ctm_obs_bias(lat_aqs, lng_aqs, 
#            do3dt2m_aqs, lat_gmi, lng_gmi, do3dt2m)
#        r_emep_bias = globalo3_calculate.ctm_obs_bias(lat_emep, lng_emep, r_emep, 
#            lat_gmi, lng_gmi, r_t2mo3)
#        do3dt2m_emep_bias = globalo3_calculate.ctm_obs_bias(lat_emep, lng_emep, 
#            do3dt2m_emep, lat_gmi, lng_gmi, do3dt2m)
        # Calculate the air quality performance index for each observational dataset
        aqpi_aqs, fac2_aqs, r_aqs, mfb_aqs = globalo3_calculate.calculate_aqpi(aqs, 
            o3_gmi, lat_gmi, lng_gmi, times_gmi, 'gridwise')
        aqpi_naps, fac2_naps, r_naps, mfb_naps = globalo3_calculate.calculate_aqpi(
            naps, o3_gmi, lat_gmi, lng_gmi, times_gmi, 'gridwise')
        aqpi_emep, fac2_emep, r_emep, mfb_emep = globalo3_calculate.calculate_aqpi(
            emep, o3_gmi, lat_gmi, lng_gmi, times_gmi, 'gridwise')
        aqpi_china, fac2_china, r_china, mfb_china = \
            globalo3_calculate.calculate_aqpi(china, o3_gmi_china, lat_gmi_china, 
            lng_gmi_china, times_gmi_china, 'gridwise')
    # Significance (alpha = 0.05) of O3-temperature-q-jet correlations
    promptsignificance = input('Calculate significance? [y/n]\n').lower()    
    if promptsignificance == 'y':
        # For O3-temperature
        significance_r_t2mo3 = \
            globalo3_calculate.calculate_r_significance(t2m_merra, o3_gmi, 
            r_t2mo3, lat_gmi, lng_gmi)
        significance_r_t2mo3_transport = \
            globalo3_calculate.calculate_r_significance(t2m_merra, 
            o3_transport_gmi, r_t2mo3_transport, lat_gmi, lng_gmi)
        significance_r_qv2mo3 = \
            globalo3_calculate.calculate_r_significance(qv2m_merra, o3_gmi, 
            r_qv2mo3, lat_gmi, lng_gmi)
        significance_r_qv2mo3_transport = \
            globalo3_calculate.calculate_r_significance(qv2m_merra, 
            o3_transport_gmi, r_qv2mo3_transport, lat_gmi, lng_gmi)
        significance_r_t2mjetdist = \
            globalo3_calculate.calculate_r_significance(t2m_merra, 
            diff_t2mjetdist, r_t2mjetdist, lat_gmi, lng_gmi)
        significance_r_o3jetdist = globalo3_calculate.calculate_r_significance(
            o3_gmi, diff_o3jetdist, r_o3jetdist, lat_gmi, lng_gmi)
        significance_r_qv2mjetdist = \
            globalo3_calculate.calculate_r_significance(o3_gmi, 
            diff_qv2mjetdist, r_qv2mjetdist, lat_gmi, lng_gmi)    
        significance_r_U10Mjetdist = \
            globalo3_calculate.calculate_r_significance(U10M, diff_U10Mjetdist, 
            r_U10Mjetdist, lat_gmi, lng_gmi)
        significance_r_V10Mjetdist = \
            globalo3_calculate.calculate_r_significance(V10M, diff_V10Mjetdist, 
            r_V10Mjetdist, lat_gmi, lng_gmi)
        # Regions where O3-temperature or O3-humidity correlation was 
        # significant in HindcastMR2 but not significant in 
        # HindcastMR2-DiurnalAvgTQ
        diff = np.where(np.logical_and(significance_r_t2mo3 != 1., 
            significance_r_t2mo3_transport == 1.))
        significance_diff_r_t2mo3 = np.empty(r_t2mo3.shape)
        significance_diff_r_t2mo3[:] = np.nan
        significance_diff_r_t2mo3[diff[0], diff[1]] = 1.    
        del diff
        diff = np.where(np.logical_and(significance_r_qv2mo3 != 1., 
            significance_r_qv2mo3_transport == 1.))
        significance_diff_r_qv2mo3 = np.empty(r_qv2mo3.shape)
        significance_diff_r_qv2mo3[:] = np.nan
        significance_diff_r_qv2mo3[diff[0], diff[1]] = 1.
    # Load tracer fields 
    prompttracer = input('Load tracer fields? [y/n]\n').lower()
    if prompttracer == 'y':
        co_50, lat_replay, lng_replay, lev_replay = globalo3_open.open_m2g_c90(
            years, 'co_50', 1000., 800., lngmin, latmax, lngmax, latmin, 
            columnmean=True)
        nh_50, lat_replay, lng_replay, lev_replay = globalo3_open.open_m2g_c90(
            years, 'nh_50', 1000., 800., lngmin, latmax, lngmax, latmin, 
            columnmean=True)
        e90, lat_replayk, lng_replay, lev_replay = globalo3_open.open_m2g_c90(
            years, 'e90', 1000., 800., lngmin, latmax, lngmax, latmin, 
            columnmean=True)
        st80_25, lat_replay, lng_replay, lev_replay = globalo3_open.open_m2g_c90(
            years, 'st80_25', 1000., 800., lngmin, latmax, lngmax, latmin, 
            columnmean=True)
        # Interpolate tracer fields
        co_50 = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
            lng_gmi, lat_replay, lng_replay, co_50, checkplot='yes')
        nh_50 = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
            lng_gmi, lat_replay, lng_replay, nh_50, checkplot='yes')
        e90 = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
            lng_gmi, lat_replay, lng_replay, e90, checkplot='yes')
        st80_25 = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
            lng_gmi, lat_replay, lng_replay, st80_25, checkplot='yes')    
        # Find tracer field in midlatitudes
        co_50_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(co_50, 
            lat_gmi, lng_gmi, 0., 360., 20., 70.)
        nh_50_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(nh_50, 
            lat_gmi, lng_gmi, 0., 360., 20., 70.)
        e90_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(e90, 
            lat_gmi, lng_gmi, 0., 360., 20., 70.)        
        st80_25_ml, lat_ml, lng_ml = globalo3_calculate.find_grid_in_bb(
            st80_25, lat_gmi, lng_gmi, 0., 360., 20., 70.)
        # Find latitude of eddy-driven jet and tracer fields at jet
        scratch, co_50_jet_ml = globalo3_calculate.find_field_atjet(co_50_ml, 
            U500_ml, lat_ml, lng_ml, 20, anom = True) 
        scratch, nh_50_jet_ml = globalo3_calculate.find_field_atjet(nh_50_ml, 
            U500_ml, lat_ml, lng_ml, 20, anom = True) 
        scratch, e90_jet_ml = globalo3_calculate.find_field_atjet(e90_ml, 
            U500_ml, lat_ml, lng_ml, 20, anom = True)         
        scratch, st80_25_jet_ml = globalo3_calculate.find_field_atjet(st80_25_ml, 
            U500_ml, lat_ml, lng_ml, 20, anom = True) 
        # Tracer-temperature correlations
        r_t2mco_50 = globalo3_calculate.calculate_r(t2m_merra, co_50, lat_gmi, 
            lng_gmi)
        r_t2mnh_50 = globalo3_calculate.calculate_r(t2m_merra, nh_50, lat_gmi, 
            lng_gmi)
        r_t2me90 = globalo3_calculate.calculate_r(t2m_merra, e90, lat_gmi, 
            lng_gmi)        
        r_t2mst80_25 = globalo3_calculate.calculate_r(t2m_merra, st80_25, 
            lat_gmi, lng_gmi)
        # Tracer-jet distance correlations and slope
        m_co_50jetdist, r_co_50jetdist, diff_co_50jetdist = \
            globalo3_calculate.calculate_fieldjet_relationship(co_50, lat_gmi, 
            lng_gmi, lat_jet_ml, lng_ml)
        m_nh_50jetdist, r_nh_50jetdist, diff_nh_50jetdist = \
            globalo3_calculate.calculate_fieldjet_relationship(nh_50, lat_gmi, 
            lng_gmi, lat_jet_ml, lng_ml)
        m_e90jetdist, r_e90jetdist, diff_e90jetdist = \
            globalo3_calculate.calculate_fieldjet_relationship(e90, lat_gmi, 
            lng_gmi, lat_jet_ml, lng_ml)            
        m_st80_25jetdist, r_st80_25jetdist, diff_st80_25jetdist = \
            globalo3_calculate.calculate_fieldjet_relationship(st80_25, 
            lat_gmi, lng_gmi, lat_jet_ml, lng_ml)
        # Tracer-jet correlations
        significance_r_co_50jetdist = \
            globalo3_calculate.calculate_r_significance(co_50, 
            diff_co_50jetdist, r_co_50jetdist, lat_gmi, lng_gmi)
        significance_r_nh_50jetdist = \
            globalo3_calculate.calculate_r_significance(nh_50, 
            diff_nh_50jetdist, r_nh_50jetdist, lat_gmi, lng_gmi)
        significance_r_st80_25jetdist = \
            globalo3_calculate.calculate_r_significance(st80_25, 
            diff_st80_25jetdist, r_st80_25jetdist, lat_gmi, lng_gmi)
        significance_r_e90jetdist = \
            globalo3_calculate.calculate_r_significance(e90, diff_e90jetdist, 
            r_e90jetdist, lat_gmi, lng_gmi)
    # Load GISS MERRA-2 cyclone database
    promptcyclones = input('Load cyclones? [y/n]\n').lower()
    if promptcyclones == 'y':    
        cyclones = globalo3_open.open_merra2_cyclones('06/01/2008', 
            '08/31/2010', months_str)
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
          eqjet_o3) = \
              globalo3_calculate.segregate_field_bylat(o3_gmi, lng_gmi, 
              lat_jet_ml, times_gmi, cyclones=cyclones)
        # Same as above but humidity on days with equatorward and poleward 
        # jet              
        (eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_qv2m, 
          eqjet_qv2m) = \
              globalo3_calculate.segregate_field_bylat(qv2m_merra, lng_gmi, 
              lat_jet_ml, times_gmi)             
        # Same as above but find O3 on days with equatorward and poleward jet 
        (eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_t2m, 
          eqjet_t2m) = \
              globalo3_calculate.segregate_field_bylat(t2m_merra, lng_gmi, 
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
        # Find O3 and O3 anomaly within 5 grid cells of cyclone centers
        o3_anom, o3_all, o3_anom_rotated, o3_rotated = \
            globalo3_calculate.o3anom_cyclone(cyclones, times_gmi, lat_gmi, 
            lng_gmi, o3_gmi)
    # Load PBL height from MERRA-2
    promptpbl = input('Load PBL height? [y/n]\n').lower()
    if promptpbl == 'y':    
        # Load daily mean PBL height
        mtime, lat_merra, lng_merra, pblh_merra = \
            globalo3_open.open_merra2pblh_specifieddomain(years, 
            np.arange(0, 24, 1), latmin, latmax, lngmin, lngmax, 
            dailyavg='yes')
#        # Load hourly PBL height 
#        mtime_hrly, lat_merra, lng_merra, pblh_merra_hrly = \
#            globalo3_open.open_merra2pblh_specifieddomain(years, 
#            np.arange(0, 24, 1), latmin, latmax, lngmin, lngmax, 
#            dailyavg='no')
        # Interpolate PBL height
        pblh_merra = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
            lng_gmi, lat_merra, lng_merra, pblh_merra, checkplot='yes')
#        pblh_merra_hrly = globalo3_open.interpolate_merra_to_ctmresolution(
#            lat_gmi, lng_gmi, lat_merra, lng_merra, pblh_merra_hrly, 
#            checkplot='yes')        
        r_pblho3 = globalo3_calculate.calculate_r(pblh_merra, o3_gmi, lat_gmi, 
            lng_gmi)
        r_pblht2m = globalo3_calculate.calculate_r(pblh_merra, t2m_merra, 
            lat_gmi, lng_gmi)
        m_pblhjetdist, r_pblhjetdist, diff_pblhjetdist = \
            globalo3_calculate.calculate_fieldjet_relationship(pblh_merra, 
            lat_gmi, lng_gmi, lat_jet_ml, lng_ml)
        # Calculate significance 
        significance_r_pblho3 = \
            globalo3_calculate.calculate_r_significance(pblh_merra, o3_gmi, 
            r_pblho3, lat_gmi, lng_gmi)
        significance_r_pblhjetdist = \
            globalo3_calculate.calculate_r_significance(pblh_merra,
            diff_pblhjetdist, r_pblhjetdist, lat_gmi, lng_gmi)
#maparea = 'nh'
## Mean O3, mean eddy-driven jet position and variability
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    np.mean(o3_gmi, axis=0), 
#    '%s O$_{\mathregular{3}}$' %season, 
#    '[ppbv]', 
#    np.linspace(25, 65, 11), 
#    'PuBu', 
#    maparea,
#    'meano3_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n=np.nanmean(lat_jet_ml, axis=0), 
#    eerr_n=np.std(lat_jet_ml, axis=0),
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])    
# # Mean 2-meter temperature
# map_hemisphere(lat_gmi, 
#     lng_gmi,
#     np.nanmean(t2m_merra, axis=0), 
#     '$\overline{\mathregular{T}}$',
#     '[K]', 
#     np.linspace(280, 310, 11), 
#     'OrRd',
#     'nh',
#     't2m_%s_%d-%d'%(season, years[0],years[-1]), 
#     oceanon='no',
#     e_lng = lng_gmi,    
#     e_n=np.nanmean(lat_jet_ml,axis=0), 
#     eerr_n=np.std(lat_jet_ml,axis=0),
#     extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#             lat_gmi.min()+1, lat_gmi.max()-5])
# # Mean 2-meter specific humidity
# map_hemisphere(lat_gmi, 
#     lng_gmi,
#     np.nanmean(qv2m_merra, axis=0), 
#     '$\overline{\mathregular{q}}$',
#     '[g kg$^{\mathregular{-1}}$]', 
#     np.linspace(0, 20, 11), 
#     'OrRd',
#     'nh',
#     'qv2m_%s_%d-%d'%(season, years[0],years[-1]), 
#     oceanon='no',
#     extend='max',
#     e_lng = lng_gmi,    
#     e_n=np.nanmean(lat_jet_ml,axis=0), 
#     eerr_n=np.std(lat_jet_ml,axis=0),
#     extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#             lat_gmi.min()+1, lat_gmi.max()-5])

# # Mean 2-meter eastward 10-meter wind
# map_hemisphere(lat_gmi, 
#     lng_gmi,
#     np.nanmean(U10M, axis=0), 
#     '$\overline{\mathregular{U}_{\mathregular{10}}}$',
#     '[m s$^{\mathregular{-1}}$]', 
#     np.linspace(-5, 5, 11), 
#     'RdBu_r',
#     'nh',
#     'U10M_%s_%d-%d'%(season, years[0],years[-1]), 
#     oceanon='no',
#     e_lng = lng_gmi,    
#     e_n=np.nanmean(lat_jet_ml,axis=0), 
#     eerr_n=np.std(lat_jet_ml,axis=0),
#     extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#             lat_gmi.min()+1, lat_gmi.max()-5])
# # Mean 2-meter northward 10-meter wind
# map_hemisphere(lat_gmi, 
#     lng_gmi,
#     np.nanmean(V10M, axis=0), 
#     '$\overline{\mathregular{V}_{\mathregular{10}}}$',
#     '[m s$^{\mathregular{-1}}$]', 
#     np.linspace(-5, 5, 11), 
#     'RdBu_r',
#     'nh',
#     'V10M_%s_%d-%d'%(season, years[0],years[-1]), 
#     oceanon='no',
#     e_lng = lng_gmi,    
#     e_n=np.nanmean(lat_jet_ml,axis=0), 
#     eerr_n=np.std(lat_jet_ml,axis=0),
#     extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#             lat_gmi.min()+1, lat_gmi.max()-5])            
## dO3/dT
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    do3dt2m, 
#    '%s dO$_{\mathregular{3}}$/dT' %season,
#    '[ppbv K$^{\mathregular{-1}}$]', 
#    np.linspace(-2, 2., 9), 
#    'bwr', 
#    maparea,
#    'do3dt2m_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch = significance_r_t2mo3,
#    hatch_freq = ['//////'],    
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])
## dO3/dq
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    do3dq, 
#    '%s dO$_{\mathregular{3}}$/dq' %season,
#    '[ppbv kg g$^{\mathregular{-1}}$]', 
#    np.linspace(-2, 2., 9), 
#    'bwr', 
#    maparea,
#    'do3dq_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch = significance_r_qo3,
#    hatch_freq = ['//////'],    
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])        
## r(T, O3)
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    r_t2mo3, 
#    r'%s $\it{r}\:$(T, O$_\mathregular{3}$)' %season,  
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rt2mo3_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    hatch = significance_r_t2mo3,
#    hatch_freq = ['//////'],    
#    extent = [lng_gmi.min()-180., lng_gmi.max()-180., 
#              lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='neither')
## r(T, O3) from HindcastMR2-DiurnalAvgTQ
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    r_t2mo3_transport,
#    r'%s $\it{r}\:$(T, O$_\mathregular{3}$)' %season,  
#    '[$\cdot$]', 
#    np.linspace(-1., 1., 9),
#    'bwr', 
#    maparea,
#    'rt2mo3_transport_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    hatch = significance_diff_r_t2mo3,
#    hatch_freq = ['////////////'],
#    extent = [lng_gmi.min()-180., lng_gmi.max()-180., 
#              lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='neither')
## r(q, O3)
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    r_qo3, 
#    r'%s $\it{r}\:$(q, O$_\mathregular{3}$)' %season,  
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rqo3_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    hatch = significance_r_qo3,
#    hatch_freq = ['//////'],    
#    extent = [lng_gmi.min()-180., lng_gmi.max()-180., 
#              lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='neither')
## r(q, O3) from HindcastMR2-DiurnalAvgTQ
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    r_qo3_transport,
#    r'%s $\it{r}\:$(q, O$_\mathregular{3}$)' %season,  
#    '[$\cdot$]', 
#    np.linspace(-1., 1., 9),
#    'bwr', 
#    maparea,
#    'rqo3_transport_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    hatch = significance_diff_r_qo3,
#    hatch_freq = ['////////////'],
#    extent = [lng_gmi.min()-180., lng_gmi.max()-180., 
#              lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='neither')
# dO3/d(jet lat - lat)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    m_o3jetdist, 
#    '%s dO$_{\mathregular{3}}$/d($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$)' %season, 
#    '[ppbv degree$^{\mathregular{-1}}$]', 
#    np.linspace(-0.3, 0.3, 7), 
#    'bwr', 
#    maparea,
#    'do3djetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])
## dT/d(jet lat - lat)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    m_t2mjetdist, 
#    '%s dT/d($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$)' %season, 
#    '[K degree$^{\mathregular{-1}}$]', 
#    np.linspace(-0.2, 0.3, 11), 
#    'gist_earth',
#    maparea,
#    'dt2mdjetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])              
## r(jet lat - lat, O3)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_o3jetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, O$_{\mathregular{3}}$)' %season, 
#    '[$\cdot$]', 
#    np.linspace(-1, 1, 9), 
#    'bwr', 
#    maparea,
#    'ro3jetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch = significance_r_o3jetdist,
#    hatch_freq = ['//////'],    
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='neither')
## r(jet lat - lat, T)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_t2mjetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, T)' %season, 
#    '[$\cdot$]', 
#    np.linspace(-1, 1, 9), 
#    'bwr',
#    maparea,
#    'rt2mjetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch = significance_r_t2mjetdist,
#    hatch_freq = ['//////'],
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='neither')
## r(jet lat - lat, q)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_qv2mjetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, q)' %season, 
#    '[$\cdot$]', 
#    np.linspace(-1, 1, 9), 
#    'bwr',
#    maparea,
#    'rqv2mjetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch = significance_r_qv2mjetdist,
#    hatch_freq = ['//////'],
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='neither')
## r(jet lat - lat, CO50)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_co_50jetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, '%season
#    +'CO$_{\mathregular{50}}$)',
#    '[$\cdot$]', 
#    np.linspace(-0.7, 0.7, 8), 
#    'bwr',
#    maparea,
#    'rco_50jetdist_jet_%s_%d-%d'%(season, years[0], years[-1]), 
#    e_lng=lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch=significance_r_co_50jetdist,
#    hatch_freq = ['//////'],
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])
## r(jet lat - lat, NH50)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_nh_50jetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, '%season
#    +'NH$_{\mathregular{50}}$)',
#    '[$\cdot$]', 
#    np.linspace(-0.7, 0.7, 8), 
#    'bwr',
#    maparea,
#    'rnh_50jetdist_jet_%s_%d-%d'%(season, years[0], years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch=significance_r_nh_50jetdist,
#    hatch_freq = ['//////'],    
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])
## r(jet lat - lat, e90)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_e90jetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, '%season
#    +'e90)',
#    '[$\cdot$]', 
#    np.linspace(-0.7, 0.7, 8), 
#    'bwr',
#    maparea,
#    're90jetdist_jet_%s_%d-%d'%(season, years[0], years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch=significance_r_e90jetdist,
#    hatch_freq = ['//////'],    
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])
## r(jet lat - lat, ST80_25)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_st80_25jetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, '%season
#    +'ST80$_{\mathregular{25}}$)',
#    '[$\cdot$]', 
#    np.linspace(-0.7, 0.7, 8), 
#    'bwr',
#    maparea,
#    'rst80_25jetdist_jet_%s_%d-%d'%(season, years[0], years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch=significance_r_st80_25jetdist,
#    hatch_freq = ['//////'],    
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])
## r(PBLH, O3)
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    r_pblho3, 
#    r'%s $\it{r}\:$(PBLH, O$_\mathregular{3}$)' %season,  
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rpblho3_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    extent = [lng_gmi.min()-180., lng_gmi.max()-180., 
#              lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='neither')
## r(T, PBLH)
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    r_pblht2m, 
#    r'%s $\it{r}\:$(PBLH, T)' %season,  
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rt2mpblh_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    extent = [lng_gmi.min()-180., lng_gmi.max()-180., 
#              lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='neither')
## r(jet lat - lat, PBLH)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_pblhjetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, '%season+
#    'PBLH)', 
#    '[$\cdot$]', 
#    np.linspace(-1., 1., 9), 
#    'bwr', 
#    maparea,
#    'rpblhjetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch=significance_r_pblhjetdist,
#    hatch_freq = ['//////'],    
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='neither')    
## Cyclone frequency and eddy-driven jet
#map_hemisphere(lat_cyclones_binned, 
#    lng_cyclones_binned,
#    cyclones_binned, 
#    '',
#    'Cyclones',
#    np.linspace(0, 50, 11), 
#    'Blues', 
#    maparea,
#    'cyclones_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#        extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='both')
## Cylone frequency when jet is poleward 
#map_hemisphere(pwjet_lat_cyclones_binned, 
#    pwjet_lng_cyclones_binned,
#    pwjet_cyclones_binned, 
#    '',
#    'Cyclones',
#    np.linspace(0, 20, 11), 
#    'Blues',
#    maparea,
#    'cyclones_pwjet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n = pwjet_lat,
#    eerr_n = pwjet_lat_var,
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='both')
## Cylone frequency when jet is equatorward 
#map_hemisphere(eqjet_lat_cyclones_binned, 
#    eqjet_lng_cyclones_binned,
#    eqjet_cyclones_binned,
#    '',
#    'Cyclones',
#    np.linspace(0, 20, 11), 
#    'Blues', 
#    maparea,
#    'cyclones_eqjet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n = eqjet_lat,
#    eerr_n = eqjet_lat_var,
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='both')
## POLEWARD/EQUATORWARD ANOMALIES
## O3 anomaly, jet position on days where the jet is extreme poleward
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    pwjet_o3-np.nanmean(o3_gmi,axis=0), 
#    'O$_{\mathregular{3,\:PW}}$ - $\overline{\mathregular{O_3}}$',
#    '[ppbv]', 
#    np.linspace(-6, 6, 13), 
#    'bwr', 
#    maparea,
#    'o3anom_pwjet_cyclones_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n = pwjet_lat,
#    eerr_n = pwjet_lat_var,
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='both')
## 2-meter temperature anomaly, jet position on days where the jet is extreme 
## poleward
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    pwjet_t2m-np.nanmean(t2m_merra,axis=0), 
#    'T$_{\mathregular{PW}}$ - $\overline{\mathregular{T}}$',
#    '[K]', 
#    np.linspace(-6, 6, 13), 
#    'bwr', 
#    maparea,
#    't2manom_pwjet_cyclones_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n = pwjet_lat,
#    eerr_n = pwjet_lat_var,
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='both')
## 2-meter specific humidity anomaly, jet position on days where the jet is 
## extreme poleward
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    pwjet_qv2m-np.nanmean(qv2m_merra,axis=0), 
#    'q$_{\mathregular{PW}}$ - $\overline{\mathregular{q}}$',
#    '[g kg$^{\mathregular{-1}}$]',
#    np.linspace(-2, 2, 9), 
#    'bwr', 
#    maparea,
#    'qv2manom_pwjet_cyclones_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n = pwjet_lat,
#    eerr_n = pwjet_lat_var,
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='both')
## O3 anomaly, jet position on days where the jet is extreme equatorward
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    eqjet_o3-np.nanmean(o3_gmi,axis=0), 
#    'O$_{\mathregular{3,\:EW}}$ - $\overline{\mathregular{O_3}}$',
#    '[ppbv]', 
#    np.linspace(-6, 6, 13), 
#    'bwr', 
#    maparea,
#    'o3anom_ewjet_cyclones_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n = eqjet_lat,
#    eerr_n = eqjet_lat_var,
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='both')
## 2-meter temperature anomaly, jet position on days where the jet is extreme 
## equatorward
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    eqjet_t2m-np.nanmean(t2m_merra,axis=0), 
#    'T$_{\mathregular{EW}}$ - $\overline{\mathregular{T}}$',
#    '[K]', 
#    np.linspace(-6, 6, 13), 
#    'bwr', 
#    maparea,
#    't2manom_ewjet_cyclones_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n = eqjet_lat,
#    eerr_n = eqjet_lat_var,
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='both')
## 2-meter specific humidity anomaly, jet position on days where the jet is 
## extreme poleward
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    eqjet_qv2m-np.nanmean(qv2m_merra,axis=0), 
#    'q$_{\mathregular{EW}}$ - $\overline{\mathregular{q}}$',
#    '[g kg$^{\mathregular{-1}}$]',
#    np.linspace(-2, 2, 9), 
#    'bwr', 
#    maparea,
#    'qv2manom_ewjet_cyclones_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n = eqjet_lat,
#    eerr_n = eqjet_lat_var,
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    extend='both')
## Composite of O3 on days with poleward minus equatorward jet 
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    pwjet_o3-eqjet_o3, 
#    'O$_\mathregular{3,\:PW}$ - O$_\mathregular{3,\:EW}$',
#    '[ppbv]', 
#    np.linspace(-10., 10., 11),
#    'bwr', 
#    maparea,
#    'pwjeto3-ewjeto3_%s_%d-%d'%(season, years[0],years[-1]), 
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='both')
## Composite of temperature on days with poleward minus equatorward jet 
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    pwjet_t2m-eqjet_t2m, 
#    'T$_\mathregular{PW}$ - T$_\mathregular{EW}$',
#    '[K]', 
#    np.linspace(-8., 8., 9),
#    'bwr', 
#    maparea,
#    'pwjett2m-ewjett2m_%s_%d-%d'%(season, years[0],years[-1]), 
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='both')
## Composite of humidity on days with poleward minus equatorward jet 
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    pwjet_qv2m-eqjet_qv2m, 
#    'q$_\mathregular{PW}$ - q$_\mathregular{EW}$',
#    '[g kg$^{\mathregular{-1}}$]', 
#    np.linspace(-4., 4., 9),
#    'bwr', 
#    maparea,
#    'pwjetqv2m-ewjetqv2m_%s_%d-%d'%(season, years[0],years[-1]), 
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='both')
## Cyclone frequency, jet position 
#map_hemisphere(lat_cyclones_binned, 
#    lng_cyclones_binned, 
#    cyclones_binned,
#    'Cyclones',
#    '[$\mathregular{\cdot}$]', 
#    np.linspace(0, 35, 8), 
#    'OrRd', 
#    maparea,
#    'cyclones_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,
#    e_n=np.nanmean(lat_jet_ml, axis=0), 
#    eerr_n=np.std(lat_jet_ml, axis=0),
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    pcolorflag=True, 
#    lbound=0, 
#    ubound=5)
# # Cyclone frequency anomaly, jet position on days where the jet is extreme 
# # poleward
# map_hemisphere(lat_cyclones_binned, 
#     lng_cyclones_binned, 
#     pwjet_cyclones_binned-cyclones_binned,
#     'Cyclones$_{\mathregular{PW}}$ - Cyclones',
#     '[$\mathregular{\cdot}$]', 
#     np.linspace(-25, 25, 11), 
#     'coolwarm', 
#     maparea,
#     'pwcyclonesanom_%s_%d-%d'%(season, years[0],years[-1]), 
#     e_lng=lng_gmi,
#     e_n=pwjet_lat, 
#     eerr_n=pwjet_lat_var, 
#     extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#             lat_gmi.min()+1, lat_gmi.max()-5], 
#     pcolorflag=True, 
#     lbound=-5, 
#     ubound=0)
## Cyclone frequency anomaly, jet position on days where the jet is extreme 
## equatorward
#map_hemisphere(lat_cyclones_binned, 
#    lng_cyclones_binned, 
#    eqjet_cyclones_binned-cyclones_binned,
#    'Cyclones$_{\mathregular{EW}}$ - Cyclones',
#    '[$\mathregular{\cdot}$]', 
#    np.linspace(-25, 25, 11), 
#    'coolwarm', 
#    maparea,
#    'ewcyclonesanom_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng=lng_gmi,
#    e_n=eqjet_lat, 
#    eerr_n=eqjet_lat_var, 
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    pcolorflag=True, 
#    lbound=-5, 
#    ubound=0)
## Difference in cyclone frequency on days with poleward minus equatorward jet 
#map_hemisphere(lat_cyclones_binned, 
#    lng_cyclones_binned, 
#    pwjet_cyclones_binned-eqjet_cyclones_binned,
#    'Cyclones$_{\mathregular{PW}}$ - Cyclones$_{\mathregular{EW}}$',
#    '[$\mathregular{\cdot}$]', 
#    np.linspace(-12, 12, 13), 
#    'coolwarm', 
#    maparea,
#    'pwcyclones-ewcyclones_%s_%d-%d'%(season, years[0],years[-1]), 
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5], 
#    pcolorflag=True, 
#    lbound=-4, 
#    ubound=4)            
## CO_50
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    np.nanmean(co_50, axis=0), 
#    r'Mean %s 1000-800 hPa CO$_\mathregular{50}$' %season,  
#    '[ppbv]', 
#    np.linspace(0., 150., 11),
#    'cubehelix_r', 
#    maparea,
#    'co_50_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    extend='max')
## NH_50
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    np.nanmean(nh_50, axis=0), 
#    r'Mean %s 1000-800 hPa NH$_\mathregular{50}$' %season,  
#    '[ppbv]', 
#    np.linspace(0., 100000., 11),
#    'cubehelix_r', 
#    maparea,
#    'nh_50_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    extend='max')
## e90
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    np.nanmean(e90, axis=0), 
#    r'Mean %s 1000-800 hPa e90' %season,  
#    '[ppbv]', 
#    np.linspace(0., 200., 11),
#    'cubehelix_r', 
#    maparea,
#    'e90_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    extend='max')
## ST80_25
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    np.nanmean(st80_25, axis=0), 
#    r'Mean %s 1000-800 hPa ST80$_{\mathregular{25}}$' %season,  
#    '[ppbv]', 
#    np.linspace(0., 1.5, 11),
#    'cubehelix_r', 
#    maparea,
#    'st80_25_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    extend='max')
## r(T, CO_50)
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    r_t2mco_50, 
#    r'%s $\it{r}\:$(T, CO$_\mathregular{50}$)' %season,  
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rt2mco_50_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    extend='neither')
## r(T, NH_50)
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    r_t2mnh_50, 
#    r'%s $\it{r}\:$(T, NH$_\mathregular{50}$)' %season,  
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rt2mnh_50_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    extend='neither')
## r(T, e90)
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    r_t2me90,
#    r'%s $\it{r}\:$(T, e90)' %season,  
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rt2me90_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    extend='neither')
## r(T, ST80_25)
#map_hemisphere(lat_gmi, 
#    lng_gmi, 
#    r_t2mst80_25, 
#    r'%s $\it{r}\:$(T, ST80$_\mathregular{25}$)' %season,  
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rt2mst80_25_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n = np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n = np.std(lat_jet_ml,axis=0),
#    extend='neither')        
## r(T, O3) using JJA mean values 
#rt2mo3_iav()       
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
## ST80_25 anomaly in vicinity of jet 
#fieldatjet(lat_ml, 
#    lng_ml, 
#    lat_jet_ml, 
#    np.nanmean(st80_25_jet_ml, axis=0),
#    'bwr', 
#    r'$\mathregular{\delta}$ ST80$_{\mathregular{25}}$ [ppbv]',  
#    np.linspace(-0.1, 0.1, 11), 
#    '%s_st80_25anom' %maparea, 
#    '%s_%d-%d'%(season, years[0],years[-1]), 
#    skiplng=6)
## NH_50 anomaly in vicinity of jet 
#fieldatjet(lat_ml, 
#    lng_ml, 
#    lat_jet_ml, 
#    np.nanmean(nh_50_jet_ml, axis=0),
#    'bwr', 
#    r'$\mathregular{\delta}$ NH$_{\mathregular{50}}$ [ppbv]',  
#    np.linspace(-6000, 6000, 9), 
#    '%s_nh_50anom' %maparea, 
#    '%s_%d-%d'%(season, years[0],years[-1]), 
#    skiplng=6)
## CO_50 anomaly in vicinity of jet 
#fieldatjet(lat_ml, 
#    lng_ml, 
#    lat_jet_ml, 
#    np.nanmean(co_50_jet_ml, axis=0),
#    'bwr', 
#    r'$\mathregular{\delta}$ CO$_{\mathregular{50}}$ [ppbv]',  
#    np.linspace(-4, 4, 9), 
#    '%s_co_50anom' %maparea, 
#    '%s_%d-%d'%(season, years[0],years[-1]), 
#    skiplng=6)

"""CORRELATE FLOW AT SURFACE WITH DISTANCE TO JET STREAM"""
## Composites of flow at surface on days with poleward or equatorward jet 
#eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_U10M, eqjet_U10M = \
#    globalo3_calculate.segregate_field_bylat(U10M, lng_gmi, lat_jet_ml, 
#    times_gmi) 
#eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_V10M, eqjet_V10M = \
#    globalo3_calculate.segregate_field_bylat(V10M, lng_gmi, lat_jet_ml, 
#    times_gmi)
#(eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_WIND10M, 
#    eqjet_WIND10M) = globalo3_calculate.segregate_field_bylat(
#    np.hypot(U10M, V10M), lng_gmi, lat_jet_ml, times_gmi)
## r(U10M, O3)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_U10Mo3, 
#    '%s r(U$_{\mathregular{10m}}$, O$_{\mathregular{3}}$)' %season, 
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rU10Mo3_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch=significance_r_U10Mo3,
#    hatch_freq = ['//////'],        
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='neither')
## r(V10M, O3)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_V10Mo3, 
#    '%s r(V$_{\mathregular{10}}$, O$_{\mathregular{3}}$)' %season, 
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rV10Mo3_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch=significance_r_V10Mo3,
#    hatch_freq = ['//////'],            
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='neither')
## r(WIND10M, O3)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_WIND10Mo3, 
#    '%s r(WIND$_{\mathregular{10m}}$, O$_{\mathregular{3}}$)' %season, 
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rWIND10Mo3_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch=significance_r_WIND10Mo3,
#    hatch_freq = ['//////'],            
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='neither')
## r(jet lat - lat, U10M)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_U10Mjetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, U$_{\mathregular{10\:m}}$)' %season, 
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rU10Mjetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='neither')
## r(jet lat - lat, V10M)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_V10Mjetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, V$_{\mathregular{10}}$)' %season, 
#    '[$\cdot$]', 
#    np.linspace(-1.,1,9),
#    'bwr', 
#    maparea,
#    'rV10Mjetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    hatch=significance_r_V10Mjetdist,
#    hatch_freq = ['//////'],    
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='neither')            
## r(jet lat - lat, WIND10M)
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    r_WIND10Mjetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, WIND$_{\mathregular{10\:m}}$)' %season, 
#    '[$\cdot$]', 
#    np.linspace(-1., 1, 9),
#    'bwr', 
#    maparea,
#    'rWIND10Mjetdist_jet_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='neither')
## Composite of 10-meter U/V on days with poleward minus equatorward jet 
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    np.nanmean(U10M, axis=0),
#    'U$_\mathregular{10}$',
#    '[m s$^{\mathregular{-1}}$]',
#    np.linspace(-3., 3., 9),
#    'bwr',
#    maparea,
#    'U10M_%s_%d-%d'%(season, years[0],years[-1]),
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0), 
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180.,
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='both')
#map_hemisphere(lat_gmi, 
#    lng_gmi,
#    pwjet_U10M-eqjet_U10M,
#    'U$_\mathregular{10m,\:PW}$ - U$_\mathregular{10m,\:EW}$',
#    '[m s$^{\mathregular{-1}}$]',
#    np.linspace(-4., 4., 9),
#    'bwr',
#    maparea,
#    'pwjetU10M-ewjetU10M_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180.,
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='both')
#map_hemisphere(lat_gmi,
#    lng_gmi,
#    np.nanmean(V10M, axis=0),
#    'V$_\mathregular{10}$',
#    '[m s$^{\mathregular{-1}}$]', 
#    np.linspace(-2., 2., 9),
#    'bwr',
#    maparea,
#    'V10M_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),    
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='both')
#map_hemisphere(lat_gmi,
#    lng_gmi,
#    pwjet_V10M-eqjet_V10M,
#    'V$_\mathregular{10,\:PW}$ - V$_\mathregular{10,\:EW}$',
#    '[m s$^{\mathregular{-1}}$]', 
#    np.linspace(-4., 4., 9),
#    'bwr',
#    maparea,
#    'pwjetV10M-ewjetV10M_%s_%d-%d'%(season, years[0],years[-1]), 
#    e_lng = lng_gmi,    
#    e_n=np.nanmean(lat_jet_ml,axis=0), 
#    eerr_n=np.std(lat_jet_ml,axis=0),    
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='both')
#map_hemisphere(lat_gmi,
#    lng_gmi,
#    pwjet_WIND10M-eqjet_WIND10M,
#    'WIND$_\mathregular{10M,\:PW}$ - WIND$_\mathregular{10M,\:EW}$',
#    '[m s$^{\mathregular{-1}}$]', 
#    np.linspace(-5., 5., 11),
#    'bwr',
#    maparea,
#    'pwjetWIND10M-ewjetWIND10M_%s_%d-%d'%(season, years[0],years[-1]), 
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5],
#    extend='both',
#    oceanon='no')

"""MODEL-OBSERVATION COMPARISON"""
## North American r(T, O3)
#map_extent([-170, -50, 24, 60], 
#    '%s $r\:$(T, O$_{\mathregular{3}}$)' %season, 
#    '[$\mathregular{\cdot}$]', 
#    'bwr', 
#    np.linspace(-1., 1, 9),
#    'northamerica',
#    'rt2mo3_gmi_aqsnaps_%s_%d-%d' %(season, years[0], years[-1]),
#    lat_ctm=lat_gmi, 
#    lng_ctm=lng_gmi, 
#    field_ctm=r_t2mo3, 
#    lat_obs=np.hstack((lat_naps, lat_aqs)), 
#    lng_obs=np.hstack((lng_naps, lng_aqs)), 
#    field_obs=np.hstack((r_naps, r_aqs)))
## North American r(T, O3) bias
#map_extent([-170, -50, 24, 60], 
#    '%s $r\:$(T, O$_{\mathregular{3}}$)' %season, 
#    '[$\mathregular{\cdot}$]', 
#    'bwr', 
#    np.linspace(-.5, .5, 6),
#    'northamerica',
#    'rt2mo3bias_gmi_aqsnaps_%s_%d-%d' %(season, years[0], years[-1]),
#    lat_obs=np.hstack((lat_naps, lat_aqs)), 
#    lng_obs=np.hstack((lng_naps, lng_aqs)), 
#    field_obs=np.hstack((r_naps_bias, r_aqs_bias)), 
#    extend='both')
## North American dO3/dT
#map_extent([-170, -50, 24, 60], 
#    '%s dO$_{\mathregular{3}}$/dT' %season,
#    '[ppbv K$^{\mathregular{-1}}$]',
#    'bwr', 
#    np.linspace(-2.5, 2.5, 11),
#    'northamerica',
#    'do3dt2m_gmi_aqsnaps_%s_%d-%d' %(season, years[0], years[-1]),
#    lat_ctm=lat_gmi, 
#    lng_ctm=lng_gmi, 
#    field_ctm=do3dt2m, 
#    lat_obs=np.hstack((lat_naps, lat_aqs)), 
#    lng_obs=np.hstack((lng_naps, lng_aqs)), 
#    field_obs=np.hstack((do3dt2m_naps, do3dt2m_aqs)), 
#    extend='both')
## North American dO3/dT bias
#map_extent([-170, -50, 24, 60], 
#    '%s dO$_{\mathregular{3}}$/dT' %season,
#    '[ppbv K$^{\mathregular{-1}}$]',
#    'bwr', 
#    np.linspace(-1., 1., 6),
#    'northamerica',
#    'do3dt2mbias_gmi_aqsnaps_%s_%d-%d' %(season, years[0], years[-1]),
#    lat_obs=np.hstack((lat_naps, lat_aqs)), 
#    lng_obs=np.hstack((lng_naps, lng_aqs)), 
#    field_obs=np.hstack((do3dt2m_naps_bias, do3dt2m_aqs_bias)), 
#    extend='both')
## European r(T, O3)
#map_extent([-14, 37, 28, 70], 
#    '%s $r\:$(T, O$_{\mathregular{3}}$)' %season, 
#    '[$\mathregular{\cdot}$]', 
#    'bwr', 
#    np.linspace(-1., 1., 9),
#    'europe',
#    'rt2mo3_gmi_emep_%s_%d-%d' %(season, years[0], years[-1]),
#    lat_ctm=lat_gmi, 
#    lng_ctm=lng_gmi, 
#    field_ctm=r_t2mo3,
#    lat_obs=lat_emep,
#    lng_obs=lng_emep, 
#    field_obs=r_emep)
## European r(T, O3) bias
#map_extent([-14, 37, 28, 70], 
#    '%s $r\:$(T, O$_{\mathregular{3}}$)' %season, 
#    '[$\mathregular{\cdot}$]', 
#    'bwr', 
#    np.linspace(-.5, .5, 6),
#    'europe',
#    'rt2mo3bias_gmi_emep_%s_%d-%d' %(season, years[0], years[-1]),
#    lat_obs=lat_emep,
#    lng_obs=lng_emep, 
#    field_obs=r_emep_bias,
#    extend='both')
## European dO3/dT
#map_extent([-14, 37, 28, 70], 
#    '%s dO$_{\mathregular{3}}$/dT' %season,
#    '[ppbv K$^{\mathregular{-1}}$]',
#    'bwr', 
#    np.linspace(-2.5, 2.5, 11),
#    'europe',
#    'do3dt2m_gmi_emep_%s_%d-%d' %(season, years[0], years[-1]),
#    lat_ctm=lat_gmi, 
#    lng_ctm=lng_gmi, 
#    field_ctm=do3dt2m, 
#    lat_obs=lat_emep,
#    lng_obs=lng_emep, 
#    field_obs=do3dt2m_emep,
#    extend='both')
## European dO3/dT bias
#map_extent([-14, 37, 28, 70], 
#    '%s dO$_{\mathregular{3}}$/dT' %season,
#    '[ppbv K$^{\mathregular{-1}}$]',
#    'bwr', 
#    np.linspace(-1., 1., 6),
#    'europe',
#    'do3dt2mbias_gmi_emep_%s_%d-%d' %(season, years[0], years[-1]),
#    lat_obs=lat_emep,
#    lng_obs=lng_emep, 
#    field_obs=do3dt2m_emep_bias,
#    extend='both')
# # Zonally-averaged modeled/observed r(T, O3)
# zonalavg_byregion(r_t2mo3, 
#     lat_gmi, 
#     lng_gmi, 
#     r_t2mo3_aqs, 
#     lat_aqs, 
#     lng_aqs, 
#     r_t2mo3_naps, 
#     lat_naps, 
#     lng_naps, 
#     r_t2mo3_emep, 
#     lat_emep, 
#     lng_emep, 
#     r_t2mo3_china, 
#     lat_china, 
#     lng_china,
#     lng_ml, 
#     lat_jet_ml, 
#     '$r\:$(T, O$_{\mathregular{3}}$)', 
#     '[$\cdot$]',
#     np.linspace(-0.5, 1, 7),
#     'r_t2mo3_%s_%d-%d' %(season, years[0], years[-1]))
# # Zonally-averaged modeled/observed r(q, O3)
# zonalavg_byregion(r_qv2mo3, 
#     lat_gmi, 
#     lng_gmi, 
#     r_qv2mo3_aqs, 
#     lat_aqs, 
#     lng_aqs, 
#     r_qv2mo3_naps, 
#     lat_naps, 
#     lng_naps, 
#     r_qv2mo3_emep, 
#     lat_emep, 
#     lng_emep, 
#     r_qv2mo3_china, 
#     lat_china, 
#     lng_china,
#     lng_ml, 
#     lat_jet_ml, 
#     '$r\:$(q, O$_{\mathregular{3}}$)', 
#     '[$\cdot$]',
#     np.linspace(-1., 1, 5),
#     'r_qv2mo3_%s_%d-%d' %(season, years[0], years[-1]))
## Zonally-averaged modeled/observed dO3/dT
#zonalavg_byregion(do3dt2m, 
#    lat_gmi, 
#    lng_gmi, 
#    do3dt2m_aqs, 
#    lat_aqs, 
#    lng_aqs, 
#    do3dt2m_naps, 
#    lat_naps, 
#    lng_naps, 
#    do3dt2m_emep, 
#    lat_emep, 
#    lng_emep, 
#    lng_ml, 
#    lat_jet_ml, 
#    'dO$_{\mathregular{3}}$/dT',
#    '[ppbv K$^{\mathregular{-1}}$]',
#    np.linspace(-1, 3, 5),
#    'do3dt2m_%s_%d-%d' %(season, years[0], years[-1]))
## Zonally-averaged modeled/observed dO3/d(jet lat - lat)
#zonalavg_byregion(m_o3jetdist, 
#    lat_gmi, 
#    lng_gmi, 
#    do3djet_aqs, 
#    lat_aqs, 
#    lng_aqs, 
#    do3djet_naps, 
#    lat_naps,
#    lng_naps,
#    do3djet_emep,
#    lat_emep,
#    lng_emep,
#    lng_ml, 
#    lat_jet_ml, 
#    'dO$_{\mathregular{3}}$/d($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$)',
#    '[ppbv degree$^{\mathregular{-1}}$]',
#    np.linspace(-0.2, 0.4, 5),
#    'do3djetdist_%s_%d-%d' %(season, years[0], years[-1]))
## Zonally-averaged modeled/observed r(jet lat - lat, O3)
#zonalavg_byregion(r_o3jetdist, 
#    lat_gmi, 
#    lng_gmi, 
#    ro3jet_aqs,
#    lat_aqs, 
#    lng_aqs,
#    ro3jet_naps,
#    lat_naps, 
#    lng_naps, 
#    ro3jet_emep,
#    lat_emep,
#    lng_emep,
#    lng_ml, 
#    lat_jet_ml, 
#    'r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, O$_{\mathregular{3}}$)',
#    '[$\cdot$]',
#    np.linspace(-0.7, 0.7, 8), 
#    'ro3jetdist_%s_%d-%d' %(season, years[0], years[-1]))
## Zonally-averaged modeled/observed dT/d(jet lat - lat)
#zonalavg_byregion(m_t2mjetdist, 
#    lat_gmi, 
#    lng_gmi, 
#    [], 
#    [], 
#    [],
#    [],
#    [],
#    [],
#    [],
#    [],
#    [],
#    lng_ml, 
#    lat_jet_ml, 
#    'dT/d($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$)',
#    '[K degree$^{\mathregular{-1}}$]',
#    np.linspace(-0.2, 0.4, 5),
#    'dt2mdjetdist_%s_%d-%d' %(season, years[0], years[-1]))
## Zonally-averaged modeled/observed r(jet lat - lat, T)
#zonalavg_byregion(r_t2mjetdist,
#    lat_gmi, 
#    lng_gmi, 
#    [], 
#    [], 
#    [],
#    [],
#    [],
#    [],
#    [],
#    [],
#    [],
#    [],
#    [],
#    [],
#    lng_ml, 
#    lat_jet_ml, 
#    'r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, T)',
#    '[$\cdot$]',
#    np.linspace(-0.7, 0.7, 8), 
#    'rt2mjetdist_%s_%d-%d' %(season, years[0], years[-1]))
## Zonally-averaged modeled/observed r(jet lat - lat, q)
#zonalavg_byregion(r_qv2mjetdist,
#    lat_gmi, 
#    lng_gmi, 
#    [], 
#    [], 
#    [],
#    [],
#    [],
#    [],
#    [],
#    [],
#    [],
#    [],
#    [],
#    [],
#    lng_ml, 
#    lat_jet_ml, 
#    'r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, q)',
#    '[$\cdot$]',
#    np.linspace(-0.7, 0.7, 8), 
#    'rqv2mjetdist_%s_%d-%d' %(season, years[0], years[-1]))

"""PLOT LAND-BASED ZONALLY AVERAGED TRACER-TEMPERATURE FIELDS"""
#import matplotlib.pyplot as plt
#land_nh = globalo3_calculate.find_grid_overland(lat_gmi[:-5], lng_gmi)  
##
#r_t2mco_50_za = r_t2mco_50[:-5]*land_nh
#r_t2mco_50_za = np.nanmean(r_t2mco_50_za, axis=1) 
##
#r_t2mnh_50_za = r_t2mnh_50[:-5]*land_nh
#r_t2mnh_50_za = np.nanmean(r_t2mnh_50_za, axis=1)
##
#r_t2mst80_25_za = r_t2mst80_25[:-5]*land_nh
#r_t2mst80_25_za = np.nanmean(r_t2mst80_25_za, axis=1)
##
#r_t2mo3_za = r_t2mo3[:-5]*land_nh
#r_t2mo3_za = np.nanmean(r_t2mo3_za, axis=1)  
## Plotting
#fig = plt.figure(figsize=(8, 3))
#ax = plt.subplot2grid((1, 1), (0, 0))
#ax.plot(r_t2mo3_za, '-', color='k', lw=3, label='O$_{\mathregular{3}}$'); 
#ax.plot(r_t2mnh_50_za, '-', color='#d95f02', label='NH$_{\mathregular{50}}$'); 
#ax.plot(r_t2mst80_25_za, '-', color = '#7570b3', label='ST80$_{\mathregular{25}}$');
#ax.plot(r_t2mco_50_za, '-b', color='#1b9e77', label='CO$_{\mathregular{50}}$');
#ax.set_xlim([0, 80])
#ax.set_xticks([0, 20, 40, 60, 80])
#ax.set_xticklabels([0, 20, 40, 60, 80], fontsize=12)
#ax.set_xlabel('Latitude [$^{\circ}$]', fontsize=14)
#ax.set_ylim([-0.25, 0.65])
#ax.set_yticks([-0.25, 0.05, 0.35, 0.65])
#ax.set_yticklabels([-0.25,  0.05,  0.35,  0.65], fontsize=12)
#ax.set_ylabel(
#    '$r\:$(T, $\mathregular{\chi}$)$_{\mathregular{\overline{\lambda}}}$', 
#    fontsize=14)
#plt.legend(loc = 1, bbox_to_anchor=(.88, 1.2), ncol = 4, fontsize = 12, 
#    frameon=False)
#plt.subplots_adjust(bottom=0.2)
#plt.savefig('/Users/ghkerr/Desktop/o3_tracer_tempcorrel.png', dpi=350)

"""PLOT VERTICALLY-INTEGRATED MERIDIONAL TRACER FLUXES"""
#co_50_column, lat_replay, lng_replay, lev_replay = globalo3_open.open_m2g_c90(
#    years, 'co_50', 1000., 800., lngmin, latmax, lngmax, latmin, 
#    columnmean=False)
#nh_50_column, lat_replay, lng_replay, lev_replay = globalo3_open.open_m2g_c90(
#    years, 'nh_50', 1000., 800., lngmin, latmax, lngmax, latmin, 
#    columnmean=False)
#st80_25_column, lat_replay, lng_replay, lev_replay = \
#    globalo3_open.open_m2g_c90(years, 'st80_25', 1000., 800., lngmin, latmax, 
#    lngmax, latmin, columnmean=False)
#e90_column, lat_replay, lng_replay, lev_replay = globalo3_open.open_m2g_c90(
#    years, 'e90', 1000., 800., lngmin, latmax, lngmax, latmin, 
#    columnmean=False)
#V_column, lat_merra, lng_merra, column = \
#    globalo3_open.open_merra2u_specifieddomain(years, 'V', 1000., 800., lngmin, 
#    latmax, lngmax, latmin, columnmean=False)
## Interpolate column CO_50 and meridional wind
#co_50_column = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
#    lng_gmi, lat_replay, lng_replay, co_50_column, checkplot='yes')
#nh_50_column = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
#    lng_gmi, lat_replay, lng_replay, nh_50_column, checkplot='yes')
#st80_25_column = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
#    lng_gmi, lat_replay, lng_replay, st80_25_column, checkplot='yes')
#e90_column = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
#    lng_gmi, lat_replay, lng_replay, e90_column, checkplot='yes')
#V_column = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
#    lng_gmi, lat_merra, lng_merra, V_column, checkplot='yes')
## Find vertically-integrated meridional tracer flux on all days
#co_50total, co_50mean, co_50eddy = \
#    globalo3_calculate.verticallyintegrated_meridional_flux(co_50_column/1e9, 
#    V_column, mtime, lat_gmi, lng_gmi, lev_replay, 950., 800., (28/28.97))
#nh_50total, nh_50mean, nh_50eddy = \
#    globalo3_calculate.verticallyintegrated_meridional_flux(nh_50_column/1e9, 
#    V_column, mtime, lat_gmi, lng_gmi, lev_replay, 950., 800., (28/28.97))
#st80_25total, st80_25mean, st80_25eddy = \
#    globalo3_calculate.verticallyintegrated_meridional_flux(st80_25_column/1e9, 
#    V_column, mtime, lat_gmi, lng_gmi, lev_replay, 950., 800., (28/28.97)) 
#e90total, e90mean, e90eddy = \
#    globalo3_calculate.verticallyintegrated_meridional_flux(e90_column/1e9, 
#    V_column, mtime, lat_gmi, lng_gmi, lev_replay, 950., 800., (28/28.97))    
## Plotting
#zonalavg_verticallyintegrated_meridional_flux(co_50total, co_50mean, co_50eddy, 
#    lat_gmi, 'CO$_{\mathregular{50}}$', '950-800_co_50')
#zonalavg_verticallyintegrated_meridional_flux(nh_50total, nh_50mean, nh_50eddy, 
#    lat_gmi, 'NH$_{\mathregular{50}}$', '950-800_nh_50')
#zonalavg_verticallyintegrated_meridional_flux(st80_25total, st80_25mean, 
#    st80_25eddy, lat_gmi, 'ST80$_{\mathregular{25}}$', '950-800_st80_25')
#zonalavg_verticallyintegrated_meridional_flux(e90total, e90mean, e90eddy, 
#    lat_gmi, 'e90', '950-800_e90')
## Sort columned fields by days with a pole- versus equatorward jet
#V_column_eqjet, V_column_pwjet, co_50_column_eqjet, co_50_column_pwjet = \
#    globalo3_calculate.sortfield_byjetlat_column(V_column, co_50_column, 
#    lat_jet_ml, lng_gmi, lat_gmi, lev_replay, psize=30)
#V_column_eqjet, V_column_pwjet, nh_50_column_eqjet, nh_50_column_pwjet = \
#    globalo3_calculate.sortfield_byjetlat_column(V_column, nh_50_column, 
#    lat_jet_ml, lng_gmi, lat_gmi, lev_replay, psize=30)
#V_column_eqjet, V_column_pwjet, st80_25_column_eqjet, st80_25_column_pwjet = \
#    globalo3_calculate.sortfield_byjetlat_column(V_column, st80_25_column, 
#    lat_jet_ml, lng_gmi, lat_gmi, lev_replay, psize=30)
#V_column_eqjet, V_column_pwjet, e90_column_eqjet, e90_column_pwjet = \
#    globalo3_calculate.sortfield_byjetlat_column(V_column, e90_column, 
#    lat_jet_ml, lng_gmi, lat_gmi, lev_replay, psize=30)
## Find vertically-integrated meridional tracer flux on days with a pole- 
## versus equatorward jet
#co_50total_eqjet, co_50mean_eqjet, co_50eddy_eqjet = \
#    globalo3_calculate.verticallyintegrated_meridional_flux(
#    co_50_column_eqjet/1e9, V_column_eqjet, mtime[:len(V_column_eqjet)], 
#    lat_gmi, lng_gmi, lev_replay, 950., 800., (28/28.97))
#co_50total_pwjet, co_50mean_pwjet, co_50eddy_pwjet = \
#    globalo3_calculate.verticallyintegrated_meridional_flux(
#    co_50_column_pwjet/1e9, V_column_pwjet, mtime[:len(V_column_pwjet)], 
#    lat_gmi, lng_gmi, lev_replay, 950., 800., (28/28.97))
#nh_50total_eqjet, nh_50mean_eqjet, nh_50eddy_eqjet = \
#    globalo3_calculate.verticallyintegrated_meridional_flux(
#    nh_50_column_eqjet/1e9, V_column_eqjet, mtime[:len(V_column_eqjet)], 
#    lat_gmi, lng_gmi, lev_replay, 950., 800., (28/28.97))
#nh_50total_pwjet, nh_50mean_pwjet, nh_50eddy_pwjet = \
#    globalo3_calculate.verticallyintegrated_meridional_flux(
#    nh_50_column_pwjet/1e9, V_column_pwjet, mtime[:len(V_column_pwjet)], 
#    lat_gmi, lng_gmi, lev_replay, 950., 800., (28/28.97))    
#st80_25total_eqjet, st80_25mean_eqjet, st80_25eddy_eqjet = \
#    globalo3_calculate.verticallyintegrated_meridional_flux(
#    st80_25_column_eqjet/1e9, V_column_eqjet, mtime[:len(V_column_eqjet)], 
#    lat_gmi, lng_gmi, lev_replay, 950., 800., (28/28.97))
#st80_25total_pwjet, st80_25mean_pwjet, st80_25eddy_pwjet = \
#    globalo3_calculate.verticallyintegrated_meridional_flux(
#    st80_25_column_pwjet/1e9, V_column_pwjet, mtime[:len(V_column_pwjet)], 
#    lat_gmi, lng_gmi, lev_replay, 950., 800., (28/28.97))        
#e90total_eqjet, e90mean_eqjet, e90eddy_eqjet = \
#    globalo3_calculate.verticallyintegrated_meridional_flux(
#    e90_column_eqjet/1e9, V_column_eqjet, mtime[:len(V_column_eqjet)], 
#    lat_gmi, lng_gmi, lev_replay, 950., 800., (28/28.97))
#e90total_pwjet, e90mean_pwjet, e90eddy_pwjet = \
#    globalo3_calculate.verticallyintegrated_meridional_flux(
#    e90_column_pwjet/1e9, V_column_pwjet, mtime[:len(V_column_pwjet)], 
#    lat_gmi, lng_gmi, lev_replay, 950., 800., (28/28.97))   
## Plotting 
#zonalavg_verticallyintegrated_meridional_flux(co_50total_eqjet, 
#    co_50mean_eqjet, co_50eddy_eqjet, lat_gmi, 
#    'CO$_{\mathregular{50,\:equatorward\;jet}}$', '950-800_co_50_eqjet')
#zonalavg_verticallyintegrated_meridional_flux(co_50total_pwjet, 
#    co_50mean_pwjet, co_50eddy_pwjet, lat_gmi, 
#    'CO$_{\mathregular{50,\:poleward\;jet}}$', '950-800_co_50_pwjet')
#zonalavg_verticallyintegrated_meridional_flux(nh_50total_eqjet, 
#    nh_50mean_eqjet, nh_50eddy_eqjet, lat_gmi, 
#    'NH$_{\mathregular{50,\:equatorward\;jet}}$', '950-800_nh_50_eqjet')
#zonalavg_verticallyintegrated_meridional_flux(nh_50total_pwjet, 
#    nh_50mean_pwjet, nh_50eddy_pwjet, lat_gmi, 
#    'NH$_{\mathregular{50,\:poleward\;jet}}$', '950-800_nh_50_pwjet')  
#zonalavg_verticallyintegrated_meridional_flux(st80_25total_eqjet, 
#    st80_25mean_eqjet, st80_25eddy_eqjet, lat_gmi, 
#    'ST80$_{\mathregular{25,\:equatorward\;jet}}$', '950-800_st80_25_eqjet')
#zonalavg_verticallyintegrated_meridional_flux(st80_25total_pwjet, 
#    st80_25mean_pwjet, st80_25eddy_pwjet, lat_gmi, 
#    'ST80$_{\mathregular{25,\:poleward\;jet}}$', '950-800_st80_25_pwjet') 
#zonalavg_verticallyintegrated_meridional_flux(e90total_eqjet, e90mean_eqjet, 
#    e90eddy_eqjet, lat_gmi, 'e90$_{\mathregular{equatorward\;jet}}$', 
#    '950-800_e90_eqjet')
#zonalavg_verticallyintegrated_meridional_flux(e90total_pwjet, e90mean_pwjet, 
#    e90eddy_pwjet, lat_gmi, 'e90$_{\mathregular{poleward\;jet}}$', 
#    '950-800_e90_pwjet') 

"""LOAD GEOSCHEM TRACERS AND DETERMINE RELATIONSHIP WITH EDDY-DRIVEN JET"""
## Load GEOSChem tracer diagnostics
#TROPIC_50, lat_gc, lng_gc, lev_gc = \
#    globalo3_open.open_geoschem_merra2_2x25_RnPbBe(years, months_str, 
#    'SpeciesConc_TROPIC_50', latmin, latmax, lngmin, lngmax, 800., 1000., 
#    'mean')
#MIDLAT_50, lat_gc, lng_gc, lev_gc = \
#    globalo3_open.open_geoschem_merra2_2x25_RnPbBe(years, months_str, 
#    'SpeciesConc_MIDLAT_50', latmin, latmax, lngmin, lngmax, 800., 1000., 
#    'mean')
#POLAR_50, lat_gc, lng_gc, lev_gc = \
#    globalo3_open.open_geoschem_merra2_2x25_RnPbBe(years, months_str, 
#    'SpeciesConc_POLAR_50', latmin, latmax, lngmin, lngmax, 800., 1000., 
#    'mean')
#GLOBAL_50, lat_gc, lng_gc, lev_gc = \
#    globalo3_open.open_geoschem_merra2_2x25_RnPbBe(years, months_str, 
#    'SpeciesConc_GLOBAL_50', latmin, latmax, lngmin, lngmax, 800., 1000., 
#    'mean')
## Load Northern Hemisphere 500 hPa MERRA-2 data for use with GEOSChem
## simulations, but inpolate 
#U500_coarse, mtime_coarse, lat_merra_coarse, lng_merra_coarse = \
#    globalo3_open.open_merra2_specifieddomain(years, months_str,
#    [0,3,6,9,12,15,18,21], 'U', 'inst3_3d_asm_Np_500hPa', lngmin, latmax, 
#    lngmax, latmin, dailyavg='yes')
## Interpolate 500 hPa u-wind to resolution of GEOSChem
#U500_coarse = globalo3_open.interpolate_merra_to_ctmresolution(lat_gc, 
#    lng_gc, lat_merra_coarse, lng_merra_coarse, U500_coarse, checkplot='yes')
## Subset fields in mid-latitudes
#U500_coarse_ml, lat_coarse_ml, lng_coarse_ml = \
#    globalo3_calculate.find_grid_in_bb(U500_coarse, lat_gc, lng_gc, 0., 360., 
#    20., 70.)
#TROPIC_50_ml, lat_coarse_ml, lng_coarse_ml = \
#    globalo3_calculate.find_grid_in_bb(TROPIC_50, lat_gc, lng_gc, 0., 360., 
#    20., 70.)
#MIDLAT_50_ml, lat_coarse_ml, lng_coarse_ml = \
#    globalo3_calculate.find_grid_in_bb(MIDLAT_50, lat_gc, lng_gc, 0., 360., 
#    20., 70.)
#POLAR_50_ml, lat_coarse_ml, lng_coarse_ml = \
#    globalo3_calculate.find_grid_in_bb(POLAR_50, lat_gc, lng_gc, 0., 360., 
#    20., 70.)
#GLOBAL_50_ml, lat_coarse_ml, lng_coarse_ml = \
#    globalo3_calculate.find_grid_in_bb(GLOBAL_50, lat_gc, lng_gc, 0., 360., 
#    20., 70.)    
## Find latitude of eddy-driven jet and fields at jet
#lat_coarse_jet_ml, TROPIC_50_jet_ml = globalo3_calculate.find_field_atjet(
#    TROPIC_50_ml, U500_coarse_ml, lat_coarse_ml, lng_coarse_ml, 20, 
#    anom = True)    
## Slope and correlation of GEOSChem tracers/jet distance and jet distance
#m_TROPIC_50jetdist, r_TROPIC_50jetdist, diff_TROPIC_50jetdist = \
#    globalo3_calculate.calculate_fieldjet_relationship(TROPIC_50, lat_gc, 
#    lng_gc, lat_coarse_jet_ml, lng_coarse_ml)
#m_MIDLAT_50jetdist, r_MIDLAT_50jetdist, diff_MIDLAT_50jetdist = \
#    globalo3_calculate.calculate_fieldjet_relationship(MIDLAT_50, lat_gc, 
#    lng_gc, lat_coarse_jet_ml, lng_coarse_ml)
#m_POLAR_50jetdist, r_POLAR_50jetdist, diff_POLAR_50jetdist = \
#    globalo3_calculate.calculate_fieldjet_relationship(POLAR_50, lat_gc, 
#    lng_gc, lat_coarse_jet_ml, lng_coarse_ml)
#m_GLOBAL_50jetdist, r_GLOBAL_50jetdist, diff_GLOBAL_50jetdist = \
#    globalo3_calculate.calculate_fieldjet_relationship(GLOBAL_50, lat_gc, 
#    lng_gc, lat_coarse_jet_ml, lng_coarse_ml)
## r(jet lat - lat, TROPIC_50)
#map_hemisphere(lat_gc, 
#    lng_gc,
#    r_TROPIC_50jetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, '%season
#    +'TROPIC$_{\mathregular{50}}$)',
#    '[$\cdot$]', 
#    np.linspace(-0.7, 0.7, 8), 
#    'bwr',
#    maparea,
#    'rTROPIC_50jetdist_jet_%s_%d-%d'%(season, years[0], years[-1]), 
#    e_lng = lng_gc,    
#    e_n=np.nanmean(lat_coarse_jet_ml,axis=0), 
#    eerr_n=np.std(lat_coarse_jet_ml,axis=0),
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])
## r(jet lat - lat, MIDLAT_50)
#map_hemisphere(lat_gc, 
#    lng_gc,
#    r_MIDLAT_50jetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, '%season
#    +'MIDLAT$_{\mathregular{50}}$)',
#    '[$\cdot$]', 
#    np.linspace(-0.7, 0.7, 8), 
#    'bwr',
#    maparea,
#    'rMIDLAT_50jetdist_jet_%s_%d-%d'%(season, years[0], years[-1]), 
#    e_lng = lng_gc,    
#    e_n=np.nanmean(lat_coarse_jet_ml,axis=0), 
#    eerr_n=np.std(lat_coarse_jet_ml,axis=0),
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])        
## r(jet lat - lat, POLAR_50)
#map_hemisphere(lat_gc, 
#    lng_gc,
#    r_POLAR_50jetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, '%season
#    +'POLAR$_{\mathregular{50}}$)',
#    '[$\cdot$]', 
#    np.linspace(-0.7, 0.7, 8), 
#    'bwr',
#    maparea,
#    'rPOLAR_50jetdist_jet_%s_%d-%d'%(season, years[0], years[-1]), 
#    e_lng = lng_gc,    
#    e_n=np.nanmean(lat_coarse_jet_ml,axis=0), 
#    eerr_n=np.std(lat_coarse_jet_ml,axis=0),
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])  
## r(jet lat - lat, GLOBAL_50)
#map_hemisphere(lat_gc, 
#    lng_gc,
#    r_GLOBAL_50jetdist, 
#    '%s r($\mathregular{\phi_{jet}}-{\mathregular{\phi}}$, '%season
#    +'GLOBAL$_{\mathregular{50}}$)',
#    '[$\cdot$]', 
#    np.linspace(-0.7, 0.7, 8), 
#    'bwr',
#    maparea,
#    'rGLOBAL_50jetdist_jet_%s_%d-%d'%(season, years[0], years[-1]), 
#    e_lng = lng_gc,    
#    e_n=np.nanmean(lat_coarse_jet_ml,axis=0), 
#    eerr_n=np.std(lat_coarse_jet_ml,axis=0),
#    extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#            lat_gmi.min()+1, lat_gmi.max()-5])

"""MEAN AND EDDY FLUXES OF OTHER TRACER GASES"""
# lat_gmi, lng_gmi, times_gmi, no2_gmi = \
#     globalo3_open.open_overpass2_specifieddomain(years, months_str, latmin, 
#     latmax, lngmin, lngmax, 'NO2', 'HindcastMR2')
# lat_gmi, lng_gmi, times_gmi, no_gmi = \
#     globalo3_open.open_overpass2_specifieddomain(years, months_str, latmin, 
#     latmax, lngmin, lngmax, 'NO', 'HindcastMR2')
# lat_gmi, lng_gmi, times_gmi, co_gmi = \
#     globalo3_open.open_overpass2_specifieddomain(years, months_str, latmin, 
#     latmax, lngmin, lngmax, 'CO', 'HindcastMR2')
# co_gmi = co_gmi*1e9
# no_gmi = no_gmi*1e9
# no2_gmi = no2_gmi*1e9
# nox_gmi = (no_gmi+no2_gmi)
# species = no2_gmi*1e9
# import numpy as np
# # Separate 10-meter meridional wind and O3 by jet position
# V10M_eqjet, V10M_pwjet, species_eqjet, species_pwjet = \
#     globalo3_calculate.sortfield_byjetlat_column(
#     np.reshape(V10M, (276,1,92,288)), np.reshape(species, (276,1,92,288)),
#     lat_jet_ml, lng_gmi, lat_gmi, np.array([1000.]), psize=30)
# species_mean, species_stationary, species_transient, species_total = \
#     globalo3_calculate.meridional_flux(V10M, species, mtime, lat_gmi, lng_gmi)
# (species_mean_pwjet, species_stationary_pwjet, species_transient_pwjet, 
#     species_total_pwjet) = globalo3_calculate.meridional_flux(V10M_pwjet[:,0], 
#     species_pwjet[:,0], mtime[:len(V10M_pwjet)], lat_gmi, lng_gmi)
# (species_mean_eqjet, species_stationary_eqjet, species_transient_eqjet, 
#     species_total_eqjet) = globalo3_calculate.meridional_flux(V10M_eqjet[:,0], 
#     species_eqjet[:,0], mtime[:len(V10M_eqjet)], lat_gmi, lng_gmi)   
# import matplotlib.pyplot as plt
# # O3 meridional flux
# fig = plt.figure(figsize=(11, 3))
# ax1 = plt.subplot2grid((1,3),(0,0))
# ax2 = plt.subplot2grid((1,3),(0,1))
# ax3 = plt.subplot2grid((1,3),(0,2))
# for ax in [ax1, ax2, ax3]:
#     ax.set_xlim([25,70])
#     ax.set_xticks([25,40,55,70])
#     ax.set_xticklabels([25, 40, 55, 70], fontsize=12)
#     ax.set_ylim([-0.20, 0.20])
#     # ax.set_yticks(np.linspace(-0.08, 0.08, 5))
#     # ax.set_yticklabels([-0.08, -0.04, 0, 0.04, 0.08], fontsize=12)
#     # ax.set_ylim([-60, 60])
#     # ax.set_yticks(np.linspace(-60, 60, 7))
#     # ax.set_yticklabels(np.linspace(-60, 60, 7), fontsize=12)
#     # Add horizontal line for field = 0
#     ax.axhline(y=0.0, color='k', lw=1., linestyle='--', zorder=1)
# ax2.set_yticklabels([]); ax3.set_yticklabels([]) 
# latweight = np.cos(np.deg2rad(lat_gmi))
# # O3 meridional flux for all days
# ax1.plot(species_total*latweight, lw=2, ls='-', color='#A2C3E0')
# ax1.plot(species_mean*latweight, lw=2, ls='--', color='#EF9802')
# ax1.plot((species_stationary+species_transient)*latweight, lw=2, ls='-', color='#3F79B7')
# # O3 meridional flux for PW jet days
# ax2.plot(species_total_pwjet*latweight, lw=2, color='#A2C3E0', label='Total')
# ax2.plot(species_mean_pwjet*latweight, ls='--', lw=2, color='#EF9802', label='Mean')
# ax2.plot((species_stationary_pwjet+species_transient_pwjet)*latweight, 
#     lw=2, color='#3F79B7', label='Eddy (Transient + Stationary)')
# # O3 meridional flux for EW jet days
# ax3.plot(species_total_eqjet*latweight, lw=2, color='#A2C3E0')
# ax3.plot(species_mean_eqjet*latweight, ls='--', lw=2, color='#EF9802')
# ax3.plot((species_stationary_eqjet+species_transient_eqjet)*latweight, lw=2, 
#     color='#3F79B7')
# # Set titles, axis labels
# ax1.text(0.04, 0.9, '(a)', ha='left', transform=ax1.transAxes, 
#     fontsize=16, zorder=20)
# ax1.set_title('All', fontsize=16)#, x=0.02, ha='left')
# ax2.set_title('PW', fontsize=16)
# ax3.set_title('EW', fontsize=16)
# ax2.text(0.04, 0.9, '(b)', ha='left', transform=ax2.transAxes, 
#     fontsize=18, zorder=20)
# ax3.text(0.04, 0.9, '(c)', ha='left', transform=ax3.transAxes, 
#     fontsize=16, zorder=20) 
# # ax1.set_ylabel('O$_{\mathregular{3}}$ flux [ppbv m s$^{\mathregular{-1}}$]', 
# #     fontsize=16)
# ax1.set_ylabel('NO$_{\mathregular{2}}$ flux [ppbv m s$^{\mathregular{-1}}$]', 
#     fontsize=16)
# ax1.get_yaxis().set_label_coords(-0.25,0.5)
# # Add legend
# ax2.legend(loc=2, bbox_to_anchor=(-0.8,-0.08), ncol=3, fontsize=16, 
#     frameon=False)
# plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#             'no2flux_pwewjet.png', dpi=300)
# # (PW - EW) composites of CO, NOx
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, co_pwjet, co_eqjet = \
#     globalo3_calculate.segregate_field_bylat(co_gmi, lng_gmi, lat_jet_ml, 
#     times_gmi)
# map_hemisphere(lat_gmi, 
#     lng_gmi,
#     co_pwjet-co_eqjet,
#     'CO$_{\mathregular{PW}}$ - CO$_{\mathregular{EW}}$',
#     '[ppbv]', 
#     np.linspace(-20, 20, 13), 
#     'bwr', 
#     maparea,
#     'pwjetco-ewjetco_%s_%d-%d'%(season, years[0],years[-1]), 
#     e_lng = lng_gmi,
#     e_n=np.nanmean(lat_jet_ml,axis=0), 
#     eerr_n=np.std(lat_jet_ml,axis=0),
#     extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#             lat_gmi.min()+1, lat_gmi.max()-5], 
#     extend='both')
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, nox_pwjet, nox_eqjet = \
#     globalo3_calculate.segregate_field_bylat(nox_gmi, lng_gmi, lat_jet_ml, 
#     times_gmi)
# map_hemisphere(lat_gmi, 
#     lng_gmi,
#     nox_pwjet-nox_eqjet,
#     'NO$_{x\mathregular{,\:PW}}$ - NO$_{x\mathregular{,\:EW}}$',
#     '[ppbv]', 
#     np.linspace(-0.1, 0.1, 11), 
#     'bwr', 
#     maparea,
#     'pwjetnox-ewjetnox_%s_%d-%d'%(season, years[0],years[-1]), 
#     e_lng = lng_gmi,
#     e_n=np.nanmean(lat_jet_ml,axis=0), 
#     eerr_n=np.std(lat_jet_ml,axis=0),
#     extent=[lng_gmi.min()-180., lng_gmi.max()-180., 
#             lat_gmi.min()+1, lat_gmi.max()-5], 
#     extend='both')



        
# import numpy as np
# import matplotlib as mpl
# mpl.rcParams['hatch.linewidth']=0.3     
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter    
# # Composites of flow at surface on days with poleward or equatorward jet 
# ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
#                                         edgecolor=None, facecolor='lightgrey')
# fig = plt.figure(figsize=(10,7))
# # PBLH(PW - EW) 
# ax1 = plt.subplot2grid((3,2), (0,0), colspan=2,
# projection=ccrs.PlateCarree(central_longitude=0.))
# ax1.set_title('(a) PBLH$_\mathregular{PW}$ $-$ PBLH$_\mathregular{EW}$', 
#               fontsize=16, x=0.02, ha='left')
# # ax1.add_feature(ocean50m, zorder=3)
# ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
# ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
#                 lat_gmi.min()+1, lat_gmi.max()-5])
# ax1.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
# lat_formatter = LatitudeFormatter()    
# ax1.yaxis.set_major_formatter(lat_formatter)
# cmap = plt.get_cmap('RdBu')
# mb = ax1.contourf(lng_gmi, lat_gmi, np.nanmean(o3_gmi, axis=0),
#                   np.linspace(30,65,7), cmap=cmap, extend='both', 
#                   transform=ccrs.PlateCarree(), zorder=1)
# # Add colorbar
# plt.gcf().subplots_adjust(left=0.15, right=0.9, hspace=0.3)
# # plt.gcf().subplots_adjust(left=0.05, right=0.86, hspace=0.3)    
# colorbar_axes = plt.gcf().add_axes([ax1.get_position().x1+0.03, 
#     ax1.get_position().y0, 0.02, (ax1.get_position().y1-
#     ax1.get_position().y0)]) 
# colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
# ticks=np.linspace(-300,300,7), extend='neither')
# colorbar.ax.tick_params(labelsize=12)
# colorbar.set_label('[m]', fontsize=16)
# ax1.outline_patch.set_zorder(20)
# # Add axis for zonal mean PBLH
# ax1l = plt.gcf().add_axes([ax1.get_position().x0-0.15, ax1.get_position().y0, 
#             0.1,(ax1.get_position().y1-ax1.get_position().y0)])
# ax1l.plot()        
# import numpy as np
# import matplotlib as mpl
# mpl.rcParams['hatch.linewidth']=0.3     
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# # Separate O3 into pole- and equatorward jet days
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_o3, eqjet_o3 = \
#     globalo3_calculate.segregate_field_bylat(o3_gmi, lng_gmi, lat_jet_ml, 
#     times_gmi)
# # Load ocean shapefiles
# ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
#     edgecolor=None, facecolor='lightgrey')
# fig = plt.figure(figsize=(9,7))
# # For wrapping around the Prime Meridian
# if lng_gmi[-1] != 360:
#     lng_gmi[-1] = 360.    
# # O3 on days with a poleward jet
# ax1 = plt.subplot2grid((2,2), (0,0), colspan=2,
#     projection=ccrs.Miller(central_longitude=0.))
# ax1.set_title('(a) O$_\mathregular{3,\:PW}$', fontsize=16, x=0.02, ha='left')
# # ax1.add_feature(ocean50m, zorder=3)
# ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
# ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
#     lat_gmi.min()+1, lat_gmi.max()-5])
# cmap = plt.get_cmap('OrRd')
# mb = ax1.contourf(lng_gmi, lat_gmi, pwjet_o3, np.linspace(30, 60, 7), 
#     cmap=cmap, extend='both', transform=ccrs.PlateCarree(), zorder=1)
# # Eddy-driven jet
# skiplng = 6
# ax1.errorbar(lng_gmi[::skiplng], np.array(pwjet_lat)[::skiplng], 
#     yerr=np.array(pwjet_lat_var)[::skiplng], color='k', markersize=3, 
#     elinewidth=1.25, ecolor='k', fmt='o', transform=ccrs.PlateCarree(), 
#     zorder=5)
# ax1.outline_patch.set_zorder(20)
# # O3 on days with an equatorward jet
# ax2 = plt.subplot2grid((2,2), (1,0), colspan=2,
#     projection=ccrs.Miller(central_longitude=0.))
# ax2.set_title('(b) O$_\mathregular{3,\:EW}$', fontsize=16, x=0.02, ha='left')
# # ax2.add_feature(ocean50m)
# ax2.coastlines(lw=0.25, resolution='50m', color='k', zorder=4)
# ax2.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
#     lat_gmi.min()+1, lat_gmi.max()-5])
# mb = ax2.contourf(lng_gmi, lat_gmi, eqjet_o3, np.linspace(30, 60, 7), 
#     cmap=cmap, extend='both', transform=ccrs.PlateCarree(), zorder=1)
# ax2.errorbar(lng_gmi[::skiplng], np.array(eqjet_lat)[::skiplng], 
#     yerr=np.array(eqjet_lat_var)[::skiplng], color='k', markersize=3, 
#     elinewidth=1.25, ecolor='k', fmt='o', transform=ccrs.PlateCarree(), 
#     zorder=5)
# ax1.outline_patch.set_zorder(20)
# ax2.outline_patch.set_zorder(20)      
# # Add colorbar
# plt.gcf().subplots_adjust(left=0.02, right=0.86, hspace=0.3)
# colorbar_axes = plt.gcf().add_axes([
#     ax1.get_position().x1+0.03, # Left
#     (ax2.get_position().y1-ax2.get_position().y0)/2.+ax2.get_position().y0, # Bottom 
#     0.02, # Width
#     ((ax1.get_position().y1-ax1.get_position().y0)/2.+ax1.get_position().y0)-
#     ((ax2.get_position().y1-ax2.get_position().y0)/2.+ax2.get_position().y0)])
# colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
#     ticks=np.linspace(30, 60, 7), extend='both')
# colorbar.ax.tick_params(labelsize=12)
# colorbar.set_label('[ppbv]', fontsize=16)
# plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#     'o3pwew_oceanoff.pdf', dpi=600)

"""FIND ZONALLY-AVERAGED 03, U10, AND V10 FOR DIFFERENT REGIONS ON PW/EW DAYS"""
# import matplotlib.pyplot as plt
# import sys
# sys.path.append('/Users/ghkerr/phd/GMI/')
# from geo_idx import geo_idx            
# # Find fields in region 
# o3_region, lat_region, lng_region = globalo3_calculate.find_grid_in_bb(
#     o3_gmi, lat_gmi, lng_gmi, 140., 180., 25., 70.)
# t2m_region, lat_region, lng_region = globalo3_calculate.find_grid_in_bb(
#     t2m_merra, lat_gmi, lng_gmi, 140., 180., 25., 70.)
# qv2m_region, lat_region, lng_region = globalo3_calculate.find_grid_in_bb(
#     qv2m_merra, lat_gmi, lng_gmi, 140., 180., 25., 70.)
# V10M_region, lat_region, lng_region = globalo3_calculate.find_grid_in_bb(
#     V10M, lat_gmi, lng_gmi, 140., 180., 25., 70.)
# ena_left = geo_idx(140., lng_ml)
# ena_right = geo_idx(180., lng_ml)
# lat_jet_region = lat_jet_ml[:, ena_left:ena_right+1]
# # Separate based on jet latitude
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_o3, eqjet_o3 = \
#     globalo3_calculate.segregate_field_bylat(o3_region, lng_region, 
#     lat_jet_region, times_gmi)
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_t2m, eqjet_t2m = \
#     globalo3_calculate.segregate_field_bylat(t2m_region, lng_region, 
#     lat_jet_region, times_gmi)
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_qv2m, eqjet_qv2m = \
#     globalo3_calculate.segregate_field_bylat(qv2m_region, lng_region, 
#     lat_jet_region, times_gmi)
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_V10M, eqjet_V10M = \
#     globalo3_calculate.segregate_field_bylat(V10M_region, lng_region, 
#     lat_jet_region, times_gmi)
# ax1 = plt.subplot2grid((2,2),(0,0))
# ax2 = plt.subplot2grid((2,2),(0,1))
# ax3 = plt.subplot2grid((2,2),(1,0))
# ax4 = plt.subplot2grid((2,2),(1,1))
# # O3
# ax1.plot(lat_region, np.nanmean(o3_region, axis=tuple((0,2))), '-k')
# ax1.plot(lat_region, np.nanmean(pwjet_o3, axis=-1), '-r')
# ax1.plot(lat_region, np.nanmean(eqjet_o3, axis=-1), '-g')
# ax1.set_title('O3')
# ax1.set_ylabel('ppbv')
# # T
# ax2.plot(lat_region, np.nanmean(t2m_region, axis=tuple((0,2))), '-k')
# ax2.plot(lat_region, np.nanmean(pwjet_t2m, axis=-1), '-r')
# ax2.plot(lat_region, np.nanmean(eqjet_t2m, axis=-1), '-g')
# ax2.set_title('T')
# ax2.set_ylabel('K')
# # q
# ax3.plot(lat_region, np.nanmean(qv2m_region, axis=tuple((0,2))), '-k')
# ax3.plot(lat_region, np.nanmean(pwjet_qv2m, axis=-1), '-r')
# ax3.plot(lat_region, np.nanmean(eqjet_qv2m, axis=-1), '-g')
# ax3.set_title('q')
# ax3.set_ylabel('g kg-1')
# ax3.set_xlabel('deg')
# # 10-meter meridional wind
# ax4.plot(lat_region, np.nanmean(V10M_region, axis=tuple((0,2))), '-k', label='all')
# ax4.plot(lat_region, np.nanmean(pwjet_V10M, axis=-1), '-r', label='PW jet')
# ax4.plot(lat_region, np.nanmean(eqjet_V10M, axis=-1), '-g', label='EW jet')
# ax4.legend(loc=2, bbox_to_anchor=(1.05, 1.))
# ax4.set_title('V10')
# ax4.set_ylabel('m s-1')
# ax4.set_xlabel('deg')
# plt.subplots_adjust(left=0.1, hspace=0.4, wspace=0.35, right=0.8)
# plt.savefig('/Users/ghkerr/Desktop/zonalavg_pwewjet_eastp.png', dpi=300)

"""LAND VERSUS OCEAN ZONALLY-AVERAGED O3, TEMPERATURE, AND HUMIDITY"""
# land = globalo3_calculate.find_grid_overland(lat_gmi[0:-5], lng_gmi)
# wherenan = np.where(land != 1)
# ocean = np.empty(shape=land.shape)
# ocean[:] = np.nan
# ocean[wherenan] = 1.
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_o3, eqjet_o3 = \
#     globalo3_calculate.segregate_field_bylat(o3_gmi, lng_gmi, lat_jet_ml, 
#     times_gmi)
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_t2m, eqjet_t2m = \
#     globalo3_calculate.segregate_field_bylat(t2m_merra, lng_gmi, lat_jet_ml, 
#     times_gmi)
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_qv2m, eqjet_qv2m = \
#     globalo3_calculate.segregate_field_bylat(qv2m_merra, lng_gmi, lat_jet_ml, 
#     times_gmi)
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_U10M, eqjet_U10M = \
#     globalo3_calculate.segregate_field_bylat(U10M, lng_gmi, lat_jet_ml, 
#     times_gmi)
# eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_V10M, eqjet_V10M = \
#     globalo3_calculate.segregate_field_bylat(V10M, lng_gmi, lat_jet_ml, 
#     times_gmi)
# fig = plt.figure(figsize=(6,8))
# ax1 = plt.subplot2grid((3,2),(0,0))
# ax2 = plt.subplot2grid((3,2),(0,1))
# ax3 = plt.subplot2grid((3,2),(1,0))
# ax4 = plt.subplot2grid((3,2),(1,1))
# ax5 = plt.subplot2grid((3,2),(2,0))
# ax6 = plt.subplot2grid((3,2),(2,1))
# latweight = np.cos(np.deg2rad(lat_gmi[:-5]))
# # Land U10M
# ax1.plot(lat_gmi[:-5], np.nanmean(U10M[:,:-5]*land, axis=tuple((0,2))),#*latweight, 
#    lw=2, color='k')
# ax1.plot(lat_gmi[:-5], np.nanmean(pwjet_U10M[:-5]*land, axis=-1),#*latweight, 
#    lw=2, color='#3F79B7')
# ax1.plot(lat_gmi[:-5], np.nanmean(eqjet_U10M[:-5]*land, axis=-1),#*latweight, 
#    lw=2, color='#EF9802')
# # Ocean U10M
# ax2.plot(lat_gmi[:-5], np.nanmean(U10M[:,:-5]*ocean, axis=tuple((0,2))),#*latweight, 
#    lw=2, color='k')
# ax2.plot(lat_gmi[:-5], np.nanmean(pwjet_U10M[:-5]*ocean, axis=-1),#*latweight, 
#    lw=2, color='#3F79B7')
# ax2.plot(lat_gmi[:-5], np.nanmean(eqjet_U10M[:-5]*ocean, axis=-1),#*latweight, 
#    lw=2, color='#EF9802')
# for ax in [ax1, ax2]:
#     ax.set_xlim([25, 70])
#     ax.set_xticks([25, 40, 55, 70])
#     ax.set_xticklabels([''])
# # ax1.set_ylim([-0.5, 1])
# # ax1.set_yticks([-0.5, -0.2,  0.1,  0.4,  0.7,  1. ])
# # ax1.set_yticklabels([-0.5, -0.2,  0.1,  0.4,  0.7,  1. ], fontsize=12)
# # ax2.set_ylim([-3, 2])    
# # ax2.set_yticks([-3, -2, -1, 0, 1, 2])
# # ax2.set_yticklabels([-3, -2, -1, 0, 1, 2], fontsize=12)
# ax1.set_ylabel('U$_{\mathregular{10}}$ [m s$^{\mathregular{-1}}$]', 
#     fontsize=16)
# ax1.set_title('Land', fontsize=16)
# ax2.set_title('Ocean', fontsize=16)
# # Land V10M
# ax3.plot(lat_gmi[:-5], np.nanmean(V10M[:,:-5]*land, axis=tuple((0,2))),#*latweight, 
#    lw=2, color='k')
# ax3.plot(lat_gmi[:-5], np.nanmean(pwjet_V10M[:-5]*land, axis=-1),#*latweight, 
#    lw=2, color='#3F79B7')
# ax3.plot(lat_gmi[:-5], np.nanmean(eqjet_V10M[:-5]*land, axis=-1),#*latweight, 
#    lw=2, color='#EF9802')
# # Ocean U10M
# ax4.plot(lat_gmi[:-5], np.nanmean(V10M[:,:-5]*ocean, axis=tuple((0,2))),#*latweight, 
#    lw=2, color='k')
# ax4.plot(lat_gmi[:-5], np.nanmean(pwjet_V10M[:-5]*ocean, axis=-1),#*latweight, 
#    lw=2, color='#3F79B7')
# ax4.plot(lat_gmi[:-5], np.nanmean(eqjet_V10M[:-5]*ocean, axis=-1),#*latweight, 
#    lw=2, color='#EF9802')
# for ax in [ax3, ax4]:
#     ax.set_xlim([25, 70])
#     ax.set_xticks([25, 40, 55, 70])
#     ax.set_xticklabels([''])
# # ax3.set_ylim([-0.7, 0.3])    
# # ax3.set_yticks([-0.7, -0.5, -0.3, -0.1,  0.1,  0.3])
# # ax3.set_yticklabels([-0.7, -0.5, -0.3, -0.1,  0.1,  0.3], fontsize=12)
# # ax4.set_ylim([-0.5, 1.5])    
# # ax4.set_yticks([-0.5, -0.1,  0.3,  0.7,  1.1,  1.5])
# # ax4.set_yticklabels([-0.5, -0.1,  0.3,  0.7,  1.1,  1.5], fontsize=12)
# ax3.set_ylabel('V$_{\mathregular{10}}$ [m s$^{\mathregular{-1}}$]', fontsize=16)
# # Land V10M
# ax5.plot(lat_gmi[:-5], np.nanmean(o3_gmi[:,:-5]*land, axis=tuple((0,2))),#*latweight, 
#    lw=2, color='k')
# ax5.plot(lat_gmi[:-5], np.nanmean(pwjet_o3[:-5]*land, axis=-1),#*latweight, 
#    lw=2, color='#3F79B7')
# ax5.plot(lat_gmi[:-5], np.nanmean(eqjet_o3[:-5]*land, axis=-1),#*latweight, 
#    lw=2, color='#EF9802')
# # Ocean U10M
# ax6.plot(lat_gmi[:-5], np.nanmean(o3_gmi[:,:-5]*ocean, axis=tuple((0,2))),#*latweight,
#     lw=2, color='k', label='All')
# ax6.plot(lat_gmi[:-5], np.nanmean(pwjet_o3[:-5]*ocean, axis=-1),#*latweight,
#     lw=2, color='#3F79B7', label='PW')
# ax6.plot(lat_gmi[:-5], np.nanmean(eqjet_o3[:-5]*ocean, axis=-1),#*latweight,
#     lw=2, color='#EF9802', label='EW')
# for ax in [ax5, ax6]:
#     ax.set_xlim([25, 70])
#     ax.set_xticks([25, 40, 55, 70])
#     ax.set_xticklabels(['25', '40', '55', '70'], fontsize=12)
#     ax.set_xlabel('Latitude [$^{\mathregular{\circ}}$N]', fontsize=16)
# # ax5.set_ylim([5, 55])    
# # ax5.set_yticks([5., 15., 25., 35., 45., 55.])
# # ax5.set_yticklabels([5, 15, 25, 35, 45, 55], fontsize=12)
# # ax6.set_ylim([6, 26])    
# # ax6.set_yticks([6., 10., 14., 18., 22., 26.])
# # ax6.set_yticklabels([6, 10, 14, 18, 22, 26], fontsize=12)    
# ax5.set_ylabel('O$_{\mathregular{3}}$ [ppbv]', fontsize=16)
# plt.subplots_adjust(wspace=0.3, top=0.95, bottom=0.15)
# # Add legend 
# leg = ax6.legend(loc=6, ncol=3, bbox_to_anchor=(-1.15, -0.5), 
#     fontsize=16, frameon=False)
# # Add labels to each subplot
# ax1.text(0.04, 0.87, '(a)', ha='left', transform=ax1.transAxes, 
#     fontsize=16, zorder=20)
# ax2.text(0.04, 0.87, '(b)', ha='left', transform=ax2.transAxes, 
#     fontsize=16, zorder=20)
# ax3.text(0.04, 0.87, '(c)', ha='left', transform=ax3.transAxes, 
#     fontsize=16, zorder=20)
# ax4.text(0.04, 0.87, '(d)', ha='left', transform=ax4.transAxes, 
#     fontsize=16, zorder=20)
# ax5.text(0.04, 0.87, '(e)', ha='left', transform=ax5.transAxes, 
#     fontsize=16, zorder=20)
# ax6.text(0.04, 0.87, '(f)', ha='left', transform=ax6.transAxes, 
#     fontsize=16, zorder=20)
# plt.savefig('/Users/ghkerr/Desktop/zonalavg_unweighted.png', dpi=300)