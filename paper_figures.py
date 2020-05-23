#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script generates the figures used in Kerr et al. (2020) "Surface ozone-
meteorology relationships: Spatial vaariations and the role of the jet stream"

Revision History
    05022020 -- initial version created
    07022020 -- edit function 'fig8' to only include the rotated O3 anomalies
                rather than also unrotated
    11022020 -- function 'fig5' edited to show difference of correlations 
                between control and transport-only simulations
    14022020 -- happy Valentine's day! Changed map projections, added lat/lon
                ticks, and made several other small edits
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
    
def fig1(lat_gmi, lng_gmi, o3_gmi, lat_jet_ml):
    """Figure 1 of Kerr et al. (2020). (a) Mean afternoon O3 from the 
    surface-level of the GMI CTM during summer (JJA) is shown with colored 
    shading. Black contours indicate O3 variability (standard deviation): thin 
    dashed contour, 8 ppbv; thick contour, 10 ppbv. (b) Annual anthropogenic 
    NOx emissions from EDGAR averaged over 2008-2010. Scatterpoints and error 
    bars in (a-b) specify the mean position and variability of the eddy-driven 
    jet, respectively.

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    o3_gmi : numpy.ndarray     
        Daily afternoon surface-level O3 from the GMI CTM, units of ppbv, 
        [time, lat, lng]        
    lat_jet_ml : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in the Northern Hemisphere mid-latitudes, units of degrees north, 
        [time, lng]

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
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    # Load EDGAR NOx
    lat_edgar, lng_edgar, nox_edgar = globalo3_open.open_edgar_specifieddomain(
        [2008, 2009, 2010], -1., 90., 0., 360., 'NOx')
    # For wrapping around the Prime Meridian
    if lng_edgar[-1] != 360:
        lng_edgar[-1] = 360.
    if lng_gmi[-1] != 360:
        lng_gmi[-1] = 360.        
    # Convert from kg NOx m-2 s-1 to kg NOx km-2 day-1
    nox_edgar = np.nanmean(nox_edgar, axis=0)*3600.*1000000.*24.
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,5))
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax1.set_title('(a) O$_{\mathregular{3}}$', fontsize=16, x=0.02, ha='left')    
    # ax1.add_feature(ocean50m)
    ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=4)
    ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax1.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax1.xaxis.set_major_formatter(lng_formatter)         
    ax1.get_xaxis().set_ticklabels([])
    ax1.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax1.yaxis.set_major_formatter(lat_formatter)
    cmap = plt.get_cmap('OrRd')
    mb = ax1.contourf(lng_gmi, lat_gmi, np.mean(o3_gmi, axis=0), 
        np.linspace(10, 60, 11), cmap=cmap, extend='both',
        transform=ccrs.PlateCarree(), zorder=1)
    csthick = ax1.contour(lng_gmi, lat_gmi, np.nanstd(o3_gmi, axis=0), [10.], 
        colors='k', linewidths=1.5, transform=ccrs.PlateCarree(), zorder=15)
    csmedium = ax1.contour(lng_gmi, lat_gmi, np.nanstd(o3_gmi, axis=0), [8.], 
        colors='k', linestyles='--', linewidths=0.75, 
        transform=ccrs.PlateCarree(), zorder=15)
    # skiplng = 6
    # ax1.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
    #     yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], zorder=10, color='k', 
    #     markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
    #     transform=ccrs.PlateCarree())
    plt.gcf().subplots_adjust(left=0.05, right=0.86, hspace=0.3)    
    colorbar_axes = plt.gcf().add_axes([ax1.get_position().x1+0.03, 
        ax1.get_position().y0, 0.02, (ax1.get_position().y1-
        ax1.get_position().y0)])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(10, 60, 6), extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[ppbv]', fontsize=16, labelpad=11)
    ax1.outline_patch.set_zorder(20)
    # NOx emissions
    ax2=plt.subplot2grid((2,2), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2.set_title('(b) NO$_{x}$', fontsize=16, x=0.02, ha='left')    
    # ax2.add_feature(ocean50m)
    ax2.coastlines(lw=0.25, resolution='50m', color='k', zorder=4)
    ax2.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax2.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax2.xaxis.set_major_formatter(lng_formatter)      
    ax2.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax2.yaxis.set_major_formatter(lat_formatter)    
    mb = ax2.pcolormesh(lng_edgar, lat_edgar, nox_edgar, 
        cmap=plt.get_cmap('Blues', len(np.linspace(0, 8, 9))-1), 
        vmin=np.linspace(0, 8, 9)[0], vmax=np.linspace(0, 8, 9)[-1], 
        norm=mpl.colors.BoundaryNorm(np.linspace(0, 8, 9), 
        ncolors=plt.get_cmap('Blues', len(np.linspace(0, 8, 9))-1).N, 
        clip=False),transform=ccrs.PlateCarree(), rasterized=True)
    skiplng = 6
    ax2.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
        yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], zorder=10, color='k', 
        markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
        transform=ccrs.PlateCarree())
    # Add colorbar
    colorbar_axes = plt.gcf().add_axes([ax2.get_position().x1+0.03, 
        ax2.get_position().y0, 0.02, (ax2.get_position().y1-
        ax2.get_position().y0)])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(0, 8, 9), extend='max')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[kg km$^{\mathregular{-2}}$ day$^{\mathregular{-1}}$]', 
        labelpad=15, fontsize=16)
    ax2.outline_patch.set_zorder(20)      
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'fig1_nojet.pdf', dpi=600)
    return

def fig2(lat_gmi, lng_gmi, r_aqs, r_naps, r_emep, r_china): 
    """Figure 2 of Kerr et al. (2020). Plot maps of the correlation coefficient 
    (a-c) calculated between modeled O3 and observed O3 from AQS, NAPS, EMEP, 
    and MEE.
    
    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]        
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    r_aqs : numpy.ndarray
        Correlation coefficients (r) calculated between O3 from AQS and GMI 
        CTM, [lat, lng]
    r_naps : numpy.ndarray
        Correlation coefficients (r) calculated between O3 from NAPS and GMI 
        CTM, [lat, lng]
    r_emep : numpy.ndarray
        Correlation coefficients (r) calculated between O3 from EMEP and GMI 
        CTM, [lat, lng]
    r_china : numpy.ndarray
        Correlation coefficients (r) calculated between O3 from China and GMI 
        CTM, [lat, lng]

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
    fig = plt.figure(figsize=(4,6))
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
    ax1.set_title('(a) r NAPS, AQS', fontsize=16, x=0.02, ha='left')  
    ax2.set_title('(b) r EMEP', fontsize=16, x=0.02, ha='left')
    ax3.set_title('(c) r MEE', fontsize=16, x=0.02, ha='left')
    # Add oceans and coasts
    for ax in [ax1, ax2, ax3]:
        ocean = cfeature.NaturalEarthFeature(category='physical', name='ocean',
                                    scale='50m', facecolor='lightgrey')
        ax.add_feature(ocean)
        ax.coastlines(lw=0.25, color='k', zorder=3, resolution='50m')
    # Set extent
    ax1.set_extent([-165, -50, 24, 55])    
    ax2.set_extent([-14, 37, 34, 70])
    ax3.set_extent([95, 139, 15, 53])
    # Define colormap 
    cmap = plt.get_cmap('coolwarm')
    # Extract all colors from the colormap
    cmaplist = cmap(np.linspace(0.5, 1, cmap.N // 2))
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist, cmap.N)
    # Correlation plot
    clevs = np.linspace(0, 1, 11)
    norm = mpl.colors.BoundaryNorm(clevs, cmap.N)
    # North America
    mb = ax1.pcolor(lng_gmi, lat_gmi, r_aqs, cmap=cmap, norm=norm, 
        transform=ccrs.PlateCarree(), zorder=20)
    ax1.pcolor(lng_gmi, lat_gmi, r_naps, cmap=cmap, norm=norm, 
        transform=ccrs.PlateCarree(), zorder=20)
    # Europe
    ax2.pcolor(lng_gmi, lat_gmi, r_emep, cmap=cmap, norm=norm, 
        transform=ccrs.PlateCarree(), zorder=20)
    # China
    ax3.pcolor(lng_gmi, lat_gmi, r_china, cmap=cmap, norm=norm, 
        transform=ccrs.PlateCarree(), zorder=20)
    # Add colorbars
    colorbar_axes = plt.gcf().add_axes([0.05+0.05,0.11,0.8,0.03])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='horizontal', 
        extend='neither', ticks=clevs)
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[$\mathregular{\cdot}$]', fontsize=15)
    plt.gcf().subplots_adjust(left=0.05, right=0.95, bottom=0.17)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'fig2.pdf', dpi=600)
    return

def fig3(lat_gmi, lng_gmi, r_t2mo3, r_qv2mo3, significance_r_t2mo3, 
    significance_r_qv2mo3, lat_jet_ml): 
    """Figure 3 of Kerr et al. (2020). (a) r(T, O3) is shown with colored 
    shading. Hatching indicates insignificance determined with moving block 
    bootstrapping. (b) Same as (a) but showing r(q, O3). Scatterpoints and 
    error bars in (a-b) specify the mean position and variability of the eddy-
    driven jet, respectively.

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    r_t2mo3 : numpy.ndarray     
        Pearson correlation coefficient between 2-meter temperature and O3, 
        [lat, lng]
    r_qv2mo3 : numpy.ndarray     
        Pearson correlation coefficient between 2-meter specific humidity and 
        O3, [lat, lng]
    significance_r_t2mo3 : numpy.ndarray        
        Significance of r(T, O3) determined with moving block bootstrapping: 1
        implies significance, NaN implies insignificance, [lat, lng]
    significance_r_qv2mo3 : numpy.ndarray
        Significance of r(q, O3) determined with moving block bootstrapping: 1
        implies significance, NaN implies insignificance, [lat, lng]    
    lat_jet_ml : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in the Northern Hemisphere mid-latitudes, units of degrees north, 
        [time, lng]

    Returns
    -------
    None    
    """
    import numpy as np
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth']=0.3     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.patches as mpatches
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter    
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,5))
    # For wrapping around the Prime Meridian
    if lng_gmi[-1] != 360:
        lng_gmi[-1] = 360.    
    # r(T, O3)
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax1.set_title(r'(a) r(T, O$_{\mathregular{3}}$)', fontsize=16, 
        x=0.02, ha='left')    
    # ax1.add_feature(ocean50m, zorder=3)
    ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax1.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax1.xaxis.set_major_formatter(lng_formatter)         
    ax1.get_xaxis().set_ticklabels([])    
    ax1.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()    
    ax1.yaxis.set_major_formatter(lat_formatter)    
    cmap = plt.get_cmap('coolwarm')
    mb = ax1.contourf(lng_gmi, lat_gmi, r_t2mo3, np.linspace(-1, 1, 9), 
        cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(T, O3)
    ax1.contourf(lng_gmi, lat_gmi, significance_r_t2mo3, hatches=['//////'], 
        colors='none', transform=ccrs.PlateCarree(), zorder=2)
    # Add boxes corresponding to regions from Figure 4
    # Western North America
    ax1.add_patch(mpatches.Rectangle(xy=[-125, 25], width=24, height=45,
        fill=None, lw=2, edgecolor='k', transform=ccrs.PlateCarree(), 
        zorder=3))
    ax1.text(-130, 15, 'Western\nNorth\nAmerica', ha='right',
        fontweight='bold', fontsize=12, transform=ccrs.PlateCarree())
    # Eastern North America    
    ax1.add_patch(mpatches.Rectangle(xy=[-99, 25], width=24, height=45,
        fill=None, lw=2, edgecolor='k', transform=ccrs.PlateCarree(), 
        zorder=3))
    ax1.text(-70, 15, 'Eastern\nNorth\nAmerica', ha='left',
        fontweight='bold', fontsize=12, transform=ccrs.PlateCarree())    
    # European Union 
    ax1.add_patch(mpatches.Rectangle(xy=[-10, 25], width=40, height=45,
        fill=None, lw=2, edgecolor='k', transform=ccrs.PlateCarree(), 
        zorder=3))
    ax1.text(-10, 74, 'Europe', ha='left', fontsize=12, 
        fontweight='bold', transform=ccrs.PlateCarree())      
    # China
    ax1.add_patch(mpatches.Rectangle(xy=[90, 25], width=35, height=45,
        fill=None, lw=2, edgecolor='k', transform=ccrs.PlateCarree(), 
        zorder=3))     
    ax1.text(128, 63, 'China', ha='left', fontsize=12, 
        fontweight='bold', transform=ccrs.PlateCarree())
    # Eddy-driven jet
    skiplng = 6
    ax1.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
        yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], color='k', markersize=3, 
        elinewidth=1.25, ecolor='k', fmt='o', transform=ccrs.PlateCarree(), 
        zorder=5)
    ax1.outline_patch.set_zorder(20)
    # r(q, O3)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2.set_title(r'(b) r(q, O$_{\mathregular{3}}$)', fontsize=16, x=0.02, 
        ha='left')
    # ax2.add_feature(ocean50m)
    ax2.coastlines(lw=0.25, resolution='50m', color='k', zorder=4)
    ax2.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax2.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax2.xaxis.set_major_formatter(lng_formatter)          
    ax2.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()    
    ax2.yaxis.set_major_formatter(lat_formatter)
    mb = ax2.contourf(lng_gmi, lat_gmi, r_qv2mo3, np.linspace(-1,1,9), 
        cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmi, lat_gmi, significance_r_qv2mo3, hatches=['//////'], 
        colors='none', transform=ccrs.PlateCarree())
    # skiplng = 6
    ax2.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
        yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], zorder=10, color='k', 
        markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
        transform=ccrs.PlateCarree())
    ax2.outline_patch.set_zorder(20)      
    # Add colorbar
    plt.gcf().subplots_adjust(left=0.02, right=0.86, hspace=0.3)
    colorbar_axes = plt.gcf().add_axes([
        ax1.get_position().x1+0.03, # Left
        (ax2.get_position().y1-ax2.get_position().y0)/2.+ax2.get_position().y0, # Bottom 
        0.02, # Width
        ((ax1.get_position().y1-ax1.get_position().y0)/2.+ax1.get_position().y0)-
        ((ax2.get_position().y1-ax2.get_position().y0)/2.+ax2.get_position().y0)])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-1, 1, 9), extend='neither')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[$\mathregular{\cdot}$]', fontsize=16, labelpad=11)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'fig3.pdf', dpi=600)

def fig4(lat_gmi, lng_gmi, r_t2mo3, r_qv2mo3, r_t2mo3_aqs, r_qv2mo3_aqs, 
    lat_aqs, lng_aqs, r_t2mo3_naps, r_qv2mo3_naps, lat_naps, lng_naps, 
    r_t2mo3_emep, r_qv2mo3_emep, lat_emep, lng_emep, r_t2mo3_china, 
    r_qv2mo3_china, lat_china, lng_china, lng_ml, lat_jet_ml):
    """Plot zonally-averaged mean values of r(T, O3) and r(q, O3) in Western 
    North America, Eastern North America, Europe, and China from both GMI CTM 
    simulations and observations and plot as a function of latitude alongside 
    the mean position and variability of the eddy-driven jet. 
    
    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    r_t2mo3 : numpy.ndarray     
        Pearson correlation coefficient between 2-meter temperature and O3, 
        [lat, lng]
    r_qv2mo3 : numpy.ndarray     
        Pearson correlation coefficient between 2-meter specific humidity and 
        O3, [lat, lng]
    r_t2mo3_aqs : list
        Pearson correlation coefficient between 2-meter temperature and O3 at
        AQS stations, [stations,]
    r_qv2mo3_aqs : list
        Pearson correlation coefficient between 2-meter specific humidity and 
        O3 at AQS stations, [stations,]
    lat_aqs : list 
        Latitude coorindate at each AQS station, units of degrees north, 
        [stations,]
    lng_aqs : list 
        Longitude coorindate at each AQS station, units of degrees east, 
        [stations,]
    r_t2mo3_naps : list
        Pearson correlation coefficient between 2-meter temperature and O3 at
        NAPS stations, [stations,]
    r_qv2mo3_naps : list
        Pearson correlation coefficient between 2-meter specific humidity and 
        O3 at NAPS stations, [stations,]
    lat_naps : list 
        Latitude coorindate at each NAPS station, units of degrees north, 
        [stations,]
    lng_naps : list 
        Longitude coorindate at each NAPS station, units of degrees east, 
        [stations,]        
    r_t2mo3_emep : list
        Pearson correlation coefficient between 2-meter temperature and O3 at
        EMEP stations, [stations,]
    r_qv2mo3_emep : list
        Pearson correlation coefficient between 2-meter specific humidity and 
        O3 at EMEP stations, [stations,]
    lat_emep : list 
        Latitude coorindate at each EMEP station, units of degrees north, 
        [stations,]
    lng_emep : list 
        Longitude coorindate at each EMEP station, units of degrees east, 
        [stations,]
    r_t2mo3_china : list
        Pearson correlation coefficient between 2-meter temperature and O3 at
        Chinese stations, [stations,]
    r_qv2mo3_china : list
        Pearson correlation coefficient between 2-meter specific humidity and 
        O3 at Chinese stations, [stations,]        
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

    Returns
    -------
    None
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx    
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth'] = 0.75
    # Find quantity (model, model error, and observations) in regions and 
    # calculate land-based zonally-averaged quantities
    # Western North America
    # Model
    field_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(
        r_t2mo3, lat_gmi, lng_gmi, 235., 260., 25., 70.) 
    # # Create grid to sub-sample model; the step to create the empty grid 
    # # occurs now (to get the correct shape of "field_wna" before the zonal 
    # # average is taken), and this grid will be filled with NaNs for grid cells
    # # without monitors and 1s for grid cells with monitors
    # subsample_wna = np.empty(shape=field_wna.shape)
    wna_left = geo_idx(235., lng_ml)
    wna_right = geo_idx(260., lng_ml)
    lat_jet_wna = lat_jet_ml[:, wna_left:wna_right+1]
    lat_jet_err_wna = np.nanstd(np.nanmean(lat_jet_wna, axis = 1))
    lat_jet_mean_wna = np.nanmean(lat_jet_wna)
    # Observations
    obs_wna = np.where((np.array(lat_naps+lat_aqs) < 70.) &
                     (np.array(lat_naps+lat_aqs) > 25.) &
                     (np.array(lng_naps+lng_aqs) > 235.) &
                     (np.array(lng_naps+lng_aqs) < 260.))[0]
    # Combine NAPS and AQS observations 
    lat_obs_wna = np.array(lat_naps+lat_aqs)[obs_wna]
    obs_wna = np.array(r_t2mo3_naps+r_t2mo3_aqs)[obs_wna]
    # # OPTIONAL: Select grid cells containing in-situ monitors (if this is not
    # # preferred, then go the route of finding land-based grid cells)
    # subsample_wna[:] = np.nan
    # for ilat, ilng in zip(lat_obs_wna, lng_obs_wna):
    #     ilat = geo_idx(ilat, lat_wna)
    #     ilng = geo_idx(ilng, lng_wna)
    #     subsample_wna[ilat, ilng] = 1.   
    land_wna = globalo3_calculate.find_grid_overland(lat_wna, lng_wna)
    field_err_wna = 2*np.nanstd(field_wna*land_wna, axis=1)
    field_wna = np.nanmean(field_wna*land_wna, axis=1)
    # Bin observations by latitude 
    lat_obs_wna, obs_mean_wna, obs_std_wna = \
        globalo3_calculate.bin_observations_bylat(lat_wna[::2], obs_wna, 
        lat_obs_wna)
    # Eastern North America
    # Model 
    field_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
        r_t2mo3, lat_gmi, lng_gmi, 260., 295., 25., 70.)
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
    obs_ena = np.array(r_t2mo3_naps+r_t2mo3_aqs)[obs_ena]
    lat_obs_ena, obs_mean_ena, obs_std_ena = \
        globalo3_calculate.bin_observations_bylat(lat_ena[::2], obs_ena, 
        lat_obs_ena)
    # Europe (since Europe crosses the prime meridian, finding grid cells over
    # 0 deg longitude will not work, so find European domain in two steps)
    # Model 
    field_eu1, lat_eu1, lng_eu1 = globalo3_calculate.find_grid_in_bb(
        r_t2mo3, lat_gmi, lng_gmi, 350., 360., 25., 70.) 
    eu1_left = geo_idx(350., lng_ml)
    eu1_right = geo_idx(360., lng_ml)
    lat_jet_eu1 = lat_jet_ml[:, eu1_left:eu1_right+1]
    land_eu1 = globalo3_calculate.find_grid_overland(lat_eu1, lng_eu1)
    field_eu2, lat_eu2, lng_eu2 = globalo3_calculate.find_grid_in_bb(
        r_t2mo3, lat_gmi, lng_gmi, 0., 30., 25., 70.) 
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
    obs_eu = np.array(r_t2mo3_emep)[emep_eu]
    lat_obs_eu = np.array(lat_emep)[emep_eu]
    lat_obs_eu, obs_mean_eu, obs_std_eu = \
        globalo3_calculate.bin_observations_bylat(lat_eu[::2], obs_eu, 
        lat_obs_eu)
    # China
    # Model 
    field_asia, lat_asia, lng_asia = globalo3_calculate.find_grid_in_bb(
    r_t2mo3, lat_gmi, lng_gmi, 90., 125., 25., 70.) 
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
    obs_asia = np.array(r_t2mo3_china)[obs_asia]
    lat_obs_china, obs_mean_china, obs_std_china = \
        globalo3_calculate.bin_observations_bylat(lat_asia[::2], obs_asia, 
        lat_obs_asia)
    # Plotting
    fig = plt.figure(figsize=(13, 6.))
    ax1 = plt.subplot2grid((4,4), (0,0), colspan=2)
    ax2 = plt.subplot2grid((4,4), (1,0), colspan=2)
    ax3 = plt.subplot2grid((4,4), (2,0), colspan=2)
    ax4 = plt.subplot2grid((4,4), (3,0), colspan=2)
    color_r = '#EF9802'
    color_jet = 'k'
    color_obs = '#3F79B7'
    for ax in [ax1, ax2, ax3, ax4]: 
        ax.set_xlim([25, 70])
        ax.set_xticks(np.linspace(25, 70, 10))
        ax.set_xticklabels([])
        ax.set_ylim([np.linspace(-0.5, 1, 7)[0], np.linspace(-0.5, 1, 7)[-1]])
        ax.set_yticks(np.linspace(-0.5, 1, 7))
        # Add horizontal line for field = 0
        ax.axhline(y = 0.0, color = 'k', lw = 1., linestyle = '--')        
    ax4.set_xticklabels([ '25', '', '35', '', '45', '', '55', '', '65', ''])    
    ax4.set_xlabel('Latitude [$^{\mathregular{\circ}}$N]', fontsize=16)
    ax2.set_ylabel('r(T, O$_{\mathregular{3}}$) [$\cdot$]', fontsize=16, 
        y=-0.15)    
    # Western North America
    ax1.text(0.98, 0.8, 'Western North America', ha = 'right', 
        transform=ax1.transAxes, fontsize=12, zorder=20)
    p1, = ax1.plot(lat_wna, field_wna, ls='-', lw=2, color=color_r, 
        zorder=6)
    p2 = ax1.fill_between(lat_wna, field_wna-field_err_wna,
        field_wna+field_err_wna, color = color_r, alpha=0.2, zorder=5)
    p3, = ax1.plot(lat_obs_wna, obs_mean_wna, '-', lw=2, color=color_obs,
    zorder=4)
    p4 = ax1.fill_between(lat_obs_wna, (obs_mean_wna-obs_std_wna),
        (obs_mean_wna+obs_std_wna), facecolor='None', hatch='////', 
        edgecolor=color_obs, linewidth=0.0, zorder=8)
    # Western North America eddy-driven jet
    p7 = ax1.errorbar(lat_jet_mean_wna, 0., xerr=[lat_jet_err_wna], 
        fmt = 'o', color = color_jet, ecolor = color_jet, elinewidth=2, 
        capsize=0, zorder=10)
    # Add legend
    leg = ax1.legend([(p3, p4), (p1, p2), p7], 
        ['Observed $\mathregular{\pm}$ 2$\mathregular{\sigma}$',
        'Modeled $\mathregular{\pm}$ 2$\mathregular{\sigma}$',
        'Jet stream'], loc=6, ncol=3, bbox_to_anchor=(0.3, 1.3), 
        fontsize=16, frameon=False)
    # Change the marker size manually for observations
    leg.legendHandles[0]._legmarker.set_markersize(5)
    # Eastern North America
    ax2.text(0.98, 0.8, 'Eastern North America', ha='right', 
        transform=ax2.transAxes, fontsize=12, zorder=20)
    ax2.plot(lat_ena, field_ena, ls='-', lw=2, color=color_r, zorder=6)
    ax2.fill_between(lat_ena, field_ena-field_err_ena, field_ena+field_err_ena, 
        color=color_r, alpha=0.2, zorder=5)
    ax2.plot(lat_obs_ena, obs_mean_ena, '-', lw=2, color=color_obs, zorder=4)
    ax2.fill_between(lat_obs_ena, (obs_mean_ena-obs_std_ena),
        (obs_mean_ena+obs_std_ena), facecolor='None', hatch='////', 
        edgecolor=color_obs, linewidth=0.0, zorder=8)
    # Eastern North America eddy-driven jet
    ax2.errorbar(lat_jet_mean_ena, 0., xerr=[lat_jet_err_ena], fmt='o', 
    color=color_jet, ecolor=color_jet, elinewidth=2, capsize=0, zorder=10)
    # Europe
    ax3.text(0.98, 0.8, 'Europe', ha='right', transform=ax3.transAxes, 
        fontsize=12, zorder=20)
    ax3.plot(lat_eu, field_eu, ls='-', lw=2, color=color_r, zorder=6)
    ax3.fill_between(lat_eu, field_eu-field_err_eu, field_eu+field_err_eu, 
        color=color_r, alpha=0.2, zorder=5)
    ax3.plot(lat_obs_eu, obs_mean_eu, '-', lw=2, color=color_obs, zorder=4)
    ax3.fill_between(lat_obs_eu, (obs_mean_eu-obs_std_eu),
        (obs_mean_eu+obs_std_eu), facecolor='None', hatch='////', 
        edgecolor=color_obs, linewidth=0.0, zorder=8)
    # Europe eddy-driven jet
    ax3.errorbar(lat_jet_mean_eu, 0., xerr=[lat_jet_err_eu], fmt='o', 
        color=color_jet, ecolor=color_jet, elinewidth=2, capsize=0, zorder=10)
    # East Asia
    ax4.text(0.98, 0.8, 'China', ha='right', transform=ax4.transAxes, 
        fontsize=12, zorder=20)
    ax4.plot(lat_asia, field_asia, ls='-', lw=2, color=color_r, zorder=6)
    ax4.fill_between(lat_asia, field_asia-field_err_asia,
        field_asia+field_err_asia, color=color_r, alpha=0.2, zorder=5)
    ax4.plot(lat_obs_china, obs_mean_china, '-', lw=2, color=color_obs, 
        zorder=4)
    ax4.fill_between(lat_obs_china, (obs_mean_china-obs_std_china),
        (obs_mean_china+obs_std_china), facecolor='None', hatch='////', 
        edgecolor=color_obs, linewidth=0.0, zorder=8)
    # East Asia eddy-driven jet
    ax4.errorbar(lat_jet_mean_asia, 0., xerr=[lat_jet_err_asia], fmt='o', 
        color=color_jet, ecolor=color_jet, elinewidth=2, capsize=0, zorder=10)
    # For O3-humidity relationship
    ax1 = plt.subplot2grid((4,4), (0,2), colspan=2)
    ax2 = plt.subplot2grid((4,4), (1,2), colspan=2)
    ax3 = plt.subplot2grid((4,4), (2,2), colspan=2)
    ax4 = plt.subplot2grid((4,4), (3,2), colspan=2)
    field_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(
        r_qv2mo3, lat_gmi, lng_gmi, 235., 260., 25., 70.) 
    wna_left = geo_idx(235., lng_ml)
    wna_right = geo_idx(260., lng_ml)
    lat_jet_wna = lat_jet_ml[:, wna_left:wna_right+1]
    lat_jet_err_wna = np.nanstd(np.nanmean(lat_jet_wna, axis = 1))
    lat_jet_mean_wna = np.nanmean(lat_jet_wna)
    land_wna = globalo3_calculate.find_grid_overland(lat_wna, lng_wna)
    field_err_wna = 2*np.nanstd(field_wna*land_wna, axis=1)
    field_wna = np.nanmean(field_wna*land_wna, axis=1)
    obs_wna = np.where((np.array(lat_naps+lat_aqs) < 70.) &
                     (np.array(lat_naps+lat_aqs) > 25.) &
                     (np.array(lng_naps+lng_aqs) > 235.) &
                     (np.array(lng_naps+lng_aqs) < 260.))[0]
    lat_obs_wna = np.array(lat_naps+lat_aqs)[obs_wna]
    obs_wna = np.array(r_qv2mo3_naps+r_qv2mo3_aqs)[obs_wna]
    lat_obs_wna, obs_mean_wna, obs_std_wna = \
        globalo3_calculate.bin_observations_bylat(lat_wna[::2], obs_wna, 
        lat_obs_wna)
    field_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
        r_qv2mo3, lat_gmi, lng_gmi, 260., 295., 25., 70.)
    ena_left = geo_idx(260., lng_ml)
    ena_right = geo_idx(295., lng_ml)
    lat_jet_ena = lat_jet_ml[:, ena_left:ena_right+1]
    lat_jet_err_ena = np.nanstd(np.nanmean(lat_jet_ena, axis = 1))
    lat_jet_mean_ena = np.nanmean(lat_jet_ena)
    land_ena = globalo3_calculate.find_grid_overland(lat_ena, lng_ena)
    field_err_ena = 2*np.nanstd(field_ena*land_ena, axis=1)
    field_ena = np.nanmean(field_ena*land_ena, axis=1)
    obs_ena = np.where((np.array(lat_naps+lat_aqs) < 70.) &
                     (np.array(lat_naps+lat_aqs) > 25.) &
                     (np.array(lng_naps+lng_aqs) > 260.) &
                     (np.array(lng_naps+lng_aqs) < 295.))[0]
    lat_obs_ena = np.array(lat_naps+lat_aqs)[obs_ena]
    obs_ena = np.array(r_qv2mo3_naps+r_qv2mo3_aqs)[obs_ena]
    lat_obs_ena, obs_mean_ena, obs_std_ena = \
        globalo3_calculate.bin_observations_bylat(lat_ena[::2], obs_ena, 
        lat_obs_ena)
    field_eu1, lat_eu1, lng_eu1 = globalo3_calculate.find_grid_in_bb(
        r_qv2mo3, lat_gmi, lng_gmi, 350., 360., 25., 70.) 
    eu1_left = geo_idx(350., lng_ml)
    eu1_right = geo_idx(360., lng_ml)
    lat_jet_eu1 = lat_jet_ml[:, eu1_left:eu1_right+1]
    land_eu1 = globalo3_calculate.find_grid_overland(lat_eu1, lng_eu1)
    field_eu2, lat_eu2, lng_eu2 = globalo3_calculate.find_grid_in_bb(
        r_qv2mo3, lat_gmi, lng_gmi, 0., 30., 25., 70.) 
    eu2_left = geo_idx(0., lng_ml)
    eu2_right = geo_idx(30., lng_ml)
    lat_jet_eu2 = lat_jet_ml[:, eu2_left:eu2_right+1]
    land_eu2 = globalo3_calculate.find_grid_overland(lat_eu2, lng_eu2)
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
    obs_eu = np.array(r_qv2mo3_emep)[emep_eu]
    lat_obs_eu = np.array(lat_emep)[emep_eu]
    lat_obs_eu, obs_mean_eu, obs_std_eu = \
        globalo3_calculate.bin_observations_bylat(lat_eu[::2], obs_eu, 
        lat_obs_eu)
    field_asia, lat_asia, lng_asia = globalo3_calculate.find_grid_in_bb(
        r_qv2mo3, lat_gmi, lng_gmi, 90., 125., 25., 70.) 
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
    obs_asia = np.array(r_qv2mo3_china)[obs_asia]
    lat_obs_china, obs_mean_china, obs_std_china = \
        globalo3_calculate.bin_observations_bylat(lat_asia[::2], obs_asia, 
    lat_obs_asia)
    for ax in [ax1, ax2, ax3, ax4]: 
        ax.set_xlim([25, 70])
        ax.set_xticks(np.linspace(25, 70, 10))
        ax.set_xticklabels([])
        ax.set_ylim([np.linspace(-1., 1, 5)[0], np.linspace(-1., 1, 5)[-1]])
        ax.set_yticks(np.linspace(-1., 1, 5))
        ax.yaxis.tick_right()
        ax.axhline(y = 0.0, color = 'k', lw = 1., linestyle = '--')        
    ax4.set_xticklabels([ '25', '', '35', '', '45', '', '55', '', '65', ''])    
    ax4.set_xlabel('Latitude [$^{\mathregular{\circ}}$N]', fontsize=16)
    ax2.set_ylabel('r(q, O$_{\mathregular{3}}$) [$\cdot$]', fontsize=16, 
        rotation=270)
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.set_label_coords(1.13,-0.15)
    ax1.text(0.98, 0.8, 'Western North America', ha = 'right', 
        transform=ax1.transAxes, fontsize=12, zorder=20)
    p1, = ax1.plot(lat_wna, field_wna, ls='-', lw=2, color=color_r, 
        zorder=6)
    p2 = ax1.fill_between(lat_wna, field_wna-field_err_wna,
        field_wna+field_err_wna, color = color_r, alpha=0.2, zorder=5)
    p3, = ax1.plot(lat_obs_wna, obs_mean_wna, '-', lw=2, color=color_obs,     
        zorder=4)
    p4 = ax1.fill_between(lat_obs_wna, (obs_mean_wna-obs_std_wna),
        (obs_mean_wna+obs_std_wna), facecolor='None', hatch='////', 
        edgecolor=color_obs, linewidth=0.0, zorder=8)
    p7 = ax1.errorbar(lat_jet_mean_wna, 0., xerr=[lat_jet_err_wna], 
        fmt='o', color=color_jet, ecolor=color_jet, elinewidth=2, 
        capsize=0, zorder=10)
    leg.legendHandles[0]._legmarker.set_markersize(5)
    ax2.text(0.98, 0.8, 'Eastern North America', ha='right', 
        transform=ax2.transAxes, fontsize=12, zorder=20)
    ax2.plot(lat_ena, field_ena, ls='-', lw=2, color=color_r, zorder=6)
    ax2.fill_between(lat_ena, field_ena-field_err_ena, field_ena+field_err_ena, 
        color=color_r, alpha=0.2, zorder=5)
    ax2.plot(lat_obs_ena, obs_mean_ena, '-', lw=2, color=color_obs, zorder=4)
    ax2.fill_between(lat_obs_ena, (obs_mean_ena-obs_std_ena),
        (obs_mean_ena+obs_std_ena), facecolor='None', hatch='////', 
        edgecolor=color_obs, linewidth=0.0, zorder=8)
    ax2.errorbar(lat_jet_mean_ena, 0., xerr=[lat_jet_err_ena], fmt='o', 
        color=color_jet, ecolor=color_jet, elinewidth=2, capsize=0, zorder=10)
    ax3.text(0.98, 0.8, 'Europe', ha='right', transform=ax3.transAxes, 
        fontsize=12, zorder=20)
    ax3.plot(lat_eu, field_eu, ls='-', lw=2, color=color_r, zorder=6)
    ax3.fill_between(lat_eu, field_eu-field_err_eu, field_eu+field_err_eu, 
        color=color_r, alpha=0.2, zorder=5)
    ax3.plot(lat_obs_eu, obs_mean_eu, '-', lw=2, color=color_obs, zorder=4)
    ax3.fill_between(lat_obs_eu, (obs_mean_eu-obs_std_eu),
        (obs_mean_eu+obs_std_eu), facecolor='None', hatch='////', 
        edgecolor=color_obs, linewidth=0.0, zorder=8)
    ax3.errorbar(lat_jet_mean_eu, 0., xerr=[lat_jet_err_eu], fmt='o', 
        color=color_jet, ecolor=color_jet, elinewidth=2, capsize=0, zorder=10)
    ax4.text(0.98, 0.8, 'China', ha='right', transform=ax4.transAxes, 
        fontsize=12, zorder=20)
    ax4.plot(lat_asia, field_asia, ls='-', lw=2, color=color_r, zorder=6)
    ax4.fill_between(lat_asia, field_asia-field_err_asia,
        field_asia+field_err_asia, color=color_r, alpha=0.2, zorder=5)
    ax4.plot(lat_obs_china, obs_mean_china, '-', lw=2, color=color_obs, 
        zorder=4)
    ax4.fill_between(lat_obs_china, (obs_mean_china-obs_std_china),
        (obs_mean_china+obs_std_china), facecolor='None', hatch='////', 
        edgecolor=color_obs, linewidth=0.0, zorder=8)
    ax4.errorbar(lat_jet_mean_asia, 0., xerr=[lat_jet_err_asia], fmt='o', 
        color=color_jet, ecolor=color_jet, elinewidth=2, capsize=0, zorder=10)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
            'fig4.pdf', dpi=600)
    return

def fig5(lat_gmi, lng_gmi, r_t2mo3, r_t2mo3_transport, r_qv2mo3, 
    r_qv2mo3_transport, significance_diff_r_t2mo3, significance_diff_r_qv2mo3, 
    lat_jet_ml): 
    """Figure 5 of Kerr et al. (2020). (a) r(T, O3) from the transport only 
    simulation is shown with colored shading. Hatching indicates regions that 
    were significant in the full simulation that became insignificant in the 
    transport only simulation shading. (b) Same as (a) but showing r(q, O3). 
    Scatterpoints and error bars in (a-b) specify the mean position and 
    variability of the eddy-driven jet, respectively.

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    r_t2mo3_transport : numpy.ndarray     
        Pearson correlation coefficient between 2-meter temperature and O3 from 
        the transport only simulation, [lat, lng]
    r_qv2mo3_transport : numpy.ndarray     
        Pearson correlation coefficient between 2-meter specific humidity and 
        O3 from the transport only simulation, [lat, lng]
    significance_diff_r_t2mo3 : numpy.ndarray        
        The difference in the significance of r(T, O3) and r(T, O3) from the 
        transport only) simulation determined with moving block bootstrapping; 
        values of 1 correspond to grid cells that were significant in the full 
        simulation that became insignificant in the transport only simulation, 
        [lat, lng]
    significance_diff_r_qv2mo3 : numpy.ndarray
        The difference in the significance of r(q, O3) and r(q, O3) from the 
        transport only) simulation determined with moving block bootstrapping; 
        values of 1 correspond to grid cells that were significant in the full 
        simulation that became insignificant in the transport only simulation, 
        [lat, lng]
    lat_jet_ml : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in the Northern Hemisphere mid-latitudes, units of degrees north, 
        [time, lng]

    Returns
    -------
    None    
    """
    import numpy as np
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth']=0.3     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,5))
    # For wrapping around the Prime Meridian
    if lng_gmi[-1] != 360:
        lng_gmi[-1] = 360.
    # r(T, O3)
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax1.set_title(r'(a) $\mathregular{\delta}$ r(T, O$_{\mathregular{3}})$', 
        fontsize=16, x=0.02, ha='left')    
    # ax1.add_feature(ocean50m, zorder=3)
    ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax1.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax1.xaxis.set_major_formatter(lng_formatter)         
    ax1.get_xaxis().set_ticklabels([])
    ax1.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()    
    ax1.yaxis.set_major_formatter(lat_formatter)    
    cmap = plt.get_cmap('coolwarm')
    mb = ax1.contourf(lng_gmi, lat_gmi, (r_t2mo3-r_t2mo3_transport), 
        np.linspace(-0.4, 0.4, 9), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(T, O3)
    ax1.contourf(lng_gmi, lat_gmi, significance_diff_r_t2mo3, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    # Eddy-driven jet
    skiplng = 6
    ax1.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
        yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], color='k', markersize=3, 
        elinewidth=1.25, ecolor='k', fmt='o', transform=ccrs.PlateCarree(), 
        zorder=5)
    ax1.outline_patch.set_zorder(20)
    # r(q, O3)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2.set_title(r'(b) $\mathregular{\delta}$ r(q, O$_{\mathregular{3}}$)', 
        fontsize=16, x=0.02, ha='left')
    # ax2.add_feature(ocean50m)
    ax2.coastlines(lw=0.25, resolution='50m', color='k', zorder=4)
    ax2.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax2.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax2.xaxis.set_major_formatter(lng_formatter)          
    ax2.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax2.yaxis.set_major_formatter(lat_formatter)    
    mb = ax2.contourf(lng_gmi, lat_gmi, r_qv2mo3-r_qv2mo3_transport, 
        np.linspace(-0.4,0.4,9), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmi, lat_gmi, significance_diff_r_qv2mo3, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree())
    skiplng = 6
    ax2.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
        yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], zorder=10, color='k', 
        markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
        transform=ccrs.PlateCarree())
    ax2.outline_patch.set_zorder(20)      
    # Add colorbar
    plt.gcf().subplots_adjust(left=0.05, right=0.86, hspace=0.3)
    colorbar_axes = plt.gcf().add_axes([
        ax1.get_position().x1+0.03, # Left
        (ax2.get_position().y1-ax2.get_position().y0)/2.+ax2.get_position().y0, # Bottom 
        0.02, # Width
        ((ax1.get_position().y1-ax1.get_position().y0)/2.+ax1.get_position().y0)-
        ((ax2.get_position().y1-ax2.get_position().y0)/2.+ax2.get_position().y0)])    
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-0.4, 0.4, 9), extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[$\mathregular{\cdot}$]', fontsize=16, labelpad=11)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'fig5.pdf', dpi=600)
    return
    
def fig6(lat_gmi, lng_gmi, o3_gmi, t2m_merra, qv2m_merra, lat_jet_ml, 
    times_gmi, significance_r_o3jetdist, significance_r_t2mjetdist, 
    significance_r_qv2mjetdist):
    """Figure 7 of Kerr et al. (2020). (a) Difference in O3 on days with a
    poleward jet versus days with an equatorward jet. Hatched grid cells 
    correspond to insignificant values of r(O3, jet latitude) determined with 
    moving block bootstrapping. (b) Same as (a) but for T. (c) Same as (a) 
    but for q. Scatterpoints and error bars in (a-c) specify the mean 
    position and variability of the eddy-driven jet, respectively, for days 
    with a poleward jet and equatorward jet. 

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    o3_gmi : numpy.ndarray     
        Daily afternoon surface-level O3 from the GMI CTM, units of ppbv, 
        [time, lat, lng]
    t2m_merra : numpy.ndarray 
        Daily maximum 2-meter temperatures interpolated to the resolution of 
        the CTM, [time, lat, lng]
    qv2m_merra : numpy.ndarray
        Daily mean 3-hourly specific humidity units of g kg-1, [time, lat, lng]
    lat_jet_ml : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in the Northern Hemisphere mid-latitudes, units of degrees north, 
        [time, lng]
    times_gmi : numpy.ndarray
        datetime.date objects corresponding to every day in measuring period, 
        [time,] 
    significance_r_o3jetdist : numpy.ndarray
        Significance of r(O3, jet latitude) determined with moving block 
        bootstrapping: 1 implies significance, NaN implies insignificance, 
        [lat, lng]
    significance_r_t2mjetdist : numpy.ndarray
        Significance of r(T, jet latitude) determined with moving block 
        bootstrapping: 1 implies significance, NaN implies insignificance, 
        [lat, lng]
    significance_r_qv2mjetdist : numpy.ndarray
        Significance of r(q, jet latitude) determined with moving block 
        bootstrapping: 1 implies significance, NaN implies insignificance, 
        [lat, lng]        

    Returns
    -------
    None    
    """
    import numpy as np
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth']=0.3     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter    
    # Composites of flow at surface on days with poleward or equatorward jet 
    eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_o3, eqjet_o3 = \
        globalo3_calculate.segregate_field_bylat(o3_gmi, lng_gmi, lat_jet_ml, 
        times_gmi)
    eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_t2m, eqjet_t2m = \
        globalo3_calculate.segregate_field_bylat(t2m_merra, lng_gmi, lat_jet_ml, 
        times_gmi)
    eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_qv2m, eqjet_qv2m = \
        globalo3_calculate.segregate_field_bylat(qv2m_merra, lng_gmi, lat_jet_ml, 
        times_gmi)
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,7.7))
    if lng_gmi[-1] != 360:
        lng_gmi[-1] = 360.    
    # O3(PW - EW) 
    ax1 = plt.subplot2grid((3,2), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax1.set_title('(a) O$_\mathregular{3,\:PW}$ $-$ O$_\mathregular{3,\:EW}$', 
        fontsize=16, x=0.02, ha='left')
    # ax1.add_feature(ocean50m, zorder=3)
    ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax1.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax1.xaxis.set_major_formatter(lng_formatter)         
    ax1.get_xaxis().set_ticklabels([])
    ax1.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()    
    ax1.yaxis.set_major_formatter(lat_formatter)    
    cmap = plt.get_cmap('coolwarm')
    mb = ax1.contourf(lng_gmi, lat_gmi, (pwjet_o3-eqjet_o3), 
        np.linspace(-8,8,9), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(O3, jetlat)
    ax1.contourf(lng_gmi, lat_gmi, significance_r_o3jetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    # Add colorbar
    plt.gcf().subplots_adjust(left=0.05, right=0.86, hspace=0.3)
    colorbar_axes = plt.gcf().add_axes([ax1.get_position().x1+0.03, 
        ax1.get_position().y0, 0.02, (ax1.get_position().y1-
        ax1.get_position().y0)]) 
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-8,8,9), extend='neither')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[ppbv]', fontsize=16, labelpad=15)
    ax1.outline_patch.set_zorder(20)
    # T(PW - EW) 
    ax2 = plt.subplot2grid((3,2), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2.set_title('(b) T$_\mathregular{PW}$ $-$ T$_\mathregular{EW}$', 
        fontsize=16, x=0.02, ha='left')
    # ax2.add_feature(ocean50m, zorder=3)
    ax2.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax2.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax2.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax2.xaxis.set_major_formatter(lng_formatter)         
    ax2.get_xaxis().set_ticklabels([])
    ax2.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()    
    ax2.yaxis.set_major_formatter(lat_formatter)    
    mb = ax2.contourf(lng_gmi, lat_gmi, (pwjet_t2m-eqjet_t2m), 
        np.linspace(-8,8,9), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmi, lat_gmi, significance_r_t2mjetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    colorbar_axes = plt.gcf().add_axes([ax2.get_position().x1+0.03, 
        ax2.get_position().y0, 0.02, (ax2.get_position().y1-
        ax2.get_position().y0)]) 
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-8,8,9), extend='neither')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[K]', fontsize=16, labelpad=15)
    ax2.outline_patch.set_zorder(20)
    # q(PW - EW) 
    ax3 = plt.subplot2grid((3,2), (2,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax3.set_title('(c) q$_\mathregular{PW}$ $-$ q$_\mathregular{EW}$', 
        fontsize=16, x=0.02, ha='left')
    # ax3.add_feature(ocean50m, zorder=3)
    ax3.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax3.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax3.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax3.xaxis.set_major_formatter(lng_formatter)          
    ax3.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()    
    ax3.yaxis.set_major_formatter(lat_formatter)    
    mb = ax3.contourf(lng_gmi, lat_gmi, (pwjet_qv2m-eqjet_qv2m), 
        np.linspace(-3,3,7), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax3.contourf(lng_gmi, lat_gmi, significance_r_qv2mjetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    colorbar_axes = plt.gcf().add_axes([ax3.get_position().x1+0.03, 
        ax3.get_position().y0, 0.02, (ax3.get_position().y1-
        ax3.get_position().y0)]) 
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-3,3,7), extend='neither')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[g kg$^{\mathregular{-1}}$]', fontsize=16, labelpad=13)
    ax3.outline_patch.set_zorder(20)
    # Mean position of the eddy-driven jet
    for ax in [ax1, ax2, ax3]:
        skiplng = 6
        ax.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml, axis=0)[::skiplng],
            yerr=np.std(lat_jet_ml,axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())       
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'fig6.pdf', dpi=600)
    return

def fig7(lat_cyclones_binned, lng_cyclones_binned, cyclones_binned,
    pwjet_cyclones_binned, eqjet_cyclones_binned, lat_gmi, lng_gmi, 
    lat_jet_ml): 
    """Figure 7 of Kerr et al. (2020). (a) Total number of cyclones detected 
    by MCMS on sub-daily (six-hourly) time scales binned to a ~4˚ x ~4˚ grid. 
    As in Figure 1, the mean latitude of the eddy-driven jet and its 
    variability are shown for reference. (b) The difference in the total 
    number of cyclones calculated between days when the jet is in a poleward 
    (PW) and equatorward (EW) position. A poleward (an equatorward) jet is 
    defined at each longitudinal band when the latitude of the jet exceeds 
    (is less than) the 70th (30th) percentile. 
    
    Parameters
    ----------
    lat_cyclones_binned : numpy.ndarray
        Binned latitude array corresponding to cyclones, units of degrees 
        north, [lat_bin,]
    lng_cyclones_binned : numpy.ndarray
        Binned longitude array corresponding to cyclones, units of degrees
        east, [lng_bin,]   
    cyclones_binned : numpy.ndarray
        Binned cyclone frequency, [lat_bin, lng_bin]
    pwjet_cyclones_binned : numpy.ndarray
        Binned cyclone frequency on the subset of days where the jet is in 
        a poleward position, [lat_bin, lng_bin]
    eqjet_cyclones_binned : numpy.ndarray
        Binned cyclone frequency on the subset of days where the jet is in 
        an equatorward position, [lat_bin, lng_bin]    
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    lat_jet_ml : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in the Northern Hemisphere mid-latitudes, units of degrees north, 
        [time, lng]

    Returns
    -------
    None
    """
    import numpy as np
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth']=0.3     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter        
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,5))
    if lng_gmi[-1] != 360:
        lng_gmi[-1] = 360.    
    # Total cyclone frequency
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax1.set_title(r'(a) Cyclone frequency', fontsize=16, x=0.02, ha='left')    
    ax1.add_feature(ocean50m, zorder=3)
    ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax1.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax1.xaxis.set_major_formatter(lng_formatter)         
    ax1.get_xaxis().set_ticklabels([])
    ax1.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()    
    ax1.yaxis.set_major_formatter(lat_formatter)    
    cmap = plt.get_cmap('OrRd')
    clevs = np.linspace(0, 35, 8)
    # Mask small frequencies of cyclones 
    cyclones_binned_mask = np.ma.masked_array(cyclones_binned, 
        ((cyclones_binned >= 0) & (cyclones_binned <= 5)))
    mb = ax1.pcolor(lng_cyclones_binned, lat_cyclones_binned, 
        cyclones_binned_mask, cmap=plt.get_cmap(cmap, len(clevs)-1),
        vmin=clevs[0], vmax=clevs[-1], norm = mpl.colors.BoundaryNorm(clevs, 
        ncolors=cmap.N, clip=False), transform=ccrs.PlateCarree(), 
        zorder=5)
    # Eddy-driven jet
    skiplng = 6
    ax1.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
        yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], color='k', markersize=3, 
        elinewidth=1.25, ecolor='k', fmt='o', transform=ccrs.PlateCarree(), 
        zorder=6)
    # Colorbar
    plt.gcf().subplots_adjust(left=0.02, right=0.86, hspace=0.3)    
    colorbar_axes = plt.gcf().add_axes([ax1.get_position().x1+0.03, 
        ax1.get_position().y0, 0.02, (ax1.get_position().y1-
        ax1.get_position().y0)])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(0, 35, 8), extend='max')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[$\mathregular{\cdot}$]', fontsize=16, labelpad=9)
    ax1.outline_patch.set_zorder(20)
    # (PW - EW) cyclone difference
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2.set_title('(b) Cyclones$_{\:\mathregular{PW}}$ - '+\
        'Cyclones$_{\: \mathregular{EW}}$', fontsize=16, x=0.02, ha='left')
    ax2.add_feature(ocean50m)
    ax2.coastlines(lw=0.25, resolution='50m', color='k', zorder=4)
    ax2.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax2.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax2.xaxis.set_major_formatter(lng_formatter)      
    ax2.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()    
    ax2.yaxis.set_major_formatter(lat_formatter)
    cmap = plt.get_cmap('coolwarm')
    clevs = np.linspace(-12, 12, 13)
    diff = pwjet_cyclones_binned-eqjet_cyclones_binned
    diff_mask = np.ma.masked_array(diff, ((diff >= -2) & (diff <= 2)))
    mb = ax2.pcolor(lng_cyclones_binned, lat_cyclones_binned, diff_mask, 
        cmap=plt.get_cmap(cmap, len(clevs)-1), vmin=clevs[0], 
        vmax=clevs[-1], norm = mpl.colors.BoundaryNorm(clevs, 
        ncolors=cmap.N, clip=False), transform=ccrs.PlateCarree(), 
        zorder=5)
    skiplng = 6
    ax2.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
        yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], zorder=10, color='k', 
        markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
        transform=ccrs.PlateCarree())
    ax2.outline_patch.set_zorder(20)      
    # Add colorbar
    colorbar_axes = plt.gcf().add_axes([ax2.get_position().x1+0.03, 
        ax2.get_position().y0, 0.02, (ax2.get_position().y1-
        ax2.get_position().y0)])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-12, 12, 7), extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[$\mathregular{\cdot}$]', fontsize=16)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'fig7.pdf', dpi=600)
    return 
    
def fig8(o3_anom_rotated):
    """Plot mean O3 anomaly and the standard deviation in the vicinity (+/- 5 
    grid cells) of cyclones for the (left) unrotated case and (right) rotated
    case, where the east side of the subplot is the direction of cyclone 
    motion. Note that finding the mean of the rotated O3 anomaly raises a 
    RuntimeWarning because the rotation process fills the outer edge with NaNs
    when applicable. 

    Parameters
    ----------
    o3_anom_rotated : list 
        The daily O3 anomaly in the vicinity of cyclone rotated such that 
        the upward direction is the direction of cyclone motion (n.b., the mean 
        rotated anomaly is further rotated 90˚ in the subplot such that the 
        direction of cyclone motion is to the right), [cyclones,]    
    
    Returns
    -------
    None
    """
    import numpy as np
    import matplotlib.pyplot as plt
    # Plotting unrotated and rotated cases
    fig = plt.figure(figsize=(5.5,4.2))
    ax = plt.subplot2grid((1,1),(0,0))
    clevs = np.linspace(-1., 1., 11)
    ax.set_aspect('equal')
    mb = ax.contourf(np.rot90(np.nanmean(np.stack(o3_anom_rotated),axis=0)),
        clevs, cmap=plt.get_cmap('coolwarm'), extend='both', zorder=1)
    CS = ax.contour(np.nanstd(np.stack(o3_anom_rotated), axis=0), 
        [6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5],
        linewidths=1.5, colors='k', zorder=3)                     
    plt.clabel(CS, fontsize=10, fmt='%.2f', inline=True)
    # Axis labels
    ax.set_xlabel('Tangential', fontsize=16)
    ax.set_ylabel('Orthogonal', fontsize=16)
    # Since rotation handles interpolation around the periphery in a non-
    # ideal way, strip off the edge grid cells in unrotated and rotated 
    # cases
    # Aesthetics (crosshairs for system's center, axes labels)
    ax.axhline(y=6, xmin=ax.get_xlim()[0], xmax=ax.get_xlim()[1], ls='--', 
        lw=0.75, color='k', zorder=2)
    ax.axvline(x=6, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], ls='--', 
        lw=0.75, color='k', zorder=2)        
    ax.set_xlim([1, o3_anom_rotated[0].shape[0]-2])
    ax.set_ylim([1, o3_anom_rotated[0].shape[0]-2])    
    ax.set_xticks(np.arange(1, o3_anom_rotated[0].shape[0]-1, 1))
    ax.set_xticklabels(['', '-4', '', '-2', '', '0', '', '2', '', '4', ''], 
        fontsize=12)
    ax.set_yticks(np.arange(1, o3_anom_rotated[0].shape[0]-1, 1))
    ax.set_yticklabels(['', '-4', '', '-2', '', '0', '', '2', '', '4', ''],
        fontsize=12)
    plt.gcf().subplots_adjust(left=0.03, right=0.86, wspace=0.35)  
    # Add colorbar
    colorbar_axes = plt.gcf().add_axes([0.8, ax.get_position().y0,
        0.04, (ax.get_position().y1-ax.get_position().y0)]) 
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-1.,1.,11), extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('$\mathregular{\delta}\:$O$_{\mathregular{3}}$ [ppbv]', 
        fontsize=16, labelpad=15)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'fig8.pdf', dpi=600)
    return 
  
def fig9(lat_gmi, lng_gmi, pblh_merra, U10M, V10M, lat_jet_ml, times_gmi, 
    significance_r_pblhjetdist, significance_r_U10Mjetdist, 
    significance_r_V10Mjetdist):
    """Figure 9 of Kerr et al. (2020). (a) Difference in PBLH on days with a
    poleward jet versus days with an equatorward jet. Hatched grid cells 
    correspond to insignificant values of r(PBLH, jet latitude) determined with 
    moving block bootstrapping. (b) Same as (a) but for U10. (c) Same as (a) 
    but for V10. Scatterpoints and error bars in (a-c) specify the mean 
    position and variability of the eddy-driven jet, respectively, for days 
    with a poleward jet and equatorward jet. 

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    pblh_merra : numpy.ndarray 
        Daily mean PBL height, units of m, [time, lat, lng]
    U10M : numpy.ndarray 
        Daily mean 10-meter zonal wind, units of m s-1 (U10), [time, lat, lng]
    V10M : numpy.ndarray 
        Daily mean 10-meter meridional wind (V10), units of m s-1, [time, lat, 
        lng]
    lat_jet_ml : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in the Northern Hemisphere mid-latitudes, units of degrees north, 
        [time, lng]
    times_gmi : numpy.ndarray
        datetime.date objects corresponding to every day in measuring period, 
        [time,] 
    significance_r_pblhjetdist : numpy.ndarray
        Significance of r(PBLH, jet latitude) determined with moving block 
        bootstrapping: 1 implies significance, NaN implies insignificance, 
        [lat, lng]
    significance_r_U10Mjetdist : numpy.ndarray
        Significance of r(U10, jet latitude) determined with moving block 
        bootstrapping: 1 implies significance, NaN implies insignificance, 
        [lat, lng]
    significance_r_V10Mjetdist : numpy.ndarray
        Significance of r(V10, jet latitude) determined with moving block 
        bootstrapping: 1 implies significance, NaN implies insignificance, 
        [lat, lng]        

    Returns
    -------
    None    
    """
    import numpy as np
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth']=0.3     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter    
    # Composites of flow at surface on days with poleward or equatorward jet 
    eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_pblh, eqjet_pblh = \
        globalo3_calculate.segregate_field_bylat(pblh_merra, lng_gmi, lat_jet_ml, 
        times_gmi)
    eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_U10M, eqjet_U10M = \
        globalo3_calculate.segregate_field_bylat(U10M, lng_gmi, lat_jet_ml, 
        times_gmi)
    eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_V10M, eqjet_V10M = \
        globalo3_calculate.segregate_field_bylat(V10M, lng_gmi, lat_jet_ml, 
        times_gmi)
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,7.7))
    if lng_gmi[-1] != 360:
        lng_gmi[-1] = 360.    
    # PBLH(PW - EW) 
    ax1 = plt.subplot2grid((3,2), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax1.set_title('(a) PBLH$_\mathregular{PW}$ $-$ PBLH$_\mathregular{EW}$', 
        fontsize=16, x=0.02, ha='left')
    # ax1.add_feature(ocean50m, zorder=3)
    ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax1.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax1.xaxis.set_major_formatter(lng_formatter)         
    ax1.get_xaxis().set_ticklabels([])    
    ax1.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()    
    ax1.yaxis.set_major_formatter(lat_formatter)
    cmap = plt.get_cmap('coolwarm')
    mb = ax1.contourf(lng_gmi, lat_gmi, (pwjet_pblh-eqjet_pblh), 
        np.linspace(-300,300,7), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(PBLH, jetlat)
    ax1.contourf(lng_gmi, lat_gmi, significance_r_pblhjetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    # Add colorbar
    plt.gcf().subplots_adjust(left=0.05, right=0.86, hspace=0.3)
    colorbar_axes = plt.gcf().add_axes([ax1.get_position().x1+0.03, 
        ax1.get_position().y0, 0.02, (ax1.get_position().y1-
        ax1.get_position().y0)]) 
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-300,300,7), extend='neither')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[m]', fontsize=16)
    ax1.outline_patch.set_zorder(20)
    # U10M-jet position relationship 
    ax2 = plt.subplot2grid((3,2), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2.set_title('(b) U$_\mathregular{10,\:PW}$ $-$ U$_\mathregular{10,\:EW}$', 
        fontsize=16, x=0.02, ha='left')
    # ax2.add_feature(ocean50m, zorder=3)
    ax2.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax2.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax2.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax2.xaxis.set_major_formatter(lng_formatter)         
    ax2.get_xaxis().set_ticklabels([])    
    ax2.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()    
    ax2.yaxis.set_major_formatter(lat_formatter)    
    mb = ax2.contourf(lng_gmi, lat_gmi, (pwjet_U10M-eqjet_U10M), 
        np.linspace(-4, 4, 9), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmi, lat_gmi, significance_r_U10Mjetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    colorbar_axes = plt.gcf().add_axes([ax2.get_position().x1+0.03, 
        ax2.get_position().y0, 0.02, (ax2.get_position().y1-
        ax2.get_position().y0)]) 
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-4,4,9), extend='neither')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[m s$^{\mathregular{-1}}$]', fontsize=16, labelpad=15)    
    ax2.outline_patch.set_zorder(20)
    # V10M-jet position relationship
    ax3 = plt.subplot2grid((3,2), (2,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax3.set_title('(c) V$_\mathregular{10,\:PW}$ $-$ V$_\mathregular{10,\:EW}$', 
        fontsize=16, x=0.02, ha='left')
    # ax3.add_feature(ocean50m, zorder=3)
    ax3.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax3.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    ax3.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax3.xaxis.set_major_formatter(lng_formatter)        
    ax3.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()    
    ax3.yaxis.set_major_formatter(lat_formatter)    
    mb = ax3.contourf(lng_gmi, lat_gmi, (pwjet_V10M-eqjet_V10M), 
        np.linspace(-3,3,7), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax3.contourf(lng_gmi, lat_gmi, significance_r_V10Mjetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    colorbar_axes = plt.gcf().add_axes([ax3.get_position().x1+0.03, 
        ax3.get_position().y0, 0.02, (ax3.get_position().y1-
        ax3.get_position().y0)]) 
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-3,3,7), extend='neither')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[m s$^{\mathregular{-1}}$]', fontsize=16, labelpad=15)
    ax3.outline_patch.set_zorder(20)    
    # Mean position of the eddy-driven jet
    for ax in [ax1, ax2, ax3]:
        skiplng = 6
        ax.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml, axis=0)[::skiplng],
            yerr=np.std(lat_jet_ml,axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'fig9.pdf', dpi=600)  

def fig10(lat_gmi, V10M, o3_gmi, t2m_merra, qv2m_merra): 
    """Figure 10 of Kerr et al. (2020). The total flux and the eddy and mean
    meridional contributions to the total flux of (a-c) O3, (d-f) temperature, 
    (g-i) and humidity. Calculations of the total flux and its components are 
    done for all summer days (first column; a, d, g), days when the jet is in 
    a PW position (second column; b, e, h), and days when the jet is in a EW 
    position (third column; c, f, i).
    
    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    V10M : numpy.ndarray 
        Daily mean 10-meter meridional wind (V10), units of m s-1, [time, lat, 
        lng]
    o3_gmi : numpy.ndarray     
        Daily afternoon surface-level O3 from the GMI CTM, units of ppbv, 
        [time, lat, lng]        
    t2m_merra : numpy.ndarray 
        Daily maximum 2-meter temperatures interpolated to the resolution of 
        the CTM, [time, lat, lng]
    qv2m_merra : numpy.ndarray
        Daily mean 3-hourly specific humidity units of g kg-1, [time, lat, lng]

    Returns
    -------
    None    
    """
    import numpy as np
    # Separate 10-meter meridional wind and O3 by jet position
    V10M_eqjet, V10M_pwjet, o3_eqjet, o3_pwjet = \
        globalo3_calculate.sortfield_byjetlat_column(
        np.reshape(V10M, (276,1,92,288)), np.reshape(o3_gmi, (276,1,92,288)),
        lat_jet_ml, lng_gmi, lat_gmi, np.array([1000.]), psize=30)
    t2m_eqjet, t2m_pwjet, t2m_eqjet, t2m_pwjet = \
        globalo3_calculate.sortfield_byjetlat_column(
        np.reshape(t2m_merra, (276,1,92,288)), np.reshape(t2m_merra, 
        (276,1,92,288)),lat_jet_ml, lng_gmi, lat_gmi, np.array([1000.]), 
        psize=30)
    qv2m_eqjet, qv2m_pwjet, qv2m_eqjet, qv2m_pwjet = \
        globalo3_calculate.sortfield_byjetlat_column(
        np.reshape(qv2m_merra, (276,1,92,288)), np.reshape(qv2m_merra, 
        (276,1,92,288)),lat_jet_ml, lng_gmi, lat_gmi, np.array([1000.]), 
        psize=30)    
    # Surface-level O3, 2-meter temperature flux for all days, pole- and 
    # equatorward jet days
    o3_mean, o3_stationary, o3_transient, o3_total = \
        globalo3_calculate.meridional_flux(V10M, o3_gmi, mtime, lat_gmi, lng_gmi)
    o3_mean_pwjet, o3_stationary_pwjet, o3_transient_pwjet, o3_total_pwjet = \
        globalo3_calculate.meridional_flux(V10M_pwjet[:,0], o3_pwjet[:,0], 
        mtime[:len(V10M_pwjet)], lat_gmi, lng_gmi)
    o3_mean_eqjet, o3_stationary_eqjet, o3_transient_eqjet, o3_total_eqjet = \
        globalo3_calculate.meridional_flux(V10M_eqjet[:,0], o3_eqjet[:,0],     
        mtime[:len(V10M_eqjet)], lat_gmi, lng_gmi)
    t2m_mean, t2m_stationary, t2m_transient, t2m_total = \
        globalo3_calculate.meridional_flux(V10M, t2m_merra, mtime, lat_gmi, 
        lng_gmi)
    t2m_mean_pwjet, t2m_stationary_pwjet, t2m_transient_pwjet, t2m_total_pwjet = \
        globalo3_calculate.meridional_flux(V10M_pwjet[:,0], t2m_pwjet[:,0], 
        mtime[:len(V10M_pwjet)], lat_gmi, lng_gmi)
    t2m_mean_eqjet, t2m_stationary_eqjet, t2m_transient_eqjet, t2m_total_eqjet = \
        globalo3_calculate.meridional_flux(V10M_eqjet[:,0], t2m_eqjet[:,0],     
        mtime[:len(V10M_eqjet)], lat_gmi, lng_gmi)
    qv2m_mean, qv2m_stationary, qv2m_transient, qv2m_total = \
        globalo3_calculate.meridional_flux(V10M, qv2m_merra, mtime, lat_gmi, 
        lng_gmi)
    qv2m_mean_pwjet, qv2m_stationary_pwjet, qv2m_transient_pwjet, qv2m_total_pwjet = \
        globalo3_calculate.meridional_flux(V10M_pwjet[:,0], qv2m_pwjet[:,0], 
        mtime[:len(V10M_pwjet)], lat_gmi, lng_gmi)
    qv2m_mean_eqjet, qv2m_stationary_eqjet, qv2m_transient_eqjet, qv2m_total_eqjet = \
        globalo3_calculate.meridional_flux(V10M_eqjet[:,0], qv2m_eqjet[:,0],     
        mtime[:len(V10M_eqjet)], lat_gmi, lng_gmi)    
    import matplotlib.pyplot as plt
    # O3 meridional flux
    fig = plt.figure(figsize=(11, 12))
    ax1 = plt.subplot2grid((3,3),(0,0))
    ax2 = plt.subplot2grid((3,3),(0,1))
    ax3 = plt.subplot2grid((3,3),(0,2))
    for ax in [ax1, ax2, ax3]:
        ax.set_xlim([25,70])
        ax.set_xticks([25,40,55,70])
        ax.set_xticklabels([])
        ax.set_ylim([-30, 35])
        ax.set_yticks(np.linspace(-30, 35, 6))
        ax.set_yticklabels([-30, -17, -4, 9, 22, 35], fontsize=12)
        # Add horizontal line for field = 0
        ax.axhline(y=0.0, color='k', lw=1., linestyle='--', zorder=1)
    ax2.set_yticklabels([]); ax3.set_yticklabels([]) 
    # O3 meridional flux for all days
    ax1.plot(o3_total, lw=2, ls='-', color='#A2C3E0')
    ax1.plot(o3_mean, lw=2, ls='--', color='#EF9802')
    ax1.plot((o3_stationary+o3_transient), lw=2, ls='-', color='#3F79B7')
    # O3 meridional flux for PW jet days
    ax2.plot(o3_total_pwjet, lw=2, color='#A2C3E0')
    ax2.plot(o3_mean_pwjet, ls='--', lw=2, color='#EF9802')
    ax2.plot((o3_stationary_pwjet+o3_transient_pwjet), 
        lw=2, color='#3F79B7')
    # O3 meridional flux for EW jet days
    ax3.plot(o3_total_eqjet, lw=2, color='#A2C3E0')
    ax3.plot(o3_mean_eqjet, ls='--', lw=2, color='#EF9802')
    ax3.plot((o3_stationary_eqjet+o3_transient_eqjet), lw=2, 
        color='#3F79B7')
    # Set titles, axis labels
    ax1.text(0.04, 0.9, '(a)', ha='left', transform=ax1.transAxes, 
        fontsize=16, zorder=20)
    ax1.set_title('All', fontsize=16)#, x=0.02, ha='left')
    ax2.set_title('PW', fontsize=16)
    ax3.set_title('EW', fontsize=16)
    ax2.text(0.04, 0.9, '(b)', ha='left', transform=ax2.transAxes, 
        fontsize=18, zorder=20)
    ax3.text(0.04, 0.9, '(c)', ha='left', transform=ax3.transAxes, 
        fontsize=16, zorder=20) 
    ax1.set_ylabel('O$_{\mathregular{3}}$ flux [ppbv m s$^{\mathregular{-1}}$]', 
        fontsize=16)
    ax1.get_yaxis().set_label_coords(-0.20,0.5)
    # 2-meter temperature meridional flux    
    ax4 = plt.subplot2grid((3,3),(1,0))
    ax5 = plt.subplot2grid((3,3),(1,1))
    ax6 = plt.subplot2grid((3,3),(1,2))
    for ax in [ax4, ax5, ax6]:
        ax.set_xlim([25,70])
        ax.set_xticks([25,40,55,70])
        ax.set_xticklabels([])
        ax.set_ylim([-200, 350])
        ax.set_yticks(np.linspace(-200, 350, 6))
        ax.set_yticklabels([-200, -90, 20, 130, 240, 350], fontsize=12)
        ax.axhline(y=0.0, color='k', lw=1., linestyle='--', zorder=1)
    ax5.set_yticklabels([]); ax6.set_yticklabels([])
    # 2-meter meridional temperature flux for all days
    ax4.plot(t2m_total, lw=2, ls='-', color='#A2C3E0')
    ax4.plot(t2m_mean, lw=2, ls='--', color='#EF9802')
    ax4.plot((t2m_stationary+t2m_transient), lw=2, ls='-', color='#3F79B7')
    # 2-meter meridional temperature for PW jet days
    ax5.plot(t2m_total_pwjet, lw=2, color='#A2C3E0')
    ax5.plot(t2m_mean_pwjet, lw=2, ls='--', color='#EF9802')
    ax5.plot((t2m_stationary_pwjet+t2m_transient_pwjet), 
        lw=2, color='#3F79B7')
    # 2-meter meridional temperature for EW jet days
    ax6.plot(t2m_total_eqjet, lw=2, color='#A2C3E0')
    ax6.plot(t2m_mean_eqjet, lw=2, ls='--', color='#EF9802')
    ax6.plot((t2m_stationary_eqjet+t2m_transient_eqjet), lw=2, 
        color='#3F79B7')
    # Set titles, axis labels
    ax4.text(0.04, 0.9, '(d)', ha='left', transform=ax4.transAxes, 
        fontsize=16, zorder=20)
    ax5.text(0.04, 0.9, '(e)', ha='left', transform=ax5.transAxes, 
        fontsize=16, zorder=20)
    ax6.text(0.04, 0.9, '(f)', ha='left', transform=ax6.transAxes, 
        fontsize=16, zorder=20)
    ax4.set_ylabel('Temperature flux [K m s$^{\mathregular{-1}}$]', 
        fontsize=16)   
    ax4.get_yaxis().set_label_coords(-0.20,0.55)
    # 2-meter specific humidity meridional flux    
    ax7 = plt.subplot2grid((3,3),(2,0))
    ax8 = plt.subplot2grid((3,3),(2,1))
    ax9 = plt.subplot2grid((3,3),(2,2))
    for ax in [ax7, ax8, ax9]:
        ax.set_xlim([25,70])
        ax.set_xticks([25,40,55,70])
        ax.set_xticklabels([25, 40, 55, 70], fontsize=12)
        ax.set_ylim([-10, 10])
        ax.set_yticks(np.linspace(-10,10,6))
        ax.set_yticklabels([-10, -6, -2, 2, 6, 10], fontsize=12)
        ax.axhline(y=0.0, color='k', lw=1., linestyle='--', zorder=1)
    ax8.set_yticklabels([]); ax9.set_yticklabels([])
    # 2-meter meridional humidity flux for all days
    ax7.plot(qv2m_total, lw=2, ls='-', color='#A2C3E0')
    ax7.plot(qv2m_mean, lw=2, ls='--', color='#EF9802')
    ax7.plot((qv2m_stationary+qv2m_transient), lw=2, ls='-', 
        color='#3F79B7')
    # 2-meter meridional humidity for PW jet days
    ax8.plot(qv2m_total_pwjet, lw=2, color='#A2C3E0', label='Total')
    ax8.plot(qv2m_mean_pwjet, lw=2, ls='--', color='#EF9802', 
        label='Mean')
    ax8.plot((qv2m_stationary_pwjet+qv2m_transient_pwjet), 
        lw=2, color='#3F79B7', label='Eddy (Transient + Stationary)')
    # 2-meter meridional humidity for EW jet days
    ax9.plot(qv2m_total_eqjet, lw=2, color='#A2C3E0')
    ax9.plot(qv2m_mean_eqjet, lw=2, ls='--', color='#EF9802')
    ax9.plot((qv2m_stationary_eqjet+qv2m_transient_eqjet), lw=2, 
        color='#3F79B7')
    # Set titles, axis labels
    ax7.text(0.04, 0.9, '(g)', ha='left', transform=ax7.transAxes, 
        fontsize=16, zorder=20)
    ax8.text(0.04, 0.9, '(h)', ha='left', transform=ax8.transAxes, 
        fontsize=16, zorder=20)
    ax9.text(0.04, 0.9, '(i)', ha='left', transform=ax9.transAxes, 
        fontsize=16, zorder=20)
    ax8.set_xlabel('Latitude [$^{\circ}$N]', fontsize=16)
    ax7.set_ylabel('Humidity flux [g m kg$^{\mathregular{-1}}$ s$^{\mathregular{-1}}$]', 
        fontsize=16)   
    ax7.get_yaxis().set_label_coords(-0.20,0.45)
    # Add legend
    ax8.legend(loc=2, bbox_to_anchor=(-0.8,-0.22), ncol=3, fontsize=16, 
        frameon=False)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'fig10.pdf', dpi=600)
    return     

def figS1(lat_gmi, lng_gmi, do3dt2m, do3dq, significance_r_t2mo3, 
    significance_r_qv2mo3, lat_jet_ml): 
    """Figure S1 of Kerr et al. (2020). (a) The slope of the ordinary least 
    squares (OLS) regression of O3 versus temperature, dO3/dT. Hatching denotes 
    regions where the correlation between O3 and temperature is insignificant,
    determined using moving block bootstrap resampling to estimate the 95% 
    confidence interval. (b) Same as (a) but for O3 versus specific humidity,
    dO3/dq, with hatching showing insignificant correlation between O3 and 
    specific humidity. Scatter points and error bars are identical in (a-b)
    and show the mean latitude of the eddy-driven jet and its variability. 

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    do3dt2m : numpy.ndarray     
        The slope of the ordinary least squares linear regression of O3 versus
        2-meter temperature, units of ppbv K-1, [lat, lng]
    do3dqv2m : numpy.ndarray     
        The slope of the ordinary least squares linear regression of O3 versus
        2-meter specific humidity, units of ppbv kg g-1, [lat, lng]
    significance_r_t2mo3 : numpy.ndarray        
        Significance of r(T, O3) determined with moving block bootstrapping: 1
        implies significance, NaN implies insignificance, [lat, lng]
    significance_r_qv2mo3 : numpy.ndarray
        Significance of r(q, O3) determined with moving block bootstrapping: 1
        implies significance, NaN implies insignificance, [lat, lng]    
    lat_jet_ml : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in the Northern Hemisphere mid-latitudes, units of degrees north, 
        [time, lng]

    Returns
    -------
    None    
    """
    import numpy as np
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth']=0.3     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter    
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,5))
    if lng_gmi[-1] != 360:
        lng_gmi[-1] = 360.    
    # dO3/dT
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax1.set_title('(a) dO$_{\mathregular{3}}$/dT', fontsize=16, 
        x=0.02, ha='left')    
    # ax1.add_feature(ocean50m, zorder=3)
    ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    cmap = plt.get_cmap('coolwarm')
    mb = ax1.contourf(lng_gmi, lat_gmi, do3dt2m, np.linspace(-2, 2, 9), 
        cmap=cmap, extend='both', transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(T, O3)
    ax1.contourf(lng_gmi, lat_gmi, significance_r_t2mo3, hatches=['//////'], 
        colors='none', transform=ccrs.PlateCarree(), zorder=2)
    ax1.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax1.xaxis.set_major_formatter(lng_formatter)         
    ax1.get_xaxis().set_ticklabels([])    
    ax1.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax1.yaxis.set_major_formatter(lat_formatter)    
    # Eddy-driven jet
    skiplng = 6
    ax1.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
        yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], color='k', markersize=3, 
        elinewidth=1.25, ecolor='k', fmt='o', transform=ccrs.PlateCarree(), 
        zorder=5)
    plt.gcf().subplots_adjust(left=0.05, right=0.86, hspace=0.3)    
    colorbar_axes = plt.gcf().add_axes([ax1.get_position().x1+0.03, 
        ax1.get_position().y0, 0.02, (ax1.get_position().y1-
        ax1.get_position().y0)])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-2, 2, 9), extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[ppbv K$^{\mathregular{-1}}$]', fontsize=16)
    ax1.outline_patch.set_zorder(20)
    # dO3/dq
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2.set_title(r'(b) dO$_{\mathregular{3}}$/dq', fontsize=16, x=0.02, 
        ha='left')
    # ax2.add_feature(ocean50m)
    ax2.coastlines(lw=0.25, resolution='50m', color='k', zorder=4)
    ax2.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    mb = ax2.contourf(lng_gmi, lat_gmi, do3dq, np.linspace(-2, 2, 9), 
        cmap=cmap, extend='both', transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmi, lat_gmi, significance_r_qv2mo3, hatches=['//////'], 
        colors='none', transform=ccrs.PlateCarree())
    ax2.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax2.xaxis.set_major_formatter(lng_formatter)      
    ax2.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax2.yaxis.set_major_formatter(lat_formatter)    
    # skiplng = 6
    ax2.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
        yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], zorder=10, color='k', 
        markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
        transform=ccrs.PlateCarree())
    colorbar_axes = plt.gcf().add_axes([ax2.get_position().x1+0.03, 
        ax2.get_position().y0, 0.02, (ax2.get_position().y1-
        ax2.get_position().y0)])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-2, 2, 9), extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[ppbv kg g$^{\mathregular{-1}}$]', fontsize=16)
    ax2.outline_patch.set_zorder(20)      
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'figS1_TLS.pdf', dpi=600)
    return    

def figS2(lat_gmi, lng_gmi, r_o3jetdist, r_t2mjetdist, r_qv2mjetdist, 
    lat_jet_ml, significance_r_o3jetdist, significance_r_t2mjetdist, 
    significance_r_qv2mjetdist): 
    """Figure S2 of Kerr et al. (2020). (a) r(O3, jet distance) is shown with 
    colored shading. Hatching indicates insignificance determined with moving 
    block bootstrapping. (b) Same as (a) but showing r(T, jet distance). (c)
    Same as (a) but showing r(q, jet distance). Scatterpoints and error bars 
    in (a-c) specify the mean position and variability of the eddy-driven jet, 
    respectively.

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    r_o3jetdist : numpy.ndarray     
        Pearson correlation coefficient between O3 and jet distance, [lat, lng]
    r_t2mjetdist : numpy.ndarray     
        Pearson correlation coefficient between 2-meter temperature and jet 
        distance, [lat, lng]
    r_qv2mjetdist : numpy.ndarray     
        Pearson correlation coefficient between 2-meter specific humidity and 
        jet distance, [lat, lng]
    lat_jet_ml : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in the Northern Hemisphere mid-latitudes, units of degrees north, 
        [time, lng]
    significance_r_o3jetdist : numpy.ndarray
        Significance of r(O3, jet latitude) determined with moving block 
        bootstrapping: 1 implies significance, NaN implies insignificance, 
        [lat, lng]
    significance_r_t2mjetdist : numpy.ndarray
        Significance of r(T, jet latitude) determined with moving block 
        bootstrapping: 1 implies significance, NaN implies insignificance, 
        [lat, lng]
    significance_r_qv2mjetdist : numpy.ndarray
        Significance of r(q, jet latitude) determined with moving block 
        bootstrapping: 1 implies significance, NaN implies insignificance, 
        [lat, lng]        

    Returns
    -------
    None        
    """
    import numpy as np
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth']=0.3     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter    
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,7.7))
    if lng_gmi[-1] != 360:
        lng_gmi[-1] = 360.    
    # r(O3, jet lat - lat)
    ax1 = plt.subplot2grid((3,2), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax1.set_title(r'(a) r(O$_{\mathregular{3}}$, '+
        '$\mathregular{\phi_{jet}}$ $-$ ${\mathregular{\phi}}$)', fontsize=16, 
        x=0.02, ha='left')    
    # ax1.add_feature(ocean50m, zorder=3)
    ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    cmap = plt.get_cmap('coolwarm')
    mb = ax1.contourf(lng_gmi, lat_gmi, r_o3jetdist, np.linspace(-1,1,9), 
        cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(O3, jet lat - lat)
    ax1.contourf(lng_gmi, lat_gmi, significance_r_o3jetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    ax1.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax1.xaxis.set_major_formatter(lng_formatter)         
    ax1.get_xaxis().set_ticklabels([])    
    ax1.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax1.yaxis.set_major_formatter(lat_formatter)    
    ax1.outline_patch.set_zorder(20)
    # r(T, jet lat - lat)
    ax2 = plt.subplot2grid((3,2), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2.set_title(r'(b) r(T, $\mathregular{\phi_{jet}}$ $-$ ' +
        '${\mathregular{\phi}}$)', fontsize=16, x=0.02, ha='left')    
    # ax2.add_feature(ocean50m, zorder=3)
    ax2.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax2.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    mb = ax2.contourf(lng_gmi, lat_gmi, r_o3jetdist,
        np.linspace(-1,1,9), cmap=cmap, extend='neither', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmi, lat_gmi, significance_r_t2mjetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    ax2.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax2.xaxis.set_major_formatter(lng_formatter)         
    ax2.get_xaxis().set_ticklabels([])    
    ax2.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax2.yaxis.set_major_formatter(lat_formatter)    
    ax2.outline_patch.set_zorder(20)
    # r(q, jet lat - lat)
    ax3 = plt.subplot2grid((3,2), (2,0), colspan=2,
            projection=ccrs.PlateCarree(central_longitude=0.))
    ax3.set_title(r'(c) r(q, $\mathregular{\phi_{jet}}$ $-$ ' +
        '${\mathregular{\phi}}$)', fontsize=16, x=0.02, ha='left')    
    # ax3.add_feature(ocean50m, zorder=3)
    ax3.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax3.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    mb = ax3.contourf(lng_gmi, lat_gmi, r_qv2mjetdist,
        np.linspace(-1,1,9), cmap=cmap, extend='neither', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax3.contourf(lng_gmi, lat_gmi, significance_r_qv2mjetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    ax3.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax3.xaxis.set_major_formatter(lng_formatter)         
    ax3.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax3.yaxis.set_major_formatter(lat_formatter)    
    ax3.outline_patch.set_zorder(20)
    for ax in [ax1, ax2, ax3]:
        skiplng = 6
        ax.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,
            axis=0)[::skiplng], yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], 
            color='k', markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree(), zorder=5)
    # Add colorbar
    plt.gcf().subplots_adjust(left=0.05, right=0.86, hspace=0.3)
    colorbar_axes = plt.gcf().add_axes([
        ax1.get_position().x1+0.03, # Left
        (ax3.get_position().y1-ax3.get_position().y0)/2.+ax3.get_position().y0, # Bottom 
        0.02, # Width
        ((ax1.get_position().y1-ax1.get_position().y0)/2.+ax1.get_position().y0)-
        ((ax3.get_position().y1-ax3.get_position().y0)/2.+ax3.get_position().y0)])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-1,1,9), extend='neither')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[$\mathregular{\cdot}$]', fontsize=16)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'figS2.pdf', dpi=600)
    return

def figS3(lat_gmi, lng_gmi, r_pblhjetdist, r_U10Mjetdist, r_V10Mjetdist, 
    lat_jet_ml, significance_r_pblhjetdist, significance_r_U10Mjetdist, 
    significance_r_V10Mjetdist): 
    """Figure S3 of Kerr et al. (2020). (a) r(PBLH, jet distance) is shown with 
    colored shading. Hatching indicates insignificance determined with moving 
    block bootstrapping. (b) Same as (a) but showing r(U10, jet distance). (c)
    Same as (a) but showing r(V10, jet distance). Scatterpoints and error bars 
    in (a-c) specify the mean position and variability of the eddy-driven jet, 
    respectively.

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    r_pblhjetdist : numpy.ndarray     
        Pearson correlation coefficient between PBLH and O3, [lat, lng]
    r_U10Mjetdist : numpy.ndarray     
        Pearson correlation coefficient between 10-meter zonal (eastward) wind 
        and jet distance, [lat, lng]
    r_V10Mjetdist : numpy.ndarray     
        Pearson correlation coefficient between 10-meter meridional 
        (northward) wind and jet distance, [lat, lng]
    lat_jet_ml : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in the Northern Hemisphere mid-latitudes, units of degrees north, 
        [time, lng]
    significance_r_pblhjetdist : numpy.ndarray
        Significance of r(PBLH, jet latitude) determined with moving block 
        bootstrapping: 1 implies significance, NaN implies insignificance, 
        [lat, lng]
    significance_r_U10Mjetdist : numpy.ndarray
        Significance of r(U10, jet latitude) determined with moving block 
        bootstrapping: 1 implies significance, NaN implies insignificance, 
        [lat, lng]
    significance_r_V10Mjetdist : numpy.ndarray
        Significance of r(V10, jet latitude) determined with moving block 
        bootstrapping: 1 implies significance, NaN implies insignificance, 
        [lat, lng]        

    Returns
    -------
    None     
    """
    import numpy as np
    import matplotlib as mpl
    mpl.rcParams['hatch.linewidth']=0.3     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter    
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,7.7))
    if lng_gmi[-1] != 360:
        lng_gmi[-1] = 360.    
    # r(PBLH, jet lat - lat)
    ax1 = plt.subplot2grid((3,2), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax1.set_title(r'(a) r(PBLH, $\mathregular{\phi_{jet}}$ $-$ '+
        '${\mathregular{\phi}}$)', fontsize=16, x=0.02, ha='left')    
    # ax1.add_feature(ocean50m, zorder=3)
    ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    cmap = plt.get_cmap('coolwarm')
    mb = ax1.contourf(lng_gmi, lat_gmi, r_pblhjetdist, np.linspace(-1,1,9), 
        cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(PBLH, jetlat)
    ax1.contourf(lng_gmi, lat_gmi, significance_r_pblhjetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    ax1.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax1.xaxis.set_major_formatter(lng_formatter)         
    ax1.get_xaxis().set_ticklabels([])    
    ax1.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax1.yaxis.set_major_formatter(lat_formatter)      
    ax1.outline_patch.set_zorder(20)
    # r(U10M, jet lat - lat)
    ax2 = plt.subplot2grid((3,2), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2.set_title(r'(b) r(U$_\mathregular{10}$, '+
        '$\mathregular{\phi_{jet}}$ $-$ ${\mathregular{\phi}}$)', fontsize=16, 
        x=0.02, ha='left')    
    # ax2.add_feature(ocean50m, zorder=3)
    ax2.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax2.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    mb = ax2.contourf(lng_gmi, lat_gmi, r_U10Mjetdist,
        np.linspace(-1, 1, 9), cmap=cmap, extend='neither', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmi, lat_gmi, significance_r_U10Mjetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    ax2.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax2.xaxis.set_major_formatter(lng_formatter)         
    ax2.get_xaxis().set_ticklabels([])    
    ax2.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax2.yaxis.set_major_formatter(lat_formatter)      
    ax2.outline_patch.set_zorder(20)
    # r(V10M, jet lat - lat)
    ax3 = plt.subplot2grid((3,2), (2,0), colspan=2,
            projection=ccrs.PlateCarree(central_longitude=0.))
    ax3.set_title(r'(c) r(V$_\mathregular{10}$, '+
        '$\mathregular{\phi_{jet}}$ $-$ ${\mathregular{\phi}}$)', fontsize=16, 
        x=0.02, ha='left')    
    # ax3.add_feature(ocean50m, zorder=3)
    ax3.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax3.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    mb = ax3.contourf(lng_gmi, lat_gmi, r_V10Mjetdist,
        np.linspace(-1, 1, 9), cmap=cmap, extend='neither', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax3.contourf(lng_gmi, lat_gmi, significance_r_V10Mjetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    ax3.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    lng_formatter = LongitudeFormatter()
    ax3.xaxis.set_major_formatter(lng_formatter)         
    ax3.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    lat_formatter = LatitudeFormatter()
    ax3.yaxis.set_major_formatter(lat_formatter)      
    ax3.outline_patch.set_zorder(20)
    for ax in [ax1, ax2, ax3]:
        skiplng = 6
        ax.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,
            axis=0)[::skiplng], yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], 
            color='k', markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree(), zorder=5)
    # Add colorbar
    plt.gcf().subplots_adjust(left=0.05, right=0.86, hspace=0.3)
    colorbar_axes = plt.gcf().add_axes([
        ax1.get_position().x1+0.03, # Left
        (ax3.get_position().y1-ax3.get_position().y0)/2.+ax3.get_position().y0, # Bottom 
        0.02, # Width
        ((ax1.get_position().y1-ax1.get_position().y0)/2.+ax1.get_position().y0)-
        ((ax3.get_position().y1-ax3.get_position().y0)/2.+ax3.get_position().y0)])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(-1,1,9), extend='neither')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[$\mathregular{\cdot}$]', fontsize=16)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'figS3.pdf', dpi=600)
    return 

import numpy as np
import pandas as pd
import sys
sys.path.append('/Users/ghkerr/phd/globalo3/')
import globalo3_open, globalo3_calculate, observations_open
sys.path.append('/Users/ghkerr/phd/transporto3/')
# Load data the first iteration only 
try:lat_gmi
except NameError:
    datapath = '/Users/ghkerr/phd/globalo3/data/parsed/'      
    import pandas as pd
    import netCDF4 as nc
    # Load data
    o3_gmi = nc.Dataset(datapath+'gmi_O3control_JJA2008-2010.nc')['O3_control'][:].data
    o3_gmi_china = nc.Dataset(datapath+'gmi_O3control_JJA2016-2017.nc')['O3_control'][:].data
    lat_gmi = nc.Dataset(datapath+'gmi_O3control_JJA2008-2010.nc')['lat'][:].data
    lng_gmi = nc.Dataset(datapath+'gmi_O3control_JJA2008-2010.nc')['lng'][:].data
    mtime = nc.Dataset(datapath+'gmi_O3control_JJA2008-2010.nc')['time'][:].data
    o3_transport_gmi = nc.Dataset(datapath+'gmi_O3transportonly_JJA2008-2010.nc')['O3_transportonly'][:].data
    t2m_merra = nc.Dataset(datapath+'merra2_t2m_JJA2008-2010.nc')['T2M'][:].data
    t2m_merra_china = nc.Dataset(datapath+'merra2_t2m_JJA2016-2017.nc')['T2M'][:].data
    qv2m_merra = nc.Dataset(datapath+'merra2_qv2m_JJA2008-2010.nc')['QV2M'][:].data
    qv2m_merra_china = nc.Dataset(datapath+'merra2_qv2m_JJA2016-2017.nc')['QV2M'][:].data    
    lat_jet_ml = nc.Dataset(datapath+'merra2_JET_JJA2008-2010.nc')['jetlatitude'][:].data
    U10M = nc.Dataset(datapath+'merra2_U10M_JJA2008-2010.nc')['U10M'][:].data
    V10M = nc.Dataset(datapath+'merra2_V10M_JJA2008-2010.nc')['V10M'][:].data
    pblh_merra = nc.Dataset(datapath+'merra2_PBLH_JJA2008-2010.nc')['PBLH'][:].data
    cyclones_binned = nc.Dataset(datapath+'mcms_CYCLONEFREQ_JJA2008-2010.nc')['cyclones'][:].data
    lat_cyclones_binned = nc.Dataset(datapath+'mcms_CYCLONEFREQ_JJA2008-2010.nc')['lat'][:].data
    lng_cyclones_binned = nc.Dataset(datapath+'mcms_CYCLONEFREQ_JJA2008-2010.nc')['lng'][:].data
    eqjet_cyclones_binned = nc.Dataset(datapath+'mcms_CYCLONEFREQ_JJA2008-2010.nc')['ewjet_cyclones'][:].data
    pwjet_cyclones_binned = nc.Dataset(datapath+'mcms_CYCLONEFREQ_JJA2008-2010.nc')['pwjet_cyclones'][:].data
    cyclones = pd.read_pickle(datapath+'cyclones.pickle')
    # Load calculations of significance
    significance_r_t2mo3 = nc.Dataset(datapath+'sig_merra2_t2m_gmi_O3control.nc')['sig_T2M_O3_control'][:].data
    significance_r_qv2mo3 = nc.Dataset(datapath+'sig_merra2_qv2m_gmi_O3control.nc')['sig_QV2M_O3_control'][:].data
    significance_r_t2mo3_transport = nc.Dataset(datapath+'sig_merra2_t2m_gmi_O3transportonly.nc')['sig_T2M_O3_transportonly'][:].data
    significance_r_qv2mo3_transport = nc.Dataset(datapath+'sig_merra2_qv2m_gmi_O3transportonly.nc')['sig_QV2M_O3_transportonly'][:].data
    significance_r_t2mjetdist = nc.Dataset(datapath+'sig_merra2_jetdist_merra2_t2m.nc')['sig_JETDIST_merra2_t2m'][:].data
    significance_r_o3jetdist = nc.Dataset(datapath+'sig_merra2_jetdist_gmi_o3control.nc')['sig_JETDIST_O3_control'][:].data
    significance_r_qv2mjetdist = nc.Dataset(datapath+'sig_merra2_jetdist_merra2_qv2m.nc')['sig_JETDIST_merra2_qv2m'][:].data 
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
    # Slope and correlation of O3/jet distance and 2-meter temperature and 
    # jet distance
    m_o3jetdist, r_o3jetdist, diff_o3jetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(o3_gmi, lat_gmi, 
        lng_gmi, lat_jet_ml, lng_gmi)
    m_t2mjetdist, r_t2mjetdist, diff_t2mjetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(t2m_merra, lat_gmi, 
        lng_gmi, lat_jet_ml, lng_gmi)
    m_qv2mjetdist, r_qv2mjetdist, diff_qv2mjetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(qv2m_merra, lat_gmi, 
        lng_gmi, lat_jet_ml, lng_gmi)
    m_U10Mjetdist, r_U10Mjetdist, diff_U10Mjetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(U10M, lat_gmi, 
        lng_gmi, lat_jet_ml, lng_gmi)
    m_V10Mjetdist, r_V10Mjetdist, diff_V10Mjetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(V10M, lat_gmi, 
        lng_gmi, lat_jet_ml, lng_gmi)  
    m_pblhjetdist, r_pblhjetdist, diff_pblhjetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(pblh_merra, 
        lat_gmi, lng_gmi, lat_jet_ml, lng_gmi)
    times_gmi = pd.date_range(start='01/01/2008', end='12/31/2010')
    where_months = np.where(np.in1d(times_gmi.month, 
        np.array([6, 7, 8]))==True)[0]
    times_gmi = times_gmi[where_months]
    times_gmi = [x.to_pydatetime().date() for x in times_gmi]        
    times_gmi = np.array(times_gmi)
    times_gmi_china = pd.date_range(start='01/01/2016', end='12/31/2017')
    where_months = np.where(np.in1d(times_gmi_china.month, 
        np.array([6, 7, 8]))==True)[0]
    times_gmi_china = times_gmi_china[where_months]
    times_gmi_china = [x.to_pydatetime().date() for x in times_gmi_china]
    times_gmi_china = np.array(times_gmi_china)
    # Cyclones 
    o3_anom, o3_all, o3_anom_rotated, o3_rotated = \
        globalo3_calculate.o3anom_cyclone(cyclones, times_gmi, lat_gmi, 
        lng_gmi, o3_gmi)    
    # Load observational datasets
    promptobs = input('Load observational datasets? [y/n]\n').lower()    
    if promptobs == 'y':    
        naps = observations_open.open_napso3([2008, 2009, 2010], [6, 7, 8], [14])
        aqs = observations_open.open_aqso3([2008, 2009, 2010], [6, 7, 8], [14])
        emep = observations_open.open_emepo3([2008, 2009, 2010], [6, 7, 8], [14])
        china = observations_open.open_chinao3([2016, 2017], [6, 7, 8], [14],
            cityavg=True)
        # Calculate O3-meteorology relationships from observational datasets
        (r_t2mo3_naps, do3dt2m_naps, ro3jet_naps, do3djet_naps, lat_naps, 
            lng_naps) = globalo3_calculate.calculate_obs_o3_temp_jet(naps, 
            t2m_merra, times_gmi, lat_gmi, lng_gmi, lat_jet_ml, lng_gmi)
        (r_t2mo3_aqs, do3dt2m_aqs, ro3jet_aqs, do3djet_aqs, lat_aqs, 
            lng_aqs) = globalo3_calculate.calculate_obs_o3_temp_jet(aqs, 
            t2m_merra, times_gmi, lat_gmi, lng_gmi, lat_jet_ml, lng_gmi)    
        (r_t2mo3_emep, do3dt2m_emep, ro3jet_emep, do3djet_emep, lat_emep, 
            lng_emep) = globalo3_calculate.calculate_obs_o3_temp_jet(emep, 
            t2m_merra, times_gmi, lat_gmi, lng_gmi, lat_jet_ml, lng_gmi)
        (r_qv2mo3_naps, do3dqv2m_naps, ro3jet_naps, do3djet_naps, lat_naps, 
            lng_naps) = globalo3_calculate.calculate_obs_o3_temp_jet(naps, 
            qv2m_merra, times_gmi, lat_gmi, lng_gmi, lat_jet_ml, lng_gmi)
        (r_qv2mo3_aqs, do3dqv2m_aqs, ro3jet_aqs, do3djet_aqs, lat_aqs, 
            lng_aqs) = globalo3_calculate.calculate_obs_o3_temp_jet(aqs, 
            qv2m_merra, times_gmi, lat_gmi, lng_gmi, lat_jet_ml, lng_gmi)    
        (r_qv2mo3_emep, do3dqv2m_emep, ro3jet_emep, do3djet_emep, lat_emep, 
            lng_emep) = globalo3_calculate.calculate_obs_o3_temp_jet(emep, 
            qv2m_merra, times_gmi, lat_gmi, lng_gmi, lat_jet_ml, lng_gmi)        
        # For China (2016-2017)
        (r_t2mo3_china, do3dt2m_china, ro3jet_china, do3djet_china, lat_china,
            lng_china) = globalo3_calculate.calculate_obs_o3_temp_jet(china, 
            t2m_merra_china, np.array(times_gmi_china), lat_gmi, lng_gmi, 
            lat_jet_ml[:len(t2m_merra_china)], lng_gmi)
        (r_qv2mo3_china, do3dqv2m_china, ro3jet_china, do3djet_china, 
            lat_china, lng_china) = \
            globalo3_calculate.calculate_obs_o3_temp_jet(china, 
            qv2m_merra_china, np.array(times_gmi_china), lat_gmi, lng_gmi, 
            lat_jet_ml[:len(t2m_merra_china)], lng_gmi)
        # Calculate the air quality performance index (AQPI), correlation 
        # coefficient (r), factor-of-2 fraction (FAC2), and mean fractional 
        # bias (MFB) for each observational dataset
        aqpi_aqs, fac2_aqs, r_aqs, mfb_aqs = globalo3_calculate.calculate_aqpi(aqs, 
            o3_gmi, lat_gmi, lng_gmi, times_gmi, 'gridwise')
        aqpi_naps, fac2_naps, r_naps, mfb_naps = globalo3_calculate.calculate_aqpi(
            naps, o3_gmi, lat_gmi, lng_gmi, times_gmi, 'gridwise')
        aqpi_emep, fac2_emep, r_emep, mfb_emep = globalo3_calculate.calculate_aqpi(
            emep, o3_gmi, lat_gmi, lng_gmi, times_gmi, 'gridwise')
        aqpi_china, fac2_china, r_china, mfb_china = \
            globalo3_calculate.calculate_aqpi(china, o3_gmi_china, lat_gmi,
            lng_gmi, times_gmi_china, 'gridwise')
    # # Significance (alpha = 0.05) of O3-temperature-q-jet correlations
    # promptsignificance = input('Calculate significance? [y/n]\n').lower()    
    # if promptsignificance == 'y':
        # significance_r_U10Mjetdist = \
        #     globalo3_calculate.calculate_r_significance(U10M, diff_U10Mjetdist, 
        #     r_U10Mjetdist, lat_gmi, lng_gmi)
        # significance_r_V10Mjetdist = \
        #     globalo3_calculate.calculate_r_significance(V10M, diff_V10Mjetdist, 
        #     r_V10Mjetdist, lat_gmi, lng_gmi)
        # significance_r_pblhjetdist = \
        #     globalo3_calculate.calculate_r_significance(pblh_merra,
        #     diff_pblhjetdist, r_pblhjetdist, lat_gmi, lng_gmi)
                        
# # FIGURE 1; Mean O3 and NOx
# fig1(lat_gmi, lng_gmi, o3_gmi, lat_jet_ml)
# FIGURE 2; model performance
# fig2(lat_gmi, lng_gmi, r_aqs, r_naps, r_emep, r_china)
# # FIGURE 3; r(T, O3) and r(q, O3)
# fig3(lat_gmi, lng_gmi, r_t2mo3, r_qv2mo3, significance_r_t2mo3, 
#     significance_r_qv2mo3, lat_jet_ml)
# # FIGURE 4; zonally-averaged O3-meteorology relationships 
# fig4(lat_gmi, lng_gmi, r_t2mo3, r_qv2mo3, r_t2mo3_aqs, r_qv2mo3_aqs, 
#     lat_aqs, lng_aqs, r_t2mo3_naps, r_qv2mo3_naps, lat_naps, lng_naps, 
#     r_t2mo3_emep, r_qv2mo3_emep, lat_emep, lng_emep, r_t2mo3_china, 
#     r_qv2mo3_china, lat_china, lng_china, lng_gmi, lat_jet_ml)
# # FIGURE 5: r(T, O3) and r(q, O3) from the transport only simulation 
# fig5(lat_gmi, lng_gmi, r_t2mo3, r_t2mo3_transport, r_qv2mo3, 
#     r_qv2mo3_transport, significance_diff_r_t2mo3, significance_diff_r_qv2mo3, 
#     lat_jet_ml)
# # FIGURE 6; difference in O3, T2M, and qv2M on days with a poleward versus
# # equatorward jet
# fig6(lat_gmi, lng_gmi, o3_gmi, t2m_merra, qv2m_merra, lat_jet_ml, 
#     times_gmi, significance_r_o3jetdist, significance_r_t2mjetdist, 
#     significance_r_qv2mjetdist)
# # FIGURE 7; cyclone frequency and poleward-equatorward jet differences
# fig7(lat_cyclones_binned, lng_cyclones_binned, cyclones_binned,
#     pwjet_cyclones_binned, eqjet_cyclones_binned, lat_gmi, lng_gmi, 
#     lat_jet_ml)
# # FIGURE 8; O3 anomaly at cyclone
# fig8(o3_anom_rotated)
# # FIGURE 9; difference in PBLH, U10, and V10 on days with a poleward versus 
# # equatorward jet
# fig9(lat_gmi, lng_gmi, pblh_merra, U10M, V10M, lat_jet_ml, times_gmi, 
#     significance_r_pblhjetdist, significance_r_U10Mjetdist, 
#     significance_r_V10Mjetdist)
# # FIGURE 10; zonally-averaged mean and eddy fluxes
# fig10(lat_gmi, V10M, o3_gmi, t2m_merra, qv2m_merra)    
# FIGURE S1; dO3/dT and dO3/dq
# figS1(lat_gmi, lng_gmi, do3dt2m, do3dq, significance_r_t2mo3, 
#     significance_r_qv2mo3, lat_jet_ml)
# # FIGURE S2: r(O3, jet distance), r(T, jet distance), and r(q, jet distance)
# figS2(lat_gmi, lng_gmi, r_o3jetdist, r_t2mjetdist, r_qv2mjetdist, 
#     lat_jet_ml, significance_r_o3jetdist, significance_r_t2mjetdist, 
#     significance_r_qv2mjetdist)
# # FIGURE S3: r(PBLH, jet distance), r(U10, jet distance), and r(V10, jet
# # distance)
# figS3(lat_gmi, lng_gmi, r_pblhjetdist, r_U10Mjetdist, r_V10Mjetdist, 
#     lat_jet_ml, significance_r_pblhjetdist, significance_r_U10Mjetdist, 
#     significance_r_V10Mjetdist)










do3dt2m_sma = globalo3_calculate.calculate_do3dt_sma(t2m_merra, o3_gmi, lat_gmi, lng_gmi)
do3dq_sma = globalo3_calculate.calculate_do3dt_sma(qv2m_merra, o3_gmi, lat_gmi, lng_gmi)





# # Very negative
# i = 27
# j = 100
# # # Very positive 
# # i = 26
# # j = 50
# lat = lat_gmi
# lng = lng_gmi
# X = t2m_merra[:, i, j]
# y = o3_gmi[:, i, j]
# XX = np.linspace(X.min(), X.max(), 1000)
# # OLS
# p = np.polyfit(X, y, 1)
# yfit = np.polyval(p,XX)
# # SMA (ratio of standard deviations)
# sx = np.std(X)
# sy = np.std(y)
# sign = np.sign(np.corrcoef(X, y)[0,1])
# sma = sign*(sy/sx)
# # the intercept of an RMA regression can simply be calculated from the 
# # equation of a line once the slope is known
# # https://www2.clarku.edu/faculty/pbergmann/Resources/Biol206%20-%20Lab04-%20Bivariate%20Regression.pdf
# sma_intercept = np.mean(y) - (sma*np.mean(X))
# fig = plt.figure()
# plt.plot(X, y, 'ko', label='GMI O3/MERRA2 T');
# plt.plot(XX, yfit, '-r', label='OLS; slope=%.2f'%p[0])
# plt.plot(XX, (sma*XX)+sma_intercept, '-b', label='SMA; slope=%.2f'%sma)

# plt.legend()
# plt.xlabel('T [K]')
# plt.ylabel('O3 [ppbv]')
# plt.title('lat = %.1f deg, lng = %.1f deg'%(lat[i], lng[j]))
# plt.savefig('/Users/ghkerr/Desktop/TLS_SMA_negative.png', dpi=300)
# plt.show()




  