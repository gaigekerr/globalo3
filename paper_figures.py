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
    30062020 -- addressing reviewer comments 
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
    from cartopy.util import add_cyclic_point
    # Load EDGAR NOx
    lat_edgar, lng_edgar, nox_edgar = globalo3_open.open_edgar_specifieddomain(
        [2008, 2009, 2010], -1., 90., 0., 360., 'NOx')
    # For wrapping around the Prime Meridian
    o3_gmi, lng_gmip = add_cyclic_point(o3_gmi, coord=lng_gmi)
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
    mb = ax1.contourf(lng_gmip, lat_gmi, np.mean(o3_gmi, axis=0), 
        np.linspace(10, 60, 11), cmap=cmap, extend='both',
        transform=ccrs.PlateCarree(), zorder=1)
    csthick = ax1.contour(lng_gmip, lat_gmi, np.nanstd(o3_gmi, axis=0), [10.], 
        colors='k', linewidths=1.5, transform=ccrs.PlateCarree(), zorder=15)
    csmedium = ax1.contour(lng_gmip, lat_gmi, np.nanstd(o3_gmi, axis=0), [8.], 
        colors='k', linestyles='--', linewidths=0.75, 
        transform=ccrs.PlateCarree(), zorder=15)
    skiplng = 6
    ax1.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
        yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], zorder=10, color='k', 
        markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
        transform=ccrs.PlateCarree())
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
        'fig1.pdf', dpi=600)
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
    from cartopy.util import add_cyclic_point
    # For wrapping around the Prime Meridian
    r_t2mo3, lng_gmip = add_cyclic_point(r_t2mo3, coord=lng_gmi)
    r_qv2mo3, lng_gmip = add_cyclic_point(r_qv2mo3, coord=lng_gmi)
    significance_r_t2mo3, lng_gmip = add_cyclic_point(significance_r_t2mo3, 
        coord=lng_gmi)
    significance_r_qv2mo3, lng_gmip = add_cyclic_point(significance_r_qv2mo3, 
        coord=lng_gmi)
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,5))
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
    mb = ax1.contourf(lng_gmip, lat_gmi, r_t2mo3, np.linspace(-1, 1, 9), 
        cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(T, O3)
    ax1.contourf(lng_gmip, lat_gmi, significance_r_t2mo3, hatches=['//////'], 
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
    mb = ax2.contourf(lng_gmip, lat_gmi, r_qv2mo3, np.linspace(-1,1,9), 
        cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmip, lat_gmi, significance_r_qv2mo3, hatches=['//////'], 
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
    from cartopy.util import add_cyclic_point
    # For wrapping around the Prime Meridian
    r_t2mo3, lng_gmip = add_cyclic_point(r_t2mo3, coord=lng_gmi)
    r_t2mo3_transport, lng_gmip = add_cyclic_point(r_t2mo3_transport, 
        coord=lng_gmi)
    significance_diff_r_t2mo3, lng_gmip = add_cyclic_point(
        significance_diff_r_t2mo3, coord=lng_gmi)
    r_qv2mo3, lng_gmip = add_cyclic_point(r_qv2mo3, coord=lng_gmi)
    r_qv2mo3_transport, lng_gmip = add_cyclic_point(r_qv2mo3_transport, 
        coord=lng_gmi)
    significance_diff_r_qv2mo3, lng_gmip = add_cyclic_point(
        significance_diff_r_qv2mo3, coord=lng_gmi)
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,5))
    # For wrapping around the Prime Meridian
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
    mb = ax1.contourf(lng_gmip, lat_gmi, (r_t2mo3-r_t2mo3_transport), 
        np.linspace(-0.4, 0.4, 9), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(T, O3)
    ax1.contourf(lng_gmip, lat_gmi, significance_diff_r_t2mo3, 
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
    mb = ax2.contourf(lng_gmip, lat_gmi, r_qv2mo3-r_qv2mo3_transport, 
        np.linspace(-0.4,0.4,9), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmip, lat_gmi, significance_diff_r_qv2mo3, 
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
    from cartopy.util import add_cyclic_point
    # For wrapping around the Prime Meridian
    pwjet_o3, lng_gmip = add_cyclic_point(pwjet_o3, coord=lng_gmi)
    eqjet_o3, lng_gmip = add_cyclic_point(eqjet_o3, coord=lng_gmi)
    pwjet_t2m, lng_gmip = add_cyclic_point(pwjet_t2m, coord=lng_gmi)
    eqjet_t2m, lng_gmip = add_cyclic_point(eqjet_t2m, coord=lng_gmi)
    pwjet_qv2m, lng_gmip = add_cyclic_point(pwjet_qv2m, coord=lng_gmi)
    eqjet_qv2m, lng_gmip = add_cyclic_point(eqjet_qv2m, coord=lng_gmi)
    significance_r_o3jetdist, lng_gmip = add_cyclic_point(
        significance_r_o3jetdist, coord=lng_gmi)
    significance_r_t2mjetdist, lng_gmip = add_cyclic_point(
        significance_r_t2mjetdist, coord=lng_gmi)
    significance_r_qv2mjetdist, lng_gmip = add_cyclic_point(
        significance_r_qv2mjetdist, coord=lng_gmi)        
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,7.7))
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
    mb = ax1.contourf(lng_gmip, lat_gmi, (pwjet_o3-eqjet_o3), 
        np.linspace(-8,8,9), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(O3, jetlat)
    ax1.contourf(lng_gmip, lat_gmi, significance_r_o3jetdist, 
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
    mb = ax2.contourf(lng_gmip, lat_gmi, (pwjet_t2m-eqjet_t2m), 
        np.linspace(-8,8,9), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmip, lat_gmi, significance_r_t2mjetdist, 
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
    mb = ax3.contourf(lng_gmip, lat_gmi, (pwjet_qv2m-eqjet_qv2m), 
        np.linspace(-3,3,7), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax3.contourf(lng_gmip, lat_gmi, significance_r_qv2mjetdist, 
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
  
def fig9(lat_gmi, lng_gmi, o3_gmi, V10M, lat_jet_ml, times_gmi, 
    significance_r_V10Mjetdist, r_V10Mo3, significance_r_V10Mo3):
    """Figure 9 of Kerr et al. (2020); (a) Difference in V10M on days with a
    poleward jet versus days with an equatorward jet. Hatched grid cells 
    correspond to non-statistically significant values of r(V10M, jet latitude) 
    determined with moving block bootstrapping. (b) r(V10M, O3) with hatching
    indicating grad cells with non-statistically significant values of 
    r(V10M, O3). Stippling corresponds to grid cells where the latitudinal  
    gradient of O3 is positive. Scatterpoints and vertical bars in (a-c) 
    specify the mean position and variability of the eddy-driven jet, 
    respectively, for days with a poleward jet and equatorward jet. 

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    o3_gmi : numpy.ndarray     
        Daily afternoon surface-level O3 from the GMI CTM, units of ppbv, 
        [time, lat, lng]  
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
    mpl.rcParams['hatch.linewidth']=0.5     
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter        
    fig = plt.figure(figsize=(9,5))
    # Composites of flow at surface on days with poleward or equatorward jet 
    eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_V10M, eqjet_V10M = \
        globalo3_calculate.segregate_field_bylat(V10M, lng_gmi, lat_jet_ml, 
        times_gmi)
    from cartopy.util import add_cyclic_point
    # For wrapping around the Prime Meridian
    pwjet_V10M, lng_gmip = add_cyclic_point(pwjet_V10M, coord=lng_gmi)
    eqjet_V10M, lng_gmip = add_cyclic_point(eqjet_V10M, coord=lng_gmi)
    significance_r_V10Mjetdist, lng_gmip = add_cyclic_point(
        significance_r_V10Mjetdist, coord=lng_gmi)
    r_V10Mo3, lng_gmip = add_cyclic_point(r_V10Mo3, coord=lng_gmi)
    significance_r_V10Mo3, lng_gmip = add_cyclic_point(significance_r_V10Mo3, 
        coord=lng_gmi)
    # (PW - EW) V10M difference
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax1.set_title(r'(a) V$_\mathregular{10,\:PW}$ $-$ V$_\mathregular{10,\:EW}$', 
        fontsize=16, x=0.02, ha='left')
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
    mb = ax1.contourf(lng_gmip, lat_gmi, (pwjet_V10M-eqjet_V10M), 
        np.linspace(-3, 3, 7), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(V10M, jetlat)
    ax1.contourf(lng_gmip, lat_gmi, significance_r_V10Mjetdist, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
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
        ticks=np.linspace(-3, 3, 7), extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[$\mathregular{\cdot}$]', fontsize=16, labelpad=24)
    ax1.outline_patch.set_zorder(20)
    # r(V10M, O3)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2.set_title(r'(b) r(V$_\mathregular{10}$, O$_\mathregular{3}$)', 
        fontsize=16, x=0.02, ha='left')    
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
    mb = ax2.contourf(lng_gmip, lat_gmi, r_V10Mo3,
        np.linspace(-1, 1, 9), cmap=cmap, extend='neither', 
        transform=ccrs.PlateCarree(), zorder=1)
    # r(V10M, O3) significance 
    ax2.contourf(lng_gmip, lat_gmi, significance_r_V10Mo3, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree(), 
        zorder=2)
    # O3 gradient (stippling where grad(O3) > 0)
    grad_o3 = np.empty(shape=o3_gmi.shape[1:])
    grad_o3[:] = np.nan
    for spine_i in np.arange(0, len(lng_gmi), 1):
        # O3 concentrations for given "spine" of latitudes at a given 
        # longitude
        spine = np.nanmean(o3_gmi, axis=0)[:,spine_i]
        grad_spine = np.gradient(spine)
        grad_o3[:,spine_i] = grad_spine
    grad_o3 = np.array(grad_o3)
    grad_o3[grad_o3 > 0] = 1.
    grad_o3[grad_o3 <= 0] = np.nan
    ax2.contourf(lng_gmi, lat_gmi, grad_o3, hatches=['..'], 
        colors='none', transform=ccrs.PlateCarree(), zorder=3)
    # Eddy-driven jet
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
        ticks=np.linspace(-1, 1, 9), extend='neither')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[$\mathregular{\cdot}$]', labelpad=8, fontsize=16)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'fig9.pdf', dpi=600)
    return

def figS1(lat_gmi, lng_gmi, do3dt2m, do3dt2m_transport, do3dq, do3dq_transport,
    significance_r_t2mo3, significance_r_qv2mo3, 
    significance_r_t2mo3_transport, significance_r_qv2mo3_transport,
    lat_jet_ml): 
    """Figure S1 of Kerr et al. (2020). (a) The slope of the ordinary least 
    squares (OLS) regression of O3 versus temperature, dO3/dT. Hatching denotes 
    regions where the correlation between O3 and temperature is insignificant,
    determined using moving block bootstrap resampling to estimate the 95% 
    confidence interval. (b) Same as (a) but for O3 versus specific humidity,
    dO3/dq, with hatching showing insignificant correlation between O3 and 
    specific humidity. (c) Same as (a) but using O3 from the transport-only
    simulation. (d) Same as (b) but using O3 from the transport-only 
    simulation. (e) and (f) show (a)-(c) and (b)-(d), respectively.     
    Scatter points and vertical bars are identical in (a-f) and show the mean 
    latitude of the eddy-driven jet and its variability. 

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    do3dt2m : numpy.ndarray     
        The slope of the ordinary least squares linear regression of O3 versus
        2-meter temperature, units of ppbv K-1, [lat, lng]
    do3dt2m_transport : numpy.ndarray     
        The slope of the ordinary least squares linear regression of O3 from 
        the transport-only simulation versus 2-meter temperature, units of ppbv 
        K-1, [lat, lng]    
    do3dqv2m : numpy.ndarray     
        The slope of the ordinary least squares linear regression of O3 versus
        2-meter specific humidity, units of ppbv kg g-1, [lat, lng]
    do3dq_transport : numpy.ndarray     
        The slope of the ordinary least squares linear regression of O3 from 
        the transport-only simulation versus 2-meter specific humidity, units 
        of ppbv kg g-1, [lat, lng]
    significance_r_t2mo3 : numpy.ndarray        
        Significance of r(T, O3) determined with moving block bootstrapping: 1
        implies significance, NaN implies insignificance, [lat, lng]
    significance_r_qv2mo3 : numpy.ndarray
        Significance of r(q, O3) determined with moving block bootstrapping: 1
        implies significance, NaN implies insignificance, [lat, lng]
    significance_r_t2mo3_transport : numpy.ndarray
        Significance of r(T, O3) from the transport-only simulation determined 
        with moving block bootstrapping: 1 implies significance, NaN implies 
        insignificance, [lat, lng]
    significance_r_qv2mo3_transport : numpy.ndarray
        Significance of r(q, O3) from the transport-only simulation determined 
        with moving block bootstrapping: 1 implies significance, NaN implies 
        insignificance, [lat, lng]
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
    from cartopy.util import add_cyclic_point
    # For wrapping around the Prime Meridian
    do3dq, lng_gmip = add_cyclic_point(do3dq, coord=lng_gmi)
    do3dq_transport, lng_gmip = add_cyclic_point(do3dq_transport, coord=lng_gmi)    
    significance_r_qv2mo3_transport, lng_gmip = add_cyclic_point(
         significance_r_qv2mo3_transport, coord=lng_gmi)        
    significance_r_qv2mo3, lng_gmip = add_cyclic_point(significance_r_qv2mo3, 
        coord=lng_gmi)    
    do3dt2m, lng_gmip = add_cyclic_point(do3dt2m, coord=lng_gmi)
    do3dt2m_transport, lng_gmip = add_cyclic_point(do3dt2m_transport, 
        coord=lng_gmi)
    significance_r_t2mo3_transport, lng_gmip = add_cyclic_point(
        significance_r_t2mo3_transport, coord=lng_gmi)    
    significance_r_t2mo3, lng_gmip = add_cyclic_point(significance_r_t2mo3, 
        coord=lng_gmi)    
    fig = plt.figure(figsize=(16,7.5))
    clevs = np.linspace(-3, 3, 9)
    clevs_diff = np.linspace(-1, 1, 9)
    ax1 = plt.subplot2grid((3,4), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2 = plt.subplot2grid((3,4), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax3 = plt.subplot2grid((3,4), (2,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax4 = plt.subplot2grid((3,4), (0,2), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax5 = plt.subplot2grid((3,4), (1,2), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax6 = plt.subplot2grid((3,4), (2,2), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
        ax.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
        ax.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
            lat_gmi.min()+1, lat_gmi.max()-5])
        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
        lng_formatter = LongitudeFormatter()
        ax.xaxis.set_major_formatter(lng_formatter)         
        ax.get_xaxis().set_ticklabels([])
        ax.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()
        ax.yaxis.set_major_formatter(lat_formatter)
        ax.get_yaxis().set_ticklabels([])    
        # Eddy-driven jet
        skiplng = 6
        ax.errorbar(lng_gmi[::skiplng], 
            np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
            yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree(), zorder=5)
        ax.outline_patch.set_zorder(20)
    for ax in [ax1, ax2, ax3]:    
        ax.set_xticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()
        ax.yaxis.set_major_formatter(lat_formatter)      
    for ax in [ax3, ax6]:    
        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
        lng_formatter = LongitudeFormatter()
        ax.xaxis.set_major_formatter(lng_formatter)      
    cmap = plt.get_cmap('coolwarm')
    # Control dO3/dT
    ax1.set_title('(a) Control dO$_{\mathregular{3}}$/dT '+
        '[ppbv K$^{\mathregular{-1}}$]', fontsize=16, x=0.02, ha='left')    
    mb = ax1.contourf(lng_gmip, lat_gmi, do3dt2m, clevs, 
        cmap=cmap, extend='both', transform=ccrs.PlateCarree(), zorder=1)
    ax1.contourf(lng_gmip, lat_gmi, significance_r_t2mo3, hatches=['//////'], 
        colors='none', transform=ccrs.PlateCarree(), zorder=2)
    # Transport-only dO3/dT
    ax2.set_title('(c) Transport-only dO$_{\mathregular{3}}$/dT '+
        '[ppbv K$^{\mathregular{-1}}$]', fontsize=16, x=0.02, ha='left')    
    mb = ax2.contourf(lng_gmip, lat_gmi, do3dt2m_transport, clevs, 
        cmap=cmap, extend='both', transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmip, lat_gmi, significance_r_t2mo3_transport, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree())
    # Difference in dO3/dT
    ax3.set_title('(e) = (a)$\mathregular{-}$(c) '+
        '[ppbv K$^{\mathregular{-1}}$]', fontsize=16, x=0.02, ha='left')    
    mb_diff = ax3.contourf(lng_gmip, lat_gmi, (do3dt2m-do3dt2m_transport), 
        clevs_diff, cmap=cmap, extend='both', transform=ccrs.PlateCarree(), 
        zorder=1)
    # Control dO3/dq
    ax4.set_title(r'(b) Control dO$_{\mathregular{3}}$/dq '+
        '[ppbv kg g$^{\mathregular{-1}}$]', fontsize=16, x=0.02, ha='left')
    mb = ax4.contourf(lng_gmip, lat_gmi, do3dq, clevs, 
        cmap=cmap, extend='both', transform=ccrs.PlateCarree(), zorder=1)
    ax4.contourf(lng_gmip, lat_gmi, significance_r_qv2mo3, hatches=['//////'], 
        colors='none', transform=ccrs.PlateCarree())
    # Transport-only dO3/dq
    ax5.set_title(r'(d) Transport-only dO$_{\mathregular{3}}$/dq '+
        '[ppbv kg g$^{\mathregular{-1}}$]', fontsize=16, x=0.02, ha='left')
    mb = ax5.contourf(lng_gmip, lat_gmi, do3dq_transport, clevs, 
        cmap=cmap, extend='both', transform=ccrs.PlateCarree(), zorder=1)
    ax5.contourf(lng_gmip, lat_gmi, significance_r_qv2mo3_transport, 
        hatches=['//////'], colors='none', transform=ccrs.PlateCarree())
    # Difference in dO3/dT
    ax6.set_title('(f) = (b)$\mathregular{-}$(d) '+
        '[ppbv kg g$^{\mathregular{-1}}$]', fontsize=16, x=0.02, ha='left')    
    mb_diff = ax6.contourf(lng_gmip, lat_gmi, (do3dq-do3dq_transport), 
        clevs_diff, cmap=cmap, extend='both', transform=ccrs.PlateCarree(), 
        zorder=1)
    # Add colorbars
    plt.gcf().subplots_adjust(left=0.05, right=0.92, hspace=0.3)    
    colorbar_axes = plt.gcf().add_axes([ax4.get_position().x1+0.015, 
        ax5.get_position().y0, 0.015, (ax4.get_position().y1-
        ax5.get_position().y0)])
    colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
        ticks=clevs, extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar_axes = plt.gcf().add_axes([ax6.get_position().x1+0.015, 
        ax6.get_position().y0, 0.015, (ax6.get_position().y1-
        ax6.get_position().y0)])
    colorbar = plt.colorbar(mb_diff, colorbar_axes, orientation='vertical', 
        ticks=clevs_diff[::2], extend='both', extendfrac='auto')
    colorbar.ax.tick_params(labelsize=12)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'figS1.pdf', dpi=600)
    return

def figS2(lat_gmi, lng_gmi,  o3_gmi, o3_transport_gmi, lat_jet_ml): 
    """Figure S2 of Kerr et al. (2020). (a) Mean O3 from the transport-only 
    simulation. (b) Difference in mean O3 between the control simulation 
    (Figure 1a) and transport-only simulation (Figure S1a). Scatter points and 
    vertical bars are identical in (a-b) and show the mean latitude of the 
    eddy-driven jet and its variability. 

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    o3_gmi : numpy.ndarray     
        Daily afternoon surface-level O3 from the GMI CTM, units of ppbv, 
        [time, lat, lng]  
    o3_transport_gmi : numpy.ndarray     
        Daily afternoon surface-level O3 from the transport-only simulation of 
        the GMI CTM, units of ppbv, [time, lat, lng]          
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
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    from cartopy.util import add_cyclic_point
    # For wrapping around the Prime Meridian
    o3_transport_gmi, lng_gmip = add_cyclic_point(o3_transport_gmi, 
        coord=lng_gmi)
    o3_gmi, lng_gmip = add_cyclic_point(o3_gmi, coord=lng_gmi)
    fig = plt.figure(figsize=(9,5))
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=2,
        projection=ccrs.PlateCarree(central_longitude=0.))
    ax1.set_title('(a) $\mathregular{O_{3}}$', fontsize=16, x=0.02, ha='left')    
    ax2.set_title('(b) $\mathregular{\delta}$ $\mathregular{O_{3}}$', 
        fontsize=16, x=0.02, ha='left')    
    for ax in [ax1, ax2]:
        ax.coastlines(lw=0.25, resolution='50m', color='k', zorder=4)
        ax.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
            lat_gmi.min()+1, lat_gmi.max()-5])
        ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
        lng_formatter = LongitudeFormatter()
        ax.xaxis.set_major_formatter(lng_formatter)         
        ax.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
        lat_formatter = LatitudeFormatter()
        ax.yaxis.set_major_formatter(lat_formatter)
        skiplng = 6
        ax.errorbar(lng_gmi[::skiplng], np.nanmean(lat_jet_ml,axis=0)[::skiplng], 
            yerr=np.nanstd(lat_jet_ml,axis=0)[::skiplng], zorder=10, color='k', 
            markersize=3, elinewidth=1.25, ecolor='k', fmt='o', 
            transform=ccrs.PlateCarree())
    ax1.get_xaxis().set_ticklabels([])        
    cmap = plt.get_cmap('OrRd')    
    # Plot mean O3/variability transport-only simulations and (control-transport-only) 
    # percentage change
    mb1 = ax1.contourf(lng_gmip, lat_gmi, np.mean(o3_transport_gmi, axis=0), 
        np.linspace(10, 60, 11), cmap=cmap, extend='both',
        transform=ccrs.PlateCarree(), zorder=1)
    csthick = ax1.contour(lng_gmip, lat_gmi, np.nanstd(o3_transport_gmi, axis=0), [10.], 
        colors='k', linewidths=1.5, transform=ccrs.PlateCarree(), zorder=15)
    csmedium = ax1.contour(lng_gmip, lat_gmi, np.nanstd(o3_transport_gmi, axis=0), [8.], 
        colors='k', linestyles='--', linewidths=0.75, 
        transform=ccrs.PlateCarree(), zorder=15)
    # Negative values imply decreases in transport-only simulation!
    pc = np.nanmean((o3_transport_gmi-o3_gmi), axis=0)
    mb2 = ax2.contourf(lng_gmip, lat_gmi, pc, np.linspace(-2, 2, 9), 
        cmap=plt.get_cmap('coolwarm'), extend='both', transform=ccrs.PlateCarree(), 
        zorder=1)
    # Add colorbars
    plt.gcf().subplots_adjust(left=0.05, right=0.86, hspace=0.3)    
    colorbar_axes = plt.gcf().add_axes([ax1.get_position().x1+0.03, 
        ax1.get_position().y0, 0.02, (ax1.get_position().y1-
        ax1.get_position().y0)])
    colorbar = plt.colorbar(mb1, colorbar_axes, orientation='vertical', 
        ticks=np.linspace(10, 60, 6), extend='both')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[ppbv]', fontsize=16, labelpad=16)
    ax1.outline_patch.set_zorder(20)
    colorbar_axes = plt.gcf().add_axes([ax2.get_position().x1+0.03, 
        ax2.get_position().y0, 0.02, (ax2.get_position().y1-
        ax2.get_position().y0)])
    colorbar = plt.colorbar(mb2, colorbar_axes, orientation='vertical', extend='max')
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label('[ppbv]', labelpad=8, fontsize=16)
    ax2.outline_patch.set_zorder(20)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'figS2.pdf', dpi=600)    
    return


def figS3(lat_gmi, lng_gmi, r_o3jetdist, r_t2mjetdist, r_qv2mjetdist, 
    lat_jet_ml, significance_r_o3jetdist, significance_r_t2mjetdist, 
    significance_r_qv2mjetdist): 
    """Figure S3 of Kerr et al. (2020). (a) r(O3, jet distance) is shown with 
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
    from cartopy.util import add_cyclic_point
    # For wrapping around the Prime Meridian
    r_o3jetdist, lng_gmip = add_cyclic_point(r_o3jetdist, 
        coord=lng_gmi)
    r_t2mjetdist, lng_gmip = add_cyclic_point(r_t2mjetdist, 
        coord=lng_gmi)    
    r_qv2mjetdist, lng_gmip = add_cyclic_point(r_qv2mjetdist, 
        coord=lng_gmi)        
    significance_r_o3jetdist, lng_gmip = add_cyclic_point(
        significance_r_o3jetdist, coord=lng_gmi)    
    significance_r_t2mjetdist, lng_gmip = add_cyclic_point(
        significance_r_t2mjetdist, coord=lng_gmi) 
    significance_r_qv2mjetdist, lng_gmip = add_cyclic_point(
        significance_r_qv2mjetdist, coord=lng_gmi)        
    # Load ocean shapefiles
    ocean50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
        edgecolor=None, facecolor='lightgrey')
    fig = plt.figure(figsize=(9,7.7))
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
    mb = ax1.contourf(lng_gmip, lat_gmi, r_o3jetdist, np.linspace(-1,1,9), 
        cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(O3, jet lat - lat)
    ax1.contourf(lng_gmip, lat_gmi, significance_r_o3jetdist, 
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
    mb = ax2.contourf(lng_gmip, lat_gmi, r_o3jetdist,
        np.linspace(-1,1,9), cmap=cmap, extend='neither', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmip, lat_gmi, significance_r_t2mjetdist, 
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
    mb = ax3.contourf(lng_gmip, lat_gmi, r_qv2mjetdist,
        np.linspace(-1,1,9), cmap=cmap, extend='neither', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax3.contourf(lng_gmip, lat_gmi, significance_r_qv2mjetdist, 
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

def figS4(lat_gmi, lng_gmi, pblh_merra, U10M, WIND10M, lat_jet_ml, times_gmi, 
    significance_r_pblhjetdist, significance_r_U10Mjetdist, 
    significance_r_WIND10Mjetdist):
    """Figure S4 of Kerr et al. (2020). (a) Difference in PBLH on days with a
    poleward jet versus days with an equatorward jet. Hatched grid cells 
    correspond to insignificant values of r(PBLH, jet latitude) determined with 
    moving block bootstrapping. (b) Same as (a) but for U10. (c) Same as (a) 
    but for WIND10. Scatterpoints and error bars in (a-c) specify the mean 
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
    significance_r_WIND10Mjetdist : numpy.ndarray
        Significance of r(WIND10, jet latitude) determined with moving block 
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
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter    
    # Composites of flow at surface on days with poleward or equatorward jet 
    eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_pblh, eqjet_pblh = \
        globalo3_calculate.segregate_field_bylat(pblh_merra, lng_gmi, lat_jet_ml, 
        times_gmi)
    eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_U10M, eqjet_U10M = \
        globalo3_calculate.segregate_field_bylat(U10M, lng_gmi, lat_jet_ml, 
        times_gmi)
    eqjet_lat, eqjet_lat_var, pwjet_lat, pwjet_lat_var, pwjet_WIND10M, eqjet_WIND10M = \
        globalo3_calculate.segregate_field_bylat(WIND10M, lng_gmi, lat_jet_ml, 
        times_gmi)
    from cartopy.util import add_cyclic_point
    # For wrapping around the Prime Meridian
    pwjet_pblh, lng_gmip = add_cyclic_point(pwjet_pblh, coord=lng_gmi)
    eqjet_pblh, lng_gmip = add_cyclic_point(eqjet_pblh, coord=lng_gmi)
    pwjet_U10M, lng_gmip = add_cyclic_point(pwjet_U10M, coord=lng_gmi)
    eqjet_U10M, lng_gmip = add_cyclic_point(eqjet_U10M, coord=lng_gmi)
    pwjet_WIND10M, lng_gmip = add_cyclic_point(pwjet_WIND10M, coord=lng_gmi)
    eqjet_WIND10M, lng_gmip = add_cyclic_point(eqjet_WIND10M, coord=lng_gmi)    
    significance_r_pblhjetdist, lng_gmip = add_cyclic_point(
        significance_r_pblhjetdist, coord=lng_gmi)      
    significance_r_U10Mjetdist, lng_gmip = add_cyclic_point(
        significance_r_U10Mjetdist, coord=lng_gmi)      
    significance_r_WIND10Mjetdist, lng_gmip = add_cyclic_point(
        significance_r_WIND10Mjetdist, coord=lng_gmi)          
    fig = plt.figure(figsize=(9,7.7))
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
    mb = ax1.contourf(lng_gmip, lat_gmi, (pwjet_pblh-eqjet_pblh), 
        np.linspace(-300,300,7), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(PBLH, jetlat)
    ax1.contourf(lng_gmip, lat_gmi, significance_r_pblhjetdist, 
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
    mb = ax2.contourf(lng_gmip, lat_gmi, (pwjet_U10M-eqjet_U10M), 
        np.linspace(-4, 4, 9), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmip, lat_gmi, significance_r_U10Mjetdist, 
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
    ax3.set_title('(c) $\overline{\mathregular{U}_\mathregular{10,\:PW}}$ $-$'+
        ' $\overline{\mathregular{U}_\mathregular{10,\:EW}}$', 
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
    mb = ax3.contourf(lng_gmip, lat_gmi, (pwjet_WIND10M-eqjet_WIND10M), 
        np.linspace(-3,3,7), cmap=cmap, extend='both', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax3.contourf(lng_gmip, lat_gmi, significance_r_WIND10Mjetdist, 
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
        'figS4.pdf', dpi=600)
    return

def figS5(lat_gmi, lng_gmi, r_pblhjetdist, r_U10Mjetdist, r_WIND10Mjetdist, 
    lat_jet_ml, significance_r_pblhjetdist, significance_r_U10Mjetdist, 
    significance_r_WIND10Mjetdist): 
    """Figure S5 of Kerr et al. (2020). (a) r(PBLH, jet distance) is shown with 
    colored shading. Hatching indicates insignificance determined with moving 
    block bootstrapping. (b) Same as (a) but showing r(U10, jet distance). (c)
    Same as (a) but showing r(WIND10, jet distance). Scatterpoints and vertical 
    bars in (a-c) specify the mean position and variability of the eddy-driven 
    jet, respectively.

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
    r_WIND10Mjetdist : numpy.ndarray     
        Pearson correlation coefficient between 10-meter total wind and jet
        distance, [lat, lng]
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
    significance_r_WIND10Mjetdist : numpy.ndarray
        Significance of r(WIND10, jet latitude) determined with moving block 
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
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    from cartopy.util import add_cyclic_point
    # For wrapping around the Prime Meridian
    r_pblhjetdist, lng_gmip = add_cyclic_point(r_pblhjetdist, coord=lng_gmi)
    r_U10Mjetdist, lng_gmip = add_cyclic_point(r_U10Mjetdist, coord=lng_gmi)
    r_WIND10Mjetdist, lng_gmip = add_cyclic_point(r_WIND10Mjetdist, 
        coord=lng_gmi)
    significance_r_pblhjetdist, lng_gmip = add_cyclic_point(
        significance_r_pblhjetdist, coord=lng_gmi)
    significance_r_U10Mjetdist, lng_gmip = add_cyclic_point(
        significance_r_U10Mjetdist, coord=lng_gmi)
    significance_r_WIND10Mjetdist, lng_gmip = add_cyclic_point(
        significance_r_WIND10Mjetdist, coord=lng_gmi)
    fig = plt.figure(figsize=(9,7.7))
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
    mb = ax1.contourf(lng_gmip, lat_gmi, r_pblhjetdist, np.linspace(-1,1,9), 
        cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
    # Hatching for significance of r(PBLH, jetlat)
    ax1.contourf(lng_gmip, lat_gmi, significance_r_pblhjetdist, 
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
    mb = ax2.contourf(lng_gmip, lat_gmi, r_U10Mjetdist,
        np.linspace(-1, 1, 9), cmap=cmap, extend='neither', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax2.contourf(lng_gmip, lat_gmi, significance_r_U10Mjetdist, 
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
    # r(WIND10M, jet lat - lat)
    ax3 = plt.subplot2grid((3,2), (2,0), colspan=2,
            projection=ccrs.PlateCarree(central_longitude=0.))
    ax3.set_title(r'(c) r($\overline{\mathregular{U}_\mathregular{10}}$, '+
        '$\mathregular{\phi_{jet}}$ $-$ ${\mathregular{\phi}}$)', fontsize=16, 
        x=0.02, ha='left')    
    # ax3.add_feature(ocean50m, zorder=3)
    ax3.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
    ax3.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
        lat_gmi.min()+1, lat_gmi.max()-5])
    mb = ax3.contourf(lng_gmip, lat_gmi, r_WIND10Mjetdist,
        np.linspace(-1, 1, 9), cmap=cmap, extend='neither', 
        transform=ccrs.PlateCarree(), zorder=1)
    ax3.contourf(lng_gmip, lat_gmi, significance_r_WIND10Mjetdist, 
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
        'figS5.pdf', dpi=600)
    return

def figS6(lat_gmi, lng_gmi, o3_gmi, V10M, lat_jet_ml): 
    """Regionally-averaged O3 from the control simulation versus 
    regionally-averaged V10M. Regional averaging is conducted for the regions 
    listed in each subplots' title and shown in Figure 3 but only within 
    +/- 5 degrees of the mean position of the jet. Red dashed lines represent 
    the OLS regression, and the inset text indicates its slope. 

    Parameters
    ----------
    lat_gmi : numpy.ndarray
        GMI CTM latitude coordinates, units of degrees north, [lat,]
    lng_gmi : numpy.ndarray
        GMI CTM longitude coordinates, units of degrees east, [lng,]
    o3_gmi : numpy.ndarray     
        Daily afternoon surface-level O3 from the GMI CTM, units of ppbv, 
        [time, lat, lng]  
    V10M : numpy.ndarray 
        Daily mean 10-meter meridional wind (V10), units of m s-1, [time, lat, 
        lng]        
    lat_jet_ml : numpy.ndarray
        The latitude of the jet, identifed by maximum zonal (U) wind at 500 hPa
        in the Northern Hemisphere mid-latitudes, units of degrees north, 
        [time, lng]
    Returns
    -------
    None    
    """
    import sys
    sys.path.append('/Users/ghkerr/phd/GMI/')
    from geo_idx import geo_idx       
    fig = plt.figure(figsize=(8,8))
    ax1 = plt.subplot2grid((2,2),(0,0))
    ax2 = plt.subplot2grid((2,2),(0,1))
    ax3 = plt.subplot2grid((2,2),(1,0))
    ax4 = plt.subplot2grid((2,2),(1,1))
    # Eastern North America
    ena_left = geo_idx(260., lng_gmi)
    ena_right = geo_idx(295., lng_gmi)
    lat_jet_ena = lat_jet_ml[:, ena_left:ena_right+1]
    o3_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
        o3_gmi, lat_gmi, lng_gmi, 260., 295., lat_jet_ena.mean()-5, 
        lat_jet_ena.mean()+5)
    V10M_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
        V10M, lat_gmi, lng_gmi, 260., 295., lat_jet_ena.mean()-5, 
        lat_jet_ena.mean()+5)
    o3_ena = np.nanmean(o3_ena, axis=tuple((1,2)))
    V10M_ena = np.nanmean(V10M_ena, axis=tuple((1,2)))
    ax2.plot(V10M_ena, o3_ena, 'ko', markersize=3)
    trend = np.polyfit(V10M_ena, o3_ena, 1)
    trendpoly = np.poly1d(trend) 
    ax2.plot(np.sort(V10M_ena), np.sort(trendpoly(V10M_ena)), lw=2, 
        color='#e41a1c', ls='--', zorder=10)
    ax2.set_title('(b) Eastern North America', fontsize=16, loc='left')
    ax2.text(0.55, 0.03, '%.2f ppbv s m$^{-1}$'%trendpoly[1], ha='left',
        transform=ax2.transAxes, fontsize=12, zorder=20)
    # Western North America
    wna_left = geo_idx(235., lng_gmi)
    wna_right = geo_idx(260., lng_gmi)
    lat_jet_wna = lat_jet_ml[:, wna_left:wna_right+1]
    o3_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(
        o3_gmi, lat_gmi, lng_gmi, 235., 260., lat_jet_wna.mean()-5, 
        lat_jet_wna.mean()-5)
    V10M_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(
        V10M, lat_gmi, lng_gmi, 235., 260., lat_jet_wna.mean()-5, 
        lat_jet_wna.mean()-5)
    o3_wna = np.nanmean(o3_wna, axis=tuple((1,2)))
    V10M_wna = np.nanmean(V10M_wna, axis=tuple((1,2)))
    ax1.plot(V10M_wna, o3_wna, 'ko', markersize=3)
    trend = np.polyfit(V10M_wna, o3_wna, 1)
    trendpoly = np.poly1d(trend) 
    ax1.plot(np.sort(V10M_wna), np.sort(trendpoly(V10M_wna)), lw=2, 
        color='#e41a1c', ls='--', zorder=10)
    ax1.set_title('(a) Western North America', fontsize=16, loc='left')
    ax1.text(0.55, 0.03, '%.2f ppbv s m$^{-1}$'%trendpoly[1], ha='left',
        transform=ax1.transAxes, fontsize=12, zorder=20)
    # Europe (since Europe crosses the prime meridian, finding grid cells over
    # 0 deg longitude will not work, so find European domain in two steps)
    eu1_left = geo_idx(350., lng_gmi)
    eu1_right = geo_idx(360., lng_gmi)
    lat_jet_eu1 = lat_jet_ml[:, eu1_left:eu1_right+1]
    eu2_left = geo_idx(0., lng_gmi)
    eu2_right = geo_idx(30., lng_gmi)
    lat_jet_eu2 = lat_jet_ml[:, eu2_left:eu2_right+1]
    lat_jet_eu = np.hstack([lat_jet_eu1, lat_jet_eu2])
    o3_eu1, lat_eu1, lng_eu1 = globalo3_calculate.find_grid_in_bb(
        o3_gmi, lat_gmi, lng_gmi, 350., 360., lat_jet_eu.mean()-5, 
        lat_jet_eu.mean()+5)
    o3_eu2, lat_eu2, lng_eu2 = globalo3_calculate.find_grid_in_bb(
        o3_gmi, lat_gmi, lng_gmi, 0., 30., lat_jet_eu.mean()-5, 
        lat_jet_eu.mean()+5) 
    o3_eu = np.dstack([o3_eu1, o3_eu2])
    V10M_eu1, lat_eu1, lng_eu1 = globalo3_calculate.find_grid_in_bb(
        V10M, lat_gmi, lng_gmi, 350., 360., lat_jet_eu.mean()-5, 
        lat_jet_eu.mean()+5)
    V10M_eu2, lat_eu2, lng_eu2 = globalo3_calculate.find_grid_in_bb(
        V10M, lat_gmi, lng_gmi, 0., 30., lat_jet_eu.mean()-5, 
        lat_jet_eu.mean()+5) 
    V10M_eu = np.dstack([V10M_eu1, V10M_eu2])
    # Add back in _eu1? (mostly ocean, though)
    o3_eu = np.nanmean(o3_eu2, axis=tuple((1,2)))
    V10M_eu = np.nanmean(V10M_eu2, axis=tuple((1,2)))
    ax3.plot(V10M_eu, o3_eu, 'ko', markersize=3)
    trend = np.polyfit(V10M_eu, o3_eu, 1)
    trendpoly = np.poly1d(trend) 
    ax3.plot(np.sort(V10M_eu), np.sort(trendpoly(V10M_eu)), lw=2, 
        color='#e41a1c', ls='--', zorder=10)
    ax3.set_title('(c) Europe', fontsize=16, loc='left')
    ax3.text(0.55, 0.03, '%.2f ppbv s m$^{-1}$'%trendpoly[1], ha='left',
        transform=ax3.transAxes, fontsize=12, zorder=20)
    # Asia
    asia_left = geo_idx(90., lng_gmi)
    asia_right = geo_idx(125., lng_gmi)
    lat_jet_asia = lat_jet_ml[:, asia_left:asia_right+1]
    o3_asia, lat_asia, lng_asia = globalo3_calculate.find_grid_in_bb(
        o3_gmi, lat_gmi, lng_gmi, 90., 125., lat_jet_asia.mean()-5, 
        lat_jet_asia.mean()+5) 
    V10M_asia, lat_asia, lng_asia = globalo3_calculate.find_grid_in_bb(
        V10M, lat_gmi, lng_gmi, 90., 125., lat_jet_asia.mean()-5, 
        lat_jet_asia.mean()+5) 
    o3_asia = np.nanmean(o3_asia, axis=tuple((1,2)))
    V10M_asia = np.nanmean(V10M_asia, axis=tuple((1,2)))
    ax4.plot(V10M_asia, o3_asia, 'ko', markersize=3)
    trend = np.polyfit(V10M_asia, o3_asia, 1)
    trendpoly = np.poly1d(trend) 
    ax4.plot(np.sort(V10M_asia), np.sort(trendpoly(V10M_asia)), lw=2, 
        color='#e41a1c', ls='--', zorder=10)
    ax4.set_title('(d) China', fontsize=16, loc='left')
    ax4.text(0.55, 0.03, '%.2f ppbv s m$^{-1}$'%trendpoly[1], ha='left',
        transform=ax4.transAxes, fontsize=12, zorder=20)
    # Aesthetics 
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_xlim([-3, 3])
        ax.set_ylim([22, 50])
    for ax in [ax1, ax3]:
        ax.set_ylabel('O$_{3}$ [ppbv]', fontsize=16)
    for ax in [ax3, ax4]:
        ax.set_xlabel('V$_{10}$ [m s$^{-1}$]', fontsize=16)
    plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
        'figS6.pdf', dpi=600)
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
    lng_gmi[-1] = 358.75
    mtime = nc.Dataset(datapath+'gmi_O3control_JJA2008-2010.nc')['time'][:].data
    o3_transport_gmi = nc.Dataset(datapath+'gmi_O3transportonly_JJA2008-2010.nc')['O3_transportonly'][:].data
    t2m_merra = nc.Dataset(datapath+'merra2_t2m_JJA2008-2010.nc')['T2M'][:].data
    t2m_merra_china = nc.Dataset(datapath+'merra2_t2m_JJA2016-2017.nc')['T2M'][:].data
    qv2m_merra = nc.Dataset(datapath+'merra2_qv2m_JJA2008-2010.nc')['QV2M'][:].data
    qv2m_merra_china = nc.Dataset(datapath+'merra2_qv2m_JJA2016-2017.nc')['QV2M'][:].data    
    lat_jet_ml = nc.Dataset(datapath+'merra2_JET_JJA2008-2010.nc')['jetlatitude'][:].data
    U10M = nc.Dataset(datapath+'merra2_U10M_JJA2008-2010.nc')['U10M'][:].data
    V10M = nc.Dataset(datapath+'merra2_V10M_JJA2008-2010.nc')['V10M'][:].data
    WIND10M = np.hypot(U10M, V10M)
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
    significance_r_pblhjetdist = nc.Dataset(datapath+'sig_merra2_jetdist_merra2_PBLH.nc')['sig_JETDIST_merra2_PBLH'][:].data 
    significance_r_U10Mjetdist = nc.Dataset(datapath+'sig_merra2_jetdist_merra2_U10M.nc')['sig_JETDIST_merra2_U10M'][:].data 
    significance_r_V10Mjetdist = nc.Dataset(datapath+'sig_merra2_jetdist_merra2_V10M.nc')['sig_JETDIST_merra2_V10M'][:].data 
    significance_r_WIND10Mjetdist = nc.Dataset(datapath+'sig_merra2_jetdist_merra2_WIND10M.nc')['sig_JETDIST_merra2_WIND10M'][:].data 
    significance_r_V10Mo3 = nc.Dataset(datapath+'sig_merra2_V10M_gmi_O3control.nc')['sig_V10M_O3_control'][:].data 
    # Calculate dO3/dT, dO3/dq, r(T, O3), and r(q, O3) from model
    do3dt2m = globalo3_calculate.calculate_do3dt(t2m_merra, o3_gmi, lat_gmi, 
        lng_gmi)
    do3dq = globalo3_calculate.calculate_do3dt(qv2m_merra, o3_gmi, lat_gmi, 
        lng_gmi)
    do3dt2m_transport = globalo3_calculate.calculate_do3dt(t2m_merra, 
        o3_transport_gmi, lat_gmi, lng_gmi)
    do3dq_transport = globalo3_calculate.calculate_do3dt(qv2m_merra, 
        o3_transport_gmi, lat_gmi, lng_gmi)    
    r_t2mo3 = globalo3_calculate.calculate_r(t2m_merra, o3_gmi, lat_gmi, 
        lng_gmi)
    r_t2mo3_transport = globalo3_calculate.calculate_r(t2m_merra, 
        o3_transport_gmi, lat_gmi, lng_gmi)
    r_qv2mo3 = globalo3_calculate.calculate_r(qv2m_merra, o3_gmi, lat_gmi, 
        lng_gmi)
    r_qv2mo3_transport = globalo3_calculate.calculate_r(qv2m_merra, 
        o3_transport_gmi, lat_gmi, lng_gmi) 
    r_V10Mo3 = globalo3_calculate.calculate_r(V10M, o3_gmi, lat_gmi, lng_gmi) 
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
    m_WIND10Mjetdist, r_WIND10Mjetdist, diff_WIND10Mjetdist = \
        globalo3_calculate.calculate_fieldjet_relationship(np.hypot(U10M, 
        V10M), lat_gmi, lng_gmi, lat_jet_ml, lng_gmi)
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
                        
# # FIGURE 1: Mean O3 and NOx
# fig1(lat_gmi, lng_gmi, o3_gmi, lat_jet_ml)
# # FIGURE 2: model performance
# fig2(lat_gmi, lng_gmi, r_aqs, r_naps, r_emep, r_china)
# # FIGURE 3: r(T, O3) and r(q, O3)
# fig3(lat_gmi, lng_gmi, r_t2mo3, r_qv2mo3, significance_r_t2mo3, 
#     significance_r_qv2mo3, lat_jet_ml)
# # FIGURE 4: zonally-averaged O3-meteorology relationships 
# fig4(lat_gmi, lng_gmi, r_t2mo3, r_qv2mo3, r_t2mo3_aqs, r_qv2mo3_aqs, 
#     lat_aqs, lng_aqs, r_t2mo3_naps, r_qv2mo3_naps, lat_naps, lng_naps, 
#     r_t2mo3_emep, r_qv2mo3_emep, lat_emep, lng_emep, r_t2mo3_china, 
#     r_qv2mo3_china, lat_china, lng_china, lng_gmi, lat_jet_ml)
# # FIGURE 5: Percentage difference in r(T, O3) and r(q, O3) between the 
# # transport-only and control simulation
# fig5(lat_gmi, lng_gmi, r_t2mo3, r_t2mo3_transport, r_qv2mo3, 
#     r_qv2mo3_transport, significance_diff_r_t2mo3, significance_diff_r_qv2mo3, 
#     lat_jet_ml)
# # FIGURE 6: difference in O3, T2M, and qv2M on days with a poleward versus
# # equatorward jet
# fig6(lat_gmi, lng_gmi, o3_gmi, t2m_merra, qv2m_merra, lat_jet_ml, 
#     times_gmi, significance_r_o3jetdist, significance_r_t2mjetdist, 
#     significance_r_qv2mjetdist)
# # FIGURE 7: cyclone frequency and poleward-equatorward jet differences
# fig7(lat_cyclones_binned, lng_cyclones_binned, cyclones_binned,
#     pwjet_cyclones_binned, eqjet_cyclones_binned, lat_gmi, lng_gmi, 
#     lat_jet_ml)
# # FIGURE 8: O3 anomaly at cyclone
# fig8(o3_anom_rotated)
# # FIGURE 9: difference in V10 on days with a poleward versus equatorward jet, 
# # r(O3, V10M) and O3 latitudinal gradient. 
# fig9(lat_gmi, lng_gmi, o3_gmi, V10M, lat_jet_ml, times_gmi, 
#     significance_r_V10Mjetdist, r_V10Mo3, significance_r_V10Mo3)
# # FIGURE S1: Mean O3 from transport-only simulation and difference in O3 
# # from control and transport-only simulations
# figS1(lat_gmi, lng_gmi, do3dt2m, do3dt2m_transport, do3dq, do3dq_transport,
#     significance_r_t2mo3, significance_r_qv2mo3, 
#     significance_r_t2mo3_transport, significance_r_qv2mo3_transport,
#     lat_jet_ml)
# # FIGURE S2; dO3/dT and dO3/dq from control and transport-only simulations
# # and differences
# figS2(lat_gmi, lng_gmi,  o3_gmi, o3_transport_gmi, lat_jet_ml)
# # FIGURE S3: r(O3, jet distance), r(T, jet distance), and r(q, jet distance)
# figS3(lat_gmi, lng_gmi, r_o3jetdist, r_t2mjetdist, r_qv2mjetdist, 
#     lat_jet_ml, significance_r_o3jetdist, significance_r_t2mjetdist, 
#     significance_r_qv2mjetdist)
# # FIGURE S4: difference in PBLH, U10M, and WIND10M on days with a poleward 
# # versus equatorward jet
# figS4(lat_gmi, lng_gmi, pblh_merra, U10M, WIND10M, lat_jet_ml, times_gmi, 
#     significance_r_pblhjetdist, significance_r_U10Mjetdist, 
#     significance_r_WIND10Mjetdist)
# # FIGURE S5: r(PBLH, jet distance), r(U10, jet distance), and r(WIND10, jet
# # distance)
# figS5(lat_gmi, lng_gmi, r_pblhjetdist, r_U10Mjetdist, r_WIND10Mjetdist, 
#     lat_jet_ml, significance_r_pblhjetdist, significance_r_U10Mjetdist, 
#     significance_r_WIND10Mjetdist)
# # FIGURE S6: Scatterplots of O3 versus V10M
# figS6(lat_gmi, lng_gmi, o3_gmi, V10M, lat_jet_ml)

# # # # # Reviewer comments 
# import numpy as np
# import matplotlib as mpl
# mpl.rcParams['hatch.linewidth']=0.3     
# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import matplotlib.patches as mpatches
# from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter   
# import sys
# sys.path.append('/Users/ghkerr/phd/globalo3/')
# import globalo3_open, globalo3_calculate, observations_open
# sys.path.append('/Users/ghkerr/phd/tracer/')
# import tracer_open, tracer_calculate
# months = [6, 7, 8]
# months_str = ['jun', 'jul', 'aug']
# # Hemisphere definition
# latmin, lngmin, latmax, lngmax = -1., 0., 90., 360.
# # Load Northern Hemisphere HindcastMR2 GMI CTM O3 
# lat_gmi, lng_gmi, times_gmi_full, o3_gmi_full = \
#     globalo3_open.open_overpass2_specifieddomain(np.arange(2000, 2011, 1), 
#     months_str, latmin, latmax, lngmin, lngmax, 'O3', 'HindcastMR2')
# o3_gmi_full = o3_gmi_full*1e9
# # Load MERRA-2 Northern Hemisphere 2-meter temperatures, specific humidity, 
# # jet
# t_merra_full, lat_merra, lng_merra, lev_merra = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(
#     np.arange(2000, 2011, 1), ['jun', 'jul', 'aug'], 'T', 0., 90., 360., -1., 
#     985., 1005.)
# qv_merra_full, lat_merra, lng_merra, lev_merra = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(
#     np.arange(2000, 2011, 1), ['jun', 'jul', 'aug'], 'QV', 0., 90., 360., -1., 
#     985., 1005.)
# edj_jja, lat_edj, lng_edj, lev_edj = \
#     tracer_open.open_merra2_inst3_3d_asm_Nv_specifieddomain(
#     np.arange(2000, 2011, 1), ['jun', 'jul', 'aug'], 'U', 0., 90., 360., -1., 
#     487., 526., operation='mean')
# # Degrade to resolution of GEOSChem    
# edj_jja = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, lng_gmi, 
#     lat_edj, lng_edj, edj_jja)
# # Subset fields in mid-latitudes
# edj_jja, lat_edj, lng_edj = globalo3_calculate.find_grid_in_bb(edj_jja, lat_gmi, 
#     lng_gmi, 0., 360., 20., 70.)
# # Determine jet latitude
# edj_jja = tracer_calculate.find_jetlatitude(edj_jja, lat_edj, lng_edj)
# # Column average
# t_merra_full = np.nanmean(t_merra_full, axis=1)
# qv_merra_full = np.nanmean(qv_merra_full, axis=1)
# # Interpolate MERRA-2 to GMI resolution 
# t_merra_full = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
#     lng_gmi, lat_merra, lng_merra, t_merra_full)
# qv_merra_full = globalo3_open.interpolate_merra_to_ctmresolution(lat_gmi, 
#     lng_gmi, lat_merra, lng_merra, qv_merra_full)
# qv_merra_full = qv_merra_full*1000. # Convert from kg kg-1 to g kg-1
# years_thisyr = [x.year for x in times_gmi_full]
# # Separate by years             
# years_thisyr = np.array(years_thisyr)
# jet_wna, jet_ena, jet_eu, jet_asia = [], [], [], []
# r_to3_wna, r_qvo3_wna = [], []
# r_to3_ena, r_qvo3_ena = [], []
# r_to3_eu, r_qvo3_eu = [], []
# r_to3_asia, r_qvo3_asia = [], []
# for year in np.arange(2000, 2011, 1):#, 2011, 1):
#     thisyr = np.where(years_thisyr==year)[0]
#     jet_thisyr = edj_jja[thisyr]
#     o3_gmi_thisyr = o3_gmi_full[thisyr]
#     t_merra_thisyr = t_merra_full[thisyr]
#     qv_merra_thisyr = qv_merra_full[thisyr]
#     # Find regionally-averaged temperature, humidity, and O3 over Western 
#     # North America 
#     jet_wna = jet_thisyr[:, np.where(lng_gmi==250.)[0][0]:
#         np.where(lng_gmi==260)[0][0]+1]
#     jet_wna = np.nanmean(jet_wna)
#     o3_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(
#         o3_gmi_thisyr[1:], lat_gmi, lng_gmi, 250., 260., jet_wna-5, jet_wna+5)
#     t_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(
#         t_merra_thisyr[1:], lat_gmi, lng_gmi, 250., 260., jet_wna-5, jet_wna+5)
#     qv_wna, lat_wna, lng_wna = globalo3_calculate.find_grid_in_bb(
#         qv_merra_thisyr[1:], lat_gmi, lng_gmi, 250., 260, jet_wna-5, jet_wna+5)
#     land_wna = globalo3_calculate.find_grid_overland(lat_wna, lng_wna)
#     o3_wna = np.nanmean(o3_wna*land_wna, axis=tuple((1,2)))
#     t_wna = np.nanmean(t_wna*land_wna, axis=tuple((1,2)))
#     qv_wna = np.nanmean(qv_wna*land_wna, axis=tuple((1,2)))
#     # Calculate r(O3, T) and r(O3, q)
#     r_to3_wna.append(np.corrcoef(t_wna, o3_wna)[0,1])
#     r_qvo3_wna.append(np.corrcoef(qv_wna, o3_wna)[0,1])
#     # Eastern North America
#     jet_ena = jet_thisyr[:, np.where(lng_gmi==270.)[0][0]:
#         np.where(lng_gmi==280.)[0][0]+1]
#     jet_ena = np.nanmean(jet_ena)
#     o3_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
#         o3_gmi_thisyr[1:], lat_gmi, lng_gmi, 270., 280., jet_ena-5, jet_wna+5)
#     t_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
#         t_merra_thisyr[1:], lat_gmi, lng_gmi, 270., 280., jet_ena-5, jet_wna+5)
#     qv_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
#         qv_merra_thisyr[1:], lat_gmi, lng_gmi, 270., 280., jet_ena-5, jet_wna+5)
#     land_ena = globalo3_calculate.find_grid_overland(lat_ena, lng_ena)
#     o3_ena = np.nanmean(o3_ena*land_ena, axis=tuple((1,2)))
#     t_ena = np.nanmean(t_ena*land_ena, axis=tuple((1,2)))
#     qv_ena = np.nanmean(qv_ena*land_ena, axis=tuple((1,2)))
#     r_to3_ena.append(np.corrcoef(t_ena, o3_ena)[0,1])
#     r_qvo3_ena.append(np.corrcoef(qv_ena, o3_ena)[0,1])
#     # Europe (since Europe crosses the prime meridian, finding grid cells over
#     # 0 deg longitude will not work, so find European domain in two steps)
#     # Model 
#     jet_eu = jet_thisyr[:, np.where(lng_gmi==5.)[0][0]:
#         np.where(lng_gmi==15.)[0][0]+1]
#     jet_eu = np.nanmean(jet_eu)    
#     o3_eu, lat_eu, lng_eu = globalo3_calculate.find_grid_in_bb(
#         o3_gmi_thisyr[1:], lat_gmi, lng_gmi, 5., 15., jet_eu-5, jet_eu+5)
#     t_eu, lat_eu, lng_eu = globalo3_calculate.find_grid_in_bb(
#         t_merra_thisyr[1:], lat_gmi, lng_gmi, 5., 15., jet_eu-5, jet_eu+5)
#     qv_eu, lat_eu, lng_eu = globalo3_calculate.find_grid_in_bb(
#         qv_merra_thisyr[1:], lat_gmi, lng_gmi, 5., 15., jet_eu-5, jet_eu+5)
#     land_eu = globalo3_calculate.find_grid_overland(lat_eu, lng_eu)
#     o3_eu = np.nanmean(o3_eu*land_eu, axis=tuple((1,2)))
#     t_eu = np.nanmean(t_eu*land_eu, axis=tuple((1,2)))
#     qv_eu = np.nanmean(qv_eu*land_eu, axis=tuple((1,2)))
#     r_to3_eu.append(np.corrcoef(t_eu, o3_eu)[0,1])
#     r_qvo3_eu.append(np.corrcoef(qv_eu, o3_eu)[0,1])
#     # China
#     jet_asia = jet_thisyr[:, np.where(lng_gmi==110.)[0][0]:
#         np.where(lng_gmi==120.)[0][0]+1]
#     jet_asia = np.nanmean(jet_asia)        
#     o3_asia, lat_asia, lng_asia = globalo3_calculate.find_grid_in_bb(
#         o3_gmi_thisyr[1:], lat_gmi, lng_gmi, 110., 120., jet_asia-5, jet_asia+5)
#     t_asia, lat_asia, lng_asia = globalo3_calculate.find_grid_in_bb(
#         t_merra_thisyr[1:], lat_gmi, lng_gmi, 110., 120., jet_asia-5, jet_asia+5)
#     qv_asia, lat_asia, lng_asia = globalo3_calculate.find_grid_in_bb(
#         qv_merra_thisyr[1:], lat_gmi, lng_gmi, 110., 120., jet_asia-5, jet_asia+5)
#     land_asia = globalo3_calculate.find_grid_overland(lat_asia, lng_asia)
#     o3_asia = np.nanmean(o3_asia*land_asia, axis=tuple((1,2)))
#     t_asia = np.nanmean(t_asia*land_asia, axis=tuple((1,2)))
#     qv_asia = np.nanmean(qv_asia*land_asia, axis=tuple((1,2)))    
#     r_to3_asia.append(np.corrcoef(t_asia, o3_asia)[0,1])
#     r_qvo3_asia.append(np.corrcoef(qv_asia, o3_asia)[0,1])
#     # Plot yearly maps of r(O3, T) and r(O3, q)
#     r_to3_thisyr = globalo3_calculate.calculate_r(t_merra_thisyr, 
#         o3_gmi_thisyr, lat_gmi, lng_gmi)
#     r_qvo3_thisyr = globalo3_calculate.calculate_r(qv_merra_thisyr, 
#         o3_gmi_thisyr, lat_gmi, lng_gmi)
#     significance_r_to3_thisyr = \
#         globalo3_calculate.calculate_r_significance(t_merra_thisyr, 
#         o3_gmi_thisyr, r_to3_thisyr, lat_gmi, lng_gmi)
#     significance_r_qvo3_thisyr = \
#         globalo3_calculate.calculate_r_significance(qv_merra_thisyr, 
#         o3_gmi_thisyr, r_qvo3_thisyr, lat_gmi, lng_gmi)        
#     fig = plt.figure(figsize=(9,5))
#     # r(T, O3)
#     ax1 = plt.subplot2grid((2,2), (0,0), colspan=2,
#         projection=ccrs.PlateCarree(central_longitude=0.))
#     ax1.set_title(r'(a) %d r(T, O$_{\mathregular{3}}$)'%year, fontsize=16, 
#         x=0.02, ha='left')    
#     ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
#     ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
#         lat_gmi.min()+1, lat_gmi.max()-5])
#     ax1.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
#     lng_formatter = LongitudeFormatter()
#     ax1.xaxis.set_major_formatter(lng_formatter)         
#     ax1.get_xaxis().set_ticklabels([])    
#     ax1.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#     lat_formatter = LatitudeFormatter()    
#     ax1.yaxis.set_major_formatter(lat_formatter)    
#     cmap = plt.get_cmap('coolwarm')
#     mb = ax1.contourf(lng_gmi, lat_gmi, r_to3_thisyr, np.linspace(-1, 1, 9), 
#         cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
#     ax1.contourf(lng_gmi, lat_gmi, significance_r_to3_thisyr, hatches=['//////'], 
#         colors='none', transform=ccrs.PlateCarree(), zorder=2)
#     # r(q, O3)
#     ax2 = plt.subplot2grid((2,2), (1,0), colspan=2,
#         projection=ccrs.PlateCarree(central_longitude=0.))
#     ax2.set_title(r'(b) %d r(q, O$_{\mathregular{3}}$)'%year, fontsize=16, 
#         x=0.02, ha='left')
#     ax2.coastlines(lw=0.25, resolution='50m', color='k', zorder=4)
#     ax2.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
#         lat_gmi.min()+1, lat_gmi.max()-5])
#     ax2.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
#     lng_formatter = LongitudeFormatter()
#     ax2.xaxis.set_major_formatter(lng_formatter)          
#     ax2.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#     lat_formatter = LatitudeFormatter()    
#     ax2.yaxis.set_major_formatter(lat_formatter)
#     mb = ax2.contourf(lng_gmi, lat_gmi, r_qvo3_thisyr, np.linspace(-1,1,9), 
#         cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
#     ax2.contourf(lng_gmi, lat_gmi, significance_r_qvo3_thisyr, 
#         hatches=['//////'], colors='none', transform=ccrs.PlateCarree())
#     ax2.outline_patch.set_zorder(20)      
#     # Add region boxes
#     for ax in [ax1, ax2]:
#         # Western U.S.        
#         ax.add_patch(mpatches.Rectangle(xy=[-125, 35], width=24, height=15,
#             fill=None, lw=2, edgecolor='k', transform=ccrs.PlateCarree(), 
#             zorder=3))
#         ax.text(-130, 25, 'Western\nU.S.', ha='right',
#             fontweight='bold', fontsize=12, transform=ccrs.PlateCarree())
#         # Eastern U.S.
#         ax.add_patch(mpatches.Rectangle(xy=[-99, 35], width=24, height=15,
#             fill=None, lw=2, edgecolor='k', transform=ccrs.PlateCarree(), 
#             zorder=3))
#         ax.text(-70, 25, 'Eastern\nU.S.', ha='left',
#             fontweight='bold', fontsize=12, transform=ccrs.PlateCarree())    
#         # European Union 
#         ax.add_patch(mpatches.Rectangle(xy=[-10, 35], width=40, height=15,
#             fill=None, lw=2, edgecolor='k', transform=ccrs.PlateCarree(), 
#             zorder=3))
#         ax.text(-10, 55, 'E.U.', ha='left', fontsize=12, 
#             fontweight='bold', transform=ccrs.PlateCarree())      
#         # China
#         ax.add_patch(mpatches.Rectangle(xy=[90, 35], width=35, height=15,
#             fill=None, lw=2, edgecolor='k', transform=ccrs.PlateCarree(), 
#             zorder=3))     
#         ax.text(120, 55, 'Asia', ha='left', fontsize=12, 
#             fontweight='bold', transform=ccrs.PlateCarree())
#         ax.outline_patch.set_zorder(20)    
#     # Add colorbar
#     plt.gcf().subplots_adjust(left=0.02, right=0.86, hspace=0.3)
#     colorbar_axes = plt.gcf().add_axes([
#         ax1.get_position().x1+0.03, # Left
#         (ax2.get_position().y1-ax2.get_position().y0)/2.+ax2.get_position().y0, # Bottom 
#         0.02, # Width
#         ((ax1.get_position().y1-ax1.get_position().y0)/2.+ax1.get_position().y0)-
#         ((ax2.get_position().y1-ax2.get_position().y0)/2.+ax2.get_position().y0)])
#     colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
#         ticks=np.linspace(-1, 1, 9), extend='neither')
#     colorbar.ax.tick_params(labelsize=12)
#     colorbar.set_label('[$\mathregular{\cdot}$]', fontsize=16, labelpad=11)
#     plt.savefig('/Users/ghkerr/phd/globalo3/figs/map_rto3_rqvo3_%d.jpg'%year, 
#         dpi=300)
#     plt.show()
# fig = plt.figure()
# ax1 = plt.subplot2grid((2,2),(0,0), colspan=2)
# ax2 = plt.subplot2grid((2,2),(1,0), colspan=2)
# # Timeseries of r(O3, T)
# ax1.plot(np.arange(2000,2011,1), r_to3_wna, lw=2, color='#a6cee3', 
#       label='Western U.S.')
# ax1.plot(np.arange(2000,2011,1), r_to3_ena, lw=2, color='#1f78b4', 
#       label='Eastern U.S.')
# ax1.plot(np.arange(2000,2011,1), r_to3_eu, lw=2, color='#b2df8a', 
#       label='E.U.')
# ax1.plot(np.arange(2000,2011,1), r_to3_asia, lw=2, color='#33a02c', 
#       label='Asia')
# ax1.set_xlim([2000,2010])
# ax1.set_xticks(np.arange(2000,2011,1))
# ax1.set_xticklabels([])
# ax1.set_ylabel('r(T, O$_{\mathregular{3}}$) [$\cdot$]', fontsize=16)    
# # Timeseries of r(O3, q)
# ax2.plot(np.arange(2000,2011,1), r_qvo3_wna, lw=2, color='#a6cee3', 
#       label='Western U.S.')
# ax2.plot(np.arange(2000,2011,1), r_qvo3_ena, lw=2, color='#1f78b4', 
#       label='Eastern U.S.')
# ax2.plot(np.arange(2000,2011,1), r_qvo3_eu, lw=2, color='#b2df8a', 
#        label='E.U.')
# ax2.plot(np.arange(2000,2011,1), r_qvo3_asia, lw=2, color='#33a02c', 
#       label='Asia')
# ax2.set_xlim([2000,2010])
# ax2.set_xticks(np.arange(2000,2011,1))
# ax2.set_xticklabels(['2000', '', '2002', '', '2004', '', '2006', '', 
#     '2008', '', '2010'])
# ax2.set_ylabel('r(q, O$_{\mathregular{3}}$) [$\cdot$]', fontsize=16)    
# ax2.legend(ncol=4, loc=3, bbox_to_anchor=[0, -0.5])
# plt.subplots_adjust(bottom=0.2, top=0.95)
# plt.savefig('/Users/ghkerr/Desktop/timeseries_rO3T_rO3q.png', dpi=350)
# years_thisyr = [x.year for x in times_gmi]
# years_thisyr = np.array(years_thisyr)
# months_thismonth = [x.month for x in times_gmi]
# months_thismonth = np.array(months_thismonth)
# # Separate by years             
# for year in np.arange(2008, 2011, 1):
#     for month in [6, 7, 8]:
#         thisyrmonth = np.where((years_thisyr==year) &
#                                 (months_thismonth==month))[0]
#         t2m_merra_thisyr = t2m_merra[thisyrmonth]
#         qv2m_merra_thisyr = qv2m_merra[thisyrmonth]
#         o3_merra_thisyr = o3_transport_gmi[thisyrmonth]
#         # Plot yearly maps of r(O3, T) and r(O3, q)
#         r_to3_thisyr = globalo3_calculate.calculate_r(t2m_merra_thisyr, 
#             o3_merra_thisyr, lat_gmi, lng_gmi)
#         r_qvo3_thisyr = globalo3_calculate.calculate_r(qv2m_merra_thisyr, 
#             o3_merra_thisyr, lat_gmi, lng_gmi)
#         significance_r_to3_thisyr = \
#             globalo3_calculate.calculate_r_significance(t2m_merra_thisyr, 
#             o3_merra_thisyr, r_to3_thisyr, lat_gmi, lng_gmi)
#         significance_r_qvo3_thisyr = \
#             globalo3_calculate.calculate_r_significance(qv2m_merra_thisyr, 
#             o3_merra_thisyr, r_qvo3_thisyr, lat_gmi, lng_gmi)        
#         fig = plt.figure(figsize=(9,5))
#         # r(T, O3)
#         lng_gmi[-1]=360.
#         ax1 = plt.subplot2grid((2,2), (0,0), colspan=2,
#             projection=ccrs.PlateCarree(central_longitude=0.))
#         ax1.set_title(r'(a) %d-%d r(T, O$_{\mathregular{3}}$)'%(month,year), 
#             fontsize=16, x=0.02, ha='left')    
#         ax1.coastlines(lw=0.25, resolution='50m', color='k', zorder=3)
#         ax1.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
#             lat_gmi.min()+1, lat_gmi.max()-5])
#         ax1.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
#         lng_formatter = LongitudeFormatter()
#         ax1.xaxis.set_major_formatter(lng_formatter)         
#         ax1.get_xaxis().set_ticklabels([])    
#         ax1.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#         lat_formatter = LatitudeFormatter()    
#         ax1.yaxis.set_major_formatter(lat_formatter)    
#         cmap = plt.get_cmap('coolwarm')
#         mb = ax1.contourf(lng_gmi, lat_gmi, r_to3_thisyr, np.linspace(-1, 1, 9), 
#             cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
#         ax1.contourf(lng_gmi, lat_gmi, significance_r_to3_thisyr, hatches=['//////'], 
#             colors='none', transform=ccrs.PlateCarree(), zorder=2)
#         # r(q, O3)
#         ax2 = plt.subplot2grid((2,2), (1,0), colspan=2,
#             projection=ccrs.PlateCarree(central_longitude=0.))
#         ax2.set_title(r'(b) %d-%d r(q, O$_{\mathregular{3}}$)'%(month,year), 
#             fontsize=16, x=0.02, ha='left')
#         ax2.coastlines(lw=0.25, resolution='50m', color='k', zorder=4)
#         ax2.set_extent([lng_gmi.min()-180., lng_gmi.max()-180., 
#             lat_gmi.min()+1, lat_gmi.max()-5])
#         ax2.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
#         lng_formatter = LongitudeFormatter()
#         ax2.xaxis.set_major_formatter(lng_formatter)          
#         ax2.set_yticks([0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#         lat_formatter = LatitudeFormatter()    
#         ax2.yaxis.set_major_formatter(lat_formatter)
#         mb = ax2.contourf(lng_gmi, lat_gmi, r_qvo3_thisyr, np.linspace(-1,1,9), 
#             cmap=cmap, extend='neither', transform=ccrs.PlateCarree(), zorder=1)
#         ax2.contourf(lng_gmi, lat_gmi, significance_r_qvo3_thisyr, 
#             hatches=['//////'], colors='none', transform=ccrs.PlateCarree())
#         ax2.outline_patch.set_zorder(20)      
#         # Add colorbar
#         plt.gcf().subplots_adjust(left=0.02, right=0.86, hspace=0.3)
#         colorbar_axes = plt.gcf().add_axes([
#             ax1.get_position().x1+0.03, # Left
#             (ax2.get_position().y1-ax2.get_position().y0)/2.+ax2.get_position().y0, # Bottom 
#             0.02, # Width
#             ((ax1.get_position().y1-ax1.get_position().y0)/2.+ax1.get_position().y0)-
#             ((ax2.get_position().y1-ax2.get_position().y0)/2.+ax2.get_position().y0)])
#         colorbar = plt.colorbar(mb, colorbar_axes, orientation='vertical', 
#             ticks=np.linspace(-1, 1, 9), extend='neither')
#         colorbar.ax.tick_params(labelsize=12)
#         colorbar.set_label('[$\mathregular{\cdot}$]', fontsize=16, labelpad=11)
#         plt.savefig('/Users/ghkerr/phd/globalo3/figs/'+
#             'map_rto3transport_rqvo3transport_%d%d_noregion.jpg'%(month,year), dpi=300)

# years_thisyr = [x.year for x in times_gmi]
# months_thismonth = [x.month for x in times_gmi]
# years_thisyr = np.array(years_thisyr)
# months_thismonth = np.array(months_thismonth)
# o3_ena_daily, t_ena_daily, qv_ena_daily = [],[],[]
# o3_ena_monthly, t_ena_monthly, qv_ena_monthly = [],[],[]
# o3_ena_jja, t_ena_jja, qv_ena_jja = [],[],[]
# month_length = []
# # Separate by year            
# for year in np.arange(2008, 2011, 1):
#     whereyr = np.where(years_thisyr==year)[0]
#     t2m_merra_yr = t2m_merra[whereyr]
#     qv2m_merra_yr = qv2m_merra[whereyr]
#     o3_merra_yr = o3_gmi[whereyr]
#     o3_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
#         o3_merra_yr[1:], lat_gmi, lng_gmi, 270., 270.5, 45., 45.5)
#     t_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
#         t2m_merra_yr[1:], lat_gmi, lng_gmi, 270., 270.5, 45., 45.5)
#     qv_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
#         qv2m_merra_yr[1:], lat_gmi, lng_gmi, 270., 270.5, 45., 45.5)
#     o3_ena_jja.append(o3_ena.mean())
#     t_ena_jja.append(t_ena.mean())
#     qv_ena_jja.append(qv_ena.mean())    
#     # Separate by year and month
#     for month in np.arange(6, 9, 1):
#         period = np.where((years_thisyr==year) & (months_thismonth==month))[0]
#         month_length.append(len(period))
#         t2m_merra_period = t2m_merra[period]
#         qv2m_merra_period = qv2m_merra[period]
#         o3_merra_period = o3_gmi[period]
#         # Eastern U.S.
#         o3_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
#             o3_merra_period, lat_gmi, lng_gmi, 270., 270.5, 45., 45.5)
#         t_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
#             t2m_merra_period, lat_gmi, lng_gmi, 270., 270.5, 45., 45.5)
#         qv_ena, lat_ena, lng_ena = globalo3_calculate.find_grid_in_bb(
#             qv2m_merra_period, lat_gmi, lng_gmi, 270., 270.5, 45., 45.5)
#         o3_ena = np.nanmean(o3_ena, axis=tuple((1,2)))
#         t_ena = np.nanmean(t_ena, axis=tuple((1,2)))
#         qv_ena = np.nanmean(qv_ena, axis=tuple((1,2)))
#         # Append
#         o3_ena_daily.append(o3_ena)
#         t_ena_daily.append(t_ena)
#         qv_ena_daily.append(qv_ena)
#         o3_ena_monthly.append(o3_ena.mean())
#         t_ena_monthly.append(t_ena.mean())
#         qv_ena_monthly.append(qv_ena.mean())
# month_length = np.hstack([[0.], month_length])
# year_length = np.array([0, 92., 92., 92.])
# # Plotting
# fig = plt.figure()
# # O3
# xlabels = []
# ax1 = plt.subplot2grid((3,3),(0,0),colspan=3)
# ax1.set_title('Example for -90$^{\circ}$E, 45$^{\circ}$N')
# ax1.plot(np.hstack(o3_ena_daily), lw=0.5, color='k', zorder=1)
# for i in np.arange(0, len(month_length)-1, 1):
#     xmin = np.sum(np.array(month_length)[:i+1])
#     xmax = xmin + np.array(month_length)[i+1]
#     xlabels.append(xmax)
#     ax1.hlines(y=o3_ena_monthly[i], xmin=xmin, xmax=xmax, ls='--', 
#         zorder=3, lw=2)
# for i in np.arange(0, len(year_length)-1, 1):
#     xmin = np.sum(np.array(year_length)[:i+1])
#     xmax = xmin + np.array(year_length)[i+1]
#     ax1.hlines(y=o3_ena_jja[i], xmin=xmin, xmax=xmax, ls='-', color='r', 
#         zorder=2, lw=2)
# xlabels = np.hstack([[0.], xlabels])
# # Temperature
# ax2 = plt.subplot2grid((3,3),(1,0),colspan=3)
# ax2.plot(np.hstack(t_ena_daily), lw=0.5, color='k', zorder=1)
# for i in np.arange(0, len(month_length)-1, 1):
#     xmin = np.sum(np.array(month_length)[:i+1])
#     xmax = xmin + np.array(month_length)[i+1]
#     ax2.hlines(y=t_ena_monthly[i], xmin=xmin, xmax=xmax, ls='--', 
#         zorder=3, lw=2)
# for i in np.arange(0, len(year_length)-1, 1):
#     xmin = np.sum(np.array(year_length)[:i+1])
#     xmax = xmin + np.array(year_length)[i+1]
#     ax2.hlines(y=t_ena_jja[i], xmin=xmin, xmax=xmax, ls='-', color='r', 
#         zorder=2, lw=2)
# # Humidity
# ax3 = plt.subplot2grid((3,3),(2,0),colspan=3)
# ax3.plot(np.hstack(qv_ena_daily), lw=0.5, color='k', zorder=1, 
#     label='Daily')
# for i in np.arange(0, len(month_length)-1, 1):
#     xmin = np.sum(np.array(month_length)[:i+1])
#     xmax = xmin + np.array(month_length)[i+1]
#     ax3.hlines(y=qv_ena_monthly[i], xmin=xmin, xmax=xmax, ls='--', 
#         zorder=3, lw=2, label="Monthly mean" if i == 0 else "")
# for i in np.arange(0, len(year_length)-1, 1):
#     xmin = np.sum(np.array(year_length)[:i+1])
#     xmax = xmin + np.array(year_length)[i+1]
#     ax3.hlines(y=qv_ena_jja[i], xmin=xmin, xmax=xmax, ls='-', color='r', 
#         zorder=2, lw=2, label="JJA mean" if i == 0 else "")
# ax1.set_ylabel('O$_{3}$ [ppbv]')
# ax2.set_ylabel('T [K]')
# ax3.set_ylabel('q [g kg$^{-1}$]')
# for ax in [ax1, ax2, ax3]:
#     ax.set_xlim([0, 276])
#     ax.set_xticks(xlabels)
#     ax.set_xticklabels(['06-08', '07-08', '08-08',
#                         '06-09', '07-09', '08-09',
#                         '06-10', '07-10', '08-10'])
# ax3.legend(ncol=3, loc=7, bbox_to_anchor=(.75, -0.5), frameon=False)
# plt.subplots_adjust(hspace=0.3)
# plt.savefig('/Users/ghkerr/Desktop/timeseries_o3tq_-90E270N.png', dpi=350)
# do3dt2m_sma = globalo3_calculate.calculate_do3dt_sma(t2m_merra, o3_gmi, 
#    lat_gmi, lng_gmi)
# do3dq_sma = globalo3_calculate.calculate_do3dt_sma(qv2m_merra, o3_gmi, 
#    lat_gmi, lng_gmi)
            

