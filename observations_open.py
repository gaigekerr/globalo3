#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Function opens PM2.5 and O3 observations from various observational networks
and cleans data to a format useable for calculating the O3-temperature 
relationship. 

Revision History
    23072019 -- initial version created
    23072019 -- functions 'open_napso3' and 'open_aqso3' added
    24072019 -- function 'open_napso3' edited to reflect different nomenclature
                for sampling hours: NAPS data begin with "hour 1" rather than 
                0 (i.e., AQS, EMEP), so this must be taken into account 
                when selecting hours of interest
    26092019 -- function 'open_chinao3' added
    21102019 -- change longitude from (-180-180) to (0-360) for AQS and NAPS
                dataset
"""
def get_merged_csv(flist, **kwargs):
    """function reads CSV files in the list comprehension loop, this list of
    DataFrames will be passed to the pd.concat() function which will return 
    single concatenated DataFrame
    From https://stackoverflow.com/questions/35973782/reading-multiple-csv-
    files-concatenate-list-of-file-names-them-into-a-singe-dat
    """
    import pandas as pd
    return pd.concat([pd.read_csv(f, **kwargs) for f in flist], 
                      ignore_index = True)
    
def open_napso3(years, months, hours):
    """function opens hourly O3 observations from the Canadian National Air
    Pollution Surveillance Program (NAPS) for the specified months, years, 
    and hours. Output data frame represents daily-averaged concentrations 
    (over the hours specified) at individual NAPS sites. Note this function 
    assumes that the reported O3 observations were taken at local time, not 
    UTC. 

    Parameters
    ----------
    years : list
        Year or range of years in measuring period
    months : list
        Months (integers with Jan = 1, Feb = 2, etc.) for the season of 
        interest
    hours : list 
        Hours over which hourly O3 observations are averaged, units of local 
        time
        
    Returns
    -------
    naps_allyr : pandas.core.frame.DataFrame
        Daily averaged O3 concentrations at individual NAPS sites for the 
        specified measuring period. DataFrame contains NAPS ID and site 
        latitude and longitude
    """
    import time
    start_time = time.time()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading NAPS O3 ...')
    import pandas as pd
    import numpy as np
    PATH_NAPS = '/Users/ghkerr/phd/observations/o3/NAPS/'
    # Open NAPS siting information 
    stations = pd.read_csv(PATH_NAPS+'NAPS_stations.csv', header = 0,
                           low_memory = False)
    stations_cols = ['NAPS_ID', 'STATION_NAME', 'STATUS', 'STREET_ADDRESS',
                     'CITY', 'PROVINCE', 'TimeZone', 'Lat_Decimal', 
                     'Long_Decimal', 'Elevation_m']
    stations = pd.DataFrame(stations, columns = stations_cols)
    # Read NAPS fixed-width files (information about the width can be found
    # under "Find Out More About NAPS Monitoring Data Download Formats" on 
    # http://maps-cartes.ec.gc.ca/rnspa-naps/data.aspx?lang=en)
    naps_col_width = {'Pollutant code' : 3, 'Station (NAPS ID)' : 6,
                      'Year' : 4, 'Month' : 2, 'Day' : 2,
                      'Average for day'	: 4, 'Minimum for day' : 4,
                      'Maximum for day' : 4, 'Hour 0' : 4, 'Hour 1' : 4,
                      'Hour 2' : 4, 'Hour 3' : 4, 'Hour 4' : 4, 'Hour 5' : 4, 
                      'Hour 6' : 4, 'Hour 7' : 4, 'Hour 8' : 4,  'Hour 9' : 4, 
                      'Hour 10' : 4, 'Hour 11' : 4, 'Hour 12' : 4, 
                      'Hour 13' : 4, 'Hour 14' : 4, 'Hour 15' : 4,
                      'Hour 16' : 4, 'Hour 17' : 4, 'Hour 18' : 4, 
                      'Hour 19' : 4, 'Hour 20' : 4, 'Hour 21' : 4, 
                      'Hour 22' : 4, 'Hour 23' : 4}
    naps_allyr = pd.DataFrame()
    for year in years:
        naps = pd.read_fwf((PATH_NAPS+'NAPS_O3_%d.csv' %year), 
                         widths=list(naps_col_width.values()),
                         names=list(naps_col_width.keys()))
        # Convert individual year, month, and day columns into single datetime
        # column and set as index
        naps.index = pd.to_datetime(naps[['Year', 'Month', 'Day']])
        # Drop columns that are not needed  
        del naps['Pollutant code'], naps['Average for day'], \
            naps['Minimum for day'], naps['Maximum for day'], naps['Year'], \
            naps['Month'], naps['Day']
        # Select months in measuring period     
        naps = naps.loc[naps.index.month.isin(months)]
        # Find latitude and longitude at each site in station database, add to 
        # O3 observations as new columns 
        naps['Latitude'] = np.nan
        naps['Longitude'] = np.nan
        naps['O3 (ppbv)'] = np.nan
        for index, row in stations.iterrows():
            station_id = row['NAPS_ID']
            station_lat = row['Lat_Decimal']
            station_lng = row['Long_Decimal']
            naps.loc[naps['Station (NAPS ID)'] == station_id, ['Latitude', 
                     'Longitude']] = station_lat, station_lng           
        # Convert missing data (-999) to NaN
        naps = naps.replace(-999, np.nan)
        # Select columns corresponding to hours of interest (local time) at 
        # station; since the NAPS observations are indicated by hour starting
        # at 1 (rather than 0), subtract 1
        hours_str = ['Hour %d'%(x-1) for x in hours]
        naps_hours = naps[hours_str]
        # Form row-wise daily average for a subset of columns with missing 
        # values and add to O3 column 
        naps['O3 (ppbv)'] = naps_hours.mean(axis=1).values
        # Name index column and convert to string
        naps.index.names = ['Date Local']
        naps.index = [x.strftime('%Y-%m-%d') for x in naps.index]        
        # Convert longitude from (-180-180) to (0-360)
        naps['Longitude'] = naps['Longitude']%360.        
        # Rename station ID column
        naps = naps.rename(columns={'Station (NAPS ID)' : 'Station ID'})
        # Drop columns with hourly values
        naps = naps.loc[:,~naps.columns.str.startswith('Hour')]
        # Append yearly (seasonal) values to multi-year DataFram
        naps_allyr = naps_allyr.append(naps)
    print('NAPS O3 for %d-%d loaded in '%(years[0], years[-1])+
          '%.2f seconds!'%(time.time() - start_time))        
    return naps_allyr

def open_aqso3(years, months, hours):
    """function loads hourly observations of AQS O3 for the specified years
    and months and produces daily averages over the averaging hours.

    Parameters
    ----------
    years : list
        Year or range of years in measuring period
    months : list
        Months (integers with Jan = 1, Feb = 2, etc.) for the season of 
        interest
    hours : list 
        Hours over which hourly O3 observations are averaged, units of local 
        time
        
    Returns
    -------
    aqs_allyr : pandas.core.frame.DataFrame
        Daily averaged O3 concentrations at individual AQS sites for the 
        specified measuring period. DataFrame contains AQS ID and site 
        latitude and longitude    
    """
    import time
    start_time = time.time()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading AQS O3 ...') 
    import pandas as pd
    import numpy as np
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
    aqs_allyr = get_merged_csv(O3files, dtype = dtype, index_col = None, 
                          usecols = list(dtype.keys())) 
    # Select months in measuring period     
    aqs_allyr = aqs_allyr.loc[pd.to_datetime(aqs_allyr['Date Local']).
                         dt.month.isin(months)]    
    aqs_allyr['Station ID'] = np.nan
    # Select hours of interest
    hours_str = ['%s:00'%str(x).zfill(2) for x in hours]
    aqs_allyr = aqs_allyr.loc[aqs_allyr['Time Local'].isin(hours_str)]
    # Combine State Code, County Code, and Site Number into single site 
    # ID column 
    aqs_allyr['Station ID'] = (aqs_allyr['State Code'] + '-' +
        aqs_allyr['County Code'] + '-' + aqs_allyr['Site Num'])
    # Drop individual code columns
    del aqs_allyr['State Code'], aqs_allyr['County Code'], \
        aqs_allyr['Site Num']
    # Find daily averaged values 
    aqs_allyr = aqs_allyr.groupby(['Date Local', 'Station ID']).mean()
    # Turn MultiIndex to column
    aqs_allyr.reset_index(inplace = True)
    # Make date index 
    aqs_allyr.set_index('Date Local', inplace=True)
    # Rename O3 column and convert from ppm to ppb
    aqs_allyr = aqs_allyr.rename(columns=
                                 {'Sample Measurement' : 'O3 (ppbv)'})
    aqs_allyr['O3 (ppbv)'] = aqs_allyr['O3 (ppbv)']*1000.
    # Convert longitude from (-180-180) to (0-360)
    aqs_allyr['Longitude'] = aqs_allyr['Longitude']%360.        
    print('AQS O3 for %d-%d loaded in '%(years[0], years[-1])+
          '%.2f seconds!'%(time.time() - start_time))          
    return aqs_allyr

def open_emepo3(years, months, hours): 
    """function opens European Monitoring and Evaluation Programme (EMEP)
    hourly O3 observations and extracts observations from the specified 
    measuring period. Note that this function assumes that observations were
    recorded at local time, not UTC. 
    
    Parameters
    ----------
    years : list
        Year or range of years in measuring period
    months : list
        Months (integers with Jan = 1, Feb = 2, etc.) for the season of 
        interest
    hours : list 
        Hours over which hourly O3 observations are averaged, units of local 
        time
        
    Returns
    -------
    emep : pandas.core.frame.DataFrame
        Daily averaged O3 concentrations at individual EMEP sites for the 
        specified measuring period. DataFrame contains EMEP ID and site 
        latitude and longitude        
    """
    import time
    start_time = time.time()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading EMEP O3 ...') 
    import re
    import glob
    import pandas as pd
    import numpy as np
    #https://stackoverflow.com/questions/33997361
    def dms2dd(s):
        # example: s = """0°51'56.29"S"""
        degrees, minutes, seconds, direction = re.split('[°\'"]+', s)
        dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60);
        if direction in ('S','W'):
            dd*= -1
        return dd
    PATH_EMEP = '/Users/ghkerr/phd/observations/o3/EMEP/'
    # Open EMEP siting information 
    stations = pd.read_csv(PATH_EMEP+'EMEP_stations.csv', header = 0, 
                           sep = '\s+')
    # Convert degree/minute/second siting to decimal degrees
    stations['latitude'] = stations['latitude'].apply(dms2dd)
    stations['longitude'] = stations['longitude'].apply(dms2dd)
    # Column names, dtypes for O3 observation files
    dtype = {'Start Hour' : np.str, 'End Hour' : np.str,
             'Concentration' : np.str, 'Dummy' : np.str, 
             'Validity' : np.str}
    emep=pd.DataFrame()
    # Loop through years and open all site files
    for year in [2008, 2009, 2010]:
        files_ty = glob.glob(PATH_EMEP + 'all_o3_%s/'%year + '*_%s.dat'%year)
        files_ty_all = []
        for filename in files_ty:
            # Find station ID, station latitude, and station longitude based 
            # on filename
            station_ID = filename.split('/')[-1]       
            station_ID = station_ID.split('_')[0]
            station_lat = stations.loc[stations['Code'] == station_ID
                                       ]['latitude'].values[0]
            station_lng = stations.loc[stations['Code'] == station_ID
                                       ]['longitude'].values[0]
            df = pd.read_csv(filename, dtype = dtype, index_col = None, 
                             header = None, sep='\s+',
                             names = list(dtype.keys()), skiprows = 3)
            # Add station ID and siting information to DataFrame
            df['Latitude'] = station_lat
            df['Longitude'] = station_lng
            df['Station ID'] = station_ID        
            files_ty_all.append(df)
        files_ty_all = pd.concat(files_ty_all, axis=0, ignore_index=True)
        emep=emep.append(files_ty_all,ignore_index=True)
    # Convert O3 from string to float
    emep['Concentration'] = emep['Concentration'].astype(float)
    # Convert missing data to NaN (n.b. most missing data is indicated
    # as 999.0 but some is 9999., etc)
    emep = emep.replace(999., np.nan)
    emep = emep.replace(999.99, np.nan)
    emep = emep.replace(9999., np.nan)    
    # Split time column into two columns (date and hour)
    emep['Date Local'], emep['Hour'] = emep['Start Hour'].str.split('T').str
    # Drop rows with NaNs (otherwise an error is raised when setting the date 
    # as the index)
    emep = emep[pd.notnull(emep['Hour'])]
    # Set date column as index, rename
    emep.index = pd.to_datetime(emep['Date Local'])
    emep.index.names = ['Date Local']
    # Convert O3 from µg/m3 to ppbv
    # n.b., this conversion assumes an ambient pressure of 1 atmosphere and a 
    # temperature of 25 degrees Celsius. If this is the case, for O3 
    # 1 ppb = 2.00 μg/m3
    # The general equation is μg/m3 = (ppb)*(12.187)*(M)/(273.15 + °C), where
    # M is the molecular weight of the gaseous pollutant. An atmospheric 
    # pressure of 1 atmosphere is assumed.
    #/www2.dmu.dk/AtmosphericEnvironment/Expost/database/docs/PPM_conversion.pdf
    emep['Concentration'] = emep['Concentration']/2.    
    # Select months in measuring period     
    emep = emep.loc[emep.index.month.isin(months)]
    # Drop columns that are not needed
    del emep['End Hour'], emep['Dummy'], emep['Start Hour'], \
        emep['Date Local'], emep['Validity']
    # Select hours of interest (n.b. if overpass2 time (13:00-14:00 local time)
    # is needed, we'd select hour 13 in the 'Hour' column (originally 'Start
    # Hour'). 
    hours_str = ['%s:00:00'%str(x).zfill(2) for x in hours]
    emep = emep.loc[emep['Hour'].isin(hours_str)]
    # Find daily averaged values 
    emep = emep.groupby(['Date Local', 'Station ID']).mean()
    # Rename O3 column 
    emep = emep.rename(columns={'Concentration' : 'O3 (ppbv)'})
    # Turn MultiIndex to column
    emep.reset_index(inplace=True)
    # Make date index 
    emep.set_index('Date Local', inplace=True)
    # Convert index to string
    emep.index = [x.strftime('%Y-%m-%d') for x in emep.index]            
    # Convert longitude from (-180-180) to (0-360)
    emep['Longitude'] = emep['Longitude']%360.     
    print('EMEP O3 for %d-%d loaded in '%(years[0], years[-1])+
          '%.2f seconds!'%(time.time() - start_time))     
    return emep

def open_chinao3(years, months, hours, cityavg=False): 
    """Hourly O3 from Chinese cities for  2016-2017 was shared by Ke Li from 
    his recent paper. In addition to be stored locally, data is also stored on 
    Dropbox. 
    Site latitudes and longitudes: 
        https://www.dropbox.com/s/643cjir13u2ppkw/sta_lon_lat.txt?dl=0
    O3 observations:     
        https://www.dropbox.com/s/yx8mntdlc77mnhk/ozone_16-17.zip?dl=0
    Site cities: 
        https://www.dropbox.com/s/s0eaxe76boywphd/staname?dl=0
    Funciton opens data and siting information and extracts O3 observations
    for the hours and months of interest. If specified, the function also 
    averages observations over Chinese cities to obtain daily average city-wide 
    O3 concentrations. 

    Parameters
    ----------
    years : list
        Year or range of years in measuring period
    months : list
        Months (integers with Jan = 1, Feb = 2, etc.) for the season of 
        interest
    hours : list 
        Hours over which hourly O3 observations are averaged, units of local 
        time
    cityavg : bool
        Default False; if True, function averages daily O3 observations over 
        cities to produce a daily, city-wide average measurement
        
    Returns
    -------
    o3df : pandas.core.frame.DataFrame
        Daily-averaged O3 concentrations at individual sites or city-wide 
        averages for the specified measuring period. DataFrame contains 
        city name and site latitude and longitude    
    
    References
    ----------
    Li, K., Jacob, D. J., Liao, H., Shen, L., Zhang, Q., and Bates, K. H., 
    (2019). Anthropogenic drivers of 2013–2017 trends in summer surface ozone 
    in China. PNAS, 116(2). doi:10.1073/pnas.1812168116
    """
    import time
    start_time = time.time()
    print('# # # # # # # # # # # # # # # # # # # # # # # # # #\n'+
          'Loading Chinese O3 ...') 
    import xarray as xr
    import pandas as pd
    import numpy as np
    PATH_CHINA = '/Users/ghkerr/phd/observations/o3/china/'
    infiles = xr.open_mfdataset(PATH_CHINA+'site/'+'china_sites_*.nc',
        concat_dim='time', combine='nested')
    o3 = infiles.o3.values
    # Data from Ke Li is for 2016-2017
    date_range = pd.date_range(start='01/01/2016', end='12/31/2017')
    date_range = np.array([x.strftime('%Y-%m-%d') for x in date_range])
    # Dimensions are (days in measuring period, hours, no. sites). From what
    # I can tell, the time is in local time; i.e., a visual inspection of 
    # the diurnal cycle of np.nanmean(o3, axis=0)[:, site] yields diurnal 
    # curves with peaks around 15 hours. 
    # Extract O3 at hour(s) of interest. This assumes that hours are in 
    # military time; i.e., 0 = 12am, 13 = 1pm, etc. 
    o3 = o3[:, np.array(hours)]
    # Compute daily average
    o3 = np.nanmean(o3, axis=1)
    o3_shape = o3.shape
    # Stack O3 observations such that first 731 elements represent the 
    # observed O3 at the first station, etc
    o3 = np.hstack(o3.T)
    # Repeat date range over the number of stations
    date_range = np.tile(date_range, o3_shape[1])
    # Open Chinese siting information
    siting = pd.read_csv(PATH_CHINA+'sta_lon_lat.txt', header=None,
                           sep = '\t')
    cities = pd.read_csv(PATH_CHINA+'staname', header=None)
    # Rename the columns
    siting.columns = ['Longitude', 'Latitude']
    cities.columns = ['Cities']
    # Extract siting information and repeat similar to the date range
    lat = siting['Latitude'].values
    lat = np.repeat(lat, o3_shape[0])
    lng = siting['Longitude'].values
    lng = np.repeat(lng, o3_shape[0])
    cities = cities['Cities'].values
    cities = np.repeat(cities, o3_shape[0])
    # Create DataFrame in the same format as DataFrame for EMEP, NAPS, AQS, 
    # etc.; n.b., in this case since Chinese sites don't have site IDs, the 
    # DataFrame column "Station ID" is the name of the city that the site
    # is located in 
    o3df = pd.DataFrame(np.array([date_range, cities, lat, lng, o3]).T, 
        columns=['Date', 'Station ID', 'Latitude', 'Longitude', 'O3 (ppbv)'])
    # Convert from object to float; set date as index
    o3df.index = pd.to_datetime(o3df['Date'])
    o3df.index.names = ['Date']
    del o3df['Date']
    o3df['O3 (ppbv)'] = o3df['O3 (ppbv)'].astype(float)
    o3df['Latitude'] = o3df['Latitude'].astype(float)
    o3df['Longitude'] = o3df['Longitude'].astype(float)
    # Select months in measuring period     
    o3df = o3df.loc[o3df.index.month.isin(months)]
    # If optional parameter "cityavg" is True, all stations located in a given 
    # city will be averaged such that there is one daily average concentration 
    # per city per day
    if cityavg==True: 
        o3df = o3df.groupby(['Station ID', o3df.index]).mean()
        # Turn MultiIndex to column
        o3df.reset_index(inplace=True)
        # Make date index 
        o3df.set_index('Date', inplace=True)
    # Convert index to string
    o3df.index = [x.strftime('%Y-%m-%d') for x in o3df.index]                    
    print('Chinese O3 for %d-%d loaded in '%(years[0], years[-1])+
          '%.2f seconds!'%(time.time() - start_time))            
    return o3df