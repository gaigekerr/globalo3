#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 16:18:43 2019

open napso3, aqs o3 on 23 july 2019

@author: ghkerr
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
        # station
        hours_str = ['Hour %d'%x for x in hours]
        naps_hours = naps[hours_str]
        # Form row-wise daily average for a subset of columns with missing 
        # values and add to O3 column 
        naps['O3 (ppbv)'] = naps_hours.mean(axis=1).values
        # Name index column and convert to string
        naps.index.names = ['Date Local']
        naps.index = [x.strftime('%Y-%m-%d') for x in naps.index]        
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
    aqs_allyr.set_index('Date Local', inplace = True)
    # Rename O3 column and convert from ppm to ppb
    aqs_allyr = aqs_allyr.rename(columns=
                                 {'Sample Measurement' : 'O3 (ppbv)'})
    aqs_allyr['O3 (ppbv)'] = aqs_allyr['O3 (ppbv)']*1000.
    print('AQS O3 for %d-%d loaded in '%(years[0], years[-1])+
          '%.2f seconds!'%(time.time() - start_time))          
    return aqs_allyr