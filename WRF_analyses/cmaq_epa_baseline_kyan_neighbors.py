import pandas as pd
import numpy as np
from netCDF4 import Dataset
from datetime import timedelta, date, datetime

# Domain and time configuration for the simulation period
domain = 'd02_WHbase'
time = 'hourly'
year = '2016'
month = '01'
stdate, enddt = 1, 31  # Start and end dates
yr = 2016
mo = 1
name = '20160101'
start_dt = datetime(yr, mo, stdate)
end_dt = datetime(yr, mo, enddt, 23)  # End at the last hour of the end date
fstart = 'COMBINE_ACONC_'

# Input and output directories for CMAQ and EPA data
dir_cmaq = '/projects/b1045/wrf-cmaq/output/CONUS4K/output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852/postprocess/'
dir_epa = '/home/ksz4578/Heat_Pump_Project/WRF_analyses/AQS_data_2016/'

# Grid file containing lat/lon information for CMAQ output
grid = "/projects/b1045/wrf-cmaq/output/CONUS4K/output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852/lat_lon_CONUS4K_d02.nc"

# EPA codes for pollutants and corresponding variable names
epa_code = ['42401', '42602', '44201', '42101', '88101', '81102']
var = ['SO2', 'NO2', 'O3', 'CO', 'PM25_TOT', 'PM10']

# Generating filenames for each EPA pollutant dataset based on the time and year
epa_files = [dir_epa + '%s_%s_%s.csv' % (time, epa_code[i], year,) for i in range(len(epa_code))]

# Load the lat/lon grid from the CMAQ netCDF file
la, lo = 'LAT', 'LON'
cmaq_lon, cmaq_lat = np.array(Dataset(grid)[lo]), np.array(Dataset(grid)[la])

# Define bounding box for cropping EPA data based on the CMAQ grid's lat/lon range
llat, ulat, llon, ulon = cmaq_lat.min(), cmaq_lat.max(), cmaq_lon.min(), cmaq_lon.max()

def find_index(stn_lon, stn_lat, wrf_lon, wrf_lat):
    # Find nearest grid points for station coordinates by minimizing the lat/lon distance
    xx, yy = [], []
    for i in range(len(stn_lat)):
        abslat = np.abs(wrf_lat - stn_lat[i])
        abslon = np.abs(wrf_lon - stn_lon[i])
        c = np.maximum(abslon, abslat)  # Handle grids with different lat/lon scales
        args = np.where(c == np.min(c))  # Find indices of the nearest grid point
        print(args)
        xx.append(args[-2])
        yy.append(args[-1])
    xx = [xx[i][0] for i in range(len(xx))]  # Flatten list of indices
    yy = [yy[i][0] for i in range(len(yy))]
    return xx, yy

def make_cmaq_fnames(fstart, year, month, stdate, enddt):
    # Create filenames for each day of CMAQ data between stdate and enddt
    return [fstart + year + month + str(i).zfill(2) + '.nc' for i in range(stdate, enddt)]

def crop_epa(file, start_dt, end_dt, llat, ulat, llon, ulon):
    # Read and crop EPA dataset to the specified lat/lon bounding box and time range
    df = pd.read_csv(file)
    df = df[(df['Latitude'] >= llat) & (df['Latitude'] <= ulat) & (df.Longitude >= llon) & (df.Longitude <= ulon)].reset_index()
    
    # Combine date and time columns into a single datetime column
    df['Datetime GMT'] = pd.to_datetime(df['Date GMT'] + ' ' + df['Time GMT'])
    
    # Filter data based on the start and end datetime
    df = df[(df['Datetime GMT'] >= pd.to_datetime(start_dt)) & (df['Datetime GMT'] <= pd.to_datetime(end_dt))]
    
    # Get unique longitude and latitude values for each station
    lon, lat = df['Longitude'].unique(), df['Latitude'].unique()
    return lon, lat, df

def get_cmaq_avg(cmaqv, x, y):
    """
    Get the average value of the grid cell at (x, y) and its neighboring cells.
    Handles boundary cases by limiting indices within the grid.
    """
    x_min, x_max = max(0, x - 1), min(cmaqv.shape[2] - 1, x + 1)
    y_min, y_max = max(0, y - 1), min(cmaqv.shape[3] - 1, y + 1)
    
    # Average the values over a 3x3 grid (center + neighbors)
    return np.nanmean(cmaqv[:, :, x_min:x_max + 1, y_min:y_max + 1], axis=(2, 3))

# Generate filenames for the CMAQ output files for the date range
fnames = make_cmaq_fnames(fstart, year, month, stdate, enddt + 1)

# Load each CMAQ file into a list
cmaq = [Dataset(dir_cmaq + fnames[i]) for i in range(len(fnames))]

# Loop through each pollutant variable
for v in range(len(var)):
    print(epa_files[v], start_dt, end_dt, llat, ulat, llon, ulon)
    
    # Crop EPA data to relevant lat/lon and time range for each pollutant
    epa_lon, epa_lat, df = crop_epa(epa_files[v], start_dt, end_dt, llat, ulat, llon, ulon)
    
    # Load the pollutant concentration from CMAQ data for the current pollutant
    cmaqv = [np.array(cmaq[d][var[v]]) for d in range(len(cmaq))]

    dataframes = []
    
    # Loop through each EPA station's longitude and latitude
    for i in range(len(epa_lon)):
        tmp = df[(df['Longitude'] == epa_lon[i]) & (df['Latitude'] == epa_lat[i])]
        
        # Extract station information
        tsn = tmp['State Name'].unique()[0]
        tcn = tmp['County Name'].unique()[0]
        um = tmp['Units of Measure'].unique()[0]
        
        # Find the corresponding grid indices for the EPA station
        x, y = find_index([epa_lon[i]], [epa_lat[i]], cmaq_lon, cmaq_lat)
        x, y = x[0], y[0]
        
        # Ensure the station is within the grid bounds and exclude edge cases
        if (x != 0) | (y != 0) | (x != cmaq_lon.shape[0]) | (y != cmaq_lat.shape[1]):
            tmp.index = tmp['Datetime GMT']
            t_index = pd.date_range(start=start_dt, end=end_dt, freq='1h')
            tmp = tmp.resample('1H').mean().reindex(t_index).fillna(np.nan)
            
            # Assign the grid indices to the dataframe
            tmp['x'], tmp['y'] = x, y
            tmp['State Name'] = tsn
            tmp['County Name'] = tcn
            tmp['Units of Measure'] = um
            
            # Calculate average pollutant values over a 3x3 grid around (x, y)
            z = [get_cmaq_avg(cmaqv[d], x, y) for d in range(len(cmaq))]
            zi = []
            for q in range(len(z)):
                zi = zi + list(z[q])

            z = np.array(zi).ravel()
            print(len(z))
            
            # Handle case where CMAQ data is shorter than EPA data
            if len(z) < len(tmp): 
                z = list(z) + [np.nan] * (len(tmp) - len(z))
            
            tmp['CMAQ'] = z
            tmp['dt'] = tmp.index
            tmp['level_0'] = tmp.index
            tmp = tmp.reset_index(drop=True)
            dataframes.append(tmp)
        else:
            print(tmp['State Name'], tmp['County Name'], epa_lon[i], epa_lat[i], "failed")

    # Concatenate all station data for the current pollutant and save to a CSV file
    final_df = pd.concat(dataframes, ignore_index=True)
    final_df.to_csv(dir_epa + var[v] + '_' + domain + '_' + str(year) + '_' + str(mo) + '_EPA_CMAQ_Combine_neighbors.csv')
    print('done ' + var[v])
