import pandas as pd
from scipy.stats import pearsonr
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from netCDF4 import Dataset

####################################################################################################################################
# functions for model evaluation
def stats(data, prediction):
    x, y = data[~np.isnan(data)], prediction[~np.isnan(data)]  # get rid of NaNs
    mu_d, mu_p = np.mean(x), np.mean(y)
    mb = np.sum(y - x) / len(x)
    ge = np.mean(np.abs(x - y))
    rmse = np.sqrt(np.mean((y - x) ** 2))
    r, p = st.pearsonr(x, y)
    return mu_d, mu_p, mb, ge, rmse, r, p

def makefig(obs_d02, d02, var): 
    d02 = d02.copy(deep = True)
    obs_d02 = obs_d02.copy(deep = True)
    # Formatting the data
    obs_d02 = obs_d02.reset_index()
    d02 = d02.reset_index()

    obs_d02['hr'] = [obs_d02.dt[i].hour for i in range(len(obs_d02))]
    d02['hr'] = [d02.dt[i].hour for i in range(len(d02))]

    d02 = d02.drop('dt', axis=1)
    obs_d02 = obs_d02.drop('dt', axis=1)

    # Plot diurnal average
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot(np.nanmean(obs_d02.groupby('hr').mean(), axis=1), linestyle='--', linewidth=1, label='NCDC')
    ax.plot(np.nanmean(d02.groupby('hr').mean()), linewidth=1, label='d02')
    ax.grid()
    ax.set_xlabel('Hour of Day (GMT)')
    ax.set_ylabel(var)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(var + '_alld02_diurnal_avg.png', dpi=300)
    plt.close()

    # Plot all hourly data
    fig, ax = plt.subplots(figsize=(6, 3))
    
    ax.plot(np.nanmean(obs_d02, axis=1), linewidth=1, linestyle='--', label='NCDC')
    ax.plot(np.nanmean(d02, axis=1), linewidth=1, label='d02')
    ax.set_xlabel('Hour of Simulation')
    ax.set_ylabel(var)
    ax.grid()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(var + '_alld02_diurnal_all.png', dpi=300)
    plt.close()

####################################################################################################################################

# Read netCDF for location information
d02 = Dataset("/projects/b1045/wrf-cmaq/output/CONUS4K/output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852/lat_lon_CONUS4K_d02.nc")
lat, lon = d02['LAT'][0][0].data, d02['LON'][0][0].data

# Directory and variables for CMAQ run
domain = 'd02'
dir_WRF = r"/home/ksz4578/Heat_Pump_Project/WRF_analyses/CMAQ_LCD/"
#dir_WRF = r"/projects/b1045/wrf-cmaq/output/CONUS4K/"
sim = ['output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852']
var = ['t2', 'rh', 'winds', 'winddir']

years = ['2016'] * 4
months = ['1'] * 4
days = ['31'] * 4
months_start = ['1'] * 4
years_start = ['2016'] * 4
ss = sim * 4

# Format file names
wrf_files_d02 = [dir_WRF + sim[i] + '/' + v + 'd02.csv' for v in var for i in range(len(sim))]
print(wrf_files_d02)

# Station data files for comparison
station_files = [
    'output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852/wrfcheck_withstations_output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852_012016.csv', 'output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852/wrfcheck_withstations_output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852_RH.csv',
'output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852/wrfcheck_withstations_output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852_Wind.csv',
'output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852/wrfcheck_withstations_output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852_WindDir.csv',
]

# Process data for `d02`
for f in range(len(wrf_files_d02)):
    stn = pd.read_csv(station_files[f], index_col=0)
    
    times = pd.read_csv(ss[f] + '/completeddata_mini_extras2.csv', index_col=0)
    times_perf = pd.date_range(start=f'{months[f]}-01-{years[f]}', end=f'{months[f]}-{days[f]} {years[f]} 23:00:00', freq='H')
    
    # Filter and process the station data
    stn_d02 = stn[stn.in_d02 == True].reset_index()
    obs_d02 = pd.DataFrame([stn_d02[str(i)].tolist() for i in range(len(times))])
    d02 = pd.read_csv(wrf_files_d02[f], index_col=0).T

    # Convert units (temperature, wind speed)
    if f == 0:
        obs_d02 = (obs_d02 - 32) * 5 / 9 + 273.15  # Convert Fahrenheit to Kelvin
    if f == 2:
        d02 = 2.236936 * d02  # Convert m/s to mph
    
    # Format time columns and align data
    times = np.array(pd.to_datetime(times['0'])[0:len(times)])
    d02['dt'] = np.array(times_perf)[0:len(d02) + 1]
    obs_d02['dt'] = times
    
    # Resample to hourly averages
    obs_d02 = obs_d02.set_index('dt').resample('h').mean()
    d02 = d02.set_index('dt')
    
    # Plot and calculate statistics
    #var_name = wrf_files_d02[f].split('/')[-2].split('4km')[0] + wrf_files_d02[f].split('/')[-1].split('d02')[0]
    var_name = var[f] + 'd02'
    makefig(obs_d02, d02, var_name)
