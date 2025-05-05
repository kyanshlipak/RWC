"""
WRF Model Evaluation and Validation Toolkit

Description:
This script performs comprehensive statistical evaluation of Weather Research and Forecasting (WRF) model outputs
against observational data. It calculates performance metrics (RMSE, bias, correlation) for multiple meteorological
variables and generates diagnostic visualizations including:
- Spatial validation maps
- Temporal profile comparisons
- Diurnal cycle analyses

Key Functionalities:
1. Statistical Analysis:
   - Computes standard metrics (RMSE, MB, GE) for continuous variables
   - Specialized circular statistics for wind direction
   - Pearson correlation (linear and circular)

2. Visualization:
   - Geospatial plots of model performance metrics
   - Time series comparisons of observed vs predicted
   - Diurnal profile analyses
   - Multi-panel diagnostic plots

3. Data Processing:
   - Handles unit conversions (Fahrenheit/Kelvin, m/s/mph)
   - Temporal alignment of model/observation data
   - Coordinate transformations (LCC to WGS84)

Input Data:
- WRF model outputs (NetCDF format)
- Observational data (CSV format)
- Station metadata including geographic coordinates

Outputs:
- Statistical summary tables (Pandas DataFrames)
- Publication-quality spatial plots (Cartopy/Matplotlib)
- Temporal analysis plots
- Diurnal cycle visualizations

Dependencies:
- xarray, pandas, numpy, scipy
- cartopy, matplotlib
- pyproj (for coordinate transformations)
- netCDF4

Kyan Shlipak May 2025
"""

import pandas as pd
import numpy as np
import scipy.stats as st
from netCDF4 import Dataset

####################################################################################################################################
# Functions for model evaluation

# Function to compute angular distance for wind direction
def angular_distance(a, b):
    # Calculate the difference
    diff = np.abs(a - b)
    
    # Ensure the difference is within the range [0, 180] (shortest angle between the two)
    return np.minimum(diff, 360 - diff)

# Function to compute circular (angular) Pearson correlation
def circular_correlation(theta1_deg, theta2_deg):
    # Convert degrees to radians
    theta1 = np.deg2rad(theta1_deg)
    theta2 = np.deg2rad(theta2_deg)
    
    # Compute mean directions
    mean1 = np.arctan2(np.mean(np.sin(theta1)), np.mean(np.cos(theta1)))
    mean2 = np.arctan2(np.mean(np.sin(theta2)), np.mean(np.cos(theta2)))
    
    # Compute numerator and denominator
    sin_diff1 = np.sin(theta1 - mean1)
    sin_diff2 = np.sin(theta2 - mean2)
    numerator = np.sum(sin_diff1 * sin_diff2)
    denominator = np.sqrt(np.sum(sin_diff1**2) * np.sum(sin_diff2**2))
    
    return numerator / denominator

# Modified stats function to account for wind direction using angular distance and circular Pearson
def stats(data, prediction, is_wind_direction=False):
	x, y = data[~np.isnan(data)], prediction[~np.isnan(data)]  # get rid of NaNs

	if x.shape[0] < 2 or y.shape[0] < 2:
		return 0, 0, 0, 0, 0, 0, 0

	# Use angular distance for wind direction or normal absolute difference
	if is_wind_direction:
		# For wind direction, use angular distance for MB and circular Pearson for r
		error = angular_distance(x, y)
		r = circular_correlation(x, y)
		# diff = (y - x + 180) % 360 - 180
		# mb = np.mean(diff)  # Mean Bias using angular distance for wind direction
		
		# plt.figure(figsize=(6, 4))
		# bins = np.arange(0, 361, 30)
		# plt.hist(y, bins=bins, alpha=0.5, label='Prediction', edgecolor='k')
		# plt.hist(x, bins=bins, alpha=0.5, label='Observation', edgecolor='k')
		# plt.xticks(bins)
		# plt.xlabel("Wind Direction (degrees)")
		# plt.ylabel("Frequency")
		# plt.title("Wind Direction Histogram")
		# plt.legend()
		# plt.grid(True)
		# plt.tight_layout()
		# plt.show()
	else:
		# For all other variables, use linear difference
		error = np.abs(x - y)
		r, _ = st.pearsonr(x, y)  # Standard Pearson correlation
		  # Mean Bias using angular distance for wind direction
		


	# Calculate each statistic
	mb = np.sum(y - x) / len(x)
	mu_d, mu_p = np.mean(x), np.mean(y)
	ge = np.mean(error)
	rmse = np.sqrt(np.mean((error) ** 2))
	p = None  # p-value for correlation is not defined here for circular correlation, can be added if needed

	return mu_d, mu_p, mb, ge, rmse, r, p

####################################################################################################################################

# Read netCDF for location information
sim_2020 = ['output_CONUS4K_d02_2020_RWC_4km_sf_rrtmg_10_10_1_v3852']
sim = ['output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852']

d02 = Dataset(f"/projects/b1045/wrf-cmaq/output/CONUS4K/output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852/lat_lon_CONUS4K_d02.nc")
lat, lon = d02['LAT'][0][0].data, d02['LON'][0][0].data

# Directory and variables for CMAQ run
domain = 'd02'
dir_WRF = r"/home/ksz4578/Heat_Pump_Project/WRF_analyses/CMAQ_LCD/"
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
	f'{sim[0]}/wrfcheck_withstations_{sim[0]}_012016.csv', 
	f'{sim[0]}/wrfcheck_withstations_{sim[0]}_RH.csv',
	f'{sim[0]}/wrfcheck_withstations_{sim[0]}_Wind.csv',
	f'{sim[0]}/wrfcheck_withstations_{sim[0]}_WindDir.csv'
]

def stats_by_station(sim):
    
    # Process data for `d02`
    results = []
    
    for f in range(len(wrf_files_d02)):
        stn = pd.read_csv(dir_WRF + station_files[f], index_col=0)
        
        times = pd.read_csv(dir_WRF + ss[f] + '/completeddata_mini_extras2.csv', index_col=0)
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
        
        is_wind = False
        if f == 3:
            is_wind = True

        # Format time columns and align data
        times = np.array(pd.to_datetime(times['0'])[0:len(times)])
        d02['dt'] = np.array(times_perf)[0:len(d02) + 1]
        obs_d02['dt'] = times
        
        # Resample to hourly averages
        obs_d02 = obs_d02.set_index('dt').resample('h').mean()
        d02 = d02.set_index('dt')
    
        # Loop through each station and calculate statistics for each row
        for col in obs_d02.columns:
            
            # Select the observed data for this station (flatten the values for comparison)
            obs_station_data = obs_d02[col].to_numpy().flatten()
            d02_station_data = d02[col][6:].to_numpy().flatten()  # Assuming d02 data starts at column 5
            
            # Calculate statistics for the station
            mu_d, mu_p, mb, ge, rmse, r, p = stats(obs_station_data, d02_station_data, is_wind_direction= is_wind)
            
            # Store results for this station
            
            results.append({
                'Station': col,
                'Variable': var[f] + 'd02',
                'Mean Observed': mu_d,
                'Mean Predicted': mu_p,
                'Mean Bias': mb,
                'Gross Error': ge,
                'RMSE': rmse,
                'Pearson r': r,
            })
    return results

# Convert results to DataFrame for summary
#results_df_by_station_2020 = pd.DataFrame(stats_by_station(sim_2020))

# Display or save the results DataFrame
results_df_by_station = pd.DataFrame(stats_by_station(sim))


# Get average station values with latitude and longitude
sim = ['output_CONUS4K_d02_2020_RWC_4km_sf_rrtmg_10_10_1_v3852']

# Station data files for comparison
station_files = [
    f'{sim[0]}/wrfcheck_withstations_{sim[0]}_012016.csv', 
    f'{sim[0]}/wrfcheck_withstations_{sim[0]}_RH.csv',
    f'{sim[0]}/wrfcheck_withstations_{sim[0]}_Wind.csv',
    f'{sim[0]}/wrfcheck_withstations_{sim[0]}_WindDir.csv'
]

# Directory and variables for CMAQ run
dir_WRF = r"/home/ksz4578/Heat_Pump_Project/WRF_analyses/CMAQ_LCD/"
var = ['t2', 'rh', 'winds', 'winddir']

# Format file names
wrf_files_d02 = [dir_WRF + sim[i] + '/' + v + 'd02.csv' for v in var for i in range(len(sim))]

# Projection parameters and origin
import pandas as pd
from pyproj import Proj, Transformer, CRS

xorig = -2292000
yorig = -1584000
cell_size = 4000

# Define the source CRS using Lambert Conformal Conic parameters
source_crs = CRS.from_proj4('+proj=lcc +lat_1=33 +lat_2=45 +lon_0=-97 +lat_0=40 +datum=WGS84')

# Create a Transformer from the source CRS (LCC) to WGS84 (EPSG:4326)
transformer = Transformer.from_crs(source_crs, "EPSG:4326", always_xy=True)

for f in range(len(wrf_files_d02[:])):
    stn = pd.read_csv(dir_WRF + station_files[f], index_col=0)
    stn_d02 = stn[stn.in_d02 == True].reset_index()
        
    # Compute Cartesian coordinates (x, y) for each grid cell
    stn_d02['x'] = xorig + stn_d02['yy_d02'] * cell_size
    stn_d02['y'] = yorig + stn_d02['xx_d02'] * cell_size
    
    # Transform Cartesian coordinates to latitude and longitude
    stn_d02['lon'], stn_d02['lat'] = transformer.transform(stn_d02['x'], stn_d02['y'])
    
    # Drop intermediate Cartesian coordinates if not needed
    stn_d02.drop(columns=['x', 'y'], inplace=True)
    

# temp_results_2020 = results_df_by_station_2020.loc[results_df_by_station_2020['Variable'] == 't2d02']
# temp_results_2020['Mean Observed'] = (temp_results_2020['Mean Observed'] - 273.15)
# temp_results_2020.loc[:, 'lat'] = stn_d02['lat'].values
# temp_results_2020.loc[:, 'lon'] = stn_d02['lon'].values

temp_results_baseline = results_df_by_station.loc[results_df_by_station['Variable'] == 't2d02']
temp_results_baseline['Mean Observed'] = (temp_results_baseline['Mean Observed'] - 273.15)
temp_results_baseline['lat'] = stn_d02['lat']
temp_results_baseline['lon'] = stn_d02['lon']

# ws_results_2020 = results_df_by_station_2020.loc[results_df_by_station_2020['Variable'] == 'windsd02'].reset_index()
# ws_results_2020['lat'] = stn_d02['lat']
# ws_results_2020['lon'] = stn_d02['lon']

ws_results_baseline = results_df_by_station.loc[results_df_by_station['Variable'] == 'windsd02'].reset_index()
ws_results_baseline.loc[:, 'lat'] = stn_d02['lat']
ws_results_baseline['lon'] = stn_d02['lon']

# rh_results_2020 = results_df_by_station_2020.loc[results_df_by_station_2020['Variable'] == 'rhd02'].reset_index()
# rh_results_2020['lat'] = stn_d02['lat']
# rh_results_2020['lon'] = stn_d02['lon']

rh_results_baseline = results_df_by_station.loc[results_df_by_station['Variable'] == 'rhd02'].reset_index()
rh_results_baseline['lat'] = stn_d02['lat']
rh_results_baseline['lon'] = stn_d02['lon']


# wd_results_2020 = results_df_by_station_2020.loc[results_df_by_station_2020['Variable'] == 'winddird02'].reset_index()
# wd_results_2020['lat'] = stn_d02['lat']
# wd_results_2020['lon'] = stn_d02['lon']

wd_results_baseline = results_df_by_station.loc[results_df_by_station['Variable'] == 'winddird02'].reset_index()
wd_results_baseline.loc[:, 'lat'] = stn_d02['lat']
wd_results_baseline['lon'] = stn_d02['lon']



import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Function to plot validation points for a given variable
def plot_validation_points(ax, data, var, title, cmap,cbar_label, difference = False):
    resol = '50m'  # use data at this scale
    state_color = 'black'
    
    bodr = cfeature.NaturalEarthFeature(category='cultural', edgecolor='black', name='admin_0_boundary_lines_land', scale='50m', facecolor='none', alpha=0.7)
    ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
    lakes = cfeature.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
    states_provinces = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', edgecolor=state_color, facecolor='none')

    # Customize plot appearance
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    # Add features to the map
    ax.add_feature(ocean, linewidth=0.2, zorder=0)
    ax.add_feature(lakes, linewidth=0.2, zorder=0)
    ax.add_feature(states_provinces, edgecolor=state_color, linewidth=0.5, zorder=0)
    ax.add_feature(bodr, edgecolor=state_color, alpha=1)

    # Title
    ax.set_title(title) 

    if not difference:
        vmin = np.quantile(data[var], 0.025)
        vmax = np.quantile(data[var], 0.975)
    else:
        vmin = -max( [abs(np.quantile(data[var], 0.025)), abs(np.quantile(data[var], 0.975))] )
        vmax = -vmin 

    
    # Scatter plot
    sc = ax.scatter(data['lon'], data['lat'],
                    c=data[var], cmap=cmap,
                    s=100, edgecolor='black', transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, zorder=3)

    # Add colorbar
    cbar = plt.colorbar(sc, ax=ax, orientation='vertical', shrink=0.8, aspect=15)
    cbar.set_label(cbar_label, fontsize=12)
    cbar.set_ticks([vmin, (vmin + vmax) / 2, vmax])

# Function to plot multiple variables on the same figure
def plot_multiple_variables(data, variables, unit, measure):
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
    
    # Create the figure and axes (2x2 grid)
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(18, 14), subplot_kw={'projection': target_proj})
    
    # Loop through the variables and plot them
    plot_validation_points(axes[0,0] , data, variables[0], title = f"{variables[0]} Observed vs. Predicted {measure}", cmap = "inferno", cbar_label = f"RMSE {unit}")
    plot_validation_points(axes[0,1] , data, variables[1], title = f"{variables[1]} Observed vs. Predicted {measure}", cmap = "inferno_r", cbar_label = "Pearson r")
    plot_validation_points(axes[1,0] , data, variables[2], title = f"{variables[2]} Observed vs. Predicted {measure}",cmap = "RdBu",  cbar_label = f"Mean Bias {unit}", difference = True)
    plot_validation_points(axes[1,1] , data, variables[3], title = f"{variables[3]} Observed vs. Predicted {measure}", cmap = "inferno", cbar_label = f"Gross Error {unit}")

    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()

# Example usage
variables = ['RMSE', 'Pearson r', 'Mean Bias', 'Gross Error']  # Replace with your actual variable names
plot_multiple_variables(temp_results_baseline, variables, unit = "C", measure = "Temperature")  # Adjust vmin, vmax as needed
plot_multiple_variables(rh_results_baseline, variables, unit = "%", measure = "Relative Humidity")  # Adjust vmin, vmax as needed
plot_multiple_variables(ws_results_baseline, variables, unit = "mph", measure = "Wind Speed")  # Adjust vmin, vmax as needed
plot_multiple_variables(wd_results_baseline, variables, unit = "degrees", measure = "Wind Direction")  # Adjust vmin, vmax as needed



#
# ====================================================================================================================================================================================================================================================================================================================================
# Plot profiles over the simulated month

import pandas as pd
import numpy as np
import scipy.stats as st
from netCDF4 import Dataset

results = []

def get_stats_by_day(wrf_files_d02):
    for f in range(len(wrf_files_d02)):
        stn = pd.read_csv(dir_WRF + station_files[f], index_col=0)
        
        times = pd.read_csv(dir_WRF + ss[f] + '/completeddata_mini_extras2.csv', index_col=0)
        times_perf = pd.date_range(start=f'{months[f]}-01-{years[f]}', end=f'{months[f]}-{days[f]} {years[f]} 23:00:00', freq='H')
        
        # Filter and process the station data
        stn_d02 = stn[stn.in_d02 == True].reset_index()
        obs_d02 = pd.DataFrame([stn_d02[str(i)].tolist() for i in range(len(times))])
        d02 = pd.read_csv(wrf_files_d02[f], index_col=0).T
        
        # Convert units (temperature, wind speed)
        #if f == 0:
           # obs_d02 = (obs_d02 - 32) * 5 / 9 + 273.15  # Convert Fahrenheit to Kelvin
        #if f == 2:
            #d02 = 2.236936 * d02  # Convert m/s to mph
        
        # Format time columns and align data
        times = np.array(pd.to_datetime(times['0'])[0:len(times)])
        d02['dt'] = np.array(times_perf)[0:len(d02) + 1]
        d02 = d02[5:]
        obs_d02['dt'] = times
        
        # Set index for both DataFrames
        obs_d02 = obs_d02.set_index('dt')
        d02 = d02.set_index('dt')
        
        # Loop through each day in the date range
        for day in pd.date_range(start=f'{years[f]}-{months[f]}-01', end=f'{years[f]}-{months[f]}-{days[f]}'):
            # Filter data for the specific day
            obs_daily = obs_d02.loc[day.strftime('%Y-%m-%d')]
            d02_daily = d02.loc[day.strftime('%Y-%m-%d')]
        
             # Ensure both datasets are aligned
            # Remove duplicate indices
            obs_daily = obs_daily[~obs_daily.index.duplicated(keep='first')]
            d02_daily = d02_daily[~d02_daily.index.duplicated(keep='first')]
        
            # Ensure both datasets are aligned
            aligned_indices = obs_daily.index.intersection(d02_daily.index)
            obs_daily = obs_daily.loc[aligned_indices]
            d02_daily = d02_daily.loc[aligned_indices]
            
            # Flatten the 2D arrays for element-wise comparison
            obs_flat = obs_daily.to_numpy().flatten()
            d02_flat = d02_daily.to_numpy().flatten()
            
            # Skip if there are no valid observations or predictions for the day
            if np.all(np.isnan(obs_flat)) or np.all(np.isnan(d02_flat)):
                continue
            
            var_name = f"{var[f]}d02 - {day.strftime('%Y-%m-%d')}"
            mu_d, mu_p, mb, ge, rmse, r, p = stats(obs_flat, d02_flat)
            
            # Store results
            results.append({
                'Variable': var_name,
                'Mean Observed': mu_d,
                'Mean Predicted': mu_p,
                'Mean Bias': mb,
                'Gross Error': ge,
                'RMSE': rmse,
                'Pearson r': r,
                'p-value': p
            })
    return results

results =  get_stats_by_day(wrf_files_d02)

# Convert results to DataFrame for summary
results_df = pd.DataFrame(results)


import matplotlib.pyplot as plt
import pandas as pd

# Assuming `results_df` is already prepared with 'Date' as the index
results_df['Date'] = pd.to_datetime(results_df['Variable'].str.split(' - ').str[1])  # Extract the date from Variable
results_df.set_index('Date', inplace=True)

# Define variable-specific y-axis labels with units
variable_ylabels = {
    "t2d02": "Temperature (°C)",
    "rhd02": "Relative Humidity (%)",
    "wsd02": "Wind Speed (m/s)",
    "wdd02": "Wind Direction (degrees)"
}

# Extract unique variable categories
variable_categories = results_df['Variable'].str.split(' - ').str[0].unique()

def plot_four_panel_diff():
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))  # Create a 2x2 grid of subplots
    axes = axes.flatten()  # Flatten the 2D array of axes for easier iteration

    for i, var in enumerate(variable_categories):
        if i >= 4:  # Ensure only four panels are plotted
            break
        
        # Filter results for the current variable category
        var_results = results_df[results_df['Variable'].str.startswith(var)]
        
        if "t2" in var:
            var_results["Mean Observed"] -= 273.15
            var_results["Mean Predicted"] -= 273.15
            
        # Plot Mean Observed and Mean Predicted
        axes[i].plot(var_results.index, var_results["Mean Observed"], label="Mean Observed", linestyle='--', marker='s', color='blue')
        axes[i].plot(var_results.index, var_results["Mean Predicted"], label="Mean Predicted", linestyle='--', marker='o', color='orange')
        
        # Set title and axis labels
        axes[i].set_title(var)
        axes[i].set_xlabel("Date")
        ylabel = variable_ylabels.get(var, "Values")  # Default to "Values" if the variable isn't in the dictionary
        axes[i].set_ylabel(ylabel)
        axes[i].tick_params(axis='x', rotation=45)
        axes[i].grid(True)
        axes[i].legend()
    
    # Remove unused subplots if there are fewer than 4 variables
    for j in range(i + 1, 4):
        fig.delaxes(axes[j])
    
    # Adjust layout
    fig.suptitle("Observed v. Predicted Meteorology Profiles Over Time")
    plt.tight_layout()
    plt.show()
    

plot_four_panel_diff()


#
# ====================================================================================================================================================================================================================================================================================================================================
# Plot diurnal profiles


def get_average_diurnal_profile(wrf_files_d02):
    diurnal_profiles = []

    for f in range(len(wrf_files_d02)):
        # Load station data and WRF output
        stn = pd.read_csv(dir_WRF + station_files[f], index_col=0)
        times = pd.read_csv(dir_WRF + ss[f] + '/completeddata_mini_extras2.csv', index_col=0)
        times_perf = pd.date_range(start=f'{months[f]}-01-{years[f]}', end=f'{months[f]}-{days[f]} {years[f]} 23:00:00', freq='H')

        # Process station data
        stn_d02 = stn[stn.in_d02 == True].reset_index()
        obs_d02 = pd.DataFrame([stn_d02[str(i)].tolist() for i in range(len(times))])
        d02 = pd.read_csv(wrf_files_d02[f], index_col=0).T

        # # Unit conversions if necessary
        # if f == 0:
        #     obs_d02 = (obs_d02 - 32) * 5 / 9 + 273.15  # Fahrenheit to Kelvin
        # if f == 2:
        #     d02 = 2.236936 * d02  # m/s to mph

        # Format time columns and align data
        times = np.array(pd.to_datetime(times['0'])[0:len(times)])
        d02['dt'] = np.array(times_perf)[0:len(d02) + 1]
        d02 = d02[5:]
        obs_d02['dt'] = times

        # Set index for both DataFrames
        obs_d02 = obs_d02.set_index('dt')
        d02 = d02.set_index('dt')

        # Add a column for the hour of the day
        obs_d02['hour'] = obs_d02.index.hour
        d02['hour'] = d02.index.hour

        # Average across stations and compute diurnal profiles
        obs_diurnal_mean = obs_d02.groupby('hour').mean().mean(axis=1)
        d02_diurnal_mean = d02.groupby('hour').mean().mean(axis=1)

        # Combine observed and predicted averages into a DataFrame
        diurnal_profile = pd.DataFrame({
            f'{var[f]}_Observed': obs_diurnal_mean,
            f'{var[f]}_Predicted': d02_diurnal_mean
        })

        # Append to the list
        diurnal_profiles.append(diurnal_profile)

    # Concatenate all profiles into a single DataFrame
    combined_profiles = pd.concat(diurnal_profiles, axis=1)
    return combined_profiles

diurnal_profiles = get_average_diurnal_profile(wrf_files_d02)
import matplotlib.pyplot as plt

def plot_four_panel_diurnal_profiles(diurnal_profiles, variable_units):
    """
    Create a 4-panel plot for average diurnal profiles of observed and predicted values using Matplotlib.

    Args:
        diurnal_profiles (pd.DataFrame): DataFrame with rows as hours and columns for observed and predicted values.
        variable_units (dict): Dictionary mapping variable names to their units (e.g., {"Temperature": "C", ...}).
    """
    # Extract variable names from DataFrame
    variables = [col.split('_')[0] for col in diurnal_profiles.columns[::2]]

    # Create a 2x2 subplot layout
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    axes = axes.flatten()  # Flatten to iterate through a single list of axes

    for i, var in enumerate(variables):
        if i >= 4:  # Ensure only four panels are plotted
            break

        # Extract observed and predicted values for the current variable
        observed = diurnal_profiles[f"{var}_Observed"]
        predicted = diurnal_profiles[f"{var}_Predicted"]

        # Plot observed and predicted values
        axes[i].plot(
            diurnal_profiles.index, observed, label="Observed", linestyle='--', marker='o', color='blue'
        )
        axes[i].plot(
            diurnal_profiles.index, predicted, label="Predicted", linestyle='--', marker='s', color='orange'
        )

        # Set title, labels, and grid
        axes[i].set_title(f"{var} Diurnal Profile")
        axes[i].set_xlabel("Hour of the Day")
        ylabel = f"{var} ({variable_units.get(var, 'Values')})"
        axes[i].set_ylabel(ylabel)
        axes[i].grid(True)
        axes[i].legend()
        axes[i].tick_params(axis='x', rotation=45)

    # Remove unused subplots if there are fewer than 4 variables
    for j in range(i + 1, 4):
        fig.delaxes(axes[j])

    # Adjust layout and add a global title
    fig.suptitle("Average Diurnal Profiles of Observed and Predicted Values")
    plt.tight_layout()
    plt.show()

# Assuming diurnal_profiles is generated
variable_units = {
    "t2": "C",
    "winds": "m/s",
    "rh": "%",
    "winddir": "degrees"
}

plot_four_panel_diurnal_profiles(diurnal_profiles, variable_units)