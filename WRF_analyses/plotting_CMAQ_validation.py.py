"""
This script performs comprehensive statistical evaluation of WRF-CMAQ model outputs against EPA Air Quality System (AQS)
monitoring data. The toolkit provides:

1. Spatial analysis of model performance across CONUS
2. Temporal evaluation (diurnal cycles, daily trends)
3. Statistical metrics calculation (NMB, NME, correlation)
4. Geographic visualization of monitoring sites and model performance

Key Functionalities:
-------------------
1. Data Processing:
   - EPA AQS data ingestion and cleaning
   - Model-observation pairing
   - Quality control (detection limits, missing data)
   - Temporal aggregation (hourly, daily)
   - Geographic coordinate handling

2. Statistical Analysis:
   - Standard metrics: NMB, NME, MB, RMSE
   - Pearson correlation (linear and circular for wind direction)
   - Time-of-day analysis
   - Geographic stratification (latitudinal swaths, timezones)

3. Visualization:
   - Spatial distribution of monitoring sites
   - Model performance maps (NMB, NME, correlation)
   - Diurnal profile comparisons
   - Time series analysis

Data Requirements:
-----------------
Inputs:
- EPA AQS Data:
  • Hourly PM2.5 measurements (88101)
  • Station metadata (lat/lon, site info)
  
- Model Outputs:
  • WRF-CMAQ predictions (NetCDF/CSV)
  • Geographic projection info

Outputs:
--------
1. Statistical Reports:
   - Site-by-site performance metrics
   - Regional aggregations
   - Temporal analysis results

2. Visualization Products:
   - Publication-quality maps
   - Time series plots
   - Diagnostic figures

Dependencies:
------------
Core:
- pandas >= 1.3.0
- numpy >= 1.21.0
- scipy >= 1.7.0

Geospatial:
- cartopy >= 0.20.0
- geopandas >= 0.10.0
- pyproj >= 3.0.0

Visualization:
- matplotlib >= 3.4.0
- plotly >= 5.0.0
"""
# %% [markdown]
# ### Plot EPA AQS Sites

# %%
import pandas as pd
pm25 = pd.read_csv('AQS_data_2016/hourly_88101_2016.csv')


# %%
import pandas as pd
import plotly.express as px

# Directory containing the EPA AQS data files
dir_epa = 'AQS_data_2016/'

# EPA codes and corresponding pollutants
epa_code = ['42401', '42602', '44201', '42101', '88101', '81102']
pollutants = ['SO2', 'NO2', 'O3', 'CO', 'PM25_TOT', 'PM10']

# Initialize an empty DataFrame to store all unique site data
all_sites_df = pd.DataFrame()

# Loop through each pollutant and read only Latitude and Longitude from each large data file
for i, code in enumerate(epa_code):
    file_path = f"{dir_epa}hourly_{code}_2016.csv"
    
    # Read only the Latitude and Longitude columns
    df = pd.read_csv(file_path, usecols=['Latitude', 'Longitude'])
    
    # Keep only the unique pairs of latitude and longitude
    df_unique = df.drop_duplicates()
    
    # Add a column to indicate the pollutant
    df_unique['Pollutant'] = pollutants[i]
    
    # Append the unique site data for this pollutant to the all_sites_df DataFrame
    all_sites_df = pd.concat([all_sites_df, df_unique], ignore_index=True)

# %%
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import geopandas as gpd
import cartopy.io.shapereader as shpreader

# Convert to GeoDataFrame
gdf = gpd.GeoDataFrame(all_sites_df, geometry=gpd.points_from_xy(all_sites_df.Longitude, all_sites_df.Latitude), crs="EPSG:4326")
conus_gdf = gdf[(gdf['Latitude'].between(24.5, 49.5)) & (gdf['Longitude'].between(-125, -66))]

# Generate map features
resol = '50m'
ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
lakes = cfeature.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])


# Define map projection
projection = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': projection})

# Customize plot appearance
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

# # Add features to the map
ax.add_feature(ocean, linewidth=0, zorder=0, facecolor="lightblue")
ax.add_feature(lakes, linewidth=0, zorder=0, facecolor="lightblue")

shapefile = shpreader.natural_earth(resolution='50m', category='cultural', name='admin_0_countries')
countries = gpd.read_file(shapefile)
mexico_canada = countries[countries['ADMIN'].isin(['Canada', 'Mexico'])]
ax.add_geometries(mexico_canada.geometry, crs=ccrs.PlateCarree(), facecolor="gray", edgecolor="none", zorder=-1)


# Define color mapping for pollutants
pollutants = conus_gdf['Pollutant'].unique()
colors = plt.cm.get_cmap('Set1', len(pollutants))
pollutant_color_map = {pollutant: colors(i) for i, pollutant in enumerate(pollutants)}

# Plot data points
for pollutant in ["PM25_TOT"]:
    subset = conus_gdf[conus_gdf['Pollutant'] == pollutant]
    ax.scatter(subset.geometry.x, subset.geometry.y, transform=ccrs.PlateCarree(),
               label=pollutant, color=pollutant_color_map[pollutant], edgecolors='k', alpha=0.7, s=50, zorder = 1)

# Add legend and title
plt.legend(title="Pollutants", loc='upper right', fontsize=10)
plt.title("EPA AQS Monitoring Sites for PM₂.₅ (January 2016)", fontsize=14)

# Show plot
plt.show()


# %% [markdown]
# ### Statistics Baseline No Neighbors

# %%
import pandas as pd
import numpy as np
#from scipy.stats import pearsonr

# Load the CSV file into a DataFrame
file_path = 'AQS_data_2016/PM25_TOT_d02_WHbase_2016_1_EPA_CMAQ_Combine.csv'
df = pd.read_csv(file_path)

# # Remove any rows with missing data in 'Sample Measurement' or 'CMAQ'
# before = len(df.index)
# df = df.dropna(subset=['Sample Measurement', 'CMAQ'])
# print(round((before - len(df.index))/ before * 100, 2), "% lost to na")

# %%
percent_under_detection_limit = ((df.loc[df['Sample Measurement'] <= 2]).shape[0])/ df.shape[0] * 100 
print(f"{round(percent_under_detection_limit,2)}% of AQS measurements under the detection limit of 2ug/m^3")


df_detection_limit = df.loc[df['Sample Measurement'] >= 2]

# %%
def pearsonr(x, y):
    """
    Compute the Pearson correlation coefficient between two arrays.

    Parameters:
    x (array-like): First dataset.
    y (array-like): Second dataset.

    Returns:
    r (float): Pearson correlation coefficient.
    p_value (float): Placeholder for compatibility (always None in this implementation).
    """
    x = np.asarray(x)
    y = np.asarray(y)

    if x.shape[0] != y.shape[0]:
        raise ValueError("Input arrays must have the same length.")

    # Compute means
    x_mean = np.mean(x)
    y_mean = np.mean(y)

    # Compute Pearson correlation
    numerator = np.sum((x - x_mean) * (y - y_mean))
    denominator = np.sqrt(np.sum((x - x_mean) ** 2) * np.sum((y - y_mean) ** 2))
    r = numerator / denominator if denominator != 0 else 0  # Avoid division by zero

    return r, None  # Returning None as a placeholder for p-value

# %%
def do_overall_stats(df, daily):
    df = df.copy(deep = True)
    if daily:
        grouped = df.groupby(['Latitude', 'Longitude'])
        groups = []
        for (lat, lon), group in grouped:
            group = group.set_index('dt')
            group.index = pd.to_datetime(group.index)

            numeric_cols = group.select_dtypes(include='number').columns
            non_numeric_cols = group.select_dtypes(exclude='number').columns
            numeric_resampled = group[['Sample Measurement', 'CMAQ']].resample('D').mean().dropna()
            non_numeric_resampled = group[non_numeric_cols].resample('D').last()

            group = pd.concat([numeric_resampled, non_numeric_resampled], axis=1).reset_index().dropna(subset = ['Sample Measurement', 'CMAQ'])
            group['Latitude'] = lat
            group['Longitude'] = lon
            groups.append(group)

        df = pd.concat(groups)
        df = pd.concat(groups).reset_index()
    
    df_all = df[['CMAQ', 'Sample Measurement']].dropna()
    observed = df_all['Sample Measurement']
    predicted = df_all['CMAQ']

    # Calculate mean observed (μd) and mean predicted (μp)
    μd = observed.mean()
    μp = predicted.mean()

    # Calculate NMB
    nmb = ((predicted - observed).sum() / observed.sum()) * 100

    # Calculate NME
    nme = (abs(predicted - observed).sum() / observed.sum()) * 100

    # Calculate Pearson correlation coefficient (r)
    if len(observed) > 1:  # Ensure there is enough data for correlation
        r_value, _ = pearsonr(observed, predicted)
    else:
        r_value = np.nan  # Not enough data to compute correlation

    # Append the calculated metrics to the DataFrame
    return pd.DataFrame({
        'μd': [μd],
        'μp': [μp],
        'NMB': [nmb],
        'NME': [nme],
        'r': [r_value]
    })

do_overall_stats(df, daily = True)

# %%
file_path = 'AQS_data_2016/PM25_TOT_d02_WHbase_2016_1_EPA_CMAQ_Combine_2020_RWC.csv'
df_2020 = pd.read_csv(file_path)

do_overall_stats(df_2020, daily = True)

# %%
def do_stats(df, daily = False):
    df_new = df.copy(deep = True)

    # Group by unique (Latitude, Longitude) pairs
    grouped = df_new.groupby(['Latitude', 'Longitude'])
    
    metrics_dfs = []
    groups = []

    # Iterate through each group and calculate the metrics
    for (lat, lon), group in grouped:
        if daily:
            group_before = group
            group = group.set_index('dt')
            group.index = pd.to_datetime(group.index)

            numeric_cols = group.select_dtypes(include='number').columns
            non_numeric_cols = group.select_dtypes(exclude='number').columns
            numeric_resampled = group[['Sample Measurement', 'CMAQ']].resample('D').mean().dropna()
            non_numeric_resampled = group[non_numeric_cols].resample('D').last()

            group = pd.concat([numeric_resampled, non_numeric_resampled], axis=1).reset_index().dropna(subset = ['Sample Measurement', 'CMAQ'])
            group['Latitude'] = lat
            group['Longitude'] = lon
            groups.append(group)
        
        observed = group['Sample Measurement']
        predicted = group['CMAQ']

        # Calculate mean observed (μd) and mean predicted (μp)
        μd = observed.mean()
        μp = predicted.mean()

        # Calculate NMB
        if observed.sum() != 0:
            nmb = ((predicted - observed).sum() / observed.sum()) * 100
    
            # Calculate NME
            nme = (abs(predicted - observed).sum() / observed.sum()) * 100
        else:
            print("observed = 0")
            nmb = 0
            nme = 0
        
        # Calculate Pearson correlation coefficient (r)
        if len(observed) > 1:  # Ensure there is enough data for correlation
            r_value, _ = pearsonr(observed, predicted)
        else:
            r_value = np.nan  # Not enough data to compute correlation

        # Append the calculated metrics to the DataFrame
        metrics_dfs.append(pd.DataFrame({
            'Latitude': [lat],
            'Longitude': [lon],
            'μd': [μd],
            'μp': [μp],
            'NMB': [nmb],
            'NME': [nme],
            'r': [r_value]
        }))
    if daily:
        return pd.concat(metrics_dfs), pd.concat(groups)
    else:
        return pd.concat(metrics_dfs)

metrics_df = do_stats(df, daily = False)
metrics_df_2020 = do_stats(df_2020, daily = False)

metrics_df_daily, df_daily  = do_stats(df, daily = True)
metrics_df_daily_2020, df_daily_2020  = do_stats(df_2020, daily = True)
metrics_df_detection_limit = do_stats(df_detection_limit)
metrics_df_detection_limit_daily = do_stats(df_detection_limit, daily = True)

# Display the metrics DataFrame

# %%
df_daily_2020 = pd.read_csv("AQS_data_2016/PM25_TOT_d02_WHbase_2016_1_EPA_CMAQ_Combine_2020_RWC_daily.csv")

# %%
# Group by 'State Name' and 'County Name' and take the first Latitude and Longitude
state_county_map = df_2020.groupby(['State Name', 'County Name'])[['Latitude', 'Longitude']].first().to_dict('index')

# Function to retrieve latitude and longitude for each row in df_daily_2020
def get_lat_lon(row):
    state_county_key = (row['State Name'], row['County Name'])
    lat_lon = state_county_map.get(state_county_key, {'Latitude': None, 'Longitude': None})
    return lat_lon['Latitude'], lat_lon['Longitude']

# Apply the function to add latitude and longitude
df_daily_2020[['Latitude', 'Longitude']] = df_daily_2020.apply(get_lat_lon, axis=1, result_type='expand')


# %%
metrics_df_daily.to_csv("AQS_data_2016/validation_baseline_daily.csv")
metrics_df_detection_limit_daily.to_csv("AQS_data_2016/validation_baseline_daily_detection_limit.csv")
metrics_df.to_csv("AQS_data_2016/validation_baseline_hourly.csv")
metrics_df_detection_limit.to_csv("AQS_data_2016/validation_baseline_hourly_detection_limit.csv")

# %% [markdown]
# ### Daytime only

# %%
df['dt'] = pd.to_datetime(df['dt'])
df['hour'] = df['dt'].dt.hour

# %%
file_path = 'AQS_data_2016/PM25_TOT_d02_WHbase_2016_1_EPA_CMAQ_Combine_2020_RWC.csv'
df_2020 = pd.read_csv(file_path)

df_2020['dt'] = pd.to_datetime(df_2020['dt'])
df_2020['hour'] = df_2020['dt'].dt.hour

# %%
for i in range(0, 12):
    day_bound = i
    night_bound = i+11
    df_day = df.loc[(df['hour'] >= day_bound) & (df['hour'] <= night_bound)]
    df_night= df.loc[(df['hour'] < day_bound) | (df['hour'] > night_bound)]
    do_overall_stats(df_day,  daily = True)
    x = do_overall_stats(df_day,  daily = True)
    print(i, round(list(x['NMB'])[0], 2))

# %%
df_2020['dt'] = pd.to_datetime(df_2020['dt'])
df_2020['hour'] = df_2020['dt'].dt.hour
day_bound = 6
night_bound = 17

df_day_2020 = df_2020.loc[(df_2020['hour'] >= day_bound) & (df_2020['hour'] <= night_bound)]
df_night_2020 = df_2020.loc[(df_2020['hour'] < day_bound) | (df_2020['hour'] > night_bound)]
df_day = df.loc[(df['hour'] >= day_bound) & (df['hour'] <= night_bound)]
df_night= df.loc[(df['hour'] < day_bound) | (df['hour'] > night_bound)]

# %%
do_overall_stats(df, daily = True)

# %%
do_overall_stats(df_2020,  daily = True)

# %%
do_overall_stats(df_day,  daily = True)

# %%
do_overall_stats(df_day_2020,  daily = True)

# %%
do_overall_stats(df_night,  daily = True)

# %%
do_overall_stats(df_night_2020,  daily = True)

# %%
out, _ = do_stats(df, daily = True)
out.mean()

# %%
out, _ = do_stats(df_day, daily = True)
out.mean()

# %%
out, _ = do_stats(df_night, daily = True)
out.mean()

# %% [markdown]
# ### timezone validation

# %%
import pandas as pd
import numpy as np
from timezonefinder import TimezoneFinder

# Initialize TimezoneFinder instance
tf = TimezoneFinder()

def get_timezone_from_coordinates(lat, lon):
    """Get timezone name from latitude and longitude."""
    return tf.timezone_at(lat=lat, lng=lon)

def map_to_standard_timezone(lat, lon):
    """Map detailed timezones to standard time zones."""
    timezone = get_timezone_from_coordinates(lat, lon)
    pst = ['America/Los_Angeles']
    mst = ['America/Denver', 'America/Boise', 'America/Phoenix']  # Exclude 'America/Phoenix'
    cst = ['America/Chicago', 'America/North_Dakota/Beulah', 'America/North_Dakota/Center']
    est = ['America/New_York', 'America/Detroit', 'America/Indiana/Indianapolis', 'America/Kentucky/Louisville']
    
    if timezone in pst:
        return 'PST'
    elif timezone in mst:
        return 'MST'
    elif timezone in cst:
        return 'CST'
    elif timezone in est:
        return 'EST'
    elif timezone in others:
        return 'Other (MST no DST)'
    else:
        return 'Other'

def do_overall_stats_per_timezone(df):
    # Ensure necessary columns exist
    if 'Latitude' not in df.columns or 'Longitude' not in df.columns:
        raise ValueError("The dataframe must contain 'Latitude' and 'Longitude' columns.")
    
    df = df.dropna(subset = ['Latitude', 'Longitude'])
    
    # Get timezone for each row based on Latitude and Longitude
    df['timezone'] = df.apply(lambda row: map_to_standard_timezone(row['Latitude'], row['Longitude']), axis=1)

    # Handle missing or invalid timezones
    df['timezone'].fillna('Unknown', inplace=True)

    # Group the data by timezone
    grouped = df.groupby('timezone')

    # Define a function to compute stats for each group
    def compute_stats(group):
        df_all = group[['CMAQ', 'Sample Measurement']].dropna()
        observed = df_all['Sample Measurement']
        predicted = df_all['CMAQ']
        
        if len(observed) == 0:  # Handle case where no data is available after dropping NaNs
            return pd.Series({'μd': np.nan, 'μp': np.nan, 'NMB': np.nan, 'NME': np.nan, 'r': np.nan})

        # Calculate mean observed (μd) and mean predicted (μp)
        μd = observed.mean()
        μp = predicted.mean()

        # Calculate NMB
        nmb = ((predicted - observed).sum() / observed.sum()) * 100

        # Calculate NME
        nme = (abs(predicted - observed).sum() / observed.sum()) * 100
        mb = predicted.mean() - observed.mean()

        # Calculate Pearson correlation coefficient (r)
        if len(observed) > 1:  # Ensure there is enough data for correlation
            r_value, _ = pearsonr(observed, predicted)
        else:
            r_value = np.nan  # Not enough data to compute correlation

        return pd.Series({'μd': μd, 'μp': μp, 'NMB': nmb, 'NME': nme, 'mb': mb,'r': r_value})
        #return pd.Series({'mb': mb,'r': r_value})

    # Apply the stats calculation function to each group and return the results as a DataFrame
    return grouped.apply(compute_stats).reset_index()

# Example usage:
overall_stats_per_timezone = do_overall_stats_per_timezone(df_daily_2020)
print("RWC 2020:")
overall_stats_per_timezone

# %%
df_daily_2020

# %% [markdown]
# ### Latidudinal Swathe

# %%
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

def get_latitudinal_swath(lat):
    """Map latitude to a specific 5-degree swath."""
    # Define the swath range
    if lat < 24.3963 or lat > 49.3845:  # Check if lat is within the continental US
        return 'Outside US'
    elif lat < 35:  # Group all latitudes under 35 into "<35"
        return '<35'
    swath = int(lat // 5) * 5  # Calculate the swath starting point

    return f"{swath}° to {swath + 5}°"

def do_overall_stats_per_latitudinal_swath(df):
    # Ensure necessary columns exist
    if 'Latitude' not in df.columns or 'Longitude' not in df.columns:
        raise ValueError("The dataframe must contain 'Latitude' and 'Longitude' columns.")
    
    df = df.dropna(subset=['Latitude', 'Longitude'])
    
    # Get swath for each row based on Latitude
    df['lat_swath'] = df['Latitude'].apply(get_latitudinal_swath)

    # Handle missing or invalid swaths
    df['lat_swath'].fillna('Unknown', inplace=True)

    # Group the data by latitudinal swath
    grouped = df.groupby('lat_swath')

    # Define a function to compute stats for each group
    def compute_stats(group):
        df_all = group[['CMAQ', 'Sample Measurement']].dropna()
        observed = df_all['Sample Measurement']
        predicted = df_all['CMAQ']
        
        if len(observed) == 0:  # Handle case where no data is available after dropping NaNs
            return pd.Series({'μᴅ': np.nan, 'μᵐ': np.nan, 'NMB (%)': np.nan, 'NME (%)': np.nan, 'r': np.nan, 'n': 0})

        # Calculate means
        mu_d = observed.mean()
        mu_p = predicted.mean()

        # Calculate Normalized Mean Bias (NMB)
        nmb = 100 * (predicted.sum() - observed.sum()) / observed.sum()

        # Calculate Normalized Mean Error (NME)
        nme = 100 * np.abs(predicted - observed).sum() / observed.sum()

        # Calculate Pearson correlation coefficient (r)
        if len(observed) > 1:  # Ensure there is enough data for correlation
            r_value, _ = pearsonr(observed, predicted)
        else:
            r_value = np.nan  # Not enough data to compute correlation
        
        # Count number of valid sensors (data points)
        n = len(group['Latitude'].unique())

        return pd.Series({'μᴅ': mu_d, 'μᵐ': mu_p, 'NMB (%)': nmb, 'NME (%)': nme, 'r': r_value, 'n': n})

    # Apply the stats calculation function to each group and return the results as a DataFrame
    return grouped.apply(compute_stats).reset_index()


# Example usage:
# Assuming df_daily_2020 is your DataFrame containing the necessary data
overall_stats_per_latitudinal_swath = do_overall_stats_per_latitudinal_swath(df_daily)
print("Overall Statistics by 5-Degree Latitudinal Swath:")
print(overall_stats_per_latitudinal_swath)


# %%
overall_stats_per_latitudinal_swath = do_overall_stats_per_latitudinal_swath(df_daily_2020)
overall_stats_per_latitudinal_swath

# %% [markdown]
# ### Validation by time of day

# %%
df

# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Convert 'dt' to datetime if not already
df['dt'] = pd.to_datetime(df['dt'])

df  = df.dropna(subset=['Sample Measurement', 'CMAQ'])


# Extract hour of the day
df['hour'] = df['dt'].dt.hour

# Define a function to calculate MB and correlation
def calculate_metrics(group):
    mb = np.mean(group['CMAQ'] - group['Sample Measurement'])  # Mean Bias
    corr, _ = st.pearsonr(group['CMAQ'], group['Sample Measurement'])  # Pearson correlation
    return pd.Series({'MB': mb, 'Correlation': corr})

# Group by hour and calculate metrics
hourly_metrics = df.groupby('hour').apply(calculate_metrics).reset_index()

# Plot MB and Correlation against hour of the day
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot Mean Bias
ax1.set_xlabel('Hour of the Day')
ax1.set_ylabel('Mean Bias (MB)', color='tab:blue')
ax1.plot(hourly_metrics['hour'], hourly_metrics['MB'], label='Mean Bias', color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue')

# Create a secondary y-axis for Correlation
ax2 = ax1.twinx()
ax2.set_ylabel('Correlation', color='tab:red')
ax2.plot(hourly_metrics['hour'], hourly_metrics['Correlation'], label='Correlation', color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')

# Titles and Legends
plt.title('Mean Bias (MB) and Correlation Grouped by Hour of the Day')
fig.tight_layout()
plt.show()


# %% [markdown]
# ## Statistics 2020 RWC

# %%
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

# Load the CSV file into a DataFrame
file_path = '/home/ksz4578/Heat_Pump_Project/WRF_analyses/AQS_data_2016/PM25_TOT_d02_WHbase_2016_1_EPA_CMAQ_Combine_2020_RWC.csv'
df = pd.read_csv(file_path)

# # Remove any rows with missing data in 'Sample Measurement' or 'CMAQ'
# before = len(df.index)
# df = df.dropna(subset=['Sample Measurement', 'CMAQ'])
# print(round((before - len(df.index))/ before * 100, 2), "% lost to na")

# %%
percent_under_detection_limit = ((df.loc[df['Sample Measurement'] <= 2]).shape[0])/ df.shape[0] * 100 
print(f"{round(percent_under_detection_limit,2)}% of AQS measurements under the detection limit of 2ug/m^3")


df_detection_limit = df.loc[df['Sample Measurement'] >= 2]

# %%
metrics_df_2020 = do_stats(df, daily = False)
metrics_df_daily_2020, df_daily_2020  = do_stats(df, daily = True)
metrics_df_detection_limit_2020 = do_stats(df_detection_limit)
metrics_df_detection_limit_daily_2020 = do_stats(df_detection_limit, daily = True)

# Display the metrics DataFrame

# %%
metrics_df_daily_2020.to_csv("AQS_data_2016/validation_RWC2020_daily.csv")
metrics_df_detection_limit_daily_2020.to_csv("AQS_data_2016/validation_RWC2020_daily_detection_limit.csv")
metrics_df_2020.to_csv("AQS_data_2016/validation_RWC2020_hourly.csv")
metrics_df_detection_limit_2020.to_csv("AQS_data_2016/validation_RWC2020_hourly_detection_limit.csv")

# %%
metrics_df_detection_limit_daily_2020.mean()

# %% [markdown]
# ## Statistics Baseline Neighbors

# %%
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

# Load the CSV file into a DataFrame
file_path = '/home/ksz4578/Heat_Pump_Project/WRF_analyses/AQS_data_2016/PM25_TOT_d02_WHbase_2016_1_EPA_CMAQ_Combine_neighbors.csv'
df = pd.read_csv(file_path)

# Remove any rows with missing data in 'Sample Measurement' or 'CMAQ'
before = len(df.index)
df = df.dropna(subset=['Sample Measurement', 'CMAQ'])
print(round((before - len(df.index))/ before * 100, 2), "% lost to na")

# %%
percent_under_detection_limit = ((df.loc[df['Sample Measurement'] <= 2]).shape[0])/ df.shape[0] * 100 
print(f"{round(percent_under_detection_limit,2)}% of AQS measurements under the detection limit of 2ug/m^3")

df_detection_limit = df.loc[df['Sample Measurement'] >= 2]

# %%
def do_stats(df):
# Initialize an empty DataFrame to store the metrics for each site
    metrics_df = pd.DataFrame(columns=['Latitude', 'Longitude', 'μd', 'μp', 'NMB', 'NME', 'r'])

    # Group by unique (Latitude, Longitude) pairs
    grouped = df.groupby(['Latitude', 'Longitude'])

    # Iterate through each group and calculate the metrics
    for (lat, lon), group in grouped:
        observed = group['Sample Measurement']
        predicted = group['CMAQ']

        # Calculate mean observed (μd) and mean predicted (μp)
        μd = observed.mean()
        μp = predicted.mean()

        # Calculate NMB
        nmb = ((predicted - observed).sum() / observed.sum()) * 100

        # Calculate NME
        nme = (abs(predicted - observed).sum() / observed.sum()) * 100

        # Calculate Pearson correlation coefficient (r)
        if len(observed) > 1:  # Ensure there is enough data for correlation
            r_value, _ = pearsonr(observed, predicted)
        else:
            r_value = np.nan  # Not enough data to compute correlation

        # Append the calculated metrics to the DataFrame
        metrics_df = metrics_df.append({
            'Latitude': lat,
            'Longitude': lon,
            'μd': μd,
            'μp': μp,
            'NMB': nmb,
            'NME': nme,
            'r': r_value
        }, ignore_index=True)
    return metrics_df
    
metrics_df_neighbors = do_stats(df)
metrics_df_detection_limit = do_stats(df_detection_limit)

# Display the metrics DataFrame

# %%
import seaborn as sns
import matplotlib.pyplot as plt

# Add a column to distinguish between the two DataFrames
metrics_df['Dataset'] = 'Baseline'
metrics_df_2020['Dataset'] = '2020 RWC'

# Concatenate the two DataFrames
df_combined = pd.concat([metrics_df[['r', 'Dataset']], metrics_df_2020[['r', 'Dataset']]])

# Create the box plot
plt.figure(figsize=(8, 6))
sns.boxplot(x='Dataset', y='r', data=df_combined)
plt.title('Box Plot of R')
plt.ylabel('R')
plt.show()

# %%
metrics_df_detection_limit.mean()

# %% [markdown]
# ## Plotting Validation

# %%
metrics_df_2020['overestimation'] = metrics_df_2020['μp'] - metrics_df_2020['μd']
metrics_df['overestimation'] = metrics_df['μp'] - metrics_df['μd']

# %%
metrics_df_2020

# %%
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import pandas as pd

def plot_chem_spatial(chosen_stat="NMB", cmap="RdBu", difference=True, stat_desc="NMB (%)", alpha=0.075):
    # Convert to GeoDataFrame
    gdf = gpd.GeoDataFrame(metrics_df_2020, geometry=gpd.points_from_xy(metrics_df_2020.Longitude, metrics_df_2020.Latitude), crs="EPSG:4326")

    # Define map projection
    projection = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(20, 12), subplot_kw={'projection': projection})

    # Add map features
    ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray', zorder=-1)
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue', zorder=-1)
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='lightblue', zorder=-1)
    ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='black', zorder=-1)

    # Normalize color scale
    if difference:
        vMAX = max(abs(gdf[chosen_stat].quantile(alpha)), abs(gdf[chosen_stat].quantile(1-alpha)))

        norm = Normalize(vmin=-vMAX, vmax=vMAX)
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        # Scatter plot with color mapping
        sc = ax.scatter(gdf.geometry.x, gdf.geometry.y, c=gdf[chosen_stat], cmap=cmap, norm=norm,
                        transform=ccrs.PlateCarree(), edgecolors='k', linewidths=0.5, alpha=1, s=80)
    else:
        vmin = gdf[chosen_stat].quantile(alpha)
        vmax = gdf[chosen_stat].quantile(1-alpha)

        norm = Normalize(vmin=vmin, vmax=vmax)
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        # Scatter plot with color mapping
        sc = ax.scatter(gdf.geometry.x, gdf.geometry.y, c=gdf[chosen_stat], cmap=cmap, norm=norm,
                        transform=ccrs.PlateCarree(), edgecolors='k', linewidths=0.5, alpha=1, s=80)

    # Add colorbar with adjusted font sizes
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', shrink=0.8)
    cbar.set_label(stat_desc, fontsize=18)
    cbar.ax.tick_params(labelsize=16)  # Adjust tick font size

    # Add title with font size
    plt.title(f"Spatial Distribution of {stat_desc.split('(')[0]}", fontsize=22)

    # Set map extent
    ax.set_extent([-120, -72, 24.5, 49.5], crs=ccrs.PlateCarree())

    # Show plot
    plt.show()

plot_chem_spatial()


# %%
metrics_df_2020

# %%
plot_chem_spatial(stat_desc = "Normalized Mean Bias (%)")
plot_chem_spatial(chosen_stat="overestimation", cmap="RdBu", difference=True, stat_desc="Mean Bias (µg/m³)", alpha=0.075)
plot_chem_spatial(chosen_stat="NME", cmap="magma_r", difference=False, stat_desc="Normalized Mean Error (%)", alpha=0.075)
plot_chem_spatial(chosen_stat="r", cmap="magma_r", difference=False, stat_desc="Pearson r", alpha=0.025)

# %%
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import pandas as pd

def plot_chem_spatial(ax, chosen_stat="NMB", cmap="RdBu", difference=True, stat_desc="NMB (%)", alpha=0.075, lower = 0.025, upper = 0.925, letter = "A"):
    # Convert to GeoDataFrame
    gdf = gpd.GeoDataFrame(metrics_df_2020, geometry=gpd.points_from_xy(metrics_df_2020.Longitude, metrics_df_2020.Latitude), crs="EPSG:4326")

    # Define map projection
    projection = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

    # Set projection for the subplot
    #ax.set_projection(projection)

    # Add map features
    ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray', zorder=-1)
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue', zorder=-1)
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='lightblue', zorder=-1)
    ax.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='black', zorder=-1)

    # Normalize color scale
    if difference:
        vMAX = max(abs(gdf[chosen_stat].quantile(lower)), abs(gdf[chosen_stat].quantile(upper)))
        norm = Normalize(vmin=-vMAX, vmax=vMAX)
    else:
        vmin = gdf[chosen_stat].quantile(lower)
        vmax = gdf[chosen_stat].quantile(upper)
        norm = Normalize(vmin=vmin, vmax=vmax)

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Scatter plot with color mapping
    sc = ax.scatter(gdf.geometry.x, gdf.geometry.y, c=gdf[chosen_stat], cmap=cmap, norm=norm,
                    transform=ccrs.PlateCarree(), edgecolors='k', linewidths=0.5, alpha=1, s=80)

    # Add colorbar with adjusted font sizes
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', shrink=0.675)
    cbar.set_label(stat_desc, fontsize=18)
    cbar.ax.tick_params(labelsize=16)  # Adjust tick font size

    # Add title with font size
    ax.set_title(f"{letter} Spatial Distribution of {stat_desc.split('(')[0]}", fontsize=22)

    # Set map extent
    ax.set_extent([-120, -72, 24.5, 49.5], crs=ccrs.PlateCarree())

# Create 2×2 subplot figure
fig, axes = plt.subplots(2, 2, figsize=(20, 16), subplot_kw={'projection': ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))})

# Generate the four maps
plot_chem_spatial(axes[0, 0], stat_desc="Normalized Mean Bias (%)")
plot_chem_spatial(axes[0, 1], chosen_stat="overestimation", cmap="RdBu", difference=True, stat_desc="Mean Bias (µg/m³)", alpha=0.075, letter = "B")
plot_chem_spatial(axes[1, 0], chosen_stat="NME", cmap="magma", difference=False, stat_desc="Normalized Mean Error (%)", lower = 0.005, upper=0.9, letter = "C")
plot_chem_spatial(axes[1, 1], chosen_stat="r", cmap="magma_r", difference=False, stat_desc="Pearson r", upper=0.999, letter = "D")

# Adjust layout
plt.tight_layout()
plt.subplots_adjust(wspace=0.03, hspace=-0.25)

# Show figure
plt.show()


# %% [markdown]
# ## Plot Over Time

# %%
df.columns

# %%
df_2020.columns

# %%
import pandas as pd
import matplotlib.pyplot as plt

# Sample DataFrame structure
# df = pd.DataFrame({
#     'timestamp': pd.date_range(start='2024-01-01', periods=5000, freq='H'),
#     'station': ['Station_A'] * 5000,
#     'measured': np.random.rand(5000) * 100,
#     'modeled': np.random.rand(5000) * 100,
# })

# Ensure 'timestamp' is a datetime object
df_copy = df_2020.copy(deep = True)
df_copy['dt'] = pd.to_datetime(df_copy['dt'])

df_copy = df_copy[df_copy['dt'] >= "2016-01-01 00:00:00"]

# Set 'timestamp' as index (optional if not already)
df_copy.set_index('dt', inplace=True)

# Extract hour and date information
df_copy['hour'] = df_copy.index.hour

# Group by station and hour, then calculate the mean for each hour of the day
daily_profile = df_copy.groupby(['Latitude', 'Longitude', 'hour']).agg({'Sample Measurement': 'mean', 'CMAQ': 'mean'}).reset_index()



def plot_station_diurnal_profile(daily_profile = daily_profile, Latitude = 36.19125, city = ""):

    # Plot for a specific station
    station_data = daily_profile[daily_profile['Latitude'] == Latitude]
    station_data['measured'] = station_data['Sample Measurement']
    station_data['modeled'] = station_data['CMAQ']

    plt.figure(figsize=(10, 6))
    plt.plot(station_data['hour'], station_data['measured'], label='Measured', marker='o')
    plt.plot(station_data['hour'], station_data['modeled'], label='Modeled', marker='x')

    plt.title(f'Average Daily Profile for {city} Station')
    plt.xlabel('Hour of the Day')
    plt.ylabel('Air Pollution (unit)')
    plt.legend()
    plt.grid(True)
    plt.xticks(range(24))  # Show all hours on the x-axis
    plt.show()
    
    
def plot_station_profile(df_copy = df_copy, Latitude = 36.19125, city = ""):
    from scipy.stats import pearsonr

    # Plot for a specific station
    station_data = df_copy[df_copy['Latitude'] == Latitude]
    station_data['measured'] = station_data['Sample Measurement']
    station_data['modeled'] = station_data['CMAQ']
    r_value, _ = pearsonr( station_data['measured'], station_data['modeled'])
    print("R:", round(r_value,3))

    plt.figure(figsize=(15, 8))
    plt.plot(station_data.index, station_data['measured'], label='Measured', marker='o')
    plt.plot(station_data.index, station_data['modeled'], label='Modeled', marker='x')

    plt.title(f'Average Daily Profile for {city} Station')
    plt.xlabel('Hour of the Day')
    plt.ylabel('Air Pollution (unit)')
    plt.legend()
    plt.grid(True)
    #plt.xticks(range(24))  # Show all hours on the x-axis
    plt.show()

plot_station_diurnal_profile(Latitude = 33.57453, city = "Phoenix")
plot_station_diurnal_profile(Latitude = 46.598056, city = "Seattle ")

plot_station_profile(Latitude = 33.57453, city = "Phoenix")
plot_station_profile(Latitude = 46.598056, city = "Seattle")




# %%
station_data  = daily_profile.groupby('hour').mean() 
station_data['measured'] = station_data['Sample Measurement']
station_data['modeled'] = station_data['CMAQ']
station_data = station_data.reset_index()
plt.figure(figsize=(10, 6))
plt.plot(station_data['hour'], station_data['measured'], label='Measured', marker='o')
plt.plot(station_data['hour'], station_data['modeled'], label='Modeled', marker='x')
plt.ylim([0,15])
plt.title(f'Average Diurnal PM₂.₅ Profile for CONUS')
plt.xlabel('Hour of the Day')
plt.ylabel('PM₂.₅ Concentration (µg/m³)')
plt.legend()
plt.grid(True)
plt.xticks(range(24))  # Show all hours on the x-axis
plt.show()

# %%
df_copy

# %%
import matplotlib.dates as mdates

full_profile = df_copy.groupby(['level_0']).agg({'Sample Measurement': 'mean', 'CMAQ': 'mean'}).reset_index()

full_profile['measured'] = full_profile['Sample Measurement']
full_profile['modeled'] = full_profile['CMAQ']
full_profile = full_profile.reset_index()
full_profile['level_0'] = pd.to_datetime(full_profile['level_0'] )
plt.figure(figsize=(10, 6))
plt.plot(full_profile['level_0'], full_profile['measured'], label='Measured', marker='o')
plt.plot(full_profile['level_0'], full_profile['modeled'], label='Modeled', marker='x')

plt.ylim([0, 30])
plt.title(f'Average January PM₂.₅ Profile for CONUS')
plt.xlabel('Date')
plt.ylabel('PM₂.₅ Concentration (µg/m³)')
plt.legend()
plt.grid(True)

# Set x-axis format to show only a few ticks
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=5))  # Every 4 hours
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))  # Format as 'Jan 01'

plt.show()
