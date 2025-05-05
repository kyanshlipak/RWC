"""
WRF Data Analysis and Visualization Script

Description:
This script processes and visualizes Weather Research and Forecasting (WRF) model output data.
It calculates derived meteorological variables (relative humidity, wind direction) and generates
geospatial plots of key atmospheric parameters.

Features:
- Computes relative humidity from temperature and moisture variables
- Calculates wind direction from U/V components
- Creates publication-quality maps of meteorological fields
- Supports visualization of both native and derived variables

Input:
- NetCDF format WRF output files (specifically '2020_wrf_avg_new.nc')

Output:
- Geographic plots of temperature, humidity, wind speed and direction
- Visualizations saved as high-resolution images

Dependencies:
- xarray, numpy, matplotlib, cartopy, netCDF4

Kyan Shlipak May. 2025
"""

# Import necessary libraries
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
#from wrf import getvar, interplevel, to_np, latlon_coords
import cartopy.crs as crs
import cartopy.feature as creature
import xarray as xr
wrf_2020 = xr.open_dataset('2020_wrf_avg_new.nc')

def compute_rh_approx(wrf_ds):
    Q2 = wrf_ds["Q2"]  # Water vapor mixing ratio (kg/kg)
    T2 = wrf_ds["T2"]  # 2-meter temperature (K)

    # Assume standard surface pressure (Pa)
    P = wrf_ds['PSFC']

    # Compute saturation vapor pressure (es in Pa)
    es = 610.94 * np.exp((17.625 * (T2 - 273.15)) / (T2 - 30.11))

    # Compute saturation mixing ratio (Qs in kg/kg)
    Qs = (0.622 * es) / (P - es)

    # Compute RH (%) and clip between 0-100%
    RH = (Q2 / Qs) * 100
    RH = RH.clip(0, 100)

    return RH


wrf_2020['RH'] = compute_rh_approx(wrf_2020)


import numpy as np

def compute_wind_direction(wrf_ds):
    """
    Computes wind direction (degrees) from U10 and V10 wind components.

    Parameters:
    - wrf_ds: xarray Dataset containing 'U10' (m/s) and 'V10' (m/s).

    Returns:
    - Wind direction (degrees, 0-360°).
    """
    U10 = wrf_ds["U10"]
    V10 = wrf_ds["V10"]

    # Compute wind direction using arctan2 (returns radians)
    wind_dir = np.arctan2(-U10, -V10) * (180 / np.pi) + 180

    return wind_dir

wrf_2020['wind_dir'] = compute_wind_direction(wrf_2020)

# Import necessary libraries
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as crs
import cartopy.feature as cfeature
import xarray as xr

# # Variables of interest
# runname = 'output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852'
# dirToWRF = f'/projects/b1045/wrf-cmaq/output/CONUS4K/{runname}/'
# picdir = f'/home/ksz4578/Heat_Pump_Project/WRF_analyses/CMAQ_LCD/{runname}/figures/'

# # Gather WRF output files for domain 2 (d02) for a specific time range or single file
# filenames_d02 = sorted(glob.glob(os.path.join(dirToWRF, "wrfout_*")))

# Load the first file as an example

def plot_variable(xr_dataset, var_name, title, label, cmap = "virids", timeidx=0, save_plot=True):
    """
    Plots a specified variable from a WRF output file.

    Parameters:
    - filename: str, path to the WRF NetCDF file.
    - var_name: str, name of the variable to plot (e.g., 'T2', 'P', 'U', etc.).
    - timeidx: int, time index to plot (default is 0).
    - save_plot: bool, whether to save the plot to file (default is True).

    Returns:
    - None. Displays and optionally saves the plot.
    """
    # Load the WRF file


    if var_name == "windspeed":
        u10 = xr_dataset["U10"]
        v10 = xr_dataset["V10"]
        variable_data = np.sqrt(u10**2 + v10**2)
        #variable_data * 3600/1609 # convert to mph
    else:
        # Extract the variable
        variable_data = xr_dataset[var_name]
        
        # Example conversion for temperature (T2) from Kelvin to Fahrenheit
        #if var_name == "T2":
        #    variable_data = (variable_data - 273.15) # Convert to Fahrenheit

    # Get latitude and longitude coordinates
    lats, lons = xr_dataset["XLAT"], xr_dataset["XLONG"]

    # Set min and max for the color scale based on data range
    min_var = np.nanmin(variable_data)
    max_var = np.nanmax(variable_data)

    # Plotting
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=crs.PlateCarree())
    var_plot = ax.contourf(lons, lats, variable_data, 
                            levels=np.linspace(min_var, max_var, 100), cmap=cmap, transform=crs.PlateCarree())

    # Add color bar with label
    cbar = plt.colorbar(var_plot, ax=ax, orientation="horizontal", pad=0.05, aspect=50)
    cbar.set_label(label, fontsize=14)  # Set colorbar label font size
    cbar.ax.tick_params(labelsize=14)  # Adjust tick size on colorbar
   # cbar.set_ticks([0, 5, 10, 15, 20, 25, 31])
    
    resol = '50m'
    ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
    lakes = cfeature.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])

    # Customize plot appearance
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    # # Add features to the map
    #ax.add_feature(ocean, linewidth=0, zorder=0, facecolor="white")
    #ax.add_feature(lakes, linewidth=0, zorder=0, facecolor="white")
    ax.coastlines('50m', linewidth=0.8, color = "white")

    ax.add_feature(cfeature.BORDERS, linestyle=':', color = "white")
    #ax.add_feature(cfeature.STATES, linestyle=':', linewidth=0.5)
    ax.axis("off")

    # Title and Labels
    plt.title(title, fontsize = 16)
    # plt.xlabel("Longitude")
    # plt.ylabel("Latitude")
    plt.gcf().set_dpi(600)  # Set DPI to 300 (or any desired value)
    plt.show()

    
# Example usage
#plot_variable(wrf_2020, var_name='T2', title = "Average Simulated CONUS January Temperatures", label = "Temperature (°C)")  # Change 'T2' to any other variable name you wish to plot.

plot_variable(wrf_2020, var_name='T2', title = "Average Simulated CONUS January Temperatures", label = "Temperature (°C)", cmap = "magma")  # Change 'T2' to any other variable name you wish to plot.
plot_variable(wrf_2020, var_name='RH', title = "Average Simulated CONUS January Relative Humidity", label = "Relative Humidity (%)", cmap = "magma")     
plot_variable(wrf_2020, var_name='windspeed', title = "Average Simulated CONUS January Wind Speed", label = "Wind Speed (m/s)", cmap = "magma")  
plot_variable(wrf_2020, var_name='wind_dir', title = "Average Simulated CONUS January Wind Direction", label = "Wind Direction (°)", cmap = "twilight")  





###
### Plotting Heating Degree Days (HDDs)
###


HDDs = xr.open_dataset('../days_below_50F.nc')
plot_variable(HDDs, var_name='T2', title = "Simulated CONUS January HDDs (minimum temperature < 50°F)", label = "Number of Heating Degree Days ", cmap = "magma_r")  

import geopandas as gpd
import pandas as pd
import xarray as xr
from shapely.geometry import Point

# Step 1: Read the CBSA shapefile
cbsa = gpd.read_file('cbsa/tl_2019_us_cbsa.shp', engine="pyogrio")

# Step 2: Assuming HDDs is an xarray Dataset with 'XLAT' and 'XLONG' as coordinates
# For example, if HDDs is a DataFrame with 'XLAT', 'XLONG', and other relevant columns:
# If you have HDDs as an xarray, you can convert it to a pandas DataFrame like so:
hdd_df = HDDs.to_dataframe().reset_index()  # Reset index to have columns 'XLAT', 'XLONG', and others

# Step 3: Create a GeoDataFrame from the HDDs DataFrame
geometry = [Point(xy) for xy in zip(hdd_df['XLONG'], hdd_df['XLAT'])]
hdd_gdf = gpd.GeoDataFrame(hdd_df, geometry=geometry)

# Ensure the coordinate reference systems (CRS) match
# Assuming the CBSA shapefile uses EPSG:4326 (WGS84)
cbsa = cbsa.to_crs(epsg=4326)
hdd_gdf = hdd_gdf.set_crs(epsg=4326, allow_override=True)

# Step 4: Perform a spatial join to merge the data
