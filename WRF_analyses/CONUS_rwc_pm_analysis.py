# ------------------------------------------------------------------------------
# Script Name: CONUS_rwc_pm_analysis.py
#
# Description:
# This script analyzes PM2.5 concentrations attributed to Residential Wood 
# Combustion (RWC) across the contiguous United States (CONUS) using 4 km 
# resolution model output data. It:
#   - Loads PM2.5 data from baseline, RWC-included, and RWC-removed scenarios.
#   - Computes absolute and percent contributions of RWC to PM2.5 for 2016 and 2020.
#   - Generates a gridded GeoDataFrame for spatial analysis using projection info.
#   - Extracts RWC impacts at specific locations and their surrounding 3×3 grid.
#   - Creates main figure plot of CONUS and CBSA RWC-related PM2.5 spatial plots
#   - Creates figures of total simulated PM and RWC-related PM2.5 for CONUS domain

#
# Inputs:
#   - avg_201601.nc: Baseline PM2.5 data (includes RWC)
#   - avg_201601_2020_RWC.nc: 2020 PM2.5 data with updated RWC
#   - no_rwc_avg_201601.nc: PM2.5 data with RWC emissions removed
#   - counties shapefile: Used for mapping and reference (not actively plotted here)
#   - cbsa shape file
#
# Requirements:
#   - xarray, numpy, pandas, geopandas, shapely, cartopy, pyproj
#   - full requirements in requirements.txt
#
# ------------------------------------------------------------------------------



#import necessary libraries
import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import pandas as pd
import geopandas as gpd

#counties shapefile
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
reader = shpreader.Reader('../SMOKE_sensitivity_analyses/counties/countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())

baseline = xr.open_dataset("avg_201601.nc")
rwc_2020 = xr.open_dataset("avg_201601_2020_RWC.nc")
no_rwc = xr.open_dataset("no_rwc_avg_201601.nc")

percent_RWC = 100*(baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:])/baseline['PM25_TOT'][0,0,:,:]
percent_RWC_2020 = 100*(rwc_2020['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:])/rwc_2020['PM25_TOT'][0,0,:,:]
diff_2020 = rwc_2020['PM25_TOT'][0,0,:,:] - baseline['PM25_TOT'][0,0,:,:]

import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np

# Projection parameters
proj_params = {'proj': 'lcc',
               'lat_1': 33,
               'lat_2': 45,
               'lon_0': -97,
               'lat_0': 40}

# Coordinates of the origin
xorig = -2292000
yorig = -1584000

# Number of grid cells in x and y directions
num_cells_x = 1155
num_cells_y = 726

# Size of each grid cell (in meters)
cell_size = 4000  # 4km

# Generate grid coordinates using NumPy
x_coords = np.linspace(xorig, xorig + cell_size * num_cells_x, num_cells_x + 1)
y_coords = np.linspace(yorig, yorig + cell_size * num_cells_y, num_cells_y + 1)

# Create vertices for all grid cells using NumPy
x1, y1 = np.meshgrid(x_coords[:-1], y_coords[:-1])
x2, y2 = np.meshgrid(x_coords[1:], y_coords[:-1])
x3, y3 = np.meshgrid(x_coords[1:], y_coords[1:])
x4, y4 = np.meshgrid(x_coords[:-1], y_coords[1:])

# Reshape to 1D arrays
x1, x2, x3, x4 = x1.ravel(), x2.ravel(), x3.ravel(), x4.ravel()
y1, y2, y3, y4 = y1.ravel(), y2.ravel(), y3.ravel(), y4.ravel()

# Calculate row and column indices for each cell
rows, cols = np.divmod(np.arange(len(x1)), num_cells_x)

# Create GeoDataFrame with polygons
polygons = [Polygon([(x1[i], y1[i]), (x2[i], y2[i]), (x3[i], y3[i]), (x4[i], y4[i])]) for i in range(len(x1))]
pmgdf = gpd.GeoDataFrame({'row': rows, 'col': cols, 'geometry': polygons}, crs=proj_params)

pmgdf['PM25'] = (baseline['PM25_TOT'][0,0,:,:]).to_numpy().flatten()
pmgdf['PM25_2020'] = (rwc_2020['PM25_TOT'][0,0,:,:]).to_numpy().flatten()
pmgdf['PM25_diff'] = (rwc_2020['PM25_TOT'][0,0,:,:] - baseline['PM25_TOT'][0,0,:,:]).to_numpy().flatten()

pmgdf['PM25_RWC'] = (baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:]).to_numpy().flatten()
pmgdf['PM25_RWC_2020'] = (rwc_2020['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:]).to_numpy().flatten()
pmgdf['PM25_RWC_percent'] =(100* (baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:])/baseline['PM25_TOT'][0,0,:,:]).to_numpy().flatten()
pmgdf['PM25_RWC_percent_2020'] = (percent_RWC_2020).to_numpy().flatten()

pmgdf['RWC_percent_diff'] = (pmgdf['PM25_RWC_2020'] - pmgdf['PM25_RWC'])/ pmgdf['PM25'] * 100 

pmgdf['centroid'] = pmgdf.geometry.centroid

# MSA data
cbsa = gpd.read_file('cbsa/tl_2019_us_cbsa.shp', engine="pyogrio")

# Convert the GeoDataFrame to the target CRS
cbsa = cbsa.to_crs("+proj=lcc +lat_1=33 +lat_2=45 +lon_0=-97 +lat_0=40")

# Create a new GeoDataFrame with centroids
pmgdf_centroids = gpd.GeoDataFrame(geometry=pmgdf['centroid'], crs=pmgdf.crs)

# Perform spatial join based on centroids within cbsa geometry
merged_gdf = gpd.sjoin(pmgdf_centroids, cbsa, predicate='within')

# add to pm geodataframe
pmgdf['CBSA'] = merged_gdf['NAME']
pmgdf['CBSAFP'] =  merged_gdf['CBSAFP']

# Read the .txt file into a pandas DataFrame
population_df = pd.read_csv("../SMOKE_sensitivity_analyses/population.txt", delimiter='\t', header=None, skiprows=25)

# # Specify coordinates to align population and emission data
population_df = population_df.rename(columns = {2:"COLS", 3:"ROWS"})
population_df["COLS"] -= 111 
population_df["ROWS"] -= 126   

PWC_df = gpd.GeoDataFrame.merge(pmgdf, population_df, left_on=['row', 'col'], right_on=['ROWS', 'COLS'])

PWC_df['population'] = PWC_df[6]

PWC_df['PWC_TOT'] = PWC_df['PM25'] * PWC_df['population']
PWC_df['PWC_RWC'] = PWC_df['PM25_RWC'] * PWC_df['population']
PWC_df['PWC_TOT_2020'] = PWC_df['PM25_2020'] * PWC_df['population']
PWC_df['PWC_RWC_2020'] = PWC_df['PM25_RWC_2020'] * PWC_df['population']


total_PWC = PWC_df['PWC_TOT'].sum()
RWC_PWC = PWC_df['PWC_RWC'].sum()

total_PWC_2020 = PWC_df['PWC_TOT_2020'].sum()
RWC_PWC_2020 = PWC_df['PWC_RWC_2020'].sum()


import pandas as pd

# most populous 20
cbsa_list = [
    "New York-Newark-Jersey City, NY-NJ-PA",
    "Los Angeles-Long Beach-Anaheim, CA",
    "Chicago-Naperville-Elgin, IL-IN-WI",
    "Dallas-Fort Worth-Arlington, TX",
    "Houston-The Woodlands-Sugar Land, TX",
    "Washington-Arlington-Alexandria, DC-VA-MD-WV",
    "Miami-Fort Lauderdale-Pompano Beach, FL",
    "Philadelphia-Camden-Wilmington, PA-NJ-DE-MD",
    "Atlanta-Sandy Springs-Alpharetta, GA",
    "Boston-Cambridge-Newton, MA-NH",
    "Phoenix-Mesa-Chandler, AZ",
    "San Francisco-Oakland-Berkeley, CA",
    "Riverside-San Bernardino-Ontario, CA",
    "Detroit-Warren-Dearborn, MI",
    "Seattle-Tacoma-Bellevue, WA",
    "Minneapolis-St. Paul-Bloomington, MN-WI",
    "San Diego-Chula Vista-Carlsbad, CA",
    "Tampa-St. Petersburg-Clearwater, FL",
    "Denver-Aurora-Lakewood, CO",
    "St. Louis, MO-IL",
    "CONUS"
]


# Initialize an empty list to hold the data
data = []

# Iterate through each CBSA in the list
for cbsa_name in cbsa_list:
    cbsa_df = PWC_df[PWC_df['CBSA'] == cbsa_name]

    if cbsa_name == "CONUS": cbsa_df = PWC_df
    
    # Calculate basic statistics for PM25_RWC_2020
    pm25_rwc_2020_mean = cbsa_df['PM25_RWC_2020'].mean()
    pm25_rwc_2020_max = cbsa_df['PM25_RWC_2020'].max()
    pm25_rwc_2020_iqr = cbsa_df['PM25_RWC_2020'].quantile(0.75) - cbsa_df['PM25_RWC_2020'].quantile(0.25)
    
    # Calculate other required statistics
    pm_2014 = cbsa_df['PM25'].mean()
    pm_2020 = cbsa_df['PM25_2020'].mean()
    pm_2014_RWC = cbsa_df['PM25_RWC'].mean()
    pm_2020_RWC = cbsa_df['PM25_RWC_2020'].mean()
    percent_2014_RWC = pm_2014_RWC / pm_2014 * 100 if pm_2014 != 0 else None
    percent_2020_RWC = pm_2020_RWC / pm_2020 * 100 if pm_2020 != 0 else None
    
    # Append the results as a dictionary
    data.append({
        "CBSA": cbsa_name,
        # "PM25_2014": pm_2014,
        # "PM25_2020": pm_2020,
        # "PM25_2014_RWC": pm_2014_RWC,
        "PM25_2020_RWC": pm_2020_RWC,
        # "Percent_2014_RWC": percent_2014_RWC,
        "Percent_2020_RWC": percent_2020_RWC,
        # "PW_PM25_2014": cbsa_df['PWC_TOT'].sum() / cbsa_df['population'].sum(),
        #"PW_PM25_2020": cbsa_df['PWC_TOT_2020'].sum()/ cbsa_df['population'].sum(),
        # "PW_PM25_2014_RWC": cbsa_df['PWC_RWC'].sum() / cbsa_df['population'].sum(),
        "PW_PM25_2020_RWC": cbsa_df['PWC_RWC_2020'].sum() / cbsa_df['population'].sum(),
        # "PW_Percent_2014_RWC": (cbsa_df['PWC_RWC'].sum() / cbsa_df['population'].sum()) / (cbsa_df['PWC_TOT'].sum() / cbsa_df['population'].sum()) * 100,
        "PW_Percent_2020_RWC": (cbsa_df['PWC_RWC_2020'].sum() / cbsa_df['population'].sum()) / (cbsa_df['PWC_TOT_2020'].sum() / cbsa_df['population'].sum()) * 100,
        #"PM25_RWC_2020_Mean": pm25_rwc_2020_mean,
        "PM25_RWC_2020_Max": round(pm25_rwc_2020_max,2),
        "PM25_RWC_2020_IQR": pm25_rwc_2020_iqr
    })

# Create a DataFrame from the data
result_df = pd.DataFrame(data)

# Display the resulting DataFrame
#result_df.to_csv('MSA_20_RWC_results.csv')

print(result_df.round(2))

# ----------------------------------------------------
###
### East vs. West Analysis
### 



from pyproj import CRS, Transformer
# Projection parameters (same as in latlon_to_grid)
proj_params = {
    'proj': 'lcc',
    'lat_1': 33,
    'lat_2': 45,
    'lon_0': -97,
    'lat_0': 40
}
transformer = Transformer.from_crs(proj_params, "EPSG:4326", always_xy=True)

def grid_to_latlon(row, col):
    """
    Convert grid row and column indices to latitude and longitude.
    
    Parameters:
        row (int): Row index in the grid.
        col (int): Column index in the grid.
        
    Returns:
        tuple: (lat, lon) coordinates in WGS84.
    """

    cell_size = 4000  # Cell size in meters
    xorig = -2292000
    yorig = -1584000

    # Compute projected coordinates from row and column indices
    x_proj = xorig + col * cell_size
    y_proj = yorig + row * cell_size
    
    # Transform projected coordinates back to latitude and longitude
    lon, lat = transformer.transform(x_proj, y_proj)
    
    return lat, lon


rwc_2020_new_surg =  xr.open_dataset('../SMOKE_sensitivity_analyses/2020_new_surrogates_avg.nc')

def add_pm25(dataset):
    # Extract variables representing PM2.5 components
    pm25_components = ['PAL', 'PCA', 'PEC', 'PFE', 'PH2O', 'PK', 'PMG', 'PMN', 'PMOTHR', 'PNCOM', 'PNH4', 'PNO3', 'POC', 'PSI', 'PSO4', 'PTI']
    pm25_concentrations = dataset[pm25_components]

    # Sum up the concentrations along the specified dimensions (if necessary)
    # Access the variables you want to sum
    #PSO4', 'PNO3', 'POC', 'PEC', 'PNH4',
    dataset['PM25_total'] = pm25_concentrations.variables[pm25_components[0]]
    for component in pm25_components[1:]:
        dataset['PM25_total'].values += pm25_concentrations.variables[component].values

            
    # Sum the emissions for each grid cell across the variables
    #pm25_total = var1+ var2 + var3 + var4 + var5
    # Optionally, you can add attributes to the new variable to provide metadata
    dataset['PM25_total'].attrs['long_name'] = 'PM25_total'
    dataset['PM25_total'].attrs['units'] = 'g/s'
    dataset['PM25_total'].attrs['var_desc'] = 'Total PM2.5 concentration'
    return dataset

rwc_2020_new_surg = add_pm25(rwc_2020_new_surg)

seconds_in_a_day = 24 * 60 * 60  # Number of seconds in a day
days_in_month = 31  # Approximate days in a month (adjust if necessary)
seconds_in_a_month = seconds_in_a_day * days_in_month
grams_to_tons = 1 / 1_000_000  # Conversion factor from grams to tons
gps_to_tons = seconds_in_a_month * grams_to_tons
rwc_2020_new_surg['PM25_total_tons'] = rwc_2020_new_surg['PM25_total'] * gps_to_tons

pmgdf['PM25_total_tons'] = rwc_2020_new_surg['PM25_total_tons'][0,0,:,:].to_numpy().flatten()
PWC_df_emissions = gpd.GeoDataFrame.merge(pmgdf, population_df, left_on=['row', 'col'], right_on=['ROWS', 'COLS'])


# Function to determine if a polygon is east or west of -100 longitude
def is_west_of_100th_meridian(geometry):
    # Get the centroid of the polygon
    x_proj, y_proj = geometry.centroid.x, geometry.centroid.y
    
    # Convert to longitude/latitude
    lon, lat = transformer.transform(x_proj, y_proj)
    
    # Return True if west of -100, False otherwise
    return lon < -100

PWC_df_emissions['west_of_100'] = PWC_df_emissions['geometry'].apply(is_west_of_100th_meridian)

# Separate the DataFrame
PWC_df_emissions_west = PWC_df_emissions[PWC_df_emissions['west_of_100']]
PWC_df_emissions_east = PWC_df_emissions[~PWC_df_emissions['west_of_100']]


print("West top decile average emissions",PWC_df_emissions_west.loc[PWC_df_emissions_west[6] >= PWC_df_emissions_west[6].quantile(0.99)]['PM25_total_tons'].mean())

print("East top decile average emissions", PWC_df_emissions_east.loc[PWC_df_emissions_east[6] >= PWC_df_emissions_east[6].quantile(0.99)]['PM25_total_tons'].mean())

import geopandas as gpd
from pyproj import Transformer

# Define projection parameters (assuming Lambert Conformal Conic, same as before)
proj_params = {
    'proj': 'lcc',
    'lat_1': 33,
    'lat_2': 45,
    'lon_0': -97,
    'lat_0': 40
}

# Create a transformer from the projected CRS to WGS84
transformer = Transformer.from_crs(proj_params, "EPSG:4326", always_xy=True)

# Function to determine if a polygon is east or west of -100 longitude
def is_west_of_100th_meridian(geometry):
    # Get the centroid of the polygon
    x_proj, y_proj = geometry.centroid.x, geometry.centroid.y
    
    # Convert to longitude/latitude
    lon, lat = transformer.transform(x_proj, y_proj)
    
    # Return True if west of -100, False otherwise
    return lon < -100

# Load your GeoDataFrame (assuming PWC_df is already loaded with geometries)
PWC_df['west_of_100'] = PWC_df['geometry'].apply(is_west_of_100th_meridian)
PWC_df['population'] = PWC_df[6]

PWC_df['population_density'] = (PWC_df['population'] / PWC_df.geometry.area) * 1_000_000 * 1.609 **2# per km  

# Separate the DataFrame
PWC_west = PWC_df[PWC_df['west_of_100']]
PWC_east = PWC_df[~PWC_df['west_of_100']]

PWC_east = PWC_east.merge(
    PWC_df_emissions_east[['row', 'col', 'PM25_total_tons']],
    on=['row', 'col'],
    how='left'
)

PWC_west = PWC_west.merge(
    PWC_df_emissions_west[['row', 'col', 'PM25_total_tons']],
    on=['row', 'col'],
    how='left'
)

print("Urban emissions West",PWC_west[PWC_west['population_density'] > 500]['PM25_total_tons'].mean())

print("Urban emissions East",PWC_east[PWC_east['population_density'] > 500]['PM25_total_tons'].mean())


print("pw rwc West",(PWC_west['PM25_RWC_2020'] * PWC_west['population']).sum()/PWC_west['population'].sum())

print("Pw rwc East",(PWC_east['PM25_RWC_2020'] * PWC_east['population']).sum()/PWC_east['population'].sum())


print("East pw rwc", PWC_east['PM25_RWC_2020'].mean())
print("West pw rwc", PWC_west['PM25_RWC_2020'].mean())


print("East % rwc over 1ug", PWC_east.loc[PWC_east['PM25_RWC_2020'] >=1 ].shape[0]/PWC_east.shape[0] * 100)


print("West % rwc over 1ug",PWC_west.loc[PWC_west['PM25_RWC_2020'] >=1 ].shape[0]/PWC_west.shape[0] * 100)


# ----------------------------------------------------
###
### Urban vs. Rural
### 

urban = PWC_df.loc[PWC_df['population_density'] >= 500]
rural = PWC_df.loc[PWC_df['population_density'] < 500]

print("RWC percent urban", urban['PM25_RWC_2020'].sum() / urban['PM25_2020'].sum() * 100)
print("RWC percent rural",rural['PM25_RWC_2020'].sum() / rural['PM25_2020'].sum() * 100)

print("RWC mean urban", urban['PM25_RWC_2020'].mean())
print("RWC mean rural",rural['PM25_RWC_2020'].mean())


# ----------------------------------------------------
###
### Plotting Main Figure
### 

from pyproj import Proj

# Define the Lambert Conformal Conic projection parameters
def map_lat_lon(proj_params, xorig, yorig, ncols, nrows):
    # Define the projection transformation
    lcc_proj = Proj(proj_params)

    # Define the grid indices
    xcell, ycell = 4000, 4000  # Size of grid cells
    
    # Generate grid points
    x = np.arange(xorig, xorig + xcell * ncols, xcell)
    y = np.arange(yorig, yorig + ycell * nrows, ycell)
    
    # Create a meshgrid from the grid points
    x_mesh, y_mesh = np.meshgrid(x, y)
    
    # Convert grid indices to latitude and longitude using the projection
    lon_4km, lat_4km = lcc_proj(x_mesh, y_mesh, inverse=True)

    # Print latitude and longitude for a sample point
    print("Latitude at index (0, 0):", lat_4km[1, 1])
    print("Longitude at index (0, 0):", lon_4km[0, 1])
    
    return lon_4km, lat_4km

lon, lat = map_lat_lon(
    proj_params = {'proj': 'lcc',
                   'lat_1': 33,
                   'lat_2': 45,
                   'lon_0': -97,
                   'lat_0': 40},
    xorig = -2292000, yorig= -1584000,  # Coordinates of the origin
    ncols = 1155, nrows = 726  # Number of columns and rows
    )

# Load the boundary shapefile
boundaries = gpd.read_file("cbsa/cb_2018_us_nation_20m.shp")

# Ensure CRS match
boundaries = boundaries.to_crs(pmgdf.crs)

# Optionally simplify USA geometry for faster processing
boundaries["geometry"] = boundaries.simplify(tolerance=0.01, preserve_topology=True)

# Create a bounding box for the USA geometry to prefilter gdf
from shapely.geometry import box
usa_bounds = box(*boundaries.total_bounds)
pmgdf = pmgdf[pmgdf.intersects(usa_bounds)]  # Prefilter to reduce unnecessary computations

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

def plot_MSA(df, cbsa, MSA, title, col="PM25", cmap="viridis", vs=False, counties=False, difference=False, state_color="white", ax=None, title_size=16):
    from shapely.geometry import box
    
    # Define the target projection
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
    
    # Generate map features
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
    
    # Add features to the map
    ax.add_feature(ocean, linewidth=0, zorder=0, facecolor="white")
    ax.add_feature(lakes, linewidth=0, zorder=0, facecolor="white")
    
    # Get the geometry for the specified MSA
    if MSA != "":
        msa_geometry = cbsa.loc[cbsa["NAME"] == MSA, "geometry"].iloc[0]
        clipped_df = df.clip(mask=msa_geometry)
    else:
        msa_geometry = boundaries.geometry
        clipped_df = df  
    
    # Set visualization limits
    if vs:
        vmin = vs[0]
        vmax = vs[1]
    elif not difference:
        vmin = np.quantile(clipped_df[col], 0.025)
        vmax = np.quantile(clipped_df[col], 0.975)
    else:
        vmin = -max(abs(np.quantile(clipped_df[col], 0.025)), abs(np.quantile(clipped_df[col], 0.975)))
        vmin = max(vmin, -20)
        vmax = -vmin
    
    # Plot clipped data
    sc = clipped_df.plot(ax=ax, column=col, cmap=cmap, legend=False, vmin=vmin, vmax=vmax, zorder=-1, edgecolor="none", antialiased=False)
    if MSA == "":
        ax.plot([-100, -100], [25, 50.5], transform=ccrs.PlateCarree(), color='black', linewidth=1, linestyle='--')
        # # Add label for the 100th meridian
        # ax.text(-100.5, 48, "100th Meridian", color='red', fontsize=14, rotation=90, transform=ccrs.PlateCarree(),
        #         ha='right', va='top', fontweight='bold')

    # Create the ScalarMappable for the color bar
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Required for colorbar
    
    # Add the title with adjusted position
    #ax.set_title(title, fontsize=title_size, loc='center', pad=10)

    # Adjust layout to ensure that titles fit
    print('done')
    ax.axis("off")
    return sc, sm


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

# Initialize the figure
target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

fig = plt.figure(figsize=(16, 14), dpi = 300)  # Larger figure to accommodate everything

# Custom positions (x1, x2, y1, y2) for each subplot, adjusted for add_axes format
positions = {
    "Seattle-Tacoma-Bellevue, WA": (0, 0.25, 0, 0.35),  # Top left
    "Minneapolis-St. Paul-Bloomington, MN-WI": (0.25, 0.53, 0, 0.3),  # Top-center-left
    "Chicago-Naperville-Elgin, IL-IN-WI": (0.53, 0.75, 0, 0.3),  # Top-center-right
    "Boston-Cambridge-Newton, MA-NH": (0.75, 1, 0, 0.4),  # Top-right
    "New York-Newark-Jersey City, NY-NJ-PA": (0.75, 1, 0.4, 0.7),  # Right-middle
    "Philadelphia-Camden-Wilmington, PA-NJ-DE-MD": (0.6, 1, 0.7, 1),  # Bottom-left
    "Washington-Arlington-Alexandria, DC-VA-MD-WV": (0.25, 0.6, 0.7, 1),  # Bottom-right
    "Los Angeles-Long Beach-Anaheim, CA": (0, 0.25, 0.68, 1),  # Right-middle
    "Denver-Aurora-Lakewood, CO": (0, 0.25, 0.35, 0.68),  # Bottom-left-center
    "": (0.25, 0.75, 0.3, 0.7),  # Bottom-center
}

# Convert positions to the format (left, bottom, width, height)
converted_positions = {
    key: (x1, 1 - y2, x2 - x1, y2 - y1)  # Calculate width and height from x1, x2, y1, y2
    for key, (x1, x2, y1, y2) in positions.items()
}

labels = list("BCDEFGHIJA")  # 10 regions, adjust as needed
i = 0
# Loop through each MSA and plot
for MSA, pos in list(converted_positions.items())[:]:
    summary_stats = result_df

    df = pmgdf  # Assuming pmgdf is the data frame

    ax = fig.add_axes(pos, projection=target_proj)
    ax.set_position([pos[0], pos[1], pos[2], pos[3] - 0.028])  # Adjust top margin (0.05) as needed

    # Dynamically determine the upper bound for the color scale
    upper_bound = 5
    if MSA != "":
        if df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(0.98) >= 5:
            bound_2020 = np.ceil(df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(0.98))
            upper_bound = bound_2020
        place = MSA.split('-')[0]
    else:
        df = PWC_df
        if df['PM25_RWC_2020'].quantile(0.98) >= 5:
            bound_2020 = np.ceil(df['PM25_RWC_2020'].quantile(0.995))
            upper_bound = bound_2020
        place = "CONUS"


    sc, sm = plot_MSA(
        df=df,
        col="PM25_RWC_2020",
        cbsa=cbsa,
        MSA=MSA,
        title=f'{place}',
        cmap="magma_r",
        vs=[0, upper_bound],
        counties=False,
        difference=False,
        ax=ax,
        title_size=18
    )
    
    # Adjust colorbar
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", shrink=0.8, extend='max', extendfrac=0.04, aspect = 15 )
    cbar.ax.tick_params(labelsize=18)
    
    # Draw a bounding box around the entire axis (including the colorbar)
    bbox = Rectangle(
        (pos[0], pos[1]),  # Starting from the bottom-left corner of the axis
        pos[2],  # Width
        pos[3],  # Height
        transform=fig.transFigure,  # Use figure-relative coordinates for the bounding box
        color="black",  # Box color
        linewidth=2,  # Line width of the box
        fill=False  # No fill for the bounding box
    )
    fig.add_artist(bbox)  # Add the bounding box to the figure

    # Adjust title position to ensure it fits within the axis
    #ax.set_title(place, fontsize=18, loc='center',pad=10)  # Adjust the pad value as necessary
    # Calculate the center x-position of the axis
    title_x = pos[0] + pos[2] / 2
    title_y = pos[1] + pos[3] - 0.028  # Slightly above the top of the bounding box

    if place == "CONUS":
        place = "CONUS RWC PM₂.₅ Concentration (µg/m³)"
        title_y = pos[1] + pos[3] - 0.048  # Slightly above the top of the bounding box


    fig.text(
        title_x, title_y,
        f"{labels[i]} {place}",
        ha='center',
        va='bottom',
        fontsize=23,
    )
    #ax.set_aspect('auto')
    #ax.set_clip_on(True)
    i+=1

# Manually adjust subplot layout to ensure everything fits well
plt.show()
#fig.savefig("msa_pm25_map.png", format="png", dpi=600, bbox_inches='tight')


# ----------------------------------------------------
###
### Plotting CONUS PM2.5 Concentrations
### 

def plot(data, label, title, lat, lon, cbar="viridis", vs=False, counties=False, bounds=False, difference=False, state_color="white"):
    if not bounds:
        bounds = [[0, lat.shape[0]], [0, lon.shape[1]]]
    
    # Create a mask for areas outside the bounds and fill with zeros
    masked_data = np.zeros_like(data)
    masked_data[bounds[0][0]:bounds[0][1], bounds[1][0]:bounds[1][1]] = data[bounds[0][0]:bounds[0][1], bounds[1][0]:bounds[1][1]]

    # Create matplotlib figure
    plt.figure(figsize=(15, 8))

    # Create subplot with the Mollweide projection
    ax = plt.subplot(111, projection=ccrs.PlateCarree(central_longitude=-97.0, globe=None))

    # Create map
    if difference:
        mm = ax.pcolormesh(
            lon,
            lat,
            masked_data,
            transform=ccrs.PlateCarree(),
            cmap=cbar
        )
    else:
        mm = ax.pcolormesh(
            lon,
            lat,
            masked_data,
            transform=ccrs.PlateCarree(),
            cmap=cbar
        )
    
    # Map coastlines on the map with resolution 110m
    #ax.coastlines(resolution='50m', color=state_color, linewidth = 0.5)
    
    colormap = plt.get_cmap(cbar)  # Use the provided colormap
    zero_color = colormap(0)  # Get the color at zero value (this is RGBA)

    # Generate features
    resol = '50m'  # use data at this scale
    bodr = cfeature.NaturalEarthFeature(category='cultural', edgecolor=state_color,
                                        name='admin_0_boundary_lines_land', scale=resol, facecolor='none', alpha=0.7)
    ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale=resol, edgecolor='white', facecolor=cfeature.COLORS['water'], zorder = 3)
    lakes = cfeature.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='white', facecolor=cfeature.COLORS['water'])
    states_provinces = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines',
                                                    scale='50m', edgecolor=state_color, facecolor='none')
    
    # Add features
    if counties:
        ax.add_feature(COUNTIES, facecolor='none', edgecolor=state_color)
    ax.add_feature(ocean, linewidth=0)
    ax.add_feature(states_provinces, edgecolor=state_color, linewidth=0.5)
    ax.add_feature(bodr, edgecolor=state_color, alpha=1, linewidth = 0.5)
    
    if vs:
        mm.set_clim(vmin=vs[0], vmax=vs[1])
    else:
        mm.set_clim(vmin=np.quantile(data, 0.025), vmax=np.quantile(data, 0.975))

    # Plot colorbar legend
    cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.7, pad=0.02)  # Shrink adjusts the height

    # Get the color corresponding to zero from the colormap

    # Generate land feature with the zero color
    #land = cfeature.NaturalEarthFeature('physical', 'land', scale='10m', edgecolor='black', facecolor=zero_color)

    # Add the land feature to the bottom layer
    #ax.add_feature(land, zorder=-10)  # Use a negative zorder to ensure it's at the bottom
    
    cbar.outline.set_visible(True)  # Make colorbar outline visible
    cbar.ax.set_ylabel(label, rotation=90, fontsize=18,  labelpad=9)
    cbar.ax.tick_params(labelsize=18)
    ax.plot([-100, -100], [25, 50.5], transform=ccrs.PlateCarree(), color='red', linewidth=2, linestyle='--')
    # Add label for the 100th meridian
    ax.text(-100.5, 48, "100th Meridian", color='red', fontsize=14, rotation=90, transform=ccrs.PlateCarree(),
            ha='right', va='top', fontweight='bold')
    
    ax.set_title(title, fontsize = 20)
    ax.set_extent([-127, -66.5, 25, 50.5], crs=ccrs.PlateCarree())
    # Remove the axis frame and borders
    ax.set_frame_on(False)
    ax.spines['geo'].set_visible(False)  # Hide geographic spines
    plt.gcf().set_dpi(600)  # Set DPI to 300 (or any desired value)
    plt.show()

plot(data = baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:],
     vs = [0,8],
     label =  "RWC PM2.5 Concentrations (µg/m³)", 
     title = "NEI 2014 Average RWC contributed PM2.5 Concentrations",
     lon = lon,
     lat = lat,
     cbar = "turbo",
     state_color = "black"
    )

plot(data = rwc_2020['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:],
     vs = [0,5],
     label =  "RWC PM₂.₅ Concentrations (µg/m³)", 
     title = "Average RWC contributed PM₂.₅ Concentrations",
     lon = lon,
     lat = lat,
     cbar = "magma_r",
     state_color = "white"
    )