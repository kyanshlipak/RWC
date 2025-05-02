# Description: 
#   This script analyzes the impact of Residential Wood Combustion (RWC) on PM2.5 
#   concentrations at different spatial resolutions (4km, 12km, 36km). It processes
#   netCDF files containing PM2.5 data, performs spatial aggregations, calculates
#   population-weighted averages, estimates health impacts, and generates visualizations.
#
# Key Functionalities:
#   1. Processes baseline, RWC, and no-RWC PM2.5 scenarios
#   2. Calculates percentage contributions of RWC to total PM2.5
#   3. Aggregates data to 12km and 36km resolutions
#   4. Merges with census tract data for demographic analysis
#   5. Estimates mortality impacts using concentration-response functions
#   6. Generates comparative visualizations at different resolutions
#
# Input Files:
#   change these paths to match yours
#   - NetCDF files: avg_201601.nc, avg_201601_2020_RWC.nc, no_rwc_avg_201601.nc
#   - Shapefiles: cbsa/tl_2019_us_cbsa.shp, US_census_shapefiles/US_tract_2018.shp
#   - Population data: population.txt
#   - Mortality data: health_data/merged_baseline_mortality_rate.csv
#
# Outputs:
#   - Shapefiles with RWC PM2.5 at census tract level (12km and 36km resolutions)
#   - Comparative visualizations of PM2.5 at different aggregation resolutions
#   - Mortality estimates and demographic breakdowns
#
# Dependencies:
#   - xarray, netCDF4, cartopy, matplotlib, numpy, pandas, geopandas, pyproj
#
# Notes:
#   - Projection used: Lambert Conformal Conic (LCC) with standard parallels at 33° and 45°
#   - Health impact calculations use relative risk values from epidemiological studies
#   - Winter-only analysis (January data) with seasonal adjustment factor of 4
#



# ====================================================================================
# Basic importing of results


#import necessary libraries
import xarray as xr
from netCDF4 import Dataset, MFDataset, num2date
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import Proj
import geopandas as gpd

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

# ====================================================================================
# Aggregation to 12km and 36km

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Polygon

# Projection parameters
proj_params = {'proj': 'lcc',
               'lat_1': 33,
               'lat_2': 45,
               'lon_0': -97,
               'lat_0': 40}

# Coordinates of the origin
xorig = -2292000
yorig = -1584000

# Original 4km grid size
base_cell_size = 4000  # 4km
num_cells_x = 1155
num_cells_y = 726

def create_grid(cell_size, group_size):
    """Generates a GeoDataFrame for a coarser grid using the given cell_size."""
    new_x_cells = num_cells_x // group_size
    new_y_cells = num_cells_y // group_size

    # Generate new grid coordinates
    x_coords = np.linspace(xorig, xorig + cell_size * new_x_cells, new_x_cells + 1)
    y_coords = np.linspace(yorig, yorig + cell_size * new_y_cells, new_y_cells + 1)

    # Create polygon vertices
    x1, y1 = np.meshgrid(x_coords[:-1], y_coords[:-1])
    x2, y2 = np.meshgrid(x_coords[1:], y_coords[:-1])
    x3, y3 = np.meshgrid(x_coords[1:], y_coords[1:])
    x4, y4 = np.meshgrid(x_coords[:-1], y_coords[1:])

    # Reshape to 1D arrays
    x1, x2, x3, x4 = x1.ravel(), x2.ravel(), x3.ravel(), x4.ravel()
    y1, y2, y3, y4 = y1.ravel(), y2.ravel(), y3.ravel(), y4.ravel()

    # Compute new row and column indices
    rows, cols = np.divmod(np.arange(len(x1)), new_x_cells)

    # Create GeoDataFrame with polygons
    polygons = [Polygon([(x1[i], y1[i]), (x2[i], y2[i]), (x3[i], y3[i]), (x4[i], y4[i])]) for i in range(len(x1))]
    grid_gdf = gpd.GeoDataFrame({'row_group': rows, 'col_group': cols, 'geometry': polygons}, crs=proj_params)

    return grid_gdf

def aggregate_pm25(pmgdf, group_size):
    """Aggregates PM25_RWC_2020 by grouping row and col in group_size x group_size blocks."""
    pmgdf["row_group"] = pmgdf["row"] // group_size
    pmgdf["col_group"] = pmgdf["col"] // group_size

    aggregated = (
        pmgdf.groupby(["row_group", "col_group"])["PM25_RWC_2020"]
        .mean()  # Change to .sum() if needed
        .reset_index()
    )
    
    return aggregated

# Aggregate PM2.5 data
pmgdf_12km = aggregate_pm25(pmgdf, 3)
pmgdf_36km = aggregate_pm25(pmgdf, 9)

# Create polygon grids
grid_12km = create_grid(12000, 3)
grid_36km = create_grid(36000, 9)

# Merge pollution data with grid polygons
pmgdf_12km = grid_12km.merge(pmgdf_12km, on=["row_group", "col_group"], how="left")
pmgdf_36km = grid_36km.merge(pmgdf_36km, on=["row_group", "col_group"], how="left")


###
#### Population weighting
###

# Read the .txt file into a pandas DataFrame
population_df = pd.read_csv("../SMOKE_sensitivity_analyses/population.txt", delimiter='\t', header=None, skiprows=25)

# # Specify coordinates to align population and emission data
population_df = population_df.rename(columns = {2:"COLS", 3:"ROWS"})
population_df["COLS"] -= 111 
population_df["ROWS"] -= 126   
population_df['population'] = population_df[6]

def aggregate_population(pmgdf, group_size):
    """Aggregates PM25_RWC_2020 by grouping row and col in group_size x group_size blocks."""
    pmgdf["row_group"] = pmgdf["ROWS"] // group_size
    pmgdf["col_group"] = pmgdf["COLS"] // group_size

    aggregated = (
        pmgdf.groupby(["row_group", "col_group"])["population"]
        .mean()  # Change to .sum() if needed
        .reset_index()
    )
    
    return aggregated

population_12km = aggregate_population(population_df, 3)
population_36km = aggregate_population(population_df, 9)



pmgdf_12km = gpd.GeoDataFrame.merge(pmgdf_12km, population_12km, left_on=['row_group', 'col_group'], right_on=['row_group', 'col_group'])
pmgdf_36km = gpd.GeoDataFrame.merge(pmgdf_36km, population_36km, left_on=['row_group', 'col_group'], right_on=['row_group', 'col_group'])


# ====================================================================================
# Plotting

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from cartopy import crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

# Define the target projection
target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

# Create the plot
fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': target_proj})

resol = '50m'

# Base map features
#ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=0.8, zorder=1)  # Country borders
ax.add_feature(cfeature.OCEAN, facecolor="white", zorder=2)
ax.add_feature(cfeature.LAKES, facecolor="white", zorder=2)

# Load shapefile for countries to mask Mexico and Canada
shapefile = shpreader.natural_earth(resolution='50m', category='cultural', name='admin_0_countries')
countries = gpd.read_file(shapefile)

# Get Mexico and Canada
mexico_canada = countries[countries['ADMIN'].isin(['Canada', 'Mexico'])]

# Add Mexico and Canada as white
ax.add_geometries(mexico_canada.geometry, crs=ccrs.PlateCarree(), facecolor="white", edgecolor="none", zorder=3)

# Plot the GeoDataFrame
vmin, vmax = 0, 5  # Define color limits
cmap = 'magma_r'

sc = pmgdf_12km.plot(
    ax=ax,
    column='PM25_RWC_2020',  
    cmap=cmap,
    legend=False,
    vmin=vmin,
    vmax=vmax,
    zorder=0  # Ensures PM2.5 data is drawn on top
)

# Customize the title and axis
ax.set_title("RWC PM₂.₅ at 12km Resolution", fontsize=16)
ax.axis("off")

# Add colorbar with a label
sm = ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
cbar = fig.colorbar(sm, ax=ax, orientation='vertical', fraction=0.03, pad=0.02)
cbar.set_label("PM₂.₅ Concentration (µg/m³)", fontsize=14)
cbar.ax.tick_params(labelsize=14)

# Show the plot
plt.show()



import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from cartopy import crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

# Define the target projection
target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

# Create the plot
fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': target_proj})

resol = '50m'

# Base map features
#ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=0.8, zorder=1)  # Country borders
ax.add_feature(cfeature.OCEAN, facecolor="white", zorder=2)
ax.add_feature(cfeature.LAKES, facecolor="white", zorder=2)

# Load shapefile for countries to mask Mexico and Canada
shapefile = shpreader.natural_earth(resolution='50m', category='cultural', name='admin_0_countries')
countries = gpd.read_file(shapefile)

# Get Mexico and Canada
mexico_canada = countries[countries['ADMIN'].isin(['Canada', 'Mexico'])]

# Add Mexico and Canada as white
ax.add_geometries(mexico_canada.geometry, crs=ccrs.PlateCarree(), facecolor="white", edgecolor="none", zorder=3)

# Plot the GeoDataFrame
vmin, vmax = 0, 5  # Define color limits
cmap = 'magma_r'

sc = pmgdf_36km.plot(
    ax=ax,
    column='PM25_RWC_2020',  
    cmap=cmap,
    legend=False,
    vmin=vmin,
    vmax=vmax,
    zorder=0  # Ensures PM2.5 data is drawn on top
)

# Customize the title and axis
ax.set_title("RWC PM₂.₅ at 36km Resolution", fontsize=16)
ax.axis("off")

# Add colorbar with a label
sm = ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
cbar = fig.colorbar(sm, ax=ax, orientation='vertical', fraction=0.03, pad=0.02)
cbar.set_label("PM₂.₅ Concentration (µg/m³)", fontsize=14)
cbar.ax.tick_params(labelsize=14)

# Show the plot
plt.show()



# ====================================================================================
# Merging 12km and 36km aggregated RWC PM estimates
# with census tract data for health impact assessment, etc.

merged_gdf = gpd.read_file("../SMOKE_sensitivity_analyses/merged_gdf/merged_gdf.shp")
merged_gdf = merged_gdf.to_crs(pmgdf.crs) # convert to the same SMOKE CRS projection 
merged_shapes = merged_gdf[['GISJOIN','geometry']]
intersection_12km = gpd.overlay(pmgdf_12km, merged_shapes, how='intersection')
intersection_36km = gpd.overlay(pmgdf_36km, merged_shapes, how='intersection')

print(intersection_36km.columns)



census_tract_df = pd.read_csv('census_tract_race.csv', header = 0, encoding='latin1')
census_tract_df.drop(['YEAR',
                      'COUSUBA',
                      'PLACEA',
                      'BLKGRPA',
                      'CONCITA',
                      'AIANHHA',
                      'RES_ONLYA',
                      'TRUSTA',
                      'AIHHTLI',
                      'ANRCA',
                      'NECTAA',
                      'CNECTAA',
                      'NECTADIVA',
                      'CDCURRA',
                    'SLDUA',
                    'SLDLA',
                    'ZCTA5A',
                    'SUBMCDA',
                    'SDELMA',
                    'SDSECA',
                    'SDUNIA',
                    'PCI',
                    'PUMAA',
                    'BTBGA'], axis=1, inplace=True)


#2. Import US Census shapefiles 
import geopandas as gpd
gdf = gpd.read_file("../SMOKE_sensitivity_analyses/US_census_shapefiles/US_tract_2018.shp")

census_tract_gdf = gdf.merge(census_tract_df, how='left', 
                       left_on='GISJOIN', right_on='GISJOIN')
census_tract_gdf = census_tract_gdf.rename(columns={
           "ALUKE001":    "Total",
           "ALUKE003":    "White alone",
           "ALUKE004":    "Black or African American alone",
           "ALUKE005":    "American Indian and Alaska Native alone",
           "ALUKE006":    "Asian alone",
           "ALUKE007":    "Native Hawaiian and Other Pacific Islander alone",
           "ALUKE008":    "Some Other Race alone",
           "ALUKE009":    "Two or More Races",
           "ALUKE010":    "Two or More Races: Two races including Some Other Race",
           "ALUKE011":    "Two or More Races: Two races excluding Some Other Race, and three or more races",
           'ALUKE012': 'Hispanic or Latino',
           'ALUKE013': 'Hispanic or Latino: White alone',
           'ALUKE014': 'Hispanic or Latino: Black or African American alone',
           'ALUKE015': 'Hispanic or Latino: American Indian and Alaska Native alone',
           'ALUKE016': 'Hispanic or Latino: Asian alone',
           'ALUKE017': 'Hispanic or Latino: Native Hawaiian and Other Pacific Islander alone',
           'ALUKE018': 'Hispanic or Latino: Some other race alone',
           'ALUKE019': 'Hispanic or Latino: Two or more races',
           'ALUKE020': 'Hispanic or Latino: Two or more races: Two races including Some other race',
           'ALUKE021': 'Hispanic or Latino: Two or more races: Two races excluding Some other race, and three or more races'
                   })

# Projection parameters
proj_params = {'proj': 'lcc',
               'lat_1': 33,
               'lat_2': 45,
               'lon_0': -97,
               'lat_0': 40}

census_tract_gdf = census_tract_gdf.to_crs(proj_params) # convert to the same SMOKE CRS projection 

from shapely.geometry import Polygon

# Define a function to calculate the area of a polygon
def calculate_polygon_area(polygon):
    return polygon.area


# get census tract level RWC PM for 12km
intersection_12km['fraction'] = intersection_12km.geometry.area/(12000 ** 2) # area of intersection / divided by area of gridcell
intersection_12km['concentrations_per_polygon'] = intersection_12km['PM25_RWC_2020'] * intersection_12km['fraction'] # concentration * area of intersection/gridcell
summed_df_12km = intersection_12km.groupby('GISJOIN')['concentrations_per_polygon'].sum().reset_index() # concentration * census tract area / gridcell area
merged_gdf_12km = merged_gdf.merge(summed_df_12km, on='GISJOIN')

merged_gdf_12km['PM25_RWC_2020'] = merged_gdf_12km['concentrations_per_polygon']/merged_gdf_12km.geometry.area * 12000 ** 2
census_tract_gdf_12km = census_tract_gdf.merge(merged_gdf_12km, on='GISJOIN')
del census_tract_gdf_12km['geometry_y']
census_tract_gdf_12km.geometry = census_tract_gdf_12km['geometry_x']
census_tract_gdf_12km.to_file("2020_rwc_census_tract_pm25_12km.shp")


# get census tract level RWC PM for 36km 
intersection_36km['fraction'] = intersection_36km.geometry.area/(36000 ** 2) # area of intersection / divided by area of gridcell
intersection_36km['concentrations_per_polygon'] = intersection_36km['PM25_RWC_2020'] * intersection_36km['fraction'] # concentration * area of intersection/gridcell
summed_df_36km = intersection_36km.groupby('GISJOIN')['concentrations_per_polygon'].sum().reset_index() # concentration * census tract area / gridcell area
merged_gdf_36km = merged_gdf.merge(summed_df_36km, on='GISJOIN')

merged_gdf_36km['PM25_RWC_2020'] = merged_gdf_36km['concentrations_per_polygon']/merged_gdf_36km.geometry.area * 36000 ** 2
census_tract_gdf_36km = census_tract_gdf.merge(merged_gdf_36km, on='GISJOIN')
del census_tract_gdf_36km['geometry_y']
census_tract_gdf_36km.geometry = census_tract_gdf_36km['geometry_x']
census_tract_gdf_36km.to_file("2020_rwc_census_tract_pm25_36km.shp")



# ====================================================================================
# Basic numerical comparisons

print("CONUS average RWC 4km", PWC_df['PM25_RWC_2020'].mean())
print("CONUS average RWC 12km", pmgdf_12km['PM25_RWC_2020'].mean())
print("CONUS average RWC 436m", pmgdf_36km['PM25_RWC_2020'].mean())

pwc_4km = (PWC_df['PM25_RWC_2020'] * PWC_df['population']).sum() /  PWC_df['population'].sum()
pwc_12km = (pmgdf_12km['PM25_RWC_2020'] * pmgdf_12km['population']).sum() /  pmgdf_12km['population'].sum()
pwc_36km = (pmgdf_36km['PM25_RWC_2020'] * pmgdf_36km['population']).sum() /  pmgdf_36km['population'].sum()

print("CONUS PW average RWC 4km", pwc_4km)
print("CONUS PW average RWC 12km", pwc_12km)
print("CONUS PW average RWC 436m", pwc_36km)

# ====================================================================================
# Mortality comparisons

# mortality data
mortality_data = pd.read_csv("health_data/merged_baseline_mortality_rate.csv")
cbsa['CBSA'] = cbsa['NAME']

def add_CBSA(rwc_census_tract_baseline, cbsa):
    # Convert centroid to GeoDataFrame to use in spatial join
    rwc_census_tract_baseline = rwc_census_tract_baseline.copy(deep = True)
    rwc_census_tract_baseline['centroid'] = rwc_census_tract_baseline.geometry.centroid
    centroids_gdf = rwc_census_tract_baseline.set_geometry('centroid')
    
    
    # Perform spatial join based on centroid within cbsa geometry
    rwc_census_tract_baseline_merged = gpd.sjoin(centroids_gdf, cbsa[['geometry', 'CBSA']], predicate='within', how = 'left')
    
    # Merge the CBSA column back to the original rwc_census_tract_baseline dataframe
    rwc_census_tract_baseline['CBSA'] = rwc_census_tract_baseline_merged['CBSA']
    return rwc_census_tract_baseline


census_tract_gdf_12km = add_CBSA(census_tract_gdf_12km, cbsa)
census_tract_gdf_36km = add_CBSA(census_tract_gdf_36km, cbsa)

pm_mortality_data_df_12km = pd.merge(census_tract_gdf_12km, mortality_data, on = "GISJOIN")
pm_mortality_data_df_36km = pd.merge(census_tract_gdf_36km, mortality_data, on = "GISJOIN")

pm_mortality_data_df_12km['Total'] = pm_mortality_data_df_12km['Total_x']
pm_mortality_data_df_36km['Total'] = pm_mortality_data_df_36km['Total_x']


def mortality_function(pm_mortality_data_df, relative_risk = 1.08, res = "4km"):
    pm_mortality_data_df = pm_mortality_data_df.copy(deep = True)
    pm_mortality_data_df['real_relative_risk'] = relative_risk ** (pm_mortality_data_df['PM25_RWC_2020'] / 10 /4) # divide by four because only winter
    pm_mortality_data_df['AF'] = (pm_mortality_data_df['real_relative_risk'] - 1) / pm_mortality_data_df['real_relative_risk']
    pm_mortality_data_df['Attributable_Mortality'] = pm_mortality_data_df['Total'] * pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']
    pm_mortality_data_df['Attributable_Mortality_Rate'] = pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']
    print(res + ": RR:", relative_risk, " - ", round(pm_mortality_data_df['Attributable_Mortality'].sum()))
    return pm_mortality_data_df

def mortality_col_only(pm_mortality_data_df, relative_risk = 1.08, res = "4km"):
    pm_mortality_data_df = pm_mortality_data_df.copy(deep = True)
    pm_mortality_data_df['real_relative_risk'] = relative_risk ** (pm_mortality_data_df['PM25_RWC_2020'] / 10 /4) # divide by four because only winter
    pm_mortality_data_df['AF'] = (pm_mortality_data_df['real_relative_risk'] - 1) / pm_mortality_data_df['real_relative_risk']
    pm_mortality_data_df['Attributable_Mortality'] = pm_mortality_data_df['Total'] * pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']
    pm_mortality_data_df['Attributable_Mortality_Rate'] = pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']
    print(res + ": RR:", relative_risk, " - ", round(pm_mortality_data_df['Attributable_Mortality'].sum()))
    return pm_mortality_data_df[['Attributable_Mortality', 'Attributable_Mortality_Rate']]

pm_mortality_data_df_12km = mortality_function(pm_mortality_data_df_12km, res = "12km", relative_risk = 1.08)
pm_mortality_data_df_12km[['Attributable_Mortality_lower', 'Attributable_Mortality_Rate_lower']] = mortality_col_only(pm_mortality_data_df_12km, res = "12km",relative_risk  = 1.06)
pm_mortality_data_df_12km[['Attributable_Mortality_upper', 'Attributable_Mortality_Rate_upper']] = mortality_col_only(pm_mortality_data_df_12km, res = "12km",relative_risk  = 1.09)
_ = mortality_function(pm_mortality_data_df_12km, res = "12km",relative_risk = 1.17)


pm_mortality_data_df_36km = mortality_function(pm_mortality_data_df_36km, res = "36km",relative_risk = 1.08)
pm_mortality_data_df_36km[['Attributable_Mortality_lower', 'Attributable_Mortality_Rate_lower']] = mortality_col_only(pm_mortality_data_df_36km, res = "36km",relative_risk  = 1.06)
pm_mortality_data_df_36km[['Attributable_Mortality_upper', 'Attributable_Mortality_Rate_upper']] = mortality_col_only(pm_mortality_data_df_36km, res = "36km",relative_risk  = 1.09)
_ = mortality_function(pm_mortality_data_df_36km, res = "36km",relative_risk = 1.17)



def urban_rural_stats(pm_mortality_data_df_2020):
	pm_mortality_data_df_2020['population_density_per_sqmile'] = (
		pm_mortality_data_df_2020['Total'] / (pm_mortality_data_df_2020.geometry.area / 2_589_988.11)
	)

	pm_mortality_data_df_2020["Urban_Rural"] = pm_mortality_data_df_2020["population_density_per_sqmile"].apply(
		lambda x: "Urban" if x >= 500 else "Rural"
	)

	pm_mortality_data_df_2020.groupby("Urban_Rural")["Attributable_Mortality"].sum().reset_index()

	return pm_mortality_data_df_2020.groupby("Urban_Rural").apply(
		lambda x: x["Attributable_Mortality"].sum() / x["Total"].sum() * 100000
	).reset_index(name="Mortality_Ratio")

print("\n\n12km urban rural stats")
print(urban_rural_stats(pm_mortality_data_df_12km))

print("\n\n36km urban rural stats")
print(urban_rural_stats(pm_mortality_data_df_36km))



result_12km = (
    pm_mortality_data_df_12km
    .groupby('CBSA')
    .agg({
        'Attributable_Mortality': 'sum',
        'Attributable_Mortality_Rate': 'mean',
        'Total': 'sum'
    })
    .sort_values('Total', ascending=False)
    .head(20)  # Top 10 rows
)

print("\n\n12km CBSA stats")
print(result_12km)



result_36km = (
    pm_mortality_data_df_36km
    .groupby('CBSA')
    .agg({
        'Attributable_Mortality': 'sum',
        'Attributable_Mortality_Rate': 'mean',
        'Total': 'sum'

    })
    .sort_values('Total', ascending=False)
    .head(20)  # Top 10 rows
)

print("\n\n36km CBSA stats")
print(result_36km)


# ====================================================================================
# MSA resolution plotting comparisons


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

    # Create the ScalarMappable for the color bar
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Required for colorbar
    
    # Add the title with adjusted position
    ax.set_title(title, fontsize=title_size, loc='center', pad=10)

    # Adjust layout to ensure that titles fit
    ax.axis("off")
    return sc, sm


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


import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

for MSA in cbsa_list[0:3]:
    df = pmgdf

    upper_bound = 5
    if df.loc[df['CBSA'] == MSA]['PM25_RWC'].quantile(0.98) >= 5 or df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(0.98) >= 5:
        bound_2014 = np.ceil(df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(0.98))
        bound_2020 = np.ceil(df.loc[df['CBSA'] == MSA]['PM25_RWC'].quantile(0.98))
        upper_bound = max(bound_2014, bound_2020)
    
    # Create a Matplotlib figure with 1 row and 3 columns
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(30, 10), subplot_kw={'projection': target_proj})  # Adjust aspect ratio
    fig.subplots_adjust(wspace=0.05)  # Decrease space between plots

    place = MSA.split('-')[0]

    sc, sm = plot_MSA(
        df=pmgdf,
        col="PM25_RWC_2020",
        cbsa=cbsa,
        MSA=MSA,
        title=f'{place} RWC PM₂.₅ 4 km',
        cmap="magma_r",
        vs=[0, upper_bound],
        counties=False,
        difference=False,
        ax=axes[0],
        title_size=24
    )

    sc, sm = plot_MSA(
        df=pmgdf_12km,
        col="PM25_RWC_2020",
        cbsa=cbsa,
        MSA=MSA,
        title=f'{place} RWC PM₂.₅ 12 km',
        cmap="magma_r",
        vs=[0, upper_bound],
        counties=False,
        difference=False,
        ax=axes[1],
        title_size=24
    )

    sc, sm = plot_MSA(
        df=pmgdf_36km,
        col="PM25_RWC_2020",
        cbsa=cbsa,
        MSA=MSA,
        title=f'{place} RWC PM₂.₅ 36 km',
        cmap="magma_r",
        vs=[0, upper_bound],
        counties=False,
        difference=False,
        ax=axes[2],
        title_size=24
    )

    # Create a separate axis for the colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # Adjust position and size
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="vertical")
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_ylabel("PM₂.₅ Concentration (µg/m³)", rotation=90, fontsize=18, labelpad=7)

    plt.show()
