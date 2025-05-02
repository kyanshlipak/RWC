# ------------------------------------------------------------------------------
# Script Name: species_analysis.py
#
# Description:
# This script analyzes the chemical composition and contributions of PM2.5 from
# Residential Wood Combustion (RWC) in the continental United States using 4 km 
# gridded CMAQ model outputs. It:
#   - Loads PM2.5 and species-specific data for baseline, RWC-included (2020), 
#     and RWC-removed scenarios.
#   - Computes RWC contributions to total PM2.5 and to individual species.
#   - Constructs a 4 km resolution GeoDataFrame (`pmgdf`) of grid cells.
#   - Merges spatial PM2.5 data with population data by grid location.
#   - Calculates and summarizes the relative importance of species in RWC PM2.5.
#
# Inputs:
#   - avg_201601.nc: Baseline PM2.5 data (includes RWC)
#   - avg_201601_2020_RWC.nc: 2020 PM2.5 data with updated RWC emissions
#   - no_rwc_avg_201601.nc: PM2.5 data with RWC emissions removed
#   - population.txt: Grid-level population data
#
# Outputs:
#   - `rwc_contribution_df`: DataFrame summarizing the species composition of 
#     RWC-attributable PM2.5 (µg/m³, % of total RWC PM2.5, and % of species from RWC)
#
# Dependencies:
#   - Python libraries: xarray, numpy, pandas, geopandas, shapely
# ------------------------------------------------------------------------------



#import necessary libraries
import xarray as xr
import numpy as np
import pandas as pd
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



# Mapping PM2.5 species shorthand to real chemical names
pm25_species_names = {
    'PM25_SO4': 'Sulfate (SO₄²⁻)',
    'PM25_NO3': 'Nitrate (NO₃⁻)',
    'PM25_NH4': 'Ammonium (NH₄⁺)',
    'PM25_OC': 'Organic Carbon (OC)',
    'PM25_OM': 'Organic Matter (OM)',
    'PM25_EC': 'Elemental Carbon (EC)',
    'PM25_SOIL': 'Soil Particles',
    'PM25_CL': 'Chloride (Cl⁻)',
    'PM25_NA': 'Sodium (Na⁺)',
    'PM25_MG': 'Magnesium (Mg²⁺)',
    'PM25_K': 'Potassium (K⁺)',
    'PM25_CA': 'Calcium (Ca²⁺)',
    'PM25_UNSPEC1': 'Unspecified PM2.5'
}

# Initialize results dictionary
results = {}
SUM = 0
for species, real_name in pm25_species_names.items():
    SUM += (rwc_2020[species][0,0,:,:]).to_numpy().flatten().mean()
    # Calculate RWC contribution for each species
    pmgdf[species] = (rwc_2020[species][0,0,:,:]).to_numpy().flatten()
    pmgdf[f'{species}_RWC_2020'] = (rwc_2020[species][0,0,:,:] - no_rwc[species][0,0,:,:]).to_numpy().flatten()
    
# Read the .txt file into a pandas DataFrame
population_df = pd.read_csv("../SMOKE_sensitivity_analyses/population.txt", delimiter='\t', header=None, skiprows=25)

# # Specify coordinates to align population and emission data
population_df = population_df.rename(columns = {2:"COLS", 3:"ROWS"})
population_df["COLS"] -= 111 
population_df["ROWS"] -= 126   

PWC_df = gpd.GeoDataFrame.merge(pmgdf, population_df, left_on=['row', 'col'], right_on=['ROWS', 'COLS'])


for species, real_name in pm25_species_names.items():

    # Find percentage of total RWC PM2.5 from this species
    species_RWC_percent_2020 = 100 * (PWC_df[f'{species}_RWC_2020'].mean() / PWC_df['PM25_RWC_2020'].mean())

    # Find percentage of total PM2.5 for each species that comes from RWC
    species_from_RWC_2020 = 100 * (PWC_df[f'{species}_RWC_2020'].mean() / PWC_df[species].mean())
    
    # Store in dictionary with real chemical names
    results[real_name] = {
        'Average  RWC PM2.5 Contribution µg/m^3':PWC_df[f'{species}_RWC_2020'].mean(),
        'Percent of RWC PM2.5': round(species_RWC_percent_2020, 1), 
        'Percent of Species from RWC': round(species_from_RWC_2020, 1)
    }

# Convert results dictionary to DataFrame
rwc_contribution_df = pd.DataFrame.from_dict(results, orient='index')

# Display summary
rwc_contribution_df
