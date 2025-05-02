# ====================================================================================
# PM2.5 Residential Wood Combustion (RWC) Impact Analysis
#
# Description:
#   This comprehensive analysis pipeline evaluates the health impacts of PM2.5 from 
#   residential wood combustion by integrating air quality modeling, demographic data,
#   and health outcome metrics. The analysis includes:
#
# 1. BASE DATA PROCESSING
#    - Imports and processes 4km resolution PM2.5 data from netCDF files:
#      * Baseline concentrations
#      * RWC scenario concentrations 
#      * No-RWC scenario concentrations
#    - Calculates RWC contribution percentages and absolute differences
#    - Creates geospatial grid with polygon geometries

# 2. POPULATION-WEIGHTED ANALYSES
#    - Integrates population data with pollution grids
#    - Computes population-weighted concentrations (PWC)
#    - Calculates total and RWC-specific PWCs
#    - Aggregates to 12km and 36km resolutions

# 3. HEALTH IMPACT ASSESSMENT
#    - Processes mortality rate data by age group
#    - Standardizes age groups to match census data
#    - Computes baseline mortality rates (all ages and 25+)
#    - Applies concentration-response functions:
#      * Core analysis (RR=1.08 per 10μg/m³)
#      * Sensitivity analyses (RR=1.06, 1.09, 1.17)
#    - Calculates attributable mortality counts and rates

# 4. SPATIAL ANALYSES
#    - Allocates pollution to census tracts
#    - Computes tract-level:
#      * Absolute concentrations
#      * Per-capita exposures
#      * Mortality impacts
#    - Aggregates results by:
#      * CBSA (metropolitan areas)
#      * State
#      * Urban/Rural classification
#      * East/West regions (divided at 100°W)

# 5. DISPARITY ANALYSES
#    - Racial/ethnic exposure differences:
#      * By concentration deciles
#      * Population-weighted averages
#    - Relative disparity calculations:
#      * PM2.5 exposure differences
#      * Baseline mortality differences
#      * Attributable mortality differences
#    - Non-white vs. white exposure comparisons

# 6. EMISSIONS-CONCENTRATION ANALYSES
#    - Processes RWC emissions data
#    - Correlates emissions with concentrations:
#      * Spatial correlations
#      * Per-capita relationships
#    - Maps emissions vs. concentration patterns

# 7. VISUAL ANALYSES
#    - Multi-panel maps showing:
#      * Emissions
#      * Concentrations
#      * Mortality impacts
#    - Decile plots of racial/ethnic exposure distribution
#    - Relative disparity bar charts
#    - MSA-specific visualizations

# 8. STATISTICAL REPORTING
#    - Correlation matrices (emissions vs. demographics)
#    - Mortality rate tables by geographic level
#    - Population-weighted average results
#    - Percentage contributions by source
#    - Seasonal adjustment factors (winter-only)

# Inputs:
#   - NetCDF files: avg_201601.nc, avg_201601_2020_RWC.nc, no_rwc_avg_201601.nc
#   - Shapefiles: Census tract boundaries, CBSA boundaries
#   - CSV files: Mortality rates, population demographics
#
# Outputs:
#   - Spatial datasets of PM2.5 concentrations and health impacts
#   - Statistical analyses of exposure disparities
#   - Publication-quality visualizations
#
# Key Functions:
#   - mortality_function(): Core health impact calculation
#   - plot_MSA(): Spatial visualization generator
#   - decile_plot(): Demographic exposure distribution
#   - MSA_relative_disparity(): Disparity analysis
#
# Dependencies:
#   - xarray, geopandas, matplotlib, numpy, pandas, cartopy, pyproj, shapely, scipy


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
from shapely.geometry import Polygon


baseline = xr.open_dataset("avg_201601.nc")
rwc_2020 = xr.open_dataset("avg_201601_2020_RWC.nc")
no_rwc = xr.open_dataset("no_rwc_avg_201601.nc")

percent_RWC = 100*(baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:])/baseline['PM25_TOT'][0,0,:,:]
percent_RWC_2020 = 100*(rwc_2020['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:])/rwc_2020['PM25_TOT'][0,0,:,:]
diff_2020 = rwc_2020['PM25_TOT'][0,0,:,:] - baseline['PM25_TOT'][0,0,:,:]



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
cbsa['CBSA'] = cbsa['NAME']
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
# Basic importing of mortality data

def import_rwc_shapefile(loc = 'census_tract_data/2016_rwc_census_tract_pm25.shp'):
    rwc_census_tract_baseline = gpd.read_file(loc)
    rwc_census_tract_baseline = rwc_census_tract_baseline.rename(columns={
               "White alon":    "White",
               "Black or A":    "Black",
               "American I":    "American Indian",
               "Asian alon":    "Asian",
               "Native Haw":    "Native Hawaiian or Pacific Islander",
               'Some Other': 'Other', 
               'Two or Mor': 'Two or More Races',
               'Hispanic o': 'Hispanic'
                       })
    
    rwc_census_tract_baseline['GEOID_x'] = rwc_census_tract_baseline['GEOID_x'].astype(int)
    rwc_census_tract_baseline['COUNTYFP'] = rwc_census_tract_baseline['COUNTYFP'].astype(int)
    rwc_census_tract_baseline['STATEFP'] = rwc_census_tract_baseline['STATEFP'].astype(float)
    return rwc_census_tract_baseline

rwc_census_tract_2020 = import_rwc_shapefile(loc = 'census_tract_data/2020_rwc_census_tract_pm25.shp')
#rwc_census_tract_baseline = import_rwc_shapefile(loc = 'census_tract_data/2016_rwc_census_tract_pm25.shp')

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


#rwc_census_tract_baseline = add_CBSA(rwc_census_tract_baseline, cbsa)
rwc_census_tract_2020 = add_CBSA(rwc_census_tract_2020, cbsa)

mortality_data = pd.read_csv("health_data/merged_baseline_mortality_rate.csv")

# pm_mortality_data_df_baseline =  pd.merge(rwc_census_tract_baseline, mortality_data, on = "GISJOIN")
# pm_mortality_data_df_baseline['2016_pm2.5'] = pm_mortality_data_df_baseline['concentrat']

pm_mortality_data_df_2020 =  pd.merge(rwc_census_tract_2020, mortality_data, on = "GISJOIN")
pm_mortality_data_df_2020['2016_pm2.5'] = pm_mortality_data_df_2020['RWC_PM25_T'] / pm_mortality_data_df_2020.geometry.area * 4000 ** 2
pm_mortality_data_df_2020['2016_pm10'] = pm_mortality_data_df_2020['RWC_PM10_w'] / pm_mortality_data_df_2020.geometry.area * 4000 ** 2

def mortality_function(pm_mortality_data_df, relative_risk = 1.08, col = '2016_pm2.5'):
    pm_mortality_data_df = pm_mortality_data_df.copy(deep = True)
    pm_mortality_data_df['real_relative_risk'] = relative_risk ** (pm_mortality_data_df[col] / 10 /4) # divide by four because only winter
    pm_mortality_data_df['AF'] = (pm_mortality_data_df['real_relative_risk'] - 1) / pm_mortality_data_df['real_relative_risk']
    pm_mortality_data_df['Attributable_Mortality'] = pm_mortality_data_df['Total'] * pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']
    pm_mortality_data_df['Attributable_Mortality_Rate'] = pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']
    print("RR:", relative_risk, " - ", round(pm_mortality_data_df['Attributable_Mortality'].sum()))
    return pm_mortality_data_df

def mortality_col_only(pm_mortality_data_df, relative_risk = 1.08, col = '2016_pm2.5'):
    pm_mortality_data_df = pm_mortality_data_df.copy(deep = True)
    pm_mortality_data_df['real_relative_risk'] = relative_risk ** (pm_mortality_data_df[col] / 10 /4) # divide by four because only winter
    pm_mortality_data_df['AF'] = (pm_mortality_data_df['real_relative_risk'] - 1) / pm_mortality_data_df['real_relative_risk']
    pm_mortality_data_df['Attributable_Mortality'] = pm_mortality_data_df['Total'] * pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']
    pm_mortality_data_df['Attributable_Mortality_Rate'] = pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']
    print("RR:", relative_risk, " - ", round(pm_mortality_data_df['Attributable_Mortality'].sum()))
    return pm_mortality_data_df[['Attributable_Mortality', 'Attributable_Mortality_Rate']]

pm_mortality_data_df_2020 = mortality_function(pm_mortality_data_df_2020, relative_risk = 1.08)
pm_mortality_data_df_2020[['Attributable_Mortality_lower', 'Attributable_Mortality_Rate_lower']] = mortality_col_only(pm_mortality_data_df_2020, relative_risk  = 1.06)
pm_mortality_data_df_2020[['Attributable_Mortality_upper', 'Attributable_Mortality_Rate_upper']] = mortality_col_only(pm_mortality_data_df_2020, relative_risk  = 1.09)
_ = mortality_function(pm_mortality_data_df_2020, relative_risk = 1.17)


# ====================================================================================
# Basic importing of emissions data


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

# 2. Import US Census shapefiles 
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

# 3. Create gridded pollutant data

# Projection parameters
def make_grid_gdf():
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
    
    # Create GeoDataFrame with polygons
    polygons = [Polygon([(x1[i], y1[i]), (x2[i], y2[i]), (x3[i], y3[i]), (x4[i], y4[i])]) for i in range(len(x1))]
    grid_gdf = gpd.GeoDataFrame(geometry=polygons, crs=proj_params)
    return grid_gdf

grid_gdf = make_grid_gdf()

grid_gdf['PM25_total_tons'] = rwc_2020_new_surg['PM25_total_tons'].to_numpy().ravel()
grid_gdf = grid_gdf.reset_index().rename(columns={'index': 'iD'})

intersection = gpd.overlay(grid_gdf, census_tract_gdf, how='intersection')

intersection['fraction'] = intersection.geometry.area/(cell_size ** 2) # area of intersection / divided by area of gridcell
intersection['PM25_total_tons_per_polygon'] = intersection['PM25_total_tons'] * intersection['fraction'] # emissions of gridcell * area of intersection/area of gridcell

summed_df = intersection.groupby('GISJOIN')['PM25_total_tons_per_polygon'].sum().reset_index()
census_tract_gdf = census_tract_gdf.merge(summed_df, on='GISJOIN')

census_tract_gdf['PM25_total_tons_per_area'] = census_tract_gdf['PM25_total_tons_per_polygon']/census_tract_gdf.geometry.area * 1_000_000 * (1.609**2)
census_tract_gdf['PM25_total_tons_per_capita'] = census_tract_gdf['PM25_total_tons_per_polygon']/census_tract_gdf['Total'] * 100_000
census_tract_gdf['White_percentage'] = census_tract_gdf['White alone']/census_tract_gdf['Total'] * 100 
census_tract_gdf['Black_percentage'] = census_tract_gdf["Black or African American alone"]/census_tract_gdf['Total'] * 100 
census_tract_gdf['American Indian_percentage'] = census_tract_gdf["American Indian and Alaska Native alone"]/census_tract_gdf['Total'] * 100 
census_tract_gdf['Asian_percentage'] = census_tract_gdf["Asian alone"]/census_tract_gdf['Total'] * 100 
census_tract_gdf["Native Hawaiian or Pacific Islander_percentage"] = census_tract_gdf["Native Hawaiian and Other Pacific Islander alone"]/census_tract_gdf['Total'] * 100 
census_tract_gdf["Hispanic_percentage"] = census_tract_gdf['Hispanic or Latino']/census_tract_gdf['Total'] * 100 
census_tract_gdf['Non-White_Fraction'] = 100 - census_tract_gdf['White_percentage']

# Perform a spatial join to add CBSA information to census_tract_gdf
census_tract_gdf = census_tract_gdf.sjoin(cbsa[['geometry', 'NAME']], how="left", predicate="intersects")

# Rename the CBSA column for clarity
census_tract_gdf.rename(columns={'NAME': 'CBSA'}, inplace=True)
census_tract_gdf.rename(columns={'NAME_right': 'CBSA'}, inplace=True)

# Emissions Race Correlation Table --------------------------------------------------------------------------------------


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
    "St. Louis, MO-IL"
]

# Create an empty list to store results
emissions_correlation_results = []

# List of demographic columns
demographic_cols = [
    "Non-White_Fraction",
    "Black_percentage",
    "American Indian_percentage",
    "Asian_percentage",
    "Native Hawaiian or Pacific Islander_percentage",
    "Hispanic_percentage",
    "White_percentage"
]

# Target variables for correlation
target_vars = ["PM25_total_tons_per_capita"]
for msa in cbsa_list:
    if msa == "CONUS":
        MSA_df = census_tract_gdf
    else:
        # Filter data for the current MSA
        MSA_df = census_tract_gdf.loc[census_tract_gdf["CBSA"] == msa]

    # Compute Pearson correlation coefficients
    correlations = MSA_df[demographic_cols + target_vars].corr(method="pearson")

    # Extract relevant correlations
    msa_correlations = {"MSA": msa.split('-')[0]}
    
    for demo in demographic_cols:
        for target in target_vars:
            msa_correlations[f"{demo.split('_')[0]}_vs_{target}"] = correlations.loc[demo, target]

    # Append to list
    emissions_correlation_results.append(msa_correlations)

# Convert results into a DataFrame
emissions_correlation_df = pd.DataFrame(emissions_correlation_results)
emissions_correlation_df

# ====================================================================================
# CONUS Health Stats

###
# CBSA LEVEL --------------------------------------------------------------------------------------
###


results = []

for cbsa in cbsa_list:
    df_cbsa = pm_mortality_data_df_2020[pm_mortality_data_df_2020["CBSA"] == cbsa]
     
    summed = df_cbsa[["Attributable_Mortality", "Attributable_Mortality_lower", "Attributable_Mortality_upper"]].sum()/df_cbsa['Total'].sum() * 100_000
    summed["CBSA"] = cbsa
    
    results.append(summed)

print("CBSA Mortality Reuslts", pd.DataFrame(results))


# State LEVEL --------------------------------------------------------------------------------------

print("STATE Mortality Results", \
      pm_mortality_data_df_2020.groupby("STATE")[["Attributable_Mortality", "Attributable_Mortality_lower", "Attributable_Mortality_upper"]].sum().reset_index())


# Urban vs. Rural  --------------------------------------------------------------------------------------
pm_mortality_data_df_2020['population_density_per_sqmile'] = (
    pm_mortality_data_df_2020['Total'] / (pm_mortality_data_df_2020.geometry.area / 2_589_988.11)
)

pm_mortality_data_df_2020["Urban_Rural"] = pm_mortality_data_df_2020["population_density_per_sqmile"].apply(
    lambda x: "Urban" if x >= 500 else "Rural"
)
print("Urban vs. rural mortality\n", \
      pm_mortality_data_df_2020.groupby("Urban_Rural")["Attributable_Mortality"].sum().reset_index()
)

print("\n\nUrban vs. rural mortality rates\n", \
pm_mortality_data_df_2020.groupby("Urban_Rural").apply(
    lambda x: x["Attributable_Mortality"].sum() / x["Total"].sum() * 100000
).reset_index(name="Mortality_Rate")
)

# East vs. West  --------------------------------------------------------------------------------------

import geopandas as gpd
from pyproj import Proj, transform

# Define projection
proj_lcc = Proj(proj='lcc', lat_1=33, lat_2=45, lon_0=-97, lat_0=40, x_0=0, y_0=0, datum='WGS84')

# Define WGS84 (lat/lon)
proj_wgs84 = Proj(proj="latlong", datum="WGS84")

# Convert 100°W longitude (with a reference latitude) to projected X coordinate
x_100w, _ = transform(proj_wgs84, proj_lcc, -100, 40)  # Using lat=40° as a reference

# Classify East vs. West
pm_mortality_data_df_2020["Region"] = pm_mortality_data_df_2020["geometry"].centroid.x.apply(
    lambda x: "West" if x < x_100w else "East"
)

print("\n\n East vs. West mortality\n", \
	pm_mortality_data_df_2020.groupby("Region")["Attributable_Mortality"].sum().reset_index()
)

print("\n\n East vs. West mortality rates\n", \
	pm_mortality_data_df_2020.groupby("Region").apply(
    	lambda x: x["Attributable_Mortality"].sum() / x["Total"].sum() * 100000
	).reset_index(name="Mortality_Ratio")
)



# ====================================================================================
# Emissions - Concentrations - Mortality Plots
from shapely.geometry import box
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


def plot_MSA(df, cbsa, MSA, title, col="PM25", cmap="viridis", vs=False, counties=False, difference=False, state_color="white", ax=None, title_size = 20):
    from shapely.geometry import box
    
    # Define the target projection
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
    
    # Generate map features
    resol = '50m'
    bodr = cfeature.NaturalEarthFeature(category='cultural', edgecolor='white', name='admin_0_boundary_lines_land', scale='50m', facecolor='none', alpha=0.7)
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
    # if counties:
    #     ax.add_feature(COUNTIES, facecolor='none', edgecolor=state_color)
    ax.add_feature(ocean, linewidth=0, zorder=0, facecolor = "white")
    ax.add_feature(lakes, linewidth=0, zorder=0, facecolor = "white")
    # ax.add_feature(states_provinces, edgecolor=state_color, linewidth=0.5, zorder=0)
    
    # Get the geometry for the specified MSA
    msa_geometry = cbsa.loc[cbsa["NAME"] == MSA, "geometry"].iloc[0]
    
    # Clip the data to the CBSA boundary
    clipped_df = df.clip(mask=msa_geometry)
    #clipped_df = gpd.overlay(df, gpd.GeoDataFrame(geometry=[msa_geometry], crs=cbsa.crs), how='intersection')

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
    legend_bool = True
    sc = clipped_df.plot(ax=ax, column=col, cmap=cmap, legend=False, vmin=vmin, vmax=vmax, zorder=-1, edgecolor="none", antialiased=False)

    # Create the ScalarMappable for the color bar
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Required for colorbar
    
    
    # Add the CBSA boundary
    ax.add_geometries([msa_geometry], crs=target_proj, edgecolor=state_color, facecolor='none', linewidth=1.5)
    
    # Add title
    ax.set_title(title, fontsize = title_size)
    ax.axis("off")
    return sc, sm



def plot_MSA_mortality(df, MSA, title, validation_df, col="PM25", cmap="viridis", vs=False, counties=False, difference=False, state_color="white", ax=None, title_size = 20):
    # Define the target projection
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
    
    # Generate map features
    resol = '50m'  # use data at this scale
    bodr = cfeature.NaturalEarthFeature(category='cultural', edgecolor='black', name='admin_0_boundary_lines_land', scale='50m', facecolor='none', alpha=0.7)
    ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
    lakes = cfeature.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
    states_provinces = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', edgecolor=state_color, facecolor='none')
    
    # # Customize plot appearance
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    # # Add features to the map
    if counties:
        ax.add_feature(COUNTIES, facecolor='none', edgecolor=state_color)
    # ax.add_feature(ocean, linewidth=0.2, zorder=0)
    # ax.add_feature(lakes, linewidth=0.2, zorder=0)
    # ax.add_feature(states_provinces, edgecolor=state_color, linewidth=0.5, zorder=0)
    # msa_bounds = df.loc[df["CBSA"] == MSA].total_bounds  # [minx, miny, maxx, maxy]

    ax.add_feature(ocean, linewidth=0, zorder=0, facecolor = "white")
    ax.add_feature(lakes, linewidth=0, zorder=0, facecolor = "white")
    # ax.add_feature(states_provinces, edgecolor=state_color, linewidth=0.5, zorder=0)
    
    # Get the geometry for the specified MSA
    msa_geometry = cbsa.loc[cbsa["NAME"] == MSA, "geometry"].iloc[0]

    sub_df = df.loc[df['CBSA'] == MSA]
    
    # Drop rows where geometry is now None (no intersection)
    try:
        clipped_df = sub_df.clip(mask=msa_geometry)
    except Exception as e:
        clipped_df = sub_df
    
    if vs:
        vmin = vs[0]
        vmax = vs[1]
    elif not difference:
        vmin = np.quantile(clipped_df[col], 0.025)
        vmax = np.quantile(clipped_df[col], 0.975)
    else:
        vmin = -max( [abs(np.quantile(clipped_df[col], 0.025)), abs(np.quantile(clipped_df[col], 0.975))] )
        vmax =-vmin
        
    # # Plot MSA area data
    sc = clipped_df.plot(ax=ax, column=col, cmap=cmap, legend=False, vmin=vmin, vmax=vmax, zorder=-1, legend_kwds={'label': col})
    
    # Create the ScalarMappable for the color bar
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(cmap=cmap, norm=norm)

    ax.set_title(title, fontsize = title_size)
    ax.axis("off")
    return sc, sm

pm_mortality_data_df_2020['Attributable_Mortality_Rate_100000'] = pm_mortality_data_df_2020['Attributable_Mortality_Rate'] * 100000


pmgdf['PM25_total_tons'] = rwc_2020_new_surg['PM25_total_tons'][0,0,:,:].to_numpy().flatten()



###
### NYC emissions vs conceontrations plot
###

ny_df = pmgdf.loc[pmgdf['CBSA'] == cbsa_list[0]]
ny_df.loc[ny_df['PM25_RWC_2020'] >= 2].shape[0]/ny_df.shape[0] * 100

import matplotlib.pyplot as plt

# Extract variables
x = ny_df['PM25_total_tons'].values
y = ny_df['PM25_RWC_2020'].values

# Fit regression line through origin: minimize ||y - mx||²
slope = np.sum(x * y) / np.sum(x**2)
reg_line = slope * x  # No intercept

# Create plot
plt.figure(figsize=(7, 6), dpi=300)
plt.scatter(x, y, alpha=0.7, edgecolors='k', s = 100)

# Plot regression line through origin
plt.plot(x, reg_line, 'r--', linewidth=2, label=f'Regression Line: y = {slope:.2f}x')

# Titles and labels
plt.title('C Emissions vs. Concentrations', fontsize=23)
plt.xlabel('PM₂.₅ Tons Emitted', fontsize=23)
plt.ylabel('PM₂.₅ Concentration (µg/m³)', fontsize=23)
plt.ylim([0, 5.5])

# Set tick font sizes
plt.xticks(fontsize=23)
plt.yticks(fontsize=23)
#plt.legend()

plt.grid(True)
plt.tight_layout()
plt.show()


### Plot cbsa emissions and concentrations on a plot

for MSA in cbsa_list[0:1]:
    df = pmgdf

    upper_bound = 5
    if df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(1) >= 5:
        upper_bound = np.ceil(df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(1))

    emissions_upper_bound = 1
    if df.loc[df['CBSA'] == MSA]['PM25_total_tons'].quantile(1) >= 1:
        emissions_upper_bound = np.ceil(df.loc[df['CBSA'] == MSA]['PM25_total_tons'].quantile(1))
    
    # Create a Matplotlib figure with 1 row and 3 columns
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 8), subplot_kw={'projection': target_proj}, dpi = 300)
    plt.subplots_adjust(wspace=0.1)

    place = MSA.split('-')[0]
    
    sc, sm = plot_MSA(
            df=df,
            col="PM25_total_tons",
            cbsa = cbsa,
            MSA=MSA,
            title=f'A {place} RWC Emissions',
            cmap="cividis_r",
            vs=[0, emissions_upper_bound],
            counties=False,
            difference=False,
            ax=axes[0],
            title_size = 20
        )
    colorbar_title_size = 20
    colorbar_tick_size = 20
    
    cbar = fig.colorbar(sm, ax=axes[0], orientation="vertical", shrink=0.6)
    cbar.ax.tick_params(labelsize=20)
    cbar.ax.set_ylabel("PM₂.₅ Tons Emitted", rotation=90, fontsize=20, labelpad = 7)


    sc, sm = plot_MSA(
            df=df,
            col="PM25_RWC_2020",
            cbsa = cbsa,
            MSA=MSA,
            title=f'B {place} RWC Concentrations',
            cmap="magma_r",
            vs=[0, upper_bound],
            counties=False,
            difference=False,
            ax=axes[1],
            title_size = 20
        )
    cbar = fig.colorbar(sm, ax=axes[1], orientation="vertical", shrink=0.6)
    cbar.ax.tick_params(labelsize=20)
    cbar.ax.set_ylabel("PM₂.₅ Concentration (µg/m³)", rotation=90, fontsize=20, labelpad = 7)   
    

# Plot CBSA emissions - concentrations - mortality plots
for MSA in cbsa_list[1:3]:
    df = pmgdf

    upper_bound = 5
    if df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(0.98) >= 5:
        upper_bound = np.ceil(df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(0.98))

    emissions_upper_bound = 1
    if df.loc[df['CBSA'] == MSA]['PM25_total_tons'].quantile(0.98) >= 1:
        emissions_upper_bound = np.ceil(df.loc[df['CBSA'] == MSA]['PM25_total_tons'].quantile(0.98))
    
    # Create a Matplotlib figure with 1 row and 3 columns
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(28, 28), subplot_kw={'projection': target_proj})
    
    place = MSA.split('-')[0]
    
    sc, sm = plot_MSA(
            df=df,
            col="PM25_total_tons",
            cbsa = cbsa,
            MSA=MSA,
            title=f'A {place} RWC Emissions',
            cmap="cividis_r",
            vs=[0, emissions_upper_bound],
            counties=False,
            difference=False,
            ax=axes[0],
            title_size = 24
        )
    colorbar_title_size = 22
    colorbar_tick_size = 22
    
    cbar = fig.colorbar(sm, ax=axes[0], orientation="vertical", shrink=0.25)
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_ylabel("PM₂.₅ Tons Emitted", rotation=90, fontsize=22, labelpad = 7)


    sc, sm = plot_MSA(
            df=df,
            col="PM25_RWC_2020",
            cbsa = cbsa,
            MSA=MSA,
            title=f'B {place} RWC Concentrations',
            cmap="magma_r",
            vs=[0, upper_bound],
            counties=False,
            difference=False,
            ax=axes[1],
            title_size = 24
        )
    cbar = fig.colorbar(sm, ax=axes[1], orientation="vertical", shrink=0.25)
    cbar.ax.tick_params(labelsize=22)
    cbar.ax.set_ylabel("PM₂.₅ Concentration (µg/m³)", rotation=90, fontsize=22, labelpad = 7)   

    sc, sm = plot_MSA_mortality(
        df=pm_mortality_data_df_2020,
        col="Attributable_Mortality_Rate_100000",
        MSA=MSA,
        title=f'C Attributable Mortality Rate',
        validation_df=pd.DataFrame(),  # No validation points for difference
        cmap="viridis_r",
        vs=False,
        counties=False,
        difference=False,
        ax=axes[2],
        state_color="black",# Pass the third axis,
        title_size = 24
    )
    
    colorbar_title_size = 22
    colorbar_tick_size = 22
    
    cbar = fig.colorbar(sm, ax=axes[2], orientation="vertical", shrink=0.25)
    cbar.ax.tick_params(labelsize=22)
    cbar.ax.set_ylabel("Deaths per 100,000", rotation=90, fontsize=22, labelpad = 7)
    plt.show()



# ====================================================================================
# Emissions vs. Concentrations Correlations Table

PWC_df_emissions = gpd.GeoDataFrame.merge(pmgdf, population_df, left_on=['row', 'col'], right_on=['ROWS', 'COLS'])
PWC_df_emissions['population'] = PWC_df_emissions[6]
PWC_df_emissions['PM25_total_tons_per_person'] = PWC_df_emissions['PM25_total_tons']/PWC_df_emissions['population']

from scipy.stats import pearsonr

# Prepare an empty list to store results
results = []

# Loop through each CBSA in the list
for MSA in cbsa_list:
    # Filter the DataFrame for the current CBSA
    cbsa_df = PWC_df_emissions.loc[PWC_df_emissions['CBSA'] == MSA]
    
    # Flatten the emissions and concentrations for this CBSA
    emissions = cbsa_df['PM25_total_tons'].values.flatten()
    concentrations = cbsa_df['PM25_RWC_2020'].values.flatten()
    
    
    # Check if both arrays have the same length to avoid errors
    if len(emissions) == len(concentrations):
        # Calculate the Pearson correlation coefficient
        correlation_coefficient, _ = pearsonr(emissions, concentrations)
        
        # Store the result in a dictionary or list
        results.append({'CBSA': MSA, 
                        'Emissions': round(emissions.sum()), 
                        'emissions_per_area': round(1000 * emissions.sum()/(cbsa.loc[cbsa['NAME'] == MSA].geometry.area.sum() / 1000000),1),
                        'emissions_per_capita' : round(emissions.sum()/cbsa_df['population'].sum() * 100000,2),
                        'spatial correlation': round(correlation_coefficient,2)})
    else:
        print(f"Warning: Data mismatch for CBSA {MSA}. Emissions and concentrations arrays have different lengths.")

# Convert results to DataFrame for better analysis
results_df = pd.DataFrame(results)
results_df = results_df.set_index('CBSA')
print("emissions vs. concentrations table", results_df)


# ====================================================================================
# Decile Plots

# Define colors matching the image
race_colors = {
    "White": "#AFCFEA",
    "Black": "#1F4E9C",
    "American Indian": "#A2CD73",
    "Asian": "#1E792C",
    "Native Hawaiian or Pacific Islander": "#F4A9A3",
    "Other": "#D7261E",
    "Two or More Races": "#F8C471",
    "Hispanic": "#E67E22"
}

def decile_plot(df_old, loc, col, version="baseline", ax=None, first = False, year = "2014", letter = ''):
    """
    Creates a figure with a horizontal decile bar plot showing the racial composition 
    within deciles of the specified data column, with tick marks between bars labeled 
    with the decile boundaries.
    
    Parameters:
    - df: DataFrame containing the data.
    - loc: Location filter for the 'CBSA Title' column. If empty, the whole DataFrame is used.
    - col: Column to create decile plots for.
    - version: Version label to include in the plot title.
    """
    
    # Copy and filter data if loc is specified
    place = letter + loc.split("-")[0]
    print(place)

    df = df_old.copy(deep=True)
    if loc == "CONUS":
        df = df
    elif loc:
        df = df[df['CBSA'] == loc]
    
    # Drop rows with NaN values in the specified column
    df = df.dropna(subset=[col])

    # Define race columns and number of deciles
    race_columns = ["White", "Black", "American Indian", "Asian", "Native Hawaiian or Pacific Islander",
                    "Other", "Two or More Races", "Hispanic"]
    num_iles = 10
    
    # Create decile bins and apply them to the DataFrame
    bins = pd.qcut(df[col], q=num_iles, retbins=True, duplicates='drop') 
    df['decile'] = pd.qcut(df[col], q=num_iles, labels=False, duplicates='drop')
    
    # Get decile ranges from the bins
    decile_ranges = bins[1]
    if col == "Attributable_Mortality_Rate":
        decile_ranges = decile_ranges * 100000

    tick_labels = [f"{value:.1f}" for value in decile_ranges]
    
    # Calculate race percentages within each decile
    for race in race_columns:
        df[race + '_percentage'] = (df[race] / df['Total']) * 100

    # Plot the horizontal decile data
    bottom = np.zeros(len(decile_ranges) - 1)
    for race in race_columns:
        race_percentages = [df[df['decile'] == decile][race + '_percentage'].mean() for decile in range(num_iles)]
        ax.barh(range(len(race_percentages)), race_percentages, left=bottom, height=0.9, label=race, color=race_colors[race])
        bottom += np.array(race_percentages)
    
    measure = ""
    if col == "Attributable_Mortality_Rate":
        measure = "Mortality"
        label = "RWC PM₂.₅ Attr. Mort. Rate Deciles\n(Deaths per 100,000)"
    if col == "2016_pm2.5":
        measure = "PM2.5"
        label = "RWC PM₂.₅ Concentration Deciles\n(µg/m³) "

    ax.tick_params(axis='x', labelsize=16)
    
    ax.set_title(f"{place}", fontsize=20)

    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_ylim(-0.45, len(race_percentages) - 0.55)  # Adjust to fit the bars tightly
    ax.set_xlim(0, 100)
    
    ax.set_xticks([0, 25, 50, 75, 100])
    #ax.set_xlabel('Percentage of Population', fontsize=14)
    
    if first:
        ax.set_ylabel(label, labelpad=45, fontsize=16)

    ax.text(0.75, -0.44, f"{tick_labels[0]}-", ha='right', va='center', fontsize=16, color="black")
        
    # Adjust the first and last bars
    for idx, label in enumerate(tick_labels):
        if idx < len(tick_labels) - 1 and idx > 0:
            ax.text(0.75, idx - 0.5, f"{label}-", ha='right', va='center', fontsize=16, color="black")
    
    # Adjust last bar (if needed)
    ax.text(0.75, len(tick_labels) - 1 - 0.55, f"{tick_labels[-1]}-", ha='right', va='center', fontsize=16, color="black")


import matplotlib.patches as patches


# Example usage for plotting on subplots:
fig, axs = plt.subplots(2, 4, figsize=(14, 12))
plt.tight_layout(pad=2.0, rect=[0, 0, 1.0, 1.0])  # Adjust layout to remove padding
plt.subplots_adjust(hspace=0.2, wspace=0.35)  # Adjust horizontal and vertical spacing
plt.subplots_adjust(bottom=0.15)  # Adjust bottom padding to 0.75 inches

#Call decile_plot for each subplot
decile_plot(pm_mortality_data_df_2020, "CONUS", "Attributable_Mortality_Rate", letter = "F ", ax=axs[1,0], year = "2020", first = True)
decile_plot(pm_mortality_data_df_2020, "Los Angeles-Long Beach-Anaheim, CA", "Attributable_Mortality_Rate", letter = "G ", ax=axs[1, 1], year = "2020")
decile_plot(pm_mortality_data_df_2020, 'Chicago-Naperville-Elgin, IL-IN-WI', "Attributable_Mortality_Rate", letter = "H ",ax=axs[1, 2], year = "2020")
decile_plot(pm_mortality_data_df_2020, "New York-Newark-Jersey City, NY-NJ-PA", "Attributable_Mortality_Rate", letter = "I ",ax=axs[1, 3], year = "2020")

decile_plot(pm_mortality_data_df_2020, "CONUS", "2016_pm2.5",  letter = "B ",ax=axs[0, 0], year = "2020", first = True)
decile_plot(pm_mortality_data_df_2020, "Los Angeles-Long Beach-Anaheim, CA", "2016_pm2.5",  letter = "C ",ax=axs[0, 1],  year = "2020")
decile_plot(pm_mortality_data_df_2020, 'Chicago-Naperville-Elgin, IL-IN-WI', "2016_pm2.5",  letter = "D ",ax=axs[0, 2], year = "2020")
decile_plot(pm_mortality_data_df_2020, "New York-Newark-Jersey City, NY-NJ-PA", "2016_pm2.5",  letter = "E ",ax=axs[0, 3], year = "2020")

fig.text(0.5, 0.04, "Racial/Ethnic Composition", ha='center', fontsize=20)
plt.legend(loc='upper center', bbox_to_anchor=(-1.5, -0.38), ncol=4, fontsize=16, borderaxespad=0.)

# Create an additional row of axes below the main plots for racial composition bars
comp_axes = [fig.add_axes([ax.get_position().x0, 0.08, ax.get_position().width, 0.035]) for ax in axs[1]]

# Create a separate axes for the label on the left
#label_ax = fig.add_axes([axs[1, 0].get_position().x0 - 0.12, 0.08, 0.1, 0.035])
label_ax = fig.add_axes([axs[1, 0].get_position().x0 - 0.015, 0.08, 0.1, 0.035])  # Shift right by ~1-2 inches

label_ax.set_xlim(0, 1)
label_ax.set_ylim(0, 1)
label_ax.text(0, 0.5, "Domain\nRacial/Ethnic\nComposition", ha='right', va='center', fontsize=16)
label_ax.set_frame_on(False)
label_ax.set_xticks([])
label_ax.set_yticks([])

# Loop through the locations and create the racial composition bars
for ax, loc in zip(comp_axes, ["CONUS", "Los Angeles-Long Beach-Anaheim, CA", "Chicago-Naperville-Elgin, IL-IN-WI", 
                                "New York-Newark-Jersey City, NY-NJ-PA"]):
    
    df = pm_mortality_data_df_2020.copy()
    if loc != "CONUS":
        df = df[df['CBSA'] == loc]
    
    # Calculate total racial composition
    race_columns = ["White", "Black", "American Indian", "Asian", "Native Hawaiian or Pacific Islander",
                    "Other", "Two or More Races", "Hispanic"]
    total_pop = df["Total"].sum()
    race_percentages = [df[race].sum() / total_pop * 100 for race in race_columns]
    
    # Plot a single stacked bar for the racial composition
    bottom = np.zeros(1)
    for race, percentage in zip(race_columns, race_percentages):
        ax.barh([0], [percentage], left=bottom, height=1, color=race_colors[race])
        bottom += np.array([percentage])
    
    #rect = patches.Rectangle((0, -0.5), 103, 1, linewidth=1, edgecolor='black', facecolor='none')

    #ax.add_patch(rect)
    ax.set_xlim(0, 100)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)  # Remove frame

#fig.savefig("main_decile_with_legend.png", format="png", dpi=600, bbox_inches='tight')


# More MSA decile plots  --------------------------------------------------------------------------------------


def plot_msa_deciles(msas, start_letter):
    fig, axs = plt.subplots(2, 4, figsize=(14, 12))
    plt.tight_layout(pad=2.0, rect=[0, 0, 1.0, 1.0])
    plt.subplots_adjust(hspace=0.2, wspace=0.35)
    plt.subplots_adjust(bottom=0.15)
    
    for i, msa in enumerate(msas):
        letter1 = chr(ord(start_letter) + i) + "1 "
        letter2 = chr(ord(start_letter) + i) + "2 "
        decile_plot(pm_mortality_data_df_2020, msa, "Attributable_Mortality_Rate", letter=letter1, ax=axs[0, i], first=(i==0), year="2020")
        decile_plot(pm_mortality_data_df_2020, msa, "2016_pm2.5", letter=letter2, ax=axs[1, i], first=(i==0), year="2020")
    
    fig.text(0.5, 0.04, "Racial/Ethnic Composition", ha='center', fontsize=20)
    
    comp_axes = [fig.add_axes([ax.get_position().x0, 0.08, ax.get_position().width, 0.035]) for ax in axs[1]]
    label_ax = fig.add_axes([axs[1, 0].get_position().x0 - 0.015, 0.08, 0.1, 0.035])
    label_ax.set_xlim(0, 1)
    label_ax.set_ylim(0, 1)
    label_ax.text(0, 0.5, "Domain\nRacial/Ethnic\nComposition", ha='right', va='center', fontsize=16)
    label_ax.set_frame_on(False)
    label_ax.set_xticks([])
    label_ax.set_yticks([])
    
    for ax, loc in zip(comp_axes, msas):
        df = pm_mortality_data_df_2020.copy()
        if loc != "CONUS":
            df = df[df['CBSA'] == loc]
        
        race_columns = ["White", "Black", "American Indian", "Asian", "Native Hawaiian or Pacific Islander", "Other", "Two or More Races", "Hispanic"]
        total_pop = df["Total"].sum()
        race_percentages = [df[race].sum() / total_pop * 100 for race in race_columns]
        
        bottom = np.zeros(1)
        for race, percentage in zip(race_columns, race_percentages):
            ax.barh([0], [percentage], left=bottom, height=1, color=race_colors[race])
            bottom += np.array([percentage])
        
        ax.set_xlim(0, 100)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)

plot_msa_deciles(cbsa_list[0:4], "A")
plot_msa_deciles(cbsa_list[4:8], "E")
plot_msa_deciles(cbsa_list[8:12], "I")
plot_msa_deciles(cbsa_list[12:16], "N")
plot_msa_deciles(cbsa_list[16:], "R")


# ====================================================================================
# Relative Disparity Plot

race_colors = {
    "White": "#AFCFEA",
    "Black": "#1F4E9C",
    "American Indian": "#A2CD73",
    "Asian": "#1E792C",
    "Native Hawaiian or Pacific Islander": "#F4A9A3",
    "Other": "#D7261E",
    "Two or More Races": "#F8C471",
    "Hispanic": "#E67E22"
}

df = pm_mortality_data_df_2020

# Define race columns and number of deciles
race_columns = ["White", "Black", "American Indian", "Asian", "Native Hawaiian or Pacific Islander",
                "Other", "Two or More Races", "Hispanic"]

race_decile_averages = {}

#The DataFrame 'df' has the necessary columns for the calculations
for race in race_columns:
    race_deaths_total = (df['Attributable_Mortality']  * df[race] / df['Total']).sum()
    race_mortality_rate = race_deaths_total / df[race].sum() * 100000

    pm_PWC = (df['2016_pm2.5']  * df[race]).sum()/ df[race].sum() # population weighted average RWC PM by race

    bmr_PW = (df['overall_mortality_rate_over_25'] * df[race]).sum()/ df[race].sum()
 
    #print(race, race_mortality_rate, pm_PWC, bmr_PW)
    print(f"{race:<40} {race_mortality_rate:.2f} {pm_PWC:.2f} {bmr_PW*100:.2f}")

    race_decile_averages[race] = {
        "PM": pm_PWC,
        "BMR": bmr_PW,
        "PM_mortality": race_mortality_rate
    }
    

# Calculate CONUS averages
CONUS_PM_avg = df['2016_pm2.5'].mean()
CONUS_BMR_avg = df['overall_mortality_rate_over_25'].mean()
CONUS_PM_mortality_avg = df['Attributable_Mortality_Rate_100000'].mean()

# Compare race-specific decile averages with CONUS averages
comparison = {}

for race, averages in race_decile_averages.items():
    comparison[race] = {
        "PM_diff": (averages["PM"] - CONUS_PM_avg) / CONUS_PM_avg * 100,
        "BMR_diff": (averages["BMR"] - CONUS_BMR_avg)/ CONUS_BMR_avg * 100,
        "PM_mortality_diff": (averages["PM_mortality"] - CONUS_PM_mortality_avg) / CONUS_PM_mortality_avg * 100
    }


race_colors = {
    "White": "#AFCFEA",
    "Black": "#1F4E9C",
    "American Indian": "#A2CD73",
    "Asian": "#1E792C",
    "Native Hawaiian or Pacific Islander": "#F4A9A3",
    "Other": "#D7261E",
    "Two or More Races": "#F8C471",
    "Hispanic": "#E67E22"
}

# Extract the data from the comparison dictionary
pm_differences = [comparison[race]["PM_diff"] for race in race_colors.keys()]
bmr_differences = [comparison[race]["BMR_diff"] for race in race_colors.keys()]
pm_mortality_differences = [comparison[race]["PM_mortality_diff"] for race in race_colors.keys()]

# Create the plot
fig, ax = plt.subplots(figsize=(6, 12))

y_pos = np.arange(len(race_colors))

# Plot PM2.5 differences as filled bars
for i, race in enumerate(race_colors.keys()):
    if abs(pm_differences[i]) > abs(bmr_differences[i]):
        ax.barh(y_pos[i], pm_differences[i], color=race_colors[race], height=0.8)
        ax.barh(y_pos[i], bmr_differences[i], edgecolor=race_colors[race], facecolor='white', height=0.8, linewidth=2)
    else:
        ax.barh(y_pos[i], bmr_differences[i], edgecolor=race_colors[race], facecolor='white', height=0.8, linewidth=2)
        ax.barh(y_pos[i], pm_differences[i], color=race_colors[race],  height=0.8)

# Plot PM mortality differences as black dots
ax.scatter(pm_mortality_differences, y_pos, color='black', zorder=3, s= 100)

# Formatting
ax.set_yticks(y_pos)
ax.set_yticklabels(race_colors.keys(), fontsize=20)
ax.axvline(0, color='black', linestyle='dashed', linewidth=1)
ax.set_xlabel("Percentage Difference Compared to CONUS", fontsize=20)
ax.set_title("A Differences in RWC PM₂.₅, Baseline Mortality Rates, \nand RWC PM₂.₅ Attributable Mortality Rates", fontsize=22, pad=20)

# Add grid lines
ax.xaxis.grid(True, linestyle='solid', alpha=0.7)
for i in range(len(race_colors.keys())):
    ax.axhline(i, color='gray', linestyle='solid', linewidth=0.5, alpha = 0.7)

# Format x-axis labels as percentages
ax.set_xticks([-25, 0, 25])
ax.set_xticklabels([f"{x}%" for x in [-25, 0, 25]], fontsize=20)

# Legend
legend_labels = ["RWC PM₂.₅", "BMR", "RWC PM₂.₅ Att. Mort. Rate"]
legend_patches = [
    plt.Rectangle((0, 0), 1, 1, fc="gray"),
    plt.Rectangle((0, 0), 1, 1, edgecolor="black", facecolor="none"),
    plt.Line2D([0], [0], marker="o", color="black", linestyle="None", markersize=15)
]

plt.legend(legend_patches, legend_labels, loc='lower center', bbox_to_anchor=(0.25, -0.2), fontsize=20, ncol = 3)
#fig.savefig("main_relative_disparity.png", format="png", dpi=600, bbox_inches='tight')


# Race Concentrations Correlation Table --------------------------------------------------------------------------------------

# Create an empty list to store results
correlation_results = []

for race in race_columns:
	pm_mortality_data_df_2020[race + '_percentage'] = (pm_mortality_data_df_2020[race] / pm_mortality_data_df_2020['Total']) * 100

# List of demographic columns
demographic_cols = [
    "Non_White_Fraction",
    "Black_percentage",
    "American Indian_percentage",
    "Asian_percentage",
    "Native Hawaiian or Pacific Islander_percentage",
    "Hispanic_percentage"
]

pm_mortality_data_df_2020['Non_White_Fraction'] = 100 - pm_mortality_data_df_2020['White_percentage']


# Target variables for correlation
target_vars = ["2016_pm2.5", "Attributable_Mortality_Rate_100000"]
cbsa_list.append("CONUS")
for msa in cbsa_list:
    if msa == "CONUS":
        MSA_df = pm_mortality_data_df_2020
    else:
        # Filter data for the current MSA
        MSA_df = pm_mortality_data_df_2020.loc[pm_mortality_data_df_2020["CBSA"] == msa]

    # Compute Pearson correlation coefficients
    correlations = MSA_df[demographic_cols + target_vars].corr(method="pearson")

    # Extract relevant correlations
    msa_correlations = {"MSA": msa.split('-')[0]}
    
    for demo in demographic_cols:
        for target in target_vars:
            msa_correlations[f"{demo.split('_')[0]}_vs_{target}"] = correlations.loc[demo, target]

    # Append to list
    correlation_results.append(msa_correlations)

# Convert results into a DataFrame
correlation_df = pd.DataFrame(correlation_results)
correlation_df


# MSA relative disparity plots --------------------------------------------------------------------------------------
def MSA_relative_disparity(df = pm_mortality_data_df_2020, MSA = "CONUS", ax = None):

	race_colors = {
		"White": "#AFCFEA",
		"Black": "#1F4E9C",
		"American Indian": "#A2CD73",
		"Asian": "#1E792C",
		"Native Hawaiian or Pacific Islander": "#F4A9A3",
		"Other": "#D7261E",
		"Two or More Races": "#F8C471",
		"Hispanic": "#E67E22"
	}

	# Define race columns and number of deciles
	race_columns = ["White", "Black", "American Indian", "Asian", "Native Hawaiian or Pacific Islander",
					"Other", "Two or More Races", "Hispanic"]

	race_decile_averages = {}

	#The DataFrame 'df' has the necessary columns for the calculations
	for race in race_columns:
		race_deaths_total = (df['Attributable_Mortality']  * df[race] / df['Total']).sum()
		race_mortality_rate = race_deaths_total / df[race].sum() * 100000

		pm_PWC = (df['2016_pm2.5']  * df[race]).sum()/ df[race].sum() # population weighted average RWC PM by race

		bmr_PW = (df['overall_mortality_rate_over_25'] * df[race]).sum()/ df[race].sum()
	
		#print(race, race_mortality_rate, pm_PWC, bmr_PW)
		#print(f"{race:<40} {race_mortality_rate:.2f} {pm_PWC:.2f} {bmr_PW*100:.2f}")

		race_decile_averages[race] = {
			"PM": pm_PWC,
			"BMR": bmr_PW,
			"PM_mortality": race_mortality_rate
		}
		
	# Calculate domain averages
	domain_PM_avg = df['2016_pm2.5'].mean()
	domain_BMR_avg = df['overall_mortality_rate_over_25'].mean()
	domain_PM_mortality_avg = df['Attributable_Mortality_Rate_100000'].mean()

	# Compare race-specific decile averages with CONUS averages
	comparison = {}

	for race, averages in race_decile_averages.items():
		comparison[race] = {
			"PM_diff": (averages["PM"] - domain_PM_avg) / domain_PM_avg * 100,
			"BMR_diff": (averages["BMR"] - domain_BMR_avg)/ domain_BMR_avg * 100,
			"PM_mortality_diff": (averages["PM_mortality"] - domain_PM_mortality_avg) / domain_PM_mortality_avg * 100
		}
	
	# Extract the data from the comparison dictionary
	pm_differences = [comparison[race]["PM_diff"] for race in race_colors.keys()]
	bmr_differences = [comparison[race]["BMR_diff"] for race in race_colors.keys()]
	pm_mortality_differences = [comparison[race]["PM_mortality_diff"] for race in race_colors.keys()]

	# Create the plot
	if ax is None:
		fig, ax = plt.subplots(figsize=(6, 12))  # Only create new figure if no ax is provided
	#fig.suptitle(MSA.split('-')[0])
	name = MSA.split('-')[0]
	y_pos = np.arange(len(race_colors))

	# Plot PM2.5 differences as filled bars
	for i, race in enumerate(race_colors.keys()):
		if abs(pm_differences[i]) > abs(bmr_differences[i]):
			ax.barh(y_pos[i], pm_differences[i], color=race_colors[race], height=0.8)
			ax.barh(y_pos[i], bmr_differences[i], edgecolor=race_colors[race], facecolor='white', height=0.8, linewidth=2)
		else:
			ax.barh(y_pos[i], bmr_differences[i], edgecolor=race_colors[race], facecolor='white', height=0.8, linewidth=2)
			ax.barh(y_pos[i], pm_differences[i], color=race_colors[race],  height=0.8)

	# Plot PM mortality differences as black dots
	ax.scatter(pm_mortality_differences, y_pos, color='black', zorder=3, s= 100)

	# Formatting
	ax.set_yticks(y_pos)
	ax.set_yticklabels(race_colors.keys(), fontsize=20)
	ax.axvline(0, color='black', linestyle='dashed', linewidth=1)
	ax.set_xlabel("Percentage Difference Compared to Average", fontsize=20)
	#ax.set_title("Differences in RWC PM₂.₅, Baseline Mortality Rates, \nand RWC PM₂.₅ Attributable Mortality Rates", fontsize=22, pad=20)
	ax.set_title(f"{name} RWC PM₂.₅ Relative Disparities", fontsize=22, pad=20)


	# Add grid lines
	ax.xaxis.grid(True, linestyle='solid', alpha=0.7)
	for i in range(len(race_colors.keys())):
		ax.axhline(i, color='gray', linestyle='solid', linewidth=0.5, alpha = 0.7)

	# Format x-axis labels as percentages
	ax.set_xticks([-25, 0, 25])
	ax.set_xticklabels([f"{x}%" for x in [-25, 0, 25]], fontsize=20)

	plt.gcf().set_dpi(300)


# CHANGE SLICE TO INCLUDE ALL PLOTS
for msa in cbsa_list[1:3]:
	MSA_df = pm_mortality_data_df_2020.loc[pm_mortality_data_df_2020['CBSA'] == msa]
	MSA_relative_disparity(df = MSA_df, MSA = msa)
     

# # Create a massinve figure
# fig, axes = plt.subplots(10, 2, figsize=(20, 100))  # Adjust size as needed
# axes = axes.flatten()  # Flatten to iterate easily

# # Loop through each CBSA and its corresponding subplot
# for i, msa in enumerate(cbsa_list):
#     MSA_df = pm_mortality_data_df_2020.loc[pm_mortality_data_df_2020['CBSA'] == msa]
    
#     # Call function but use the current subplot
#     MSA_relative_disparity(df=MSA_df, MSA=msa, ax=axes[i])  

# # Adjust layout for better spacing
# plt.tight_layout()
# plt.subplots_adjust(hspace=0.25)  # Adjust vertical spacing

# # Legend
# legend_labels = ["RWC PM₂.₅", "BMR", "RWC PM₂.₅ Att. Mort. Rate"]
# legend_patches = [
#     plt.Rectangle((0, 0), 1, 1, fc="gray"),
#     plt.Rectangle((0, 0), 1, 1, edgecolor="black", facecolor="none"),
#     plt.Line2D([0], [0], marker="o", color="black", linestyle="None", markersize=20)
# ]

# plt.legend(legend_patches, legend_labels, loc='lower center', bbox_to_anchor=(-0.85, -0.3), fontsize=25, ncol = 3)
# #plt.savefig("../full_MSA_disparity.png", dpi=300, bbox_inches="tight")
