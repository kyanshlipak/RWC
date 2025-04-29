# %% [markdown]
# ## Backend

# %%
#import necessary libraries
import xarray as xr
from netCDF4 import Dataset, MFDataset, num2date
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import Proj
import geopandas as gpd

# %%
#counties shapefile
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
reader = shpreader.Reader('../SMOKE_sensitivity_analyses/counties/countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())

# %%
baseline = xr.open_dataset("avg_201601.nc")
rwc_2020 = xr.open_dataset("avg_201601_2020_RWC.nc")
no_rwc = xr.open_dataset("no_rwc_avg_201601.nc")

percent_RWC = 100*(baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:])/baseline['PM25_TOT'][0,0,:,:]
percent_RWC_2020 = 100*(rwc_2020['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:])/rwc_2020['PM25_TOT'][0,0,:,:]
diff_2020 = rwc_2020['PM25_TOT'][0,0,:,:] - baseline['PM25_TOT'][0,0,:,:]

# %%
# validation info

file_path = 'AQS_data_2016/validation_baseline_hourly.csv'
validation_baseline = pd.read_csv(file_path)

file_path = 'AQS_data_2016/validation_RWC2020_hourly.csv'
validation_2020_RWC = pd.read_csv(file_path)

# %%
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



# %%
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

# %%
# MSA data
cbsa = gpd.read_file('cbsa/tl_2019_us_cbsa.shp', engine="pyogrio")

# Convert the GeoDataFrame to the target CRS
cbsa = cbsa.to_crs("+proj=lcc +lat_1=33 +lat_2=45 +lon_0=-97 +lat_0=40")

# %%
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

# %%
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Assuming `pmgdf` is your GeoDataFrame
# Add a dummy value column if it doesn't already exist
# Define the target projection
from cartopy import crs as ccrs
import cartopy.feature as cfeature

target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

# Create the plot
fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': target_proj})

# Define visualization limits
vmin = PWC_df['population'].quantile(0.025)
vmax = PWC_df['population'].quantile(0.975)


resol = '50m'
ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
lakes = cfeature.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
states_provinces = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', edgecolor='white', facecolor='none')

ax.add_feature(ocean, linewidth=0, facecolor="white")
ax.add_feature(lakes, linewidth=0, facecolor="white")
ax.add_feature(states_provinces, edgecolor='white', linewidth=0.5)

# Plot the GeoDataFrame
sc = PWC_df.plot(
    ax=ax,
    column='population',  # Replace 'value' with your desired column
    cmap='viridis',
    vmin=vmin,
    vmax=vmax,
    legend=False
)

# Add a color bar
norm = Normalize(vmin=vmin, vmax=vmax)
sm = ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])
fig.colorbar(sm, ax=ax, orientation='horizontal', label='PM25')

# Customize the title and axis
ax.set_title("Spatial Heatmap of Variable", fontsize=16)
ax.axis("off")

# Show the plot
plt.show()


# %%
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

result_df.round(2)

# %%
summary_df = result_df

# %%


# %% [markdown]
# ## Species Analyses

# %%
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


# %%
1.9475588 - 0.381070 - 0.500711

# %%
PWC_df['PM25_RWC_2020'].mean()

# %%
rwc_contribution_df['Average  RWC PM2.5 Contribution µg/m^3'].sum()

# %% [markdown]
# ## Domain Map

# %%
import numpy as np
import pyproj

# Projection parameters
proj_params = {'proj': 'lcc',
               'lat_1': 33,
               'lat_2': 45,
               'lon_0': -97,
               'lat_0': 40,
               'a': 6370000,  # Earth radius approximation
               'b': 6370000}

# Define projection transformations
proj_lcc = pyproj.Proj(proj_params)
proj_latlon = pyproj.Proj(proj="latlong", datum="WGS84")

# Coordinates of the origin
xorig = -2292000
yorig = -1584000

# Number of grid cells in x and y directions
num_cells_x = 1155
num_cells_y = 726

# Size of each grid cell (in meters)
cell_size = 4000  # 4km

# Generate boundary points along the edges of the grid
x_range = np.arange(xorig, xorig + (num_cells_x + 1) * cell_size, cell_size)
y_range = np.arange(yorig, yorig + (num_cells_y + 1) * cell_size, cell_size)

# Bottom edge (left to right)
x_boundary = list(x_range)
y_boundary = [yorig] * len(x_boundary)

# Right edge (bottom to top)
x_boundary += [x_range[-1]] * len(y_range)
y_boundary += list(y_range)

# Top edge (right to left)
x_boundary += list(reversed(x_range))
y_boundary += [y_range[-1]] * len(x_range)

# Left edge (top to bottom)
x_boundary += [xorig] * len(y_range)
y_boundary += list(reversed(y_range))

# Convert to latitude and longitude
lon_boundary, lat_boundary = pyproj.transform(proj_lcc, proj_latlon, x_boundary, y_boundary)

# Create a list of boundary coordinates
boundary_coordinates = list(zip(lat_boundary, lon_boundary))

# %%
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cfeature

# Create the figure and axis with PlateCarree projection
plt.figure(figsize=(12, 8))
ax = plt.axes(projection=crs.PlateCarree())

#ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
#lakes = cfeature.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])


# Add land feature for background coloring
ax.stock_img()  # Uses Cartopy's built-in satellite-style imagery

#ax.add_feature(cfeature.LAND, color='lightgray')
ax.add_feature(cfeature.OCEAN, color='lightblue')
ax.add_feature(cfeature.LAKES, color='lightblue')
#ax.add_feature(cfeature.STATES, linewidth=0.5, color = "gray")
ax.add_feature(cfeature.BORDERS, linestyle=":", color = "gray")

# Extract latitudes and longitudes from boundary coordinates
lat_boundary, lon_boundary = zip(*boundary_coordinates)

# Add coastlines and borders for context
ax.coastlines()

# Set axis limits based on boundary coordinates
ax.set_xlim(min(lon_boundary) - 5, max(lon_boundary) + 5 )
ax.set_ylim(min(lat_boundary) - 5, max(lat_boundary) + 5 )

# Plot boundary points as black dots
ax.scatter(lon_boundary, lat_boundary, color='black', s=1, transform=crs.PlateCarree())

# Display the plot
plt.show()


# %%
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cfeature
import rasterio
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter

# Create the figure and axis with PlateCarree projection
plt.figure(figsize=(12, 8))
ax = plt.axes(projection=crs.PlateCarree())

# Load and display the raster image (NE2_50M_SR.tif)
raster_path = "NE2_50M_SR.tif"  
with rasterio.open(raster_path) as src:
    # Read the image data
    img_data = src.read(1)  # Read the first band (if the image is single-band)

    # Transform the image coordinates to match the map projection
    extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

# Display the image on the map
ax.imshow(img_data, origin='upper', extent=extent, transform=crs.PlateCarree(), alpha=1, zorder = 0)

ax.add_feature(cfeature.OCEAN, color='lightblue', zorder = 1)
ax.add_feature(cfeature.LAKES, color='lightblue', zorder = 1)
ax.add_feature(cfeature.BORDERS, color="black", zorder = 2, linewidth = 1)

# Extract latitudes and longitudes from boundary coordinates
lat_boundary, lon_boundary = zip(*boundary_coordinates)

# Add coastlines and borders for context
ax.coastlines()

# Set axis limits based on boundary coordinates
ax.set_xlim(min(lon_boundary) - 5, max(lon_boundary) + 5)
ax.set_ylim(min(lat_boundary) - 5, max(lat_boundary) + 5)

# Plot boundary points as black dots
ax.scatter(lon_boundary, lat_boundary, color='black', s=1, transform=crs.PlateCarree())

# Add gridlines and labels for latitude and longitude
gl = ax.gridlines(crs=crs.PlateCarree(), draw_labels=True, color='gray', linestyle='--', linewidth=0.5)
gl.xlabels_top = False  # Optional: Disable top longitude labels
gl.ylabels_right = False  # Optional: Disable right latitude labels

# Set the format for latitude and longitude labels
gl.xformatter = LongitudeFormatter()
gl.yformatter = LatitudeFormatter()


# Display the plot
plt.show()


# %% [markdown]
# ## Averaging to 12km & 36km

# %%
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


# %%
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

# PWC_df['population'] = PWC_df[6]

# PWC_df['PWC_TOT'] = PWC_df['PM25'] * PWC_df['population']
# PWC_df['PWC_RWC'] = PWC_df['PM25_RWC'] * PWC_df['population']
# PWC_df['PWC_TOT_2020'] = PWC_df['PM25_2020'] * PWC_df['population']
# PWC_df['PWC_RWC_2020'] = PWC_df['PM25_RWC_2020'] * PWC_df['population']


# total_PWC = PWC_df['PWC_TOT'].sum()
# RWC_PWC = PWC_df['PWC_RWC'].sum()

# total_PWC_2020 = PWC_df['PWC_TOT_2020'].sum()
# RWC_PWC_2020 = PWC_df['PWC_RWC_2020'].sum()

# %%
pmgdf_12km['PM25_RWC_2020'].mean()

# %%
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


# %%
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


# %% [markdown]
# ### Merge with census tract data

# %%
merged_gdf = gpd.read_file("../SMOKE_sensitivity_analyses/merged_gdf/merged_gdf.shp")
merged_gdf = merged_gdf.to_crs(pmgdf.crs) # convert to the same SMOKE CRS projection 
merged_shapes = merged_gdf[['GISJOIN','geometry']]
intersection_12km = gpd.overlay(pmgdf_12km, merged_shapes, how='intersection')
intersection_36km = gpd.overlay(pmgdf_36km, merged_shapes, how='intersection')

print(intersection_36km.columns)

# %%
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

# %%
from shapely.geometry import Polygon

# Define a function to calculate the area of a polygon
def calculate_polygon_area(polygon):
    return polygon.area

intersection_12km['fraction'] = intersection_12km.geometry.area/(12000 ** 2) # area of intersection / divided by area of gridcell
intersection_12km['concentrations_per_polygon'] = intersection_12km['PM25_RWC_2020'] * intersection_12km['fraction'] # concentration * area of intersection/gridcell
summed_df_12km = intersection_12km.groupby('GISJOIN')['concentrations_per_polygon'].sum().reset_index() # concentration * census tract area / gridcell area
merged_gdf_12km = merged_gdf.merge(summed_df_12km, on='GISJOIN')

merged_gdf_12km['PM25_RWC_2020'] = merged_gdf_12km['concentrations_per_polygon']/merged_gdf_12km.geometry.area * 12000 ** 2
census_tract_gdf_12km = census_tract_gdf.merge(merged_gdf_12km, on='GISJOIN')
del census_tract_gdf_12km['geometry_y']
census_tract_gdf_12km.geometry = census_tract_gdf_12km['geometry_x']
census_tract_gdf_12km.to_file("2020_rwc_census_tract_pm25_12km.shp")

# %%
from shapely.geometry import Polygon

# Define a function to calculate the area of a polygon
def calculate_polygon_area(polygon):
    return polygon.area

intersection_36km['fraction'] = intersection_36km.geometry.area/(36000 ** 2) # area of intersection / divided by area of gridcell
intersection_36km['concentrations_per_polygon'] = intersection_36km['PM25_RWC_2020'] * intersection_36km['fraction'] # concentration * area of intersection/gridcell
summed_df_36km = intersection_36km.groupby('GISJOIN')['concentrations_per_polygon'].sum().reset_index() # concentration * census tract area / gridcell area
merged_gdf_36km = merged_gdf.merge(summed_df_36km, on='GISJOIN')

merged_gdf_36km['PM25_RWC_2020'] = merged_gdf_36km['concentrations_per_polygon']/merged_gdf_36km.geometry.area * 36000 ** 2
census_tract_gdf_36km = census_tract_gdf.merge(merged_gdf_36km, on='GISJOIN')
del census_tract_gdf_36km['geometry_y']
census_tract_gdf_36km.geometry = census_tract_gdf_36km['geometry_x']
census_tract_gdf_36km.to_file("2020_rwc_census_tract_pm25_36km.shp")

# %%
### 

# %% [markdown]
# ### Basic Numerical Comparisons

# %%
PWC_df['PM25_RWC_2020'].mean()

# %%
pmgdf_12km['PM25_RWC_2020'].mean()

# %%
pmgdf_36km['PM25_RWC_2020'].mean()

# %%
(PWC_df['PM25_RWC_2020'] * PWC_df['population']).sum() /  PWC_df['population'].sum()

# %%
(pmgdf_12km['PM25_RWC_2020'] * pmgdf_12km['population']).sum() /  pmgdf_12km['population'].sum()

# %%
(pmgdf_36km['PM25_RWC_2020'] * pmgdf_36km['population']).sum() /  pmgdf_36km['population'].sum()

# %% [markdown]
# ### Mortality Comparisons

# %%
mortality_data = pd.read_csv("health_data/merged_baseline_mortality_rate.csv")

# %%
# MSA data
cbsa = gpd.read_file('cbsa/tl_2019_us_cbsa.shp', engine="pyogrio")

# Convert the GeoDataFrame to the target CRS
cbsa = cbsa.to_crs(census_tract_gdf_36km.crs)
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

# %%
census_tract_gdf_12km = add_CBSA(census_tract_gdf_12km, cbsa)
census_tract_gdf_36km = add_CBSA(census_tract_gdf_36km, cbsa)

# %%
pm_mortality_data_df_12km = pd.merge(census_tract_gdf_12km, mortality_data, on = "GISJOIN")
pm_mortality_data_df_36km = pd.merge(census_tract_gdf_36km, mortality_data, on = "GISJOIN")

# %%
pm_mortality_data_df_12km['Total'] = pm_mortality_data_df_12km['Total_x']
pm_mortality_data_df_36km['Total'] = pm_mortality_data_df_36km['Total_x']

# %%
def mortality_function(pm_mortality_data_df, relative_risk = 1.08):
    pm_mortality_data_df = pm_mortality_data_df.copy(deep = True)
    pm_mortality_data_df['real_relative_risk'] = relative_risk ** (pm_mortality_data_df['PM25_RWC_2020'] / 10 /4) # divide by four because only winter
    pm_mortality_data_df['AF'] = (pm_mortality_data_df['real_relative_risk'] - 1) / pm_mortality_data_df['real_relative_risk']
    pm_mortality_data_df['Attributable_Mortality'] = pm_mortality_data_df['Total'] * pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']
    pm_mortality_data_df['Attributable_Mortality_Rate'] = pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']
    print("RR:", relative_risk, " - ", round(pm_mortality_data_df['Attributable_Mortality'].sum()))
    return pm_mortality_data_df

def mortality_col_only(pm_mortality_data_df, relative_risk = 1.08):
    pm_mortality_data_df = pm_mortality_data_df.copy(deep = True)
    pm_mortality_data_df['real_relative_risk'] = relative_risk ** (pm_mortality_data_df['PM25_RWC_2020'] / 10 /4) # divide by four because only winter
    pm_mortality_data_df['AF'] = (pm_mortality_data_df['real_relative_risk'] - 1) / pm_mortality_data_df['real_relative_risk']
    pm_mortality_data_df['Attributable_Mortality'] = pm_mortality_data_df['Total'] * pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']
    pm_mortality_data_df['Attributable_Mortality_Rate'] = pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']
    print("RR:", relative_risk, " - ", round(pm_mortality_data_df['Attributable_Mortality'].sum()))
    return pm_mortality_data_df[['Attributable_Mortality', 'Attributable_Mortality_Rate']]


# %%
pm_mortality_data_df_12km['PM25_RWC_2020']

# %%
pm_mortality_data_df_12km = mortality_function(pm_mortality_data_df_12km, relative_risk = 1.08)
pm_mortality_data_df_12km[['Attributable_Mortality_lower', 'Attributable_Mortality_Rate_lower']] = mortality_col_only(pm_mortality_data_df_12km, relative_risk  = 1.06)
pm_mortality_data_df_12km[['Attributable_Mortality_upper', 'Attributable_Mortality_Rate_upper']] = mortality_col_only(pm_mortality_data_df_12km, relative_risk  = 1.09)
_ = mortality_function(pm_mortality_data_df_12km, relative_risk = 1.17)

# %%
pm_mortality_data_df_36km = mortality_function(pm_mortality_data_df_36km, relative_risk = 1.08)
pm_mortality_data_df_36km[['Attributable_Mortality_lower', 'Attributable_Mortality_Rate_lower']] = mortality_col_only(pm_mortality_data_df_36km, relative_risk  = 1.06)
pm_mortality_data_df_36km[['Attributable_Mortality_upper', 'Attributable_Mortality_Rate_upper']] = mortality_col_only(pm_mortality_data_df_36km, relative_risk  = 1.09)
_ = mortality_function(pm_mortality_data_df_36km, relative_risk = 1.17)

# %%
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


urban_rural_stats(pm_mortality_data_df_12km)


# %%
urban_rural_stats(pm_mortality_data_df_36km)


# %%
result = (
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

result

# %%
result = (
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

result

# %% [markdown]
# ### MSA Plot with NYC and others

# %%
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
    print('done')
    ax.axis("off")
    return sc, sm


# %%
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

# %%
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


# %% [markdown]
# # Quick Numerical Analyses

# %%
def quick_numerical(var):
    # total and mean for baseline
    total_baseline = baseline[var].to_numpy().sum()
    mean_baseline = baseline[var].to_numpy().mean()

    total_rwc_2020 = rwc_2020[var].to_numpy().sum()
    mean_rwc_2020 = rwc_2020[var].to_numpy().mean()

    # total and mean for no_rwc
    total_no_rwc = no_rwc[var].to_numpy().sum()
    mean_no_rwc = no_rwc[var].to_numpy().mean()

    print(f"Average {var} concentration for baseline: {str(round(mean_baseline - mean_no_rwc,4))} micrograms per cubic meter")
    print(f"Average {var} concentration for 2020 RWC: {str(round(mean_rwc_2020 - mean_no_rwc,4))} micrograms per cubic meter")
    # print(f"Summed {var} concentration for baseline: {str(round(total_baseline,4))} micrograms per cubic meter")
    # print(f"Summed {var} concentration for no RWC  : {str(round(total_no_rwc,4))} micrograms per cubic meter")
    # print("\n")
    percent_diff = 100 * (total_baseline - total_no_rwc)/total_baseline
    print(f"Percent of RWC contribution for base {var}: {str( round(percent_diff,3) )} %")
    percent_diff_2020 = 100 * (total_rwc_2020 - total_no_rwc)/total_baseline
    print(f"Percent of RWC contribution for 2020 {var}: {str( round(percent_diff_2020,3) )} %")

# %%
quick_numerical("PM25_TOT")

# %%
quick_numerical("PM25_EC")

# %%
quick_numerical("PM25_OC")

# %%
quick_numerical("NOX") # doesn't change much

# %%
quick_numerical("O3") ## doesn't change much

# %%
quick_numerical("CO") # changes a little

# %%
#create dataframe of RWC relative contributions to baseline
variables = list(baseline.variables)
relative_change = pd.DataFrame(columns = ["var", "percent RWC contribution"])
percent_diffs = []

for var in variables:
    percent_diff = 100 * (baseline[var].to_numpy().sum() - no_rwc[var].to_numpy().sum())/baseline[var].to_numpy().sum()
    percent_diffs.append(percent_diff)
    
relative_change['var'] = variables
relative_change['percent RWC contribution'] = percent_diffs

# %%
pd.set_option('display.max_rows', None)
relative_change = relative_change.sort_values("percent RWC contribution")
relative_change.head(124)


# %%
pd.set_option('display.max_rows', 20)

# %%
baseline

# %%
PWC_df

# %% [markdown]
# # Emissions NEI Input Plotting

# %%
my_dir = '../SMOKE_sensitivity_analyses/'
rwc_2016 = pd.read_csv(my_dir +'2016ff.csv',skiprows=29)
rwc_2020 = pd.read_csv(my_dir +'rwc_2020.csv',skiprows=18)

# %%
#dictionary mapping state integer FIPS codes to their 2 letter abbreviations
state_mapping = {
        1: 'AL', 2: 'AK', 4: 'AZ', 5: 'AR', 6: 'CA',
        8: 'CO', 9: 'CT', 10: 'DE', 11:'DC',12: 'FL', 13: 'GA',
        15: 'HI', 16: 'ID', 17: 'IL', 18: 'IN', 19: 'IA',
        20: 'KS', 21: 'KY', 22: 'LA', 23: 'ME', 24: 'MD',
        25: 'MA', 26: 'MI', 27: 'MN', 28: 'MS', 29: 'MO',
        30: 'MT', 31: 'NE', 32: 'NV', 33: 'NH', 34: 'NJ',
        35: 'NM', 36: 'NY', 37: 'NC', 38: 'ND', 39: 'OH',
        40: 'OK', 41: 'OR', 42: 'PA', 44: 'RI', 45: 'SC',
        46: 'SD', 47: 'TN', 48: 'TX', 49: 'UT', 50: 'VT',
        51: 'VA', 53: 'WA', 54: 'WV', 55: 'WI', 56: 'WY'
}

#for a smoke input file, get FIPS codes, state codes
def process_SMOKE_file(rwc):
    #rename column
    rwc = rwc.rename(columns={"poll": "emission"})

    #make sure all FIPS are five digits
    
    rwc.loc[rwc["region_cd"] < 10000, 'region_cd'] = '0' + rwc.loc[rwc["region_cd"] < 10000]['region_cd'].astype(str)

    # Converting the column of five-digit integers to strings and extracting the first two digits
    rwc['state_int'] = rwc['region_cd'].astype(str).str[:2]

    # Converting the extracted first two digits back to integers
    rwc['state_int'] = rwc['state_int'].astype(int)

    # Dictionary mapping integers to state abbreviations
    state_mapping = {
        1: 'AL', 2: 'AK', 4: 'AZ', 5: 'AR', 6: 'CA',
        8: 'CO', 9: 'CT', 10: 'DE', 11:'DC',12: 'FL', 13: 'GA',
        15: 'HI', 16: 'ID', 17: 'IL', 18: 'IN', 19: 'IA',
        20: 'KS', 21: 'KY', 22: 'LA', 23: 'ME', 24: 'MD',
        25: 'MA', 26: 'MI', 27: 'MN', 28: 'MS', 29: 'MO',
        30: 'MT', 31: 'NE', 32: 'NV', 33: 'NH', 34: 'NJ',
        35: 'NM', 36: 'NY', 37: 'NC', 38: 'ND', 39: 'OH',
        40: 'OK', 41: 'OR', 42: 'PA', 44: 'RI', 45: 'SC',
        46: 'SD', 47: 'TN', 48: 'TX', 49: 'UT', 50: 'VT',
        51: 'VA', 53: 'WA', 54: 'WV', 55: 'WI', 56: 'WY'
    }

    # Mapping integers to state abbreviations
    rwc['state'] = rwc['state_int'].map(state_mapping)
    rwc['ann_value'].fillna(0.001, inplace=True)
    rwc = rwc[['country_cd', 'region_cd', 'scc', 'emission', 'ann_value', 'state_int', 'state']]
    return rwc

rwc_2016 = process_SMOKE_file(rwc_2016)
rwc_2020 = process_SMOKE_file(rwc_2020)

# %%
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Polygon

# Load US county shapefile
shapefile_path = my_dir + "counties/tl_2024_us_county.shp"  # Update this path
gdf_counties = gpd.read_file(shapefile_path)
target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))


def plot_counties_matplotlib(dataframe, scc=False, emission="PM25-PRI", key="ann_value", comparison=False, endpoints=False, title = f"Annual Total RWC Primary PM₂.₅ Emissions"):
    # Filter data
    if not scc:
        dataframe = dataframe.loc[dataframe['emission'] == emission].groupby('region_cd')[key].sum().reset_index()
    else:
        dataframe = dataframe.loc[(dataframe["scc"] == scc) & (dataframe["emission"] == emission)]

    dataframe['region_cd'] = dataframe['region_cd'].astype(str)

    # Merge emissions data with county geometries
    merged_gdf = gdf_counties.merge(dataframe, left_on="GEOID", right_on="region_cd", how="left").fillna(0)
    merged_gdf = merged_gdf[merged_gdf.geometry.within(Polygon([[-125, 24], [-66, 24], [-66, 50], [-125, 50]]))]

    print(merged_gdf.shape)
    # Define colormap
    cmap = "coolwarm" if comparison else "viridis"

    # Define color range
    if not endpoints:
        if comparison:
            ninety = max(abs(dataframe[key].quantile(0.025)), abs(dataframe[key].quantile(0.975)))
            endpoints = [-ninety, ninety]
        else:
            endpoints = [0, dataframe[key].quantile(0.9875)]

    # Create figure and axis with projection
    fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={"projection": target_proj})
    resol = '10m'
    lakes = cfeature.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])

    # Add basemap features (transform lakes and ocean to match target_proj)
    ax.add_feature(cfeature.OCEAN, facecolor="white", edgecolor="none", zorder = 1)
    ax.add_feature(lakes, facecolor = 'white', linewidth = 0)

    # Transform merged_gdf to target projection CRS before plotting
    merged_gdf = merged_gdf.to_crs(target_proj.proj4_init)

    # Plot county emissions
    im = merged_gdf.plot(column=key, cmap=cmap, linewidth=0, edgecolor="black", legend=False, vmin=endpoints[0], vmax=endpoints[1], ax=ax, zorder=0)

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=endpoints[0], vmax=endpoints[1]))
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, aspect=20)
    cbar.set_label(f"Primary PM₂.₅ Emissions", fontsize=16)
    cbar.ax.tick_params(labelsize=16)  # Increase tick label size

    # Title and display
    ax.set_title(title, fontsize=18)
    ax.set_extent([-120, -74, 23, 50])
    ax.axis("off")

    plt.show()
    return dataframe, gdf_counties, merged_gdf


# %%
dataframe, gdf_counties, merged_gdf = plot_counties_matplotlib(rwc_2020)

# %%
outdoor_wb = rwc_2020.loc[(rwc_2020['scc'] ==2104008700) & (rwc_2020['emission'] == "PM25-PRI")]['ann_value'].sum()

# %%
sum = rwc_2020.loc[(rwc_2020['emission'] == "PM25-PRI")]['ann_value'].sum()
print(f"Percent total emissions from outdoor woodburning {outdoor_wb/sum * 100}")

# %%
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Polygon

# Define appliance categories and their relevant SCC codes
relevant_sccs = {
    2104009000: "Firelog",
    2104008100: "Fireplace: General",
    2104008210: "Woodstove: fireplace inserts; non-EPA certified ",
    2104008220: "Woodstove: fireplace inserts; EPA certified; non-catalytic ",
    2104008230: "Woodstove: fireplace inserts; EPA certified; catalytic",
    2104008310: "Woodstove: freestanding, non-EPA certified ",
    2104008320: "Woodstove: freestanding, EPA certified, non-catalytic",
    2104008330: "Woodstove: freestanding, EPA certified, catalytic",
    2104008400: "Woodstove: pellet-fired, general (freestanding or FP insert) ",
    2104008510: "Furnace: Indoor, cordwood-fired, non-EPA certified",
    2104008610: "Hydronic heater: outdoor",  
    2104008700: "Outdoor wood burning device",
}

woodstoves = [
    2104008210, 2104008220, 2104008230, 2104008300, 2104008310, 
    2104008320, 2104008330, 2104008400
]

hydronics = [2104008610, 2104008630, 2104008620]
furnaces = [2104008510, 2104008530]
firelogs = [2104009000]
fireplaces = [2104008100]
outdoors = [2104008700]

# Function to aggregate the data based on appliance category
def aggregate(df, scc_array):
    return df.loc[df['scc'].isin(scc_array)].groupby(['region_cd', 'emission'])['ann_value'].sum().reset_index()


# Appliance categories
appliance_categories = [
    (woodstoves, "Woodstove"),
    (hydronics, "Hydronic"),
    (furnaces, "Furnace"),
    (firelogs, "Firelog"),
    (fireplaces, "Fireplace"),
    (outdoors, "Outdoor Wood Burning Device")
]

for (appliance_sccs, appliance_string) in appliance_categories[-1:]:
    agg_data = aggregate(rwc_2020, appliance_sccs)
    agg_data['region_cd'] = agg_data['region_cd'].astype(str)
    _, _, _ = plot_counties_matplotlib(agg_data, title = f"Annual Total {appliance_string} Primary PM₂.₅ Emission")
    

# %%
rwc_2020.loc[rwc_2020['emission'] == "PM25-PRI"]['ann_value'].sum()

# %%
rwc_2020.loc[rwc_2020['emission'] == "PM25-PRI"]['ann_value'].sum()

# %%
emis_state = rwc_2020.loc[rwc_2020['emission'] == "PM25-PRI"].groupby('state')['ann_value'].sum().reset_index().set_index('state').sort_values(by='ann_value', ascending = False)
emis_state['ann_value'] = (emis_state['ann_value'].astype(int))
emis_state

# %% [markdown]
# # PM2.5 Concentrations Plotting

# %%
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

# %%
#Challenge #2 change projection type and 
def plot(data, label, title, lat, lon,  cbar = "viridis", vs = False, counties = False, bounds = False, difference = False, state_color = "white" ):
    if not bounds:
        bounds = [[0, lat.shape[0]], [0, lon.shape[1]]]
    
    #create matplotlib figure
    plt.figure(figsize=(15,8))

    #create subplot with the an Mollweide projection
    ax = plt.subplot(111, projection=ccrs.PlateCarree(central_longitude=-97.0, globe=None))

    #create map with 
    if difference:
            mm = ax.pcolormesh(lon[ bounds[0][0]:bounds[0][1], bounds[1][0]:bounds[1][1] ],
                       lat[ bounds[0][0]:bounds[0][1], bounds[1][0]:bounds[1][1] ],
                       data[ bounds[0][0]:bounds[0][1], bounds[1][0]:bounds[1][1] ],
                       transform=ccrs.PlateCarree(),
                       cmap = "RdBu" #use inferno colormap
                      )
    else:
            mm = ax.pcolormesh(lon[ bounds[0][0]:bounds[0][1], bounds[1][0]:bounds[1][1] ],
                       lat[ bounds[0][0]:bounds[0][1], bounds[1][0]:bounds[1][1] ],
                       data[ bounds[0][0]:bounds[0][1], bounds[1][0]:bounds[1][1] ],
                       transform=ccrs.PlateCarree(),
                       cmap = cbar #use inferno colormap
                      )
    
    #map coastlines on the map with resolution 100m
    ax.coastlines(resolution='110m')

    #generate features
    resol = '10m'  # use data at this scale
    bodr = cfeature.NaturalEarthFeature(category='cultural',  edgecolor = state_color,name='admin_0_boundary_lines_land', scale=resol, facecolor='none', alpha=0.7)
    ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
    lakes = cfeature.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
    #rivers = cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', scale=resol, edgecolor='black', facecolor='none')
    states_provinces = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m',edgecolor = state_color, facecolor='none')
    
    #add features
    if counties:
        ax.add_feature(COUNTIES, facecolor='none', edgecolor=state_color)
   # ax.add_feature(rivers, linewidth=0.5,edgecolor='blue')
    ax.add_feature(ocean, linewidth=0.2 )
    ax.add_feature(lakes, linewidth = 0)
    ax.add_feature(states_provinces, edgecolor=state_color, linewidth=0.5)
    ax.add_feature(bodr, linestyle='--', edgecolor=state_color, alpha=1)
    
    #ax.stock_img(); #land images
    if vs:
        mm.set_clim(vmin=vs[0], vmax=vs[1])
    else:
        mm.set_clim(vmin = np.quantile(data, 0.025), vmax = np.quantile(data, 0.975))

    #plot colorbar legend
    cbar = plt.colorbar(mm)
    cbar.ax.set_ylabel(label,rotation = 90)

    ax.set_title(title)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

# %% [markdown]
# ## Plotting Overall PM2.5 Figure

# %%
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

# %%
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


# %%
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
    summary_stats = summary_df

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


# %%
pos

# %%
fig.savefig("msa_pm25_map.png", format="png", dpi=600, bbox_inches='tight')

# %%
summary_df.loc[summary_df['CBSA'] == MSA]

# %%
for MSA, pos in list(converted_positions.items())[0:2]:
	print(df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].mean())

# %%
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

# Initialize the figure
target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

fig = plt.figure(figsize=(16, 14))  # Larger figure to accommodate everything

# Custom positions (x1, x2, y1, y2) for each subplot, adjusted for add_axes format
positions = {
    "Seattle-Tacoma-Bellevue, WA": (0, 0.2, 0, 0.3),  # Top left
    "Minneapolis-St. Paul-Bloomington, MN-WI": (0.2, 0.52, 0, 0.3),  # Top-center-left
    "Chicago-Naperville-Elgin, IL-IN-WI": (0.52, 0.75, 0, 0.3),  # Top-center-right
    "Boston-Cambridge-Newton, MA-NH": (0.75, 1, 0, 0.4),  # Top-right
    "New York-Newark-Jersey City, NY-NJ-PA": (0.75, 1, 0.4, 0.7),  # Right-middle
    "Philadelphia-Camden-Wilmington, PA-NJ-DE-MD": (0.6, 1, 0.7, 1),  # Bottom-left
    "Washington-Arlington-Alexandria, DC-VA-MD-WV": (0.25, 0.6, 0.7, 1),  # Bottom-right
    "Los Angeles-Long Beach-Anaheim, CA": (0, 0.25, 0.65, 1),  # Right-middle
    "Denver-Aurora-Lakewood, CO": (0, 0.25, 0.3, 0.65),  # Bottom-left-center
    "": (0.25, 0.75, 0.3, 0.7),  # Bottom-center
}

# Convert positions to the format (left, bottom, width, height)
converted_positions = {
    key: (x1, 1 - y2, x2 - x1, y2 - y1)  # Calculate width and height from x1, x2, y1, y2
    for key, (x1, x2, y1, y2) in positions.items()
}

# Loop through each MSA and plot
for MSA, pos in list(converted_positions.items())[0:]:
    df = pmgdf  # Assuming pmgdf is the data frame

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': target_proj})

    # Dynamically determine the upper bound for the color scale
    upper_bound = 5
    if MSA != "":
        if df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(0.99) >= 5:
            bound_2020 = np.ceil(df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(0.99))
            upper_bound = bound_2020
        place = MSA
    else:
        df = PWC_df
        if df['PM25_RWC_2020'].quantile(0.995) >= 5:
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
        title_size=15
    )
    
    # Adjust colorbar
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", shrink=0.7)
    cbar.ax.tick_params(labelsize=14)
    
    # Adjust title position to ensure it fits within the axis
    ax.set_title(MSA, fontsize=18, pad=10)  # Adjust the pad value as necessary
    
    #ax.set_aspect('auto')
    #ax.set_clip_on(True)

# Manually adjust subplot layout to ensure everything fits well
plt.show()


# %% [markdown]
# ## Poster Plotting

# %%
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
    "St. Louis, MO-IL"
]


# Initialize an empty list to hold the data
data = []

# Iterate through each CBSA in the list
for cbsa in cbsa_list:
    cbsa_df = pmgdf[pmgdf['CBSA'] == cbsa]
    pm_2014 = cbsa_df['PM25'].mean()
    pm_2020 = cbsa_df['PM25_2020'].mean()
    pm_2014_RWC = cbsa_df['PM25_RWC'].mean()
    pm_2020_RWC = cbsa_df['PM25_RWC_2020'].mean()
    percent_2014_RWC = pm_2014_RWC / pm_2014 * 100 if pm_2014 != 0 else None
    percent_2020_RWC = pm_2020_RWC / pm_2020 * 100 if pm_2020 != 0 else None
    
    # Append the results as a dictionary
    data.append({
        "CBSA": cbsa,
        "PM25_2014": pm_2014,
        "PM25_2020": pm_2020,
        "PM25_2014_RWC": pm_2014_RWC,
        "PM25_2020_RWC": pm_2020_RWC,
        "Percent_2014_RWC": percent_2014_RWC,
        "Percent_2020_RWC": percent_2020_RWC
    })

# Create a DataFrame from the data
result_df = pd.DataFrame(data)

# Display the resulting DataFrame
result_df

# %%
import matplotlib.pyplot as plt
from pyproj import CRS, Transformer
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import box
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.colors import TwoSlopeNorm


def plot_MSA(df, cbsa, MSA, title, col="PM25", cmap="viridis", vs=False, counties=False, difference=False, state_color="white", ax=None, title_size = 16):
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
    legend_bool = True
    sc = clipped_df.plot(ax=ax, column=col, cmap=cmap, legend=False, vmin=vmin, vmax=vmax, zorder=-1, edgecolor="none", antialiased=False)

    # Create the ScalarMappable for the color bar
    norm = Normalize(vmin=vmin, vmax=vmax)
    #if difference:
    #    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    
        
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Required for colorbar
    
    
    # Add the CBSA boundary
    #ax.add_geometries([msa_geometry], crs=target_proj, edgecolor=state_color, facecolor='none', linewidth=1.5)
    ax.set_title(title, fontsize=title_size)
    # Add title
    #set_centered_title(ax, title, title_size)
    #ax.set_title(title, fontsize=title_size, loc='center', pad=5, x=0.5, y=-1.05)
    #ax.text(0.5, ax.get_ylim()[1] + 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0]), title, ha='center', va='bottom', fontsize=title_size, transform=ax.transAxes)

    ax.axis("off")
    return sc, sm




# %%
# MSA data
cbsa = gpd.read_file('cbsa/tl_2019_us_cbsa.shp', engine="pyogrio")

# Convert the GeoDataFrame to the target CRS
cbsa = cbsa.to_crs("+proj=lcc +lat_1=33 +lat_2=45 +lon_0=-97 +lat_0=40")

for MSA in cbsa_list[1:2]:
    df = pmgdf

    upper_bound = 5
    if df.loc[df['CBSA'] == MSA]['PM25_RWC'].quantile(0.98) >= 5 or  df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(0.98) >= 5:
        bound_2014 = np.ceil(df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(0.98))
        bound_2020 = np.ceil(df.loc[df['CBSA'] == MSA]['PM25_RWC'].quantile(0.98))
        upper_bound = max(bound_2014, bound_2020)
    
    # Create a Matplotlib figure with 1 row and 3 columns
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(30, 30), subplot_kw={'projection': target_proj})
    
    place = MSA.split('-')[0]
    
    sc, sm = plot_MSA(
            df=df,
            col="PM25_RWC",
            cbsa = cbsa,
            MSA=MSA,
            title=f'A RWC PM₂.₅ (NEI 2014)',
            cmap="magma_r",
            vs=[0, upper_bound],
            counties=False,
            difference=False,
            ax=axes[0],
            title_size = 24
        )
    colorbar_title_size = 18
    colorbar_tick_size = 18
    
    cbar = fig.colorbar(sm, ax=axes[0], orientation="vertical", shrink=0.25)
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_ylabel("PM₂.₅ Concentration (µg/m³)", rotation=90, fontsize=18, labelpad = 7)


    sc, sm = plot_MSA(
            df=df,
            col="PM25_RWC_2020",
            cbsa = cbsa,
            MSA=MSA,
            title=f'B RWC PM₂.₅ (NEI 2020)',
            cmap="magma_r",
            vs=[0, upper_bound],
            counties=False,
            difference=False,
            ax=axes[1],
            title_size = 24
        )
    cbar = fig.colorbar(sm, ax=axes[1], orientation="vertical", shrink=0.25)
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_ylabel("PM2.5 Concentration (µg/m³)", rotation=90, fontsize=18, labelpad = 7)   
    
    sc, sm = plot_MSA(
            df=df,
            col="PM25_diff",
            cbsa = cbsa,
            MSA=MSA,
            title=f'C PM₂.₅ Difference (2020 - 2014)',
            cmap="RdBu",
            vs=[-1, 3],
            counties=False,
            difference=True,
            ax=axes[2],
            title_size = 24
        )
    colorbar_title_size = 18
    colorbar_tick_size = 18

    sm.set_clim(-1, 3)  # Set the color limits from -1 to 3

    # Adjust the colorbar with the correct ticks and labeling
    cbar = fig.colorbar(sm, ax=axes[2], orientation="vertical", shrink=0.25)
        # Set the limits of the colorbar to be from -1 to 3 (cut off)

    cbar.ax.tick_params(labelsize=colorbar_tick_size)
    cbar.ax.set_ylabel("PM2.5 Concentration (µg/m³)", rotation=90, fontsize=colorbar_title_size, labelpad=7)

    # cbar = fig.colorbar(sm, ax=axes[2], orientation="vertical", shrink=0.25)
    # cbar.set_ticks([-1, 0, 1, 2, 3])  # Colorbar ticks from -1 to 3
    # cbar.ax.tick_params(labelsize=18)
    # cbar.ax.set_ylabel("PM2.5 Concentration (µg/m³)", rotation=90, fontsize=18, labelpad = 7)
    # plt.show()


# %%
df = pmgdf
MSA='Los Angeles-Long Beach-Anaheim, CA'

# Create a Matplotlib figure with 1 row and 3 columns
target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), subplot_kw={'projection': target_proj})

place = MSA.split('-')[0]

sc = plot_MSA(
        df=df,
        col="PM25_RWC",
        cbsa = cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 from RWC (NEI 2014)',
        cmap="magma",
        vs=[0, 5],
        counties=False,
        difference=False,
        ax=axes,
        title_size = 18
    )
colorbar_title_size = 18
colorbar_tick_size = 18

cbar = axes.get_figure().axes[-1]  # The colorbar is typically the last axis
cbar.set_ylabel(r"$\mu g/m^3$", fontsize=colorbar_title_size)  # Add label
cbar.tick_params(labelsize=colorbar_tick_size)
plt.show()


fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), subplot_kw={'projection': target_proj})

sc = plot_MSA(
        df=df,
        col="PM25_RWC_2020",
        cbsa = cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 from RWC (NEI 2020)',
        cmap="magma",
        vs=[0, 5],
        counties=False,
        difference=False,
        ax=axes,
        title_size = 18
    )

cbar = axes.get_figure().axes[-1]  # The colorbar is typically the last axis
cbar.set_ylabel(r"$\mu g/m^3$", fontsize=colorbar_title_size)  # Add label
cbar.tick_params(labelsize=colorbar_tick_size)
plt.show()

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), subplot_kw={'projection': target_proj})

sc = plot_MSA(
        df=df,
        col="PM25_diff",
        cbsa = cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 Difference (2020 - 2014)',
        cmap="RdBu",
        vs=[-3, 3],
        counties=False,
        difference=True,
        ax=axes,
        title_size = 18
    )
colorbar_title_size = 18
colorbar_tick_size = 18

cbar = axes.get_figure().axes[-1]  # The colorbar is typically the last axis
cbar.set_ylabel(r"$\mu g/m^3$", fontsize=colorbar_title_size)  # Add label
cbar.tick_params(labelsize=colorbar_tick_size)
plt.show()


# %%
df = pmgdf
MSA='Chicago-Naperville-Elgin, IL-IN-WI'

# Create a Matplotlib figure with 1 row and 3 columns
target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), subplot_kw={'projection': target_proj})

place = MSA.split('-')[0]

sc = plot_MSA(
        df=df,
        col="PM25_RWC",
        cbsa = cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 from RWC (NEI 2014)',
        cmap="magma",
        vs=[0, 5],
        counties=False,
        difference=False,
        ax=axes,
        title_size = 18
    )
colorbar_title_size = 18
colorbar_tick_size = 18

cbar = axes.get_figure().axes[-1]  # The colorbar is typically the last axis
cbar.set_ylabel(r"$\mu g/m^3$", fontsize=colorbar_title_size)  # Add label
cbar.tick_params(labelsize=colorbar_tick_size)
plt.show()


fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), subplot_kw={'projection': target_proj})

sc = plot_MSA(
        df=df,
        col="PM25_RWC_2020",
        cbsa = cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 from RWC (NEI 2020)',
        cmap="magma",
        vs=[0, 5],
        counties=False,
        difference=False,
        ax=axes,
        title_size = 18
    )

cbar = axes.get_figure().axes[-1]  # The colorbar is typically the last axis
cbar.set_ylabel(r"$\mu g/m^3$", fontsize=colorbar_title_size)  # Add label
cbar.tick_params(labelsize=colorbar_tick_size)
plt.show()

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), subplot_kw={'projection': target_proj})

sc = plot_MSA(
        df=df,
        col="PM25_diff",
        cbsa = cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 Difference (2020 - 2014)',
        cmap="RdBu",
        vs=[-3, 3],
        counties=False,
        difference=True,
        ax=axes,
        title_size = 18
    )
colorbar_title_size = 18
colorbar_tick_size = 18

cbar = axes.get_figure().axes[-1]  # The colorbar is typically the last axis
cbar.set_ylabel(r"$\mu g/m^3$", fontsize=colorbar_title_size)  # Add label
cbar.tick_params(labelsize=colorbar_tick_size)
plt.show()


# %%
MSA = 'New York-Newark-Jersey City, NY-NJ-PA'
# Create a Matplotlib figure with 1 row and 3 columns
target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), subplot_kw={'projection': target_proj})

place = MSA.split('-')[0]

sc = plot_MSA(
        df=df,
        col="PM25_RWC",
        cbsa = cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 from RWC (NEI 2014)',
        cmap="magma",
        vs=[0, 7],
        counties=False,
        difference=False,
        ax=axes,
        title_size = 18
    )
colorbar_title_size = 18
colorbar_tick_size = 18

cbar = axes.get_figure().axes[-1]  # The colorbar is typically the last axis
cbar.set_ylabel(r"$\mu g/m^3$", fontsize=colorbar_title_size)  # Add label
cbar.tick_params(labelsize=colorbar_tick_size)
plt.show()


fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), subplot_kw={'projection': target_proj})

sc = plot_MSA(
        df=df,
        col="PM25_RWC_2020",
        cbsa = cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 from RWC (NEI 2020)',
        cmap="magma",
        vs=[0, 7],
        counties=False,
        difference=False,
        ax=axes,
        title_size = 18
    )

cbar = axes.get_figure().axes[-1]  # The colorbar is typically the last axis
cbar.set_ylabel(r"$\mu g/m^3$", fontsize=colorbar_title_size)  # Add label
cbar.tick_params(labelsize=colorbar_tick_size)
plt.show()

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), subplot_kw={'projection': target_proj})

sc = plot_MSA(
        df=df,
        col="PM25_diff",
        cbsa = cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 Difference (2020 - 2014)',
        cmap="RdBu",
        vs=[-3, 3],
        counties=False,
        difference=True,
        ax=axes,
        title_size = 18
    )
colorbar_title_size = 18
colorbar_tick_size = 18

cbar = axes.get_figure().axes[-1]  # The colorbar is typically the last axis
cbar.set_ylabel(r"$\mu g/m^3$", fontsize=colorbar_title_size)  # Add label
cbar.tick_params(labelsize=colorbar_tick_size)
plt.show()


# %%
# MSA data
cbsa = gpd.read_file('cbsa/tl_2019_us_cbsa.shp', engine="pyogrio")

# Convert the GeoDataFrame to the target CRS
cbsa = cbsa.to_crs("+proj=lcc +lat_1=33 +lat_2=45 +lon_0=-97 +lat_0=40")

# %%
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# List of MSAs to plot
MSAs = [
    'Los Angeles-Long Beach-Anaheim, CA',
    'New York-Newark-Jersey City, NY-NJ-PA',
    'Chicago-Naperville-Elgin, IL-IN-WI'
]

# Create a Matplotlib figure with 3 rows and 3 columns
target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(30, 30), subplot_kw={'projection': target_proj})

# Adjust spacing between subplots
fig.subplots_adjust(hspace=-0.05, wspace=0.05)

df = pmgdf

# Iterate through MSAs and plot data
for i, MSA in enumerate(MSAs):
    place = MSA.split('-')[0]

    if place == "Los Angeles": place = "LA"
    if place == "New York": place = "NYC"

    # Plot for 2014 PM2.5 data
    sc, sm = plot_MSA(
        df=df,
        col="PM25_RWC",
        cbsa=cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 from RWC (NEI 2014)',
        cmap="magma_r",
        vs=[0, 5],
        counties=False,
        difference=False,
        ax=axes[i, 0],  # Use the first column for 2014 data
        title_size=24
    )
    #Add color bar
    cbar = fig.colorbar(sm, ax=axes[i, 0], orientation="vertical", shrink=0.8)
    cbar.ax.tick_params(labelsize=22)
    cbar.ax.set_ylabel("PM2.5 Concentration (µg/m³)", rotation=90, fontsize=22, labelpad = 7)
    
    # Plot for 2020 PM2.5 data
    sc, sm = plot_MSA(
        df=df,
        col="PM25_RWC_2020",
        cbsa=cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 from RWC (NEI 2020)',
        cmap="magma_r",
        vs=[0, 5],
        counties=False,
        difference=False,
        ax=axes[i, 1],  # Use the second column for 2020 data
        title_size=24
    )
    # Add color bar
    cbar = fig.colorbar(sm, ax=axes[i, 1], orientation="vertical", shrink=0.8)
    cbar.ax.tick_params(labelsize=22)
    cbar.ax.set_ylabel("PM2.5 Concentration (µg/m³)", rotation=90, fontsize=22, labelpad = 7)

    # Plot for PM2.5 difference data
    sc, sm  = plot_MSA(
        df=df,
        col="PM25_diff",
        cbsa=cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 Difference (2020 - 2014)',
        cmap="RdBu",
        vs=[-3, 3],
        counties=False,
        difference=True,
        ax=axes[i, 2],  # Use the third column for difference data
        title_size=24
    )
    # Add color bar
    cbar = fig.colorbar(sm, ax=axes[i, 2], orientation="vertical", shrink=0.8)
    cbar.ax.tick_params(labelsize=22)
    cbar.ax.set_ylabel("PM2.5 Concentration (µg/m³)", rotation=90, fontsize=22, labelpad = 7)
# Show the plot
plt.show()


# %% [markdown]
# ## City Plotting With Validation

# %%
import matplotlib.pyplot as plt
from pyproj import CRS, Transformer
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import box


def plot_MSA(df, MSA, title, validation_df, col="PM25", cmap="viridis", vs=False, counties=False, difference=False, state_color="white", ax=None):
    # Define the target projection
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
    ax.add_feature(ocean, linewidth=0.2, zorder=0)
    ax.add_feature(lakes, linewidth=0.2, zorder=0)
    ax.add_feature(states_provinces, edgecolor=state_color, linewidth=0.5, zorder=0)
    #ax.add_feature(bodr, linestyle='--', edgecolor=state_color, alpha=1)
    msa_bounds = df.loc[df["CBSA"] == MSA].total_bounds  # [minx, miny, maxx, maxy]
    
    #ax.set_extent([msa_bounds[0], msa_bounds[2], msa_bounds[1], msa_bounds[3]], crs=target_proj)

    buffer_size = 20000  # Example buffer size in the unit of the projection (e.g., meters)
    
    # Apply the buffer to expand the bounding box
    buffered_bounds = [
        msa_bounds[0] - buffer_size,  # minx
        msa_bounds[1] - buffer_size,  # miny
        msa_bounds[2] + buffer_size,  # maxx
        msa_bounds[3] + buffer_size   # maxy
    ]
    
    # Create a box polygon with the buffered bounds
    plot_area = box(*buffered_bounds)
    ax.set_extent([msa_bounds[0] -buffer_size , msa_bounds[2] + buffer_size, msa_bounds[1] - buffer_size, msa_bounds[3] + buffer_size], crs=target_proj)

    
    # Filter the main dataframe for data within the buffered plot area
    sub_df = df[df.geometry.within(plot_area)]
    
    if vs:
        vmin = vs[0]
        vmax = vs[1]
    elif not difference:
        vmin = np.quantile(sub_df[col], 0.025)
        vmax = np.quantile(sub_df[col], 0.975)
    else:
        vmin = -max( [abs(np.quantile(sub_df[col], 0.025)), abs(np.quantile(sub_df[col], 0.975))] )
        if vmin > -20:
            vmin = -20
        vmax = -vmin
    
    # # Plot MSA area data
    legend_bool = (validation_df.empty or difference == True)
    sc = sub_df.plot(ax=ax, column=col, cmap=cmap, legend=legend_bool, vmin=vmin, vmax=vmax, zorder=-1)
    
    msa_geometry = df.loc[df["CBSA"] == MSA, "geometry"].unary_union  # Get combined geometry for the MSA
    # Simplify the geometry to reduce jagged edges
    tolerance = 4000  # Adjust this value as needed for smoother lines (in coordinate units)
    simplified_msa_geometry = msa_geometry.simplify(tolerance, preserve_topology=True)
    
    # Plot the simplified geometry
    ax.add_geometries([simplified_msa_geometry], crs=target_proj, edgecolor=state_color, facecolor='none', linewidth=1.5)
    
    # Plot validation points
    if not validation_df.empty:
        # MSA bounds in Lambert Conformal coordinates (convert from kilometers to degrees)
        #bounds = sub_df.total_bounds
    
        # Define the projection used for the MSA bounds (Lambert Conformal)
        msa_proj = CRS(proj='lcc', lat_1=33, lat_2=45, lon_0=-97, lat_0=40, x_0=0, y_0=0)
    
        # Define the transformer for converting from Lambert Conformal to geographic (EPSG:4326)
        transformer = Transformer.from_crs(msa_proj, CRS.from_epsg(4326), always_xy=True)
    
        # Transform bounding box to geographic coordinates (longitude/latitude)
        lon_min, lat_min = transformer.transform(msa_bounds[0], msa_bounds[1])
        lon_max, lat_max = transformer.transform(msa_bounds[2], msa_bounds[3])
    
        # Subselect validation points that fall within the MSA's transformed bounding box
        validation_within_MSA = validation_df[
            (validation_df['Longitude'] >= lon_min) & (validation_df['Longitude'] <= lon_max) &
            (validation_df['Latitude'] >= lat_min) & (validation_df['Latitude'] <= lat_max)
        ]
    
        # Plot validation points as large dots
        sc = ax.scatter(validation_within_MSA['Longitude'], validation_within_MSA['Latitude'],
                        c=validation_within_MSA['μd'], cmap=cmap,
                        s=100, edgecolor='black', transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, zorder=3)
    
    
    # Add title and return the subset of the data
    ax.set_title(title)
    return sc

def plot_multiple_MSA(df, validation_df, MSA, vs1, vs2, vs3 = [0,5]):
    # Create a Matplotlib figure with 1 row and 3 columns
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(18, 14), subplot_kw={'projection': target_proj})

    place = MSA.split('-')[0]
    
    # Plot for PM2.5
    sc1 = plot_MSA(
        df=df,
        col="PM25",
        MSA=MSA,
        title=f'{place} PM2.5',
        validation_df=validation_df,
        cmap="magma",
        vs=[vs1[0], vs1[1]],
        counties=False,
        difference=False,
        ax=axes[0][0]  # Pass the first axis
    )

    # Plot for PM2.5 2020
    sc2 = plot_MSA(
        df=df,
        col="PM25_2020",
        MSA=MSA,
        title=f'{place} PM2.5 2020',
        validation_df=validation_df,
        cmap="magma",
        vs=[vs1[0], vs1[1]],
        counties=False,
        difference=False,
        ax=axes[0][1]  # Pass the second axis
    )

    # Plot for PM2.5 difference
    sc3 = plot_MSA(
        df=df,
        col="PM25_diff",
        MSA=MSA,
        title=f'{place} PM2.5 2020 - baseline',
        validation_df=pd.DataFrame(),  # No validation points for difference
        cmap="RdBu",
        vs=[vs2[0], vs2[1]],
        counties=False,
        difference=True,
        ax=axes[0][2],
        state_color="black"# Pass the third axis
    )

    #rwc_bounds = max

        # Plot for PM2.5 difference
    sc4 = plot_MSA(
        df=df,
        col="PM25_RWC",
        MSA=MSA,
        title=f'{place} PM2.5 RWC baseline',
        validation_df=pd.DataFrame(),  # No validation points for difference
        cmap="magma",
        vs=[vs3[0], vs3[1]],
        counties=False,
        difference=False,
        ax=axes[1][0],
    )

        # Plot for PM2.5 difference
    sc5 = plot_MSA(
        df=df,
        col="PM25_RWC_2020",
        MSA=MSA,
        title=f'{place} PM2.5 RWC 2020',
        validation_df=pd.DataFrame(),  # No validation points for difference
        cmap="magma",
        vs=[vs3[0], vs3[1]],
        counties=False,
        difference=False,
        ax=axes[1][1],
    )

    sc6 = plot_MSA(
        df=df,
        col="RWC_percent_diff",
        MSA=MSA,
        title=f'{place} PM2.5 RWC % Difference',
        validation_df=pd.DataFrame(),  # No validation points for difference
        cmap="RdBu",
        vs=False,
        counties=False,
        difference=True,
        ax=axes[1][2],
        state_color="black"# Pass the third axis
    )

    # Add colorbars for the first two subplots
    cbar1 = fig.colorbar(sc1, ax=axes[0][0], orientation='vertical')
    cbar1.set_label('Mean Concentration')

    cbar2 = fig.colorbar(sc2, ax=axes[0][1], orientation='vertical')
    cbar2.set_label('Mean Concentration')

    # cbar4 = fig.colorbar(sc4, ax=axes[0][0], orientation='vertical')
    # cbar4.set_label('Mean Concentration')

    # cbar5 = fig.colorbar(sc5, ax=axes[0][1], orientation='vertical')
    # cbar5.set_label('Mean Concentration')

    plt.show()

# %%

def plot_multiple_MSA_RWC(df, validation_df, MSA, vs3 = [0,5]):
    # Create a Matplotlib figure with 1 row and 3 columns
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 20), subplot_kw={'projection': target_proj})

    place = MSA.split('-')[0]
    
    #rwc_bounds = max

        # Plot for PM2.5 difference
    sc4 = plot_MSA(
        df=df,
        col="PM25_RWC",
        MSA=MSA,
        title=f'{place} PM2.5 RWC baseline',
        validation_df=pd.DataFrame(),  # No validation points for difference
        cmap="turbo",
        vs=[vs3[0], vs3[1]],
        counties=False,
        difference=False,
        ax=axes[0],
    )

        # Plot for PM2.5 difference
    sc5 = plot_MSA(
        df=df,
        col="PM25_RWC_2020",
        MSA=MSA,
        title=f'{place} PM2.5 RWC 2020',
        validation_df=pd.DataFrame(),  # No validation points for difference
        cmap="turbo",
        vs=[vs3[0], vs3[1]],
        counties=False,
        difference=False,
        ax=axes[1],
    )

    sc6 = plot_MSA(
        df=df,
        col="RWC_percent_diff",
        MSA=MSA,
        title=f'{place} PM2.5 RWC % Difference',
        validation_df=pd.DataFrame(),  # No validation points for difference
        cmap="RdBu",
        vs=False,
        counties=False,
        difference=True,
        ax=axes[2],
        state_color="black"# Pass the third axis
    )

    plt.show()

# %%


# %%
plot_multiple_MSA_RWC(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Los Angeles-Long Beach-Anaheim, CA',
    #vs2=[-5, 5]
)

# %%
plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Los Angeles-Long Beach-Anaheim, CA',
    vs1=[0, 25],
    vs2=[-5, 5]
)

# %%

plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Chicago-Naperville-Elgin, IL-IN-WI',
    vs1=[0, 25],
    vs2=[-5, 5]
)


# %%

plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Minneapolis-St. Paul-Bloomington, MN-WI',
    vs1=[0, 25],
    vs2=[-10, 10],
    vs3 = [0,10]
)

# %%
plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Seattle-Tacoma-Bellevue, WA',
    vs1=[0, 25],
    vs2=[-10, 10],
    vs3 = [0, 10]
)

# %%
plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='New York-Newark-Jersey City, NY-NJ-PA',
    vs1=[0, 25],
    vs2=[-5, 5],
    vs3 = [0,8]
)

# %%
plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Houston-The Woodlands-Sugar Land, TX',
    vs1=[0, 25],
    vs2=[-5, 5]
)

# %%
plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Phoenix-Mesa-Chandler, AZ',
    vs1=[0, 25],
    vs2=[-5, 5]
)

# %%
plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Philadelphia-Camden-Wilmington, PA-NJ-DE-MD',
    vs1=[0, 25],
    vs2=[-5, 5]
)

# %%

plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Washington-Arlington-Alexandria, DC-VA-MD-WV',
    vs1=[0, 25],
    vs2=[-5, 5]
)

# %%
plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Atlanta-Sandy Springs-Alpharetta, GA',
    vs1=[0, 25],
    vs2=[-5, 5]
)

# %%
plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Denver-Aurora-Lakewood, CO',
    vs1=[0, 25],
    vs2=[-5, 5]
)

# %%
plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Salt Lake City, UT',
    vs1=[0, 25],
    vs2=[-5, 5]
)

# %%

plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Kansas City, MO-KS',
    vs1=[0, 25],
    vs2=[-5, 5]
)

# %%

plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='St. Louis, MO-IL',
    vs1=[0, 25],
    vs2=[-5, 5]
)

# %%

plot_multiple_MSA(
    df=pmgdf,
    validation_df=validation_baseline,
    MSA='Dallas-Fort Worth-Arlington, TX',
    vs1=[0, 25],
    vs2=[-5, 5]
)

# %% [markdown]
# ## CONUS Plotting

# %%
plot(data = baseline['PM25_TOT'][0,0,:,:],
     vs = [0,20],
     label =  "PM2.5 Concentrations (u/m^3)", 
     title = "January 2016 PM2.5 Concentrations",
     lon = lon,
     lat = lat
    )

# %%
plot(data = rwc_2020['PM25_TOT'][0,0,:,:],
     vs = [0,20],
     label =  "PM2.5 Concentrations (u/m^3)", 
     title = "January 2016 PM2.5 Concentrations",
     lon = lon,
     lat = lat
    )

# %%
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


# %%
plot(data = baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:],
     vs = [0,8],
     label =  "RWC PM2.5 Concentrations (µg/m³)", 
     title = "NEI 2014 Average RWC contributed PM2.5 Concentrations",
     lon = lon,
     lat = lat,
     cbar = "turbo",
     state_color = "black"
    )

# %%
plot(data = rwc_2020['PM25_TOT'][0,0,:,:],
     vs = [0,20],
     label =  "PM₂.₅ Concentrations (µg/m³)", 
     title = "CONUS January Average PM₂.₅ Concentrations",
     lon = lon,
     lat = lat,
     cbar = "magma_r",
     state_color = "white"
    )

# %%


plot(data = rwc_2020['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:],
     vs = [0,5],
     label =  "RWC PM₂.₅ Concentrations (µg/m³)", 
     title = "Average RWC contributed PM₂.₅ Concentrations",
     lon = lon,
     lat = lat,
     cbar = "magma_r",
     state_color = "white"
    )

# %%
plot(data = rwc_2020['PM25_TOT'][0,0,:,:] - baseline['PM25_TOT'][0,0,:,:],
     vs = [-5,5],
     label =  "PM2.5 Concentration (µg/m³)", 
     title = "NEI 2020 - NEI 2014 RWC contributed PM2.5 Concentrations",
     lon = lon,
     lat = lat,
     difference = True,
     state_color = "black",
     cbar = 'RdBu'
    )

# %%
plot(data = (rwc_2020['PM25_TOT'][0,0,:,:] - baseline['PM25_TOT'][0,0,:,:])/(baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:]) * 100,
     vs = [-100,100],
     label =  "Percent Difference in RWC PM2.5 Concentrations", 
     title = "Percent Difference from 2020_RWC to Baseline January 2016 RWC Contributed PM2.5 Concentrations",
     lon = lon,
     lat = lat,
     difference = True,
    )

# %%
plot(data = (rwc_2020['PM25_TOT'][0,0,:,:] - no_w- baseline['PM25_TOT'][0,0,:,:])/baseline['PM25_TOT'][0,0,:,:] * 100,
     vs = [-50,50],
     label =  "Percent Difference in PM2.5 Concentrations", 
     title = "Percent Difference from 2020_RWC to Baseline January 2016 PM2.5 Concentrations",
     lon = lon,
     lat = lat,
     difference = True,
    )

# %%
rwc_2020['PM25_EC'][0,0,:,:].sum()

# %%
rwc_2020['PM25_OC'][0,0,:,:].sum()

# %%
rwc_2020['PM25_TOT'][0,0,:,:].sum()

# %%
baseline['PM25_TOT'][0,0,:,:].sum()

# %% [markdown]
# ## RWC City Plotting --> Comparisons to SA study Data

# %%
# plot MSA data
def plot_MSA(df, MSA, title, col = "PM25",cmap = "viridis", vs = False, counties = False, difference = False,  state_color = "white" ):
    # Define the target projection
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

    #create matplotlib figure
    fig, ax = plt.subplots(figsize=(12,7), subplot_kw={'projection': target_proj}) #subplot_kw={'projection': ccrs.PlateCarree(central_longitude=-97.0, globe=None)})

    # subselect df
    sub_df = df.loc[df["CBSA"] == MSA]
    
    if vs:
        vmin=vs[0]
        vmax=vs[1]
    else:
        vmin = np.quantile(sub_df[col], 0.025)
        vmax = np.quantile(sub_df[col], 0.975)

    #plot
    sub_df.plot(ax = ax, column = col, cmap=cmap, legend=True, vmin = vmin, vmax = vmax)


    #generate features
    resol = '10m'  # use data at this scale
    bodr = cfeature.NaturalEarthFeature(category='cultural',  edgecolor = 'black',name='admin_0_boundary_lines_land', scale=resol, facecolor='none', alpha=0.7)
    ocean = cfeature.NaturalEarthFeature('physical', 'ocean', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
    lakes = cfeature.NaturalEarthFeature('physical', 'lakes', scale=resol, edgecolor='black', facecolor=cfeature.COLORS['water'])
    #rivers = cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', scale=resol, edgecolor='black', facecolor='none')
    states_provinces = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m',edgecolor = state_color, facecolor='none')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])

    #add features
    if counties:
        ax.add_feature(COUNTIES, facecolor='none', edgecolor=state_color)
    #ax.add_feature(rivers, linewidth=0.5,edgecolor='blue')
    ax.add_feature(ocean, linewidth=0.2 )
    ax.add_feature(lakes, linewidth = 0)
    ax.add_feature(states_provinces, edgecolor=state_color, linewidth=0.5)
    ax.add_feature(bodr, linestyle='--', edgecolor=state_color, alpha=1)

    #plot colorbar legend
    #cbar = plt.colorbar(mm)
    #cbar.ax.set_ylabel(label,rotation = 90)

    ax.set_title(title)
    return sub_df

# %%

la_df = plot_MSA(
    df = pmgdf,
    MSA = 'Los Angeles-Long Beach-Anaheim, CA',
    title = 'LA PM2.5 from RWC',
    cmap = "magma",
    vs = [0,5],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA = 'Los Angeles-Long Beach-Anaheim, CA',
    title = 'LA PM2.5 % from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [0,45],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA = 'Los Angeles-Long Beach-Anaheim, CA',
    title = 'LA PM2.5 % from RWC',
    col = "PM25_percent_2020",
    cmap = "magma",
    vs = [0,45],
    counties = True,
    difference = False
)

la_df['PM25_percent'].mean()

# %%
# San Jose, USA 42% heatingseason PM10 (SA) – – Chow et al. (1995)
# https://trid.trb.org/View/425534

#very much underestimated here ~20%, but less understimated in the 2020 version :)

sf_df = plot_MSA(
    df = pmgdf,
    MSA = 'San Jose-Sunnyvale-Santa Clara, CA',
    title = 'San Jose PM2.5 from RWC',
    cmap = "magma",
    vs = [0,5],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA = 'San Jose-Sunnyvale-Santa Clara, CA',
    title = 'San Jose PM2.5 % from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [0,45],
    counties = True,
    difference = False
)

sf_df['PM25_percent'].mean()

# %%
print(np.quantile(sf_df['PM25'], 0.9))
print(np.quantile(sf_df['PM25'], 1))
print(np.quantile(sf_df['PM25_percent'], 0.9))
print(np.quantile(sf_df['PM25_percent'], 1))

# %%
chi_df = plot_MSA(
    df = pmgdf,
    MSA =  'Chicago-Naperville-Elgin, IL-IN-WI',
    title = 'Chicago PM2.5 from RWC',
    cmap = "magma",
    vs = [0,5],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA =  'Chicago-Naperville-Elgin, IL-IN-WI',
    title = 'Chicago PM2.5 % from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [0,45],
    counties = True,
    difference = False
)

chi_df['PM25_percent'].mean()

# %%
pmgdf['CBSA'].unique()

# %%
print(np.quantile(chi_df['PM25'], 0.8))
print(np.quantile(chi_df['PM25'], 1))
print(np.quantile(chi_df['PM25_percent'], 0.8))
print(np.quantile(chi_df['PM25_percent'], 1))

# %%
# Atlanta, USA: 11% winter PM2.5 ––Polissar et al. (2001), Sarnat et al., 2008
#  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2290994/
atlanta_df = plot_MSA(
    df = pmgdf,
    MSA = 'Atlanta-Sandy Springs-Alpharetta, GA',
    title = 'Atlanta PM2.5 from RWC',
    cmap = "magma",
    vs = [0,3],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA = 'Atlanta-Sandy Springs-Alpharetta, GA',
    title = 'Atlanta PM2.5 from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [10,30],
    counties = True,
    difference = False
)

atlanta_df['PM25_percent'].mean()

# %%
print(np.quantile(atlanta_df['PM25'], 0.9))
print(np.quantile(atlanta_df['PM25'], 1))
print(np.quantile(atlanta_df['PM25_percent'], 0.9))
print(np.quantile(atlanta_df['PM25_percent'], 1))

# %%


burlington_df = plot_MSA(
    df = pmgdf,
    MSA = 'Burlington-South Burlington, VT',
    title = 'Burlington PM2.5 from RWC',
    cmap = "magma",
    vs = [0,10],
    counties = True,
    difference = False
)

burlington_df['PM25_percent'].mean()

# %%
# Montana (5 communities), USA 55–77% heating-season PM2.5 (SA) 7–11 – Ward et al. (2010)
# measured in Libby Montana, the northwest corner city (the tiny spot)
# https://www.tandfonline.com/doi/pdf/10.3155/1047-3289.60.6.688

# Substantially increased Montana emissions in 2020 version would improve on this

plot(data = baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:],
     vs = [0,3],
     label =  "RWC PM2.5 Contribution (u/m^3)", 
     title = "Montana January 2016 RWC contributed PM2.5 Concentrations",
     lon = lon,
     lat = lat,
     cbar = "magma",
     counties = False,
     bounds = [ [660,666], [227,233] ]
    )

plot(data = percent_RWC,
     vs = [0,50],
     label =  "RWC PM2.5 Contribution (%)", 
     title = "Montana % January 2016 RWC contributed PM2.5 Concentrations",
     lon = lon,
     lat = lat,
     cbar = "magma",
     counties = False,
     bounds = [ [660,666], [227,233] ]
    )
print(f"Bounds Max: {int(percent_RWC[662:665,229:231].max())}")
print(f"Bounds Max: {int(percent_RWC[662:665,229:231].mean())}")

# %%
# Seattle: 11 ug/m3, ~30% in heating season Allen et al. (2008), Kim & Hopke 2008, Wu et al. (2007)
# This actually matches up very well

seattle_df = plot_MSA(
    df = pmgdf,
    MSA = 'Seattle-Tacoma-Bellevue, WA',
    title = 'Seattle PM2.5 from RWC',
    cmap = "magma",
    vs = [0,15],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA = 'Seattle-Tacoma-Bellevue, WA',
    title = 'Seattle PM2.5 % from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [0,80],
    counties = True,
    difference = False
)

seattle_df['PM25_percent'].mean()

# %%
print(np.quantile(seattle_df['PM25_percent'], 0.7))
print(np.quantile(seattle_df['PM25_percent'], 1))

# %%
print(np.quantile(seattle_df['PM25'], 0.85))
print(np.quantile(seattle_df['PM25'], 1))

# %%
# Portlan: ~27% in heating season Kim & Hopke 2008
# This  matches up pretty well
portland_df = plot_MSA(
    df = pmgdf,
    MSA = 'Portland-Vancouver-Hillsboro, OR-WA',
    title = 'Portland PM2.5 from RWC',
    cmap = "magma",
    vs = [0,15],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA = 'Portland-Vancouver-Hillsboro, OR-WA',
    title = 'Portland PM2.5 % from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [0,80],
    counties = True,
    difference = False
)

portland_df['PM25_percent'].mean()

# %%
print(np.quantile(portland_df['PM25'], 0.95))
print(np.quantile(portland_df['PM25'], 1))
print(np.quantile(portland_df['PM25_percent'], 0.75))
print(np.quantile(portland_df['PM25_percent'], 0.95))

# %%
#Las Vegas 11–21% annual PM2.5 (SA) – Individual sites: 11.3 ±9.8%, 15.9 ±12.9%, 11.1 ±8.0, 20.8 ±12.5% contributions to annual PM2.5
# OC: 8–16% contribution from residential wood combustion; EC: 3–7% contribution from residential wood combustion Green et al. (2013)"""

LV_df = plot_MSA(
    df = pmgdf,
    MSA = 'Las Vegas-Henderson-Paradise, NV',
    title = 'Las Vegas PM2.5 from RWC',
    cmap = "magma",
    vs = [0,5],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA = 'Las Vegas-Henderson-Paradise, NV',
    title = 'Las Vegas PM2.5 % from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [0,30],
    counties = True,
    difference = False
)

LV_df['PM25_percent'].mean()

# %%
print(np.quantile(LV_df['PM25_percent'], 0.2))
print(np.quantile(LV_df['PM25_percent'], 0.9))

# %%
print(np.quantile(LV_df['PM25'], 0.96))
print(np.quantile(LV_df['PM25'], 1))

# %%
# San Joaquin Valley 24% winter PM2.5 https://pubmed.ncbi.nlm.nih.gov/17533844/
# 2014 RWC underesimates this, but it would do a little better in the 2020 version

plot(data = baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:],
     vs = [0,4],
     label =  "RWC PM2.5 Concentrations (u/m^3)", 
     title = "San Joaquin Valley January 2016 RWC contributed PM2.5 Concentrations",
     lon = lon,
     lat = lat,
     cbar = "magma",
     counties = True,
     bounds = [ [300,440], [30,100] ]
    )

plot(data = 100*(baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:])/baseline['PM25_TOT'][0,0,:,:],
     vs = [0,25],
     label =  "RWC PM2.5 Contribution (%)", 
     title = "San Joaquin Valley % January 2016 RWC contributed PM2.5 Concentrations",
     lon = lon,
     lat = lat,
     cbar = "magma",
     counties = True,
     bounds = [ [300,440], [30,100] ]
    )

# print(f"Bounds Mean: {int(percent_RWC[300:340,150:190].mean())}")

# %%
# Truckee, USA 11–15% winter PM2.5 (SA) – – Chen et al. (2012) 
# Truckee has nothing on the map from RWC for 2014 but much more for 2020 RWC.

truckee_df = plot_MSA(
    df = pmgdf,
    MSA = 'Truckee-Grass Valley, CA',
    title = 'Truckee PM2.5 from RWC',
    cmap = "magma",
    vs = [0,2],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA = 'Truckee-Grass Valley, CA',
    title = 'Truckee PM2.5 % from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [0,50],
    counties = True,
    difference = False
)


# %%
print(np.quantile(truckee_df['PM25_percent'], 0.2))
print(np.quantile(truckee_df['PM25_percent'], 0.95))

# %%
# Rochester, New York, USA 17% winter PM2.5 (SA) 3.2 ug/m3 
# Wood smoke contribution to PM2.5 increased to 27% when the corresponding hourly PM2.5 concentrations were greater than 15 μg/m3 
# Wang et al. (2011)

# measures about 3.8 ug/m3 but more like 20-25% winter pm25
rochester_df = plot_MSA(
    df = pmgdf,
    MSA = 'Rochester, NY',
    title = 'Rochester PM2.5 from RWC',
    cmap = "magma",
    vs = [0,5],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA = 'Rochester, NY',
    title = 'Rochester PM2.5 % from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [0,50],
    counties = True,
    difference = False
)

print(rochester_df["PM25"].mean())

# %%
print(np.quantile(rochester_df['PM25_percent'], 0.2))
print(np.quantile(rochester_df['PM25_percent'], 0.9))

print('\n',np.quantile(rochester_df['PM25'], 0.75))
print(np.quantile(rochester_df['PM25'], 0.975))
print(np.quantile(rochester_df['PM25'], 1))

# %%
pmgdf['CBSA'].unique()

# %%

mn_df = plot_MSA(
    df = pmgdf,
    MSA = 'Minneapolis-St. Paul-Bloomington, MN-WI',
    title = 'Minneapolis PM2.5 from RWC',
    cmap = "magma",
    vs = [0,12],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA = 'Minneapolis-St. Paul-Bloomington, MN-WI',
    title = 'Minneapolis PM2.5 % from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [0,80],
    counties = True,
    difference = False
)

mn_df['PM25_percent'].mean()

# %%
print(np.quantile(mn_df['PM25'], 0.975))
print(np.quantile(mn_df['PM25'], 1))

# %%
phil_df = plot_MSA(
    df = pmgdf,
    MSA =  'Philadelphia-Camden-Wilmington, PA-NJ-DE-MD',
    title = 'Philadelphia PM2.5 from RWC',
    cmap = "magma",
    vs = [0,7],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA =  'Philadelphia-Camden-Wilmington, PA-NJ-DE-MD',
    title = 'Philadelphia PM2.5 % from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [0,50],
    counties = True,
    difference = False
)

phil_df['PM25_percent'].mean()

# %%

nyc_df = plot_MSA(
    df = pmgdf,
    MSA =  'New York-Newark-Jersey City, NY-NJ-PA',
    title = 'NYC PM2.5 from RWC',
    cmap = "magma",
    vs = [0,7],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA =  'New York-Newark-Jersey City, NY-NJ-PA',
    title = 'NYC PM2.5 % from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [0,50],
    counties = True,
    difference = False
)

nyc_df['PM25_percent'].mean()

# %%
print(np.quantile(mn_df['PM25'], 0.7))
print(np.quantile(mn_df['PM25'], 0.96))

# %%
print(np.quantile(mn_df['PM25_percent'], 0.01))
print(np.quantile(mn_df['PM25_percent'], 0.9))

# %%

phoenix_df = plot_MSA(
    df = pmgdf,
    MSA =  'Phoenix-Mesa-Chandler, AZ',
    title = 'Phoenix PM2.5 from RWC',
    cmap = "magma",
    vs = [0,7],
    counties = True,
    difference = False
)

plot_MSA(
    df = pmgdf,
    MSA =  'Phoenix-Mesa-Chandler, AZ',
    title = 'Phoenix PM2.5 % from RWC',
    col = "PM25_percent",
    cmap = "magma",
    vs = [0,50],
    counties = True,
    difference = False
)

phoenix_df['PM25_percent'].mean()

# %%
print(np.quantile(mn_df['PM25'], 0.7))
print(np.quantile(mn_df['PM25'], 0.96))

# %% [markdown]
# # Population weighting

# %%
# Read the .txt file into a pandas DataFrame
population_df = pd.read_csv("../SMOKE_sensitivity_analyses/population.txt", delimiter='\t', header=None, skiprows=25)

# Specify coordinates to align population and emission data
population_df = population_df.rename(columns = {2:"COLS", 3:"ROWS"})
population_df["COLS"] -= 111 
population_df["ROWS"] -= 126   

# %%
len(population_df['COLS'].unique())

# %%
len(pmgdf['col'].unique())

# %%
pm_df = baseline['PM25_TOT'][0,0,:,:].to_dataframe()
pm_df['PM25_2020'] = (rwc_2020['PM25_TOT'][0,0,:,:]).to_dataframe()
pm_df['PM25_diff'] = (rwc_2020['PM25_TOT'][0,0,:,:] - baseline['PM25_TOT'][0,0,:,:]).to_dataframe()

pm_df['PM25_RWC'] = (baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:]).to_dataframe()
pm_df['PM25_RWC_2020'] = (rwc_2020['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:]).to_dataframe()

pm_df['PM25_RWC_percent'] =(100* (baseline['PM25_TOT'][0,0,:,:] - no_rwc['PM25_TOT'][0,0,:,:])/baseline['PM25_TOT'][0,0,:,:]).to_dataframe()
pm_df['PM25_RWC_percent_2020'] = (percent_RWC_2020).to_dataframe()

pm_df['RWC_percent_diff'] = (pm_df['PM25_RWC_2020'] - pm_df['PM25_RWC'])/ pm_df['PM25_TOT'] * 100 
pm_df = pm_df.reset_index()

# %%
PWC_df = pd.merge(pm_df, population_df, left_on=['ROW', 'COL'], right_on=['ROWS', 'COLS'])
PWC_df['population'] = PWC_df[6]

PWC_df['PWC_TOT'] = PWC_df['PM25_TOT'] * PWC_df['population']
PWC_df['PWC_RWC'] = PWC_df['PM25_RWC'] * PWC_df['population']
PWC_df['PWC_TOT_2020'] = PWC_df['PM25_2020'] * PWC_df['population']
PWC_df['PWC_RWC_2020'] = PWC_df['PM25_RWC_2020'] * PWC_df['population']


total_PWC = PWC_df['PWC_TOT'].sum()
RWC_PWC = PWC_df['PWC_RWC'].sum()

total_PWC_2020 = PWC_df['PWC_TOT_2020'].sum()
RWC_PWC_2020 = PWC_df['PWC_RWC_2020'].sum()

# %%
print(f"Percent of PWC from RWC: {round(RWC_PWC/total_PWC * 100,3)}%")

# %%
print(f"population weighted mean concentration {round(RWC_PWC/PWC_df['population'].sum(), 2)} ug/m^3") 
#population weighted mean concentration

# %%
print(f"population weighted mean concentration 2020 {round(RWC_PWC_2020/PWC_df['population'].sum(), 2)} ug/m^3") 

# %%
print(f"Percent of PWC from RWC 2020: {round(PWC_df['PWC_RWC_2020'].sum()/total_PWC_2020 * 100,3)}%")

# %%
print(f"mean concentration {str(round(PWC_df['PM25_RWC'].mean(), 2))} ug/m^3") 

# %%
print(f"mean concentration {str(round(PWC_df['PM25_RWC_2020'].mean(), 2))} ug/m^3") 

# %%
100 * np.mean(PWC_df['PM25_RWC'] / PWC_df['PM25_TOT'])

# %%
100 * np.mean(PWC_df['PM25_RWC_2020'] / PWC_df['PM25_2020'])

# %%
PWC_df.shape

# %%
pm_df.loc[pm_df['PM25_TOT'] == pm_df['PM25_TOT'].max()] #southwestern Kansas

# %%
pm_df.loc[pm_df['PM25_2020'] == pm_df['PM25_2020'].max()]

# %%
from pyproj import Proj, Transformer

# Grid indices (row 328, column 466)
row, col = 328, 466
xorig = -2292000
yorig = -1584000

# Calculate the center of the grid cell in the projection coordinates
x_center = xorig + (col) * cell_size
y_center = yorig + (row) * cell_size

# Set up a transformer to convert from LCC to latitude/longitude (WGS84)
transformer = Transformer.from_crs(proj_params, "EPSG:4326", always_xy=True)

# Transform to latitude and longitude
lon, lat = transformer.transform(x_center, y_center)

# Print the results
print(f"Latitude: {lat}, Longitude: {lon}")

# %%
pm_df.loc[pm_df['PM25_TOT'] == pm_df['PM25_TOT'].min()] #southwestern Kansas

# %%
pm_df.loc[pm_df['PM25_2020'] == pm_df['PM25_2020'].min()]

# %%
from pyproj import Proj, Transformer

# Grid indices (row 328, column 466)
row, col = 650, 106

# Calculate the center of the grid cell in the projection coordinates
x_center = xorig + (col) * cell_size
y_center = yorig + (row) * cell_size

# Set up a transformer to convert from LCC to latitude/longitude (WGS84)
transformer = Transformer.from_crs(proj_params, "EPSG:4326", always_xy=True)

# Transform to latitude and longitude
lon, lat = transformer.transform(x_center, y_center)

# Print the results
print(f"Latitude: {lat}, Longitude: {lon}")

# %%
pm_df.loc[pm_df['PM25_RWC'] == pm_df['PM25_RWC'].max()] #southwestern Kansas

# %%
x_center = xorig + (610) * cell_size
y_center = yorig + (503) * cell_size
lon, lat = transformer.transform(x_center, y_center)
print(f"Latitude: {lat}, Longitude: {lon}")

# %%
pm_df.loc[pm_df['PM25_RWC_2020'] == pm_df['PM25_RWC_2020'].max()]

# %%
x_center = xorig + (182) * cell_size
y_center = yorig + (537) * cell_size
lon, lat = transformer.transform(x_center, y_center)
print(f"Latitude: {lat}, Longitude: {lon}")

# %% [markdown]
# ## East vs. West Analysis

# %%
def grid_to_latlon(row, col):
    """
    Convert grid row and column indices to latitude and longitude.
    
    Parameters:
        row (int): Row index in the grid.
        col (int): Column index in the grid.
        
    Returns:
        tuple: (lat, lon) coordinates in WGS84.
    """
    # Projection parameters (same as in latlon_to_grid)
    proj_params = {
        'proj': 'lcc',
        'lat_1': 33,
        'lat_2': 45,
        'lon_0': -97,
        'lat_0': 40
    }
    transformer = Transformer.from_crs(proj_params, "EPSG:4326", always_xy=True)
    cell_size = 4000  # Cell size in meters
    xorig = -2292000
    yorig = -1584000

    # Compute projected coordinates from row and column indices
    x_proj = xorig + col * cell_size
    y_proj = yorig + row * cell_size
    
    # Transform projected coordinates back to latitude and longitude
    lon, lat = transformer.transform(x_proj, y_proj)
    
    return lat, lon


# %%
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

# %%
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

# %%
PWC_df_emissions_east.loc[PWC_df_emissions_east['PM25_total_tons']>0.02].shape[0]/PWC_df_emissions_east.shape[0] * 100

# %%
PWC_df_emissions_east[6].sum()

# %%
np.corrcoef(PWC_df_emissions_east['PM25_total_tons'], PWC_df_emissions_east[6])

# %%
np.corrcoef(PWC_df_emissions_west['PM25_total_tons'], PWC_df_emissions_west[6])

# %%
PWC_df_emissions_west.loc[PWC_df_emissions_west[6] >= PWC_df_emissions_west[6].quantile(0.99)]['PM25_total_tons'].mean()

# %%
PWC_df_emissions_east.loc[PWC_df_emissions_east[6] >= PWC_df_emissions_east[6].quantile(0.99)]['PM25_total_tons'].mean()

# %%
PWC_df_emissions_east['PM25_total_tons'].mean()

# %%
PWC_df_emissions_west['PM25_total_tons'].mean()

# %%
PWC_df_emissions_west.loc[PWC_df_emissions_west['PM25_total_tons']>0.02].shape[0]/PWC_df_emissions_west.shape[0] * 100

# %%
PWC_df_emissions_east['PM25_total_tons'].quantile(0.99)

# %%
PWC_df_emissions_west['PM25_total_tons'].quantile(0.99)

# %%
PWC_df_emissions_west.loc[PWC_df_emissions_west['pop']]

# %%
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

# Separate the DataFrame
PWC_west = PWC_df[PWC_df['west_of_100']]
PWC_east = PWC_df[~PWC_df['west_of_100']]


# %%
PWC_east = PWC_east.merge(
    PWC_df_emissions_east[['row', 'col', 'PM25_total_tons']],
    on=['row', 'col'],
    how='left'
)

# %%
PWC_west = PWC_west.merge(
    PWC_df_emissions_west[['row', 'col', 'PM25_total_tons']],
    on=['row', 'col'],
    how='left'
)

# %%
PWC_west[PWC_west['population_density'] > 193.134]['PM25_total_tons_y'].mean()

# %%
PWC_east[PWC_east['population_density'] > 193.134]['PM25_total_tons_y'].mean()

# %%
(PWC_east['PM25_RWC_2020'] * PWC_east['population']).sum()/PWC_east['population'].sum()

# %%
(PWC_west['PM25_RWC_2020'] * PWC_west['population']).sum()/PWC_west['population'].sum()

# %%
PWC_east['PM25_RWC_2020'].mean()

# %%
PWC_west['PM25_RWC_2020'].mean()

# %%
PWC_east.loc[PWC_east['PM25_RWC_2020'] >=1 ].shape[0]/PWC_east.shape[0] * 100

# %%
PWC_west.loc[PWC_west['PM25_RWC_2020'] >=1 ].shape[0]/PWC_west.shape[0] * 100

# %%
PWC_east.shape

# %%
PWC_west['PM25_RWC_2020'].quantile(0.999)

# %%
PWC_east['PM25_RWC_2020'].quantile(0.999)

# %% [markdown]
# ## Urban vs. Rural Cells

# %%
PWC_df['population'] = PWC_df[6]

PWC_df['population_density'] = (PWC_df['population'] / PWC_df.geometry.area) * 1_000_000 * 1.609 **2# per km  

# %%
urban = PWC_df.loc[PWC_df['population_density'] >= 500]
rural = PWC_df.loc[PWC_df['population_density'] < 500]

# %%
PWC_df['PM25_RWC_2020'].median()

# %%
urban['PM25_RWC_2020'].sum() / urban['PM25_2020'].sum() * 100

# %%
rural['PM25_RWC_2020'].sum() / rural['PM25_2020'].sum() * 100

# %%
urban.shape[0]

# %%
rural.shape[0]

# %%
urban['PM25_RWC_2020'].median()

# %%
rural['PM25_RWC_2020'].median()

# %% [markdown]
# ## Checking SA Studies

# %%
from pyproj import Transformer

def latlon_to_grid(lat, lon):
    """
    Convert latitude and longitude to grid row and column indices.
    
    Parameters:
        lat (float): Latitude in WGS84.
        lon (float): Longitude in WGS84.
        cell_size (float): Grid cell size in projection units.
        proj_params (str or dict): Projection parameters for the grid.
        xorig (float): X-coordinate origin of the grid.
        yorig (float): Y-coordinate origin of the grid.
        
    Returns:
        tuple: (row, col) indices corresponding to the latitude and longitude.
    """
    # Set up a transformer to convert from latitude/longitude (WGS84) to the grid projection
    proj_params = {'proj': 'lcc',
               'lat_1': 33,
               'lat_2': 45,
               'lon_0': -97,
               'lat_0': 40}
    transformer = Transformer.from_crs("EPSG:4326", proj_params, always_xy=True)
    cell_size = 4000  # Example cell size in meters
    xorig = -2292000
    yorig = -1584000

    # Transform latitude and longitude to projected coordinates
    x_proj, y_proj = transformer.transform(lon, lat)
    
    # Compute the column and row indices
    col = round((x_proj - xorig) / cell_size)
    row = round((y_proj - yorig) / cell_size)
    
    return row, col

#Phillips neighborhood
lat, lon = 44.9543561,-93.2691361  # Replace with actual coordinates
row, col = latlon_to_grid(lat, lon)
print(f"Row: {row}, Col: {col}")

pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]

# %%
## NYC 
lat, lon = 40.73708, -73.82158
row, col = latlon_to_grid(lat, lon)

# Select the 3×3 neighborhood
neighbors = pmgdf.loc[
    (pmgdf['row'].between(row - 1, row + 1)) & 
    (pmgdf['col'].between(col - 1, col + 1))
]

# Compute the average PM2.5
print(round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_2020']),1), \
    round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_percent_2020']),1) )

print(round(neighbors['PM25_RWC_2020'].min(), 1), 
      round(neighbors['PM25_RWC_2020'].max(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].min(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].max(), 1))


# %%
lat, lon  = 33.8128, -112.2381
#phoenx - neighborhood
row, col = latlon_to_grid(lat, lon)

# Select the 3×3 neighborhood
neighbors = pmgdf.loc[
    (pmgdf['row'].between(row - 1, row + 1)) & 
    (pmgdf['col'].between(col - 1, col + 1))
]

# Compute the average PM2.5
print(round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_2020']),1), \
    round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_percent_2020']),1) )

print(round(neighbors['PM25_RWC_2020'].min(), 1), 
      round(neighbors['PM25_RWC_2020'].max(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].min(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].max(), 1))


# %%
lat, lon  = 33.7750752, -84.417395
#jefferson st. atlanta - neighborhood
row, col = latlon_to_grid(lat, lon)
print(f"Row: {row}, Col: {col}")

pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]

# %%
# north birmingham
lat, lon = 33.5516709,-86.8399714
row, col = latlon_to_grid(lat, lon)
print(f"Row: {row}, Col: {col}")

pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]

# %%
# centreville
lat, lon = 32.9649163,-87.1426895

row, col = latlon_to_grid(lat, lon)
print(f"Row: {row}, Col: {col}")

pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]

# %%
#ccheaka peak 

lat, lon = 48.3000793,-124.636648


row, col = latlon_to_grid(lat, lon)
print(f"Row: {row}, Col: {col}")

pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]


# %%
# seattle
lat, lon = 47.5857665,-122.3331929
row, col = latlon_to_grid(lat, lon)
pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]


# %%
# brownseville
lat, lon = 26.06973, -97.16222

row, col = latlon_to_grid(lat, lon)

# Select the 3×3 neighborhood
neighbors = pmgdf.loc[
    (pmgdf['row'].between(row - 1, row + 1)) & 
    (pmgdf['col'].between(col - 1, col + 1))
]

# Compute the average PM2.5
print(round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_2020']),1), \
    round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_percent_2020']),1) )

print(round(neighbors['PM25_RWC_2020'].min(), 1), 
      round(neighbors['PM25_RWC_2020'].max(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].min(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].max(), 1))

# %%
lat, lon = 34.067186, -118.227172


row, col = latlon_to_grid(lat, lon)
pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]

# %%
import numpy as np

#Los angeles

# Define target location
lat, lon = 34.067186, -118.227172

# Convert to grid coordinates
row, col = latlon_to_grid(lat, lon)

# Select the 3×3 neighborhood
neighbors = pmgdf.loc[
    (pmgdf['row'].between(row - 1, row + 1)) & 
    (pmgdf['col'].between(col - 1, col + 1))
]

# Compute the average PM2.5
print(round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_2020']),1), \
    round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_percent_2020']),1) )

print(round(neighbors['PM25_RWC_2020'].min(), 1), 
      round(neighbors['PM25_RWC_2020'].max(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].min(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].max(), 1))

# %%
# roubidoux
lat, lon = 34.002296, -117.417796

# Convert to grid coordinates
row, col = latlon_to_grid(lat, lon)

# Select the 3×3 neighborhood
neighbors = pmgdf.loc[
    (pmgdf['row'].between(row - 1, row + 1)) & 
    (pmgdf['col'].between(col - 1, col + 1))
]

# Compute the average PM2.5
print(round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_2020']),1), \
    round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_percent_2020']),1) )

print(round(neighbors['PM25_RWC_2020'].min(), 1), 
      round(neighbors['PM25_RWC_2020'].max(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].min(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].max(), 1))


# %%
# Detroit
lat, lon = 42.228611, -83.20833

# Convert to grid coordinates
row, col = latlon_to_grid(lat, lon)

# Select the 3×3 neighborhood
neighbors = pmgdf.loc[
    (pmgdf['row'].between(row - 1, row + 1)) & 
    (pmgdf['col'].between(col - 1, col + 1))
]

# Compute the average PM2.5
print(round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_2020']),1), \
    round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_percent_2020']),1) )

print(round(neighbors['PM25_RWC_2020'].min(), 1), 
      round(neighbors['PM25_RWC_2020'].max(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].min(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].max(), 1))

# %%
# Chicago
lat, lon = 41.7514, -87.713488

# Convert to grid coordinates
row, col = latlon_to_grid(lat, lon)

# Select the 3×3 neighborhood
neighbors = pmgdf.loc[
    (pmgdf['row'].between(row - 1, row + 1)) & 
    (pmgdf['col'].between(col - 1, col + 1))
]

# Compute the average PM2.5
print(round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_2020']),1), \
    round(float(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]['PM25_RWC_percent_2020']),1) )

print(round(neighbors['PM25_RWC_2020'].min(), 1), 
      round(neighbors['PM25_RWC_2020'].max(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].min(), 1), 
      round(neighbors['PM25_RWC_percent_2020'].max(), 1))

# %% [markdown]
# ### New York Study

# %%
pmgdf['PM25_RWC_percent_2020']

# %%
dictionary = {
	'Rochester':[43.1566, -77.6088],
	'Buffalo':[42.8864, -78.8784],
	'Bronx':[40.8448, -73.8648],
	'Manhattan':[40.7209, -74.0007],
	'Queens':[40.737, -73.8203], 
	'Albany':[42.6526, -73.7562],
	'Whiteface':[44.2814426,-74.1392044],
	'Pinnacle':[42.0978901,-77.2096275],
}

for key, value in dictionary.items():
	lat, lon = value[0], value[1]

	row, col = latlon_to_grid(lat, lon)

	loc_df = pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]
	print(key, float(loc_df['PM25_RWC_2020']), float(loc_df['PM25_RWC_percent_2020']))

# %%
import pandas as pd

# Given dictionary of locations with latitudes and longitudes
locations = {
    "Chattanooga": [35.0456, -85.3097],
    "Murray": [34.7928, -84.7389],
    "Rome": [34.257, -85.1647],
    "Yorkville (YRK)": [33.924, -85.038],
    "Jefferson St. (JST), Atlanta": [33.7451237, -84.4197378],
    "South DeKalb (SDK)": [33.698, -84.234],
    "Athens": [33.9519, -83.3576],
    "Augusta": [33.4735, -81.9751],
    "Macon": [32.8407, -83.6324],
    "Columbus": [32.4600, -84.9877],
    "Douglas": [31.508, -82.8507],
    "Charlton": [30.7915, -82.0029],
    "Tallahassee": [30.4383, -84.2807]
}

results = []

for key, (lat, lon) in locations.items():
    row, col = latlon_to_grid(lat, lon)
    
    # Filter the DataFrame
    loc_df = pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]
    
    if not loc_df.empty:
        pm25_rwc_2020 = float(loc_df['PM25_RWC_2020'].values[0])
        pm25_rwc_percent_2020 = float(loc_df['PM25_RWC_percent_2020'].values[0])
        
        # Append to results list
        results.append([key,pm25_rwc_2020, pm25_rwc_percent_2020])

# Create DataFrame from results
df_results = pd.DataFrame(results, columns=["Location","PM25_RWC_2020", "PM25_RWC_percent_2020"])

# Display DataFrame
df_results.sort_values(by = 'PM25_RWC_2020', ascending = False)


# %%
locations = {
   # "Fairbanks": [64.8407, -147.7225],
    "Bakersfield": [35.3561, -119.0412],
    "Klamath Falls": [42.1889, -121.7225],
    "Lakeview": [42.1889, -120.3519],
    "Oakridge": [43.7444, -122.4805],
    "Portland": [45.4965, -122.6034],
    "Bountiful": [40.9030, -111.8845],
    "Lindon": [40.3414, -111.7136],
    "Seattle (Duwamish)": [47.5632, -122.3405],
    "Tacoma (Alexander Ave.)": [47.2656, -122.3858]
}

results = []

for key, (lat, lon) in locations.items():
    row, col = latlon_to_grid(lat, lon)
    
    # Filter the DataFrame
    loc_df = pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col)]
    
    if not loc_df.empty:
        pm25_rwc_2020 = round(float(loc_df['PM25_RWC_2020'].values[0]),1)
        pm25_rwc_percent_2020 = round(float(loc_df['PM25_RWC_percent_2020'].values[0]), 1)
        # Append to results list
        results.append([key,pm25_rwc_2020, pm25_rwc_percent_2020])

# Create DataFrame from results
df_results = pd.DataFrame(results, columns=["Location","PM25_RWC_2020", "PM25_RWC_percent_2020"])

# Display DataFrame
df_results#[['PM25_RWC_2020', 'PM25_RWC_percent_2020']].set_index("PM25_RWC_2020")
        #.sort_values(by = 'PM25_RWC_2020', ascending = False)

# %%
import pandas as pd

locations = {
    "Bakersfield": [35.3561, -119.0412],
    "Klamath Falls": [42.1889, -121.7225],
    "Lakeview": [42.1889, -120.3519],
    "Oakridge": [43.7444, -122.4805],
    "Portland": [45.4965, -122.6034],
    "Bountiful": [40.9030, -111.8845],
    "Lindon": [40.3414, -111.7136],
    "Seattle (Duwamish)": [47.5632, -122.3405],
    "Tacoma (Alexander Ave.)": [47.2656, -122.3858]
}

results = []

for key, (lat, lon) in locations.items():
    row, col = latlon_to_grid(lat, lon)
    
    # Select the 3×3 neighborhood
    neighbors = pmgdf.loc[
        (pmgdf['row'].between(row - 1, row + 1)) & 
        (pmgdf['col'].between(col - 1, col + 1))
    ]
    
    if not neighbors.empty:
        pm25_rwc_2020_min = round(neighbors['PM25_RWC_2020'].min(), 1)
        pm25_rwc_2020_max = round(neighbors['PM25_RWC_2020'].max(), 1)
        
        pm25_rwc_percent_2020_min = round(neighbors['PM25_RWC_percent_2020'].min(), 1)
        pm25_rwc_percent_2020_max = round(neighbors['PM25_RWC_percent_2020'].max(), 1)
        
        results.append([
            key, 
            pm25_rwc_2020_min, pm25_rwc_2020_max, 
            pm25_rwc_percent_2020_min, pm25_rwc_percent_2020_max, 
        ])

# Create DataFrame from results
df_results = pd.DataFrame(results, columns=[
    "Location", 
    "PM25_RWC_2020_Min", "PM25_RWC_2020_Max", 
    "PM25_RWC_percent_2020_Min", "PM25_RWC_percent_2020_Max", 
])

# Display DataFrame
df_results


# %%
df_results

# %% [markdown]
# # Population Characteristics

# %% [markdown]
# ## Backend to generate pm25 per census tract

# %%
import pandas as pd
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

# %%
census_tract_df[['GISJOIN', 'STUSAB', 'REGIONA', 'DIVISIONA', 'STATE', 'STATEA',
       'COUNTY', 'COUNTYA', 'TRACTA', 'AITSCEA', 'CBSAA', 'CSAA', 'METDIVA',
       'UAA', 'GEOID', 'BTTRA']]

# %%
import geopandas as gpd
gdf = gpd.read_file("../SMOKE_sensitivity_analyses/US_census_shapefiles/US_tract_2018.shp")

# %%
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

# Create GeoDataFrame with polygons
polygons = [Polygon([(x1[i], y1[i]), (x2[i], y2[i]), (x3[i], y3[i]), (x4[i], y4[i])]) for i in range(len(x1))]
grid_gdf = gpd.GeoDataFrame(geometry=polygons, crs=proj_params)

# %%
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

# %%
var = 'PM25_TOT'
RWC_var = (baseline[var] - no_rwc[var])[0,0,:,:].to_numpy().ravel()
grid_gdf["RWC_"+var] = RWC_var
grid_gdf = grid_gdf.reset_index().rename(columns={'index': 'iD'})

# %%
grid_gdf = grid_gdf.rename(columns = {"RWC_PM25_TOT":"PM25_conc"})

# %%
# # Create intersections shapefile
grid_gdf = grid_gdf.reset_index().rename(columns={'index': 'iD'})
merged_shapes = merged_gdf[['GISJOIN','geometry']]
intersection = gpd.overlay(grid_gdf, merged_shapes, how='intersection')
intersection.to_file("intersection.shp")

# %% [markdown]
# 1) merge together
# 2) get the fraction of each census tract made up of each intersection
# 3) Sum the concentration * the fraction of the area to get a total averaged concentration
# 4) Save the geodataframe

# %%
import geopandas as gpd
intersectin = gpd.read_file('intersection_pm25/intersection.shp')

# %%
# intersectin = intersectin.merge(grid_gdf[['iD', 'PM25_conc']], on='iD')
# # merge on GISJOIN
# intersectin = intersectin.merge(census_tract_df[['GISJOIN', 'STUSAB', 'REGIONA', 'DIVISIONA', 'STATE', 'STATEA',
#        'COUNTY', 'COUNTYA', 'TRACTA', 'AITSCEA', 'CBSAA', 'CSAA', 'METDIVA',
#        'UAA', 'GEOID', 'BTTRA']], on='GISJOIN')
# intersectin = intersectin.drop('RWC_PM25_T', axis = 1)
# intersectin.to_pickle("grid_intersections.pkl")

# %%
from shapely.geometry import Polygon

# Assuming df is your DataFrame with a column of polygons
# Assuming the polygon column is named 'geometry'

# Define a function to calculate the area of a polygon
def calculate_polygon_area(polygon):
    return polygon.area

# Apply the function to calculate the area for each polygon
intersectin['census_tract_area'] = intersectin['geometry_y'].apply(calculate_polygon_area)
intersectin['intersection_area'] = intersectin['geometry_x'].apply(calculate_polygon_area)
intersectin['area_fraction'] = intersectin['intersection_area']/intersectin['census_tract_area']

# %%
# concentration * fraction column
# sum concentration fractions and save
intersectin['concentration_fraction'] = intersectin['RWC_PM25_TOT'] * intersectin['area_fraction']

# %%
summed_df = intersectin.groupby('GISJOIN')["concentration_fraction"].sum().reset_index()

# %%
census_tract_gdf = census_tract_gdf.merge(summed_df, on='GISJOIN')
census_tract_gdf['concentration_fraction'].mean()

# %%
census_tract_gdf.to_file("2016_rwc_census_tract_pm25.shp")

# %% [markdown]
# ## Population Data Science

# %% [markdown]
# ### Mortality Calculations

# %% [markdown]
# #### Import Census Tract Data and Add CBSA

# %%
import pandas as pd
EJ_index = pd.read_csv('EJ_index.csv')

# %%
import geopandas as gpd
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
    rwc_census_tract_baseline = rwc_census_tract_baseline.merge(EJ_index[['GEOID', 'SPL_EJI']], left_on='GEOID_x', right_on='GEOID', how='left')
    rwc_census_tract_baseline['COUNTYFP'] = rwc_census_tract_baseline['COUNTYFP'].astype(int)
    rwc_census_tract_baseline['STATEFP'] = rwc_census_tract_baseline['STATEFP'].astype(float)
    return rwc_census_tract_baseline

# %%
rwc_census_tract_2020 = import_rwc_shapefile(loc = 'census_tract_data/2020_rwc_census_tract_pm25.shp')
#rwc_census_tract_baseline = import_rwc_shapefile(loc = 'census_tract_data/2016_rwc_census_tract_pm25.shp')

# %%
tot_census_tract_2020 = import_rwc_shapefile(loc = 'census_tract_data/2020_tot_census_tract_pm25.shp')

# %%
tot_census_tract_2020['PM25_TOT'] = tot_census_tract_2020['emissions_']/tot_census_tract_2020.geometry.area * 4000 ** 2

# %%
rwc_census_tract_2020['PM25_TOT'] = tot_census_tract_2020['PM25_TOT']

# %%
rwc_census_tract_2020['PM25_TOT'].mean()

# %%
# MSA data
cbsa = gpd.read_file('cbsa/tl_2019_us_cbsa.shp', engine="pyogrio")

# Convert the GeoDataFrame to the target CRS
cbsa = cbsa.to_crs(rwc_census_tract_2020.crs)
cbsa['CBSA'] = cbsa['NAME']

# %%
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

# %%
#rwc_census_tract_baseline = add_CBSA(rwc_census_tract_baseline, cbsa)
rwc_census_tract_2020 = add_CBSA(rwc_census_tract_2020, cbsa)

# %% [markdown]
# #### Create baseline mortality rate statistics

# %%
import pandas as pd
mortality_data = pd.read_csv("BenMAP_mortality.csv")

# Ensure the column is of string type
mortality_data['Row'] = mortality_data['Row'].astype(str)
mortality_data['Column'] = mortality_data['Column'].astype(str)

# Step 1: Format 'Row' to be a 5-digit County FIPS code
mortality_data['County_FIPS'] = mortality_data['Row'].astype(str).str.zfill(5)

# Step 2: Format 'Column' to be a 6-digit Tract FIPS code
mortality_data['Tract_FIPS'] = mortality_data['Column'].astype(str).str.zfill(6)

# Step 3: Derive State FIPS Code (first 2 digits of the Tract FIPS code)
mortality_data['State_FIPS'] = mortality_data['County_FIPS'].str[:2]
mortality_data['County_FIPS'] = mortality_data['County_FIPS'].str[2:]

# Step 4: Create GISJOIN by concatenating 'G' + State FIPS + County FIPS + Tract FIPS + Suffix
mortality_data['GISJOIN'] = 'G' + mortality_data['State_FIPS'] + '0' +  mortality_data['County_FIPS'] + '0' + mortality_data['Tract_FIPS']

mortality_data["Age Range"] = mortality_data["Start Age"].astype(str) + "-" + mortality_data["End Age"].astype(str)

# Step 2: Pivot the table
pivoted_df = mortality_data.pivot_table(index=["GISJOIN", "Endpoint"],
                            columns="Age Range", values="Value").reset_index()

age_mapping = {
    "0-0": ["0-5"],
    "1-4": ["0-5"],
    "5-14": ["5-9", "10-14"],
    "15-24": ["15-17", "18-19", "20-21", "22-24"],
    "25-34": ["25-29", "30-34"],
    "35-44": ["35-39", "40-44"],
    "45-54": ["45-49", "50-54"],
    "55-64": ["55-59", "60-64"],
    "65-74": ["65-69", "70-74"],
    "75-84": ["75-79", "80-84"],
    "85-99": ["85+"]
}


# Function to sum columns based on the mapping
# Create a new DataFrame for the mapped age ranges
mapped_df = pd.DataFrame()
mapped_df['GISJOIN'] = pivoted_df['GISJOIN']

for mapped_range, age_ranges in age_mapping.items():
    # Sum the columns corresponding to the age ranges in the mapping
    for col in age_ranges:
        mapped_df[col] = pivoted_df[mapped_range]

# Add non-age columns (GISJOIN, Endpoint) to the new DataFrame
mapped_df['Endpoint'] = pivoted_df['Endpoint']
mortality_data = mapped_df

# %%
age_data = pd.read_csv("nhgis_ages/nhgis0006_ds244_20195_tract.csv")

# Define age range groups for summing, with column ranges for males and females
age_ranges = {
    '0-5': ['ALT0M003', 'ALT0M027'],
    '5-9': ['ALT0M004', 'ALT0M028'],
    '10-14': ['ALT0M005', 'ALT0M029'],
    '15-17': ['ALT0M006', 'ALT0M030'],
    '18-19': ['ALT0M007', 'ALT0M031'],
    '20-21': ['ALT0M008', 'ALT0M009', 'ALT0M032', 'ALT0M033'],
    '22-24': ['ALT0M010', 'ALT0M034'],
    '25-29': ['ALT0M011', 'ALT0M035'],
    '30-34': ['ALT0M012', 'ALT0M036'],
    '35-39': ['ALT0M013', 'ALT0M037'],
    '40-44': ['ALT0M014', 'ALT0M038'],
    '45-49': ['ALT0M015', 'ALT0M039'],
    '50-54': ['ALT0M016', 'ALT0M040'],
    '55-59': ['ALT0M017', 'ALT0M041'],
    '60-64': ['ALT0M018', 'ALT0M019', 'ALT0M042', 'ALT0M043'],
    '65-69': ['ALT0M020', 'ALT0M021', 'ALT0M044', 'ALT0M045'],
    '70-74': ['ALT0M022', 'ALT0M046'],
    '75-79': ['ALT0M023', 'ALT0M047'],
    '80-84': ['ALT0M024', 'ALT0M048'],
    '85+': ['ALT0M025', 'ALT0M049']
}

# Calculate total population for each age range
for age_range, columns in age_ranges.items():
    age_data[age_range] = age_data[columns].sum(axis=1)

# Drop the original columns
age_data = age_data.drop(columns=[col for cols in age_ranges.values() for col in cols])


for col in [
    "YEAR", "STUSAB", "REGIONA", "DIVISIONA", "STATE", "STATEA", "COUNTY", 
    "COUNTYA", "COUSUBA", "PLACEA", "TRACTA", "BLKGRPA", "CONCITA", "AIANHHA", 
    "RES_ONLYA", "TRUSTA", "AIHHTLI", "AITSCEA", "ANRCA", "CBSAA", "CSAA", 
    "METDIVA", "NECTAA", "CNECTAA", "NECTADIVA", "UAA", "CDCURRA", "SLDUA", 
    "SLDLA", "ZCTA5A", "SUBMCDA", "SDELMA", "SDSECA", "SDUNIA", "PCI", 
    "PUMAA", "GEOID", "BTTRA", "BTBGA", "NAME_E", "ALT0E001", "ALT0E002", 
    "ALT0E003", "ALT0E004", "ALT0E005", "ALT0E006", "ALT0E007", "ALT0E008", 
    "ALT0E009", "ALT0E010", "ALT0E011", "ALT0E012", "ALT0E013", "ALT0E014", 
    "ALT0E015", "ALT0E016", "ALT0E017", "ALT0E018", "ALT0E019", "ALT0E020", 
    "ALT0E021", "ALT0E022", "ALT0E023", "ALT0E024", "ALT0E025", "ALT0E026", 
    "ALT0E027", "ALT0E028", "ALT0E029", "ALT0E030", "ALT0E031", "ALT0E032", 
    "ALT0E033", "ALT0E034", "ALT0E035", "ALT0E036", "ALT0E037", "ALT0E038", 
    "ALT0E039", "ALT0E040", "ALT0E041", "ALT0E042", "ALT0E043", "ALT0E044", 
    "ALT0E045", "ALT0E046", "ALT0E047", "ALT0E048", "ALT0E049", "NAME_M", 
    "ALT0M001", "ALT0M002", "ALT0M026"
]:
    del age_data[col]


# %%
old_data = pd.read_csv("mortality_data.csv")
old_data['Tract ID'].unique().shape # missing lots of tracts

# %%
merged_data = pd.merge(age_data, mortality_data, on='GISJOIN', suffixes=('_population', '_mortality'))
# Compute mortality * population for each age bin
age_bins = ['0-5', '5-9', '10-14', '15-17', '18-19', '20-21', '22-24',
            '25-29', '30-34', '35-39', '40-44', '45-49', '50-54', '55-59',
            '60-64', '65-69', '70-74', '75-79', '80-84']#, '85+']

# Initialize a list to store the results
mortality_rates = []

# Iterate over each GISJOIN
for _, row in merged_data.iterrows():
    total_population = row[[col + '_population' for col in age_bins]].sum()  # Sum of population for all age bins
    total_mortality = sum(row[age_bin + '_mortality'] * row[age_bin + '_population'] for age_bin in age_bins)  # Sum of mortality * population for each bin
    overall_mortality_rate = total_mortality / total_population if total_population != 0 else 0  # Avoid division by zero
    mortality_rates.append(overall_mortality_rate)

# Add the calculated mortality rates to the DataFrame
merged_data['overall_mortality_rate'] = mortality_rates

# %%
# Age bins to consider for age > 25
age_bins_over_25 = ['25-29', '30-34', '35-39', '40-44', '45-49', '50-54', '55-59',
                    '60-64', '65-69', '70-74', '75-79', '80-84']#,'85+']

# Initialize a list to store the results
mortality_rates_over_25 = []

# Iterate over each GISJOIN
for _, row in merged_data.iterrows():
    total_population_over_25 = row[[col + '_population' for col in age_bins_over_25]].sum()  # Sum of population for age bins > 25
    total_mortality_over_25 = sum(row[age_bin + '_mortality'] * row[age_bin + '_population'] for age_bin in age_bins_over_25)  # Sum of mortality * population for age bins > 25
    overall_mortality_rate_over_25 = total_mortality_over_25 / total_population_over_25 if total_population_over_25 != 0 else 0  # Avoid division by zero
    mortality_rates_over_25.append(overall_mortality_rate_over_25)

# Add the calculated mortality rates for age > 25 to the DataFrame
merged_data['overall_mortality_rate_over_25'] = mortality_rates_over_25

# %%
merged_data.to_csv("merged_baseline_mortality_rate.csv")

# %% [markdown]
# #### Calculate Old Mortality Rates

# %%
mortality_data = pd.read_csv("mortality_data.csv")

# %%
# Define a function to convert Census Tract ID to GISJOIN
def tract_to_gisjoin(tract_id):
    # Extract state, county, and tract components
    tract_id = f'{tract_id:011}'
    state_code = tract_id[:2]
    county_code = tract_id[2:5]
    tract_code = tract_id[5:]
    
    # Format county code to meet GISJOIN requirements
    gisjoin_county_code = f'{int(county_code):03}'
    
    # Format tract code to meet GISJOIN requirements
    gisjoin_tract_code = f'{int(tract_code):07}'
    
    # Combine components to create GISJOIN
    gisjoin = f'G{state_code}0{gisjoin_county_code}{gisjoin_tract_code}'
    return gisjoin

# Apply function to the DataFrame
mortality_data['GISJOIN'] = mortality_data['Tract ID'].apply(tract_to_gisjoin)

# %%
mortality_data['Age Group'] = mortality_data['Age Group'].replace({'4-Jan': '1-4', '14-May': '5-14'})

# First, filter the rows with the age group '25-34'
split_25_34 = mortality_data[mortality_data['Age Group'] == '25-34']

# Create two new age groups: '25-30' and '30-34'
split_25_30 = split_25_34.copy()
split_30_34 = split_25_34.copy()

# Update the age group labels
split_25_30['Age Group'] = '25-30'
split_30_34['Age Group'] = '30-34'

# Now, combine these two back with the rest of the data
mortality_data = pd.concat([mortality_data[mortality_data['Age Group'] != '25-34'], split_25_30, split_30_34])

# Define a function to map age groups to their duration in years
def age_group_duration(age_group):
    if age_group == 'Under 1':
        return 1
    elif '-' in age_group:
        # Split the range, convert to integers, and find the difference + 1 for inclusivity
        start, end = map(int, age_group.split('-'))
        return end - start + 1
    else:
        # Handle "85 and older" or similar cases
        return 100  # Assuming a large number as an open-ended interval (you may adjust this)

# Apply the function to create the 'Total Years' column
mortality_data['Total Years'] = mortality_data['Age Group'].apply(age_group_duration)
mortality_data['Expected Deaths per Year'] = mortality_data['nd(x)'] / mortality_data['Total Years']
mortality_data['Mortality Rate'] = mortality_data['Expected Deaths per Year'] / mortality_data['l(x)']

# %%
# Step 2: Calculate weighted baseline mortality rate for each GISJOIN across all age groups
# Group by GISJOIN, then calculate a weighted average of the mortality rate using l(x) as weights
df_baseline_mortality = mortality_data.groupby('GISJOIN').apply(
    lambda x: x['Expected Deaths per Year'].sum() / x['l(x)'].sum()
).reset_index(name='Overall Baseline Mortality Rate')

# %%
df_baseline_mortality

# %% [markdown]
# #### Attributable Mortality Calculations

# %%
mortality_data = pd.read_csv("health_data/merged_baseline_mortality_rate.csv")

# %%
# pm_mortality_data_df_baseline =  pd.merge(rwc_census_tract_baseline, mortality_data, on = "GISJOIN")
# pm_mortality_data_df_baseline['2016_pm2.5'] = pm_mortality_data_df_baseline['concentrat']

pm_mortality_data_df_2020 =  pd.merge(rwc_census_tract_2020, mortality_data, on = "GISJOIN")
pm_mortality_data_df_2020['2016_pm2.5'] = pm_mortality_data_df_2020['RWC_PM25_T'] / pm_mortality_data_df_2020.geometry.area * 4000 ** 2
pm_mortality_data_df_2020['2016_pm10'] = pm_mortality_data_df_2020['RWC_PM10_w'] / pm_mortality_data_df_2020.geometry.area * 4000 ** 2

# %%
pm_mortality_data_df_2020['2016_pm2.5'].mean()

# %%
pm_mortality_data_df_2020['PM25_TOT'].mean()

# %%
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

# pm_mortality_data_df_baseline = mortality_function(pm_mortality_data_df_baseline, relative_risk = 1.08)
# mortality_function(pm_mortality_data_df_baseline, relative_risk  = 1.06)
# mortality_function(pm_mortality_data_df_baseline, relative_risk = 1.09)
# _ = mortality_function(pm_mortality_data_df_baseline, relative_risk = 1.17)

# %%
pm_mortality_data_df_2020 = mortality_function(pm_mortality_data_df_2020, relative_risk = 1.08)
pm_mortality_data_df_2020[['Attributable_Mortality_lower', 'Attributable_Mortality_Rate_lower']] = mortality_col_only(pm_mortality_data_df_2020, relative_risk  = 1.06)
pm_mortality_data_df_2020[['Attributable_Mortality_upper', 'Attributable_Mortality_Rate_upper']] = mortality_col_only(pm_mortality_data_df_2020, relative_risk  = 1.09)
_ = mortality_function(pm_mortality_data_df_2020, relative_risk = 1.17)

# %%
pm_mortality_data_df_2020['Attributable_Mortality'].sum()/pm_mortality_data_df_2020['Total'].sum() * 100_000

# %%
pm_mortality_data_df_2020['Attributable_Mortality_lower'].sum()/pm_mortality_data_df_2020['Total'].sum() * 100_000

# %%
pm_mortality_data_df_2020['Attributable_Mortality_upper'].sum()/pm_mortality_data_df_2020['Total'].sum() * 100_000

# %%
# pm_mortality_data_df_2020 = mortality_function(pm_mortality_data_df_2020, relative_risk = 1.08)
# pm_mortality_data_df_2020[['Attributable_Mortality_lower', 'Attributable_Mortality_Rate_lower']] = mortality_col_only(pm_mortality_data_df_2020, relative_risk  = 1.06)
# pm_mortality_data_df_2020[['Attributable_Mortality_upper', 'Attributable_Mortality_Rate_upper']] = mortality_col_only(pm_mortality_data_df_2020, relative_risk  = 1.09)
_ = mortality_function(pm_mortality_data_df_2020, relative_risk = 1.08, col = "PM25_TOT")
_ = mortality_function(pm_mortality_data_df_2020, relative_risk = 1.06, col = "PM25_TOT")
_ = mortality_function(pm_mortality_data_df_2020, relative_risk = 1.09, col = "PM25_TOT")

# %%
15529/70636

# %%
def short_term_mortality_function(pm_mortality_data_df, relative_risk=1.0132):
    # Convert yearly overall mortality rate to a daily rate
    pm_mortality_data_df['daily_mortality_rate'] = pm_mortality_data_df['overall_mortality_rate'] / 365
    
    # Calculate real relative risk for one day for 1 ug/m3
    pm_mortality_data_df['real_relative_risk'] = relative_risk ** (pm_mortality_data_df['2016_pm2.5'] / 10)
    
    # Calculate attributable fraction (AF) for short-term exposure
    pm_mortality_data_df['AF'] = (pm_mortality_data_df['real_relative_risk'] - 1) / pm_mortality_data_df['real_relative_risk']
    
    # Calculate daily attributable mortality
    pm_mortality_data_df['Daily_Attributable_Mortality'] = pm_mortality_data_df['Total'] * pm_mortality_data_df['daily_mortality_rate'] * pm_mortality_data_df['AF']
    
    # Scale daily mortality to the winter season (122 days)
    pm_mortality_data_df['Winter_Attributable_Mortality'] = pm_mortality_data_df['Daily_Attributable_Mortality'] * 122
    
    # Print summary of results
    print("Daily RR:", relative_risk, " - Winter Attributable Mortality 2020:", round(pm_mortality_data_df['Winter_Attributable_Mortality'].sum()))
    
    return pm_mortality_data_df

# %%
_ = short_term_mortality_function(pm_mortality_data_df_baseline, relative_risk = 1.0132)
_ = short_term_mortality_function(pm_mortality_data_df_baseline, relative_risk  = 1.0093)
_ = short_term_mortality_function(pm_mortality_data_df_baseline, relative_risk = 1.0172)

# %%
_ = short_term_mortality_function(pm_mortality_data_df_2020, relative_risk = 1.0132)
_ = short_term_mortality_function(pm_mortality_data_df_2020, relative_risk  = 1.0093)
_ = short_term_mortality_function(pm_mortality_data_df_2020, relative_risk = 1.0172)

# %%
pm_mortality_data_df_2020

# %% [markdown]
# ### Quick CONUS health stats

# %%
df_cbsa['Attributable_Mortality_Rate']

# %%
results = []

for cbsa in cbsa_list:
    df_cbsa = pm_mortality_data_df_2020[pm_mortality_data_df_2020["CBSA"] == cbsa]
     
    summed = df_cbsa[["Attributable_Mortality", "Attributable_Mortality_lower", "Attributable_Mortality_upper"]].sum()/df_cbsa['Total'].sum() * 100_000
    summed["CBSA"] = cbsa
    
    results.append(summed)

pd.DataFrame(results)

# %%
pm_mortality_data_df_2020.groupby("STATE")[["Attributable_Mortality", "Attributable_Mortality_lower", "Attributable_Mortality_upper"]].sum().reset_index()

# %%
# First calculation: Mortality rate per 100k
state_mortality_rate_df = pm_mortality_data_df_2020.groupby("STATE").apply(
    lambda x: x["Attributable_Mortality"].sum() / x["Total"].sum() * 100000
).reset_index(name="Mortality_Rate_per_100k")

# Second calculation: Summed mortality values
state_mortality_sums_df = pm_mortality_data_df_2020.groupby("STATE")[
    ["Attributable_Mortality", "Attributable_Mortality_lower", "Attributable_Mortality_upper"]
].sum().reset_index()

# Merge both DataFrames on 'STATE'
state_overall_df = pd.merge(state_mortality_sums_df, state_mortality_rate_df, on="STATE")

# Sort by Mortality Rate per 100k
state_overall_df = state_overall_df.sort_values(by="Mortality_Rate_per_100k", ascending=False)
state_overall_df.to_csv("State_health_results.csv")

# %%
pm_mortality_data_df_2020['population_density_per_sqmile'] = (
    pm_mortality_data_df_2020['Total'] / (pm_mortality_data_df_2020.geometry.area / 2_589_988.11)
)

pm_mortality_data_df_2020["Urban_Rural"] = pm_mortality_data_df_2020["population_density_per_sqmile"].apply(
    lambda x: "Urban" if x >= 500 else "Rural"
)

pm_mortality_data_df_2020.groupby("Urban_Rural")["Attributable_Mortality"].sum().reset_index()

# %%
pm_mortality_data_df_2020.groupby("Urban_Rural").apply(
    lambda x: x["Attributable_Mortality"].sum() / x["Total"].sum() * 100000
).reset_index(name="Mortality_Ratio")


# %%
pm_mortality_data_df_2020.loc[pm_mortality_data_df_2020["Urban_Rural"] == "Urban"]['Total'].sum()

# %%
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

pm_mortality_data_df_2020.groupby("Region")["Attributable_Mortality"].sum().reset_index()

# %%
pm_mortality_data_df_2020.groupby("Region")["Attributable_Mortality"].sum().reset_index()

# %%
pm_mortality_data_df_2020.groupby("Region").apply(
    lambda x: x["Attributable_Mortality"].sum() / x["Total"].sum() * 100000
).reset_index(name="Mortality_Ratio")


# %% [markdown]
# #### Economic Data Input

# %%
income_data = pd.read_csv("income_data/income_data.csv")
income_data = income_data[income_data['ALW1E001'] != '.']
income_data['median_household_income'] = income_data['ALW1E001'].astype(int)

# %%
pm_mortality_data_df_baseline = pd.merge(pm_mortality_data_df_baseline, income_data[['GISJOIN','median_household_income']], on = "GISJOIN")
pm_mortality_data_df_2020 = pd.merge(pm_mortality_data_df_2020, income_data[['GISJOIN','median_household_income']], on = "GISJOIN")

# %%
import matplotlib.pyplot as plt
def income_decile_comparison(df, loc = "", columns = ["2016_pm2.5", "Attributable_Mortality_Rate", "overall_mortality_rate_over_25"], version="baseline"):
    """
    Creates a figure with three panels comparing:
    - PM2.5 concentration
    - Attributable mortality rate
    - Overall mortality rate
    across income deciles.
    
    Parameters:
    - df: DataFrame containing the data.
    - version: Version label to include in the plot title.
    """
    # Define the columns for comparison
    #columns = ["2016_pm2.5", "Attributable_Mortality_Rate", "overall_mortality_rate_over_25"]
    
    labels = ["PM2.5 Concentration (µg/m³)", 
              "Attributable Mortality Rate (per 100,000)", 
              "Overall Mortality Rate (per 100,000)"]
    # Create income deciles
    # Copy and filter data if loc is specified
    df = df.copy(deep=True)
    if loc:
        df = df[df['CBSA'] == loc]
    
    df['income_decile'] = pd.qcut(df['median_household_income'], q=10, labels=range(1, 11))
    
    # Initialize the figure
    fig, axs = plt.subplots(1, 3, figsize=(18, 6), sharey=False)
    fig.suptitle(f"Comparison Across Income Deciles ({version})")

    for i, (col, label) in enumerate(zip(columns, labels)):
        # Calculate mean values of the column within each income decile
        decile_means = df.groupby('income_decile', observed = True)[col].mean()

        # Plot the results
        axs[i].bar(decile_means.index, decile_means.values, color='skyblue')
        axs[i].set_title(label)
        axs[i].set_xlabel("Income Deciles")
        axs[i].set_ylabel(f"Mean {label}")
        axs[i].set_xticks(range(1, 11))
        axs[i].set_xticklabels([str(d) for d in range(1, 11)])

    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

income_decile_comparison(pm_mortality_data_df_baseline, version="baseline")
income_decile_comparison(pm_mortality_data_df_2020, version="2020")


# %%
import pandas as pd
import matplotlib.pyplot as plt

def income_decile_comparison(df_baseline, df_2020, loc="", 
                              columns=["2016_pm2.5", "Attributable_Mortality_Rate", "overall_mortality_rate_over_25"]):
    """
    Creates a figure with three panels comparing:
    - PM2.5 concentration
    - Attributable mortality rate
    - Overall mortality rate
    across income deciles, showing baseline and 2020 data side by side.
    
    Parameters:
    - df_baseline: DataFrame containing the baseline data.
    - df_2020: DataFrame containing the 2020 data.
    - loc: Optional location filter for the data.
    - columns: Columns to compare across income deciles.
    """
    labels = ["PM2.5 Concentration (µg/m³)", 
              "Attributable Mortality Rate (per 100,000)", 
              "Overall Mortality Rate (per 100,000)"]
    
    # Filter data by location if specified
    df_baseline = df_baseline.copy()
    df_2020 = df_2020.copy()
    if loc:
        df_baseline = df_baseline[df_baseline['CBSA'] == loc]
        df_2020 = df_2020[df_2020['CBSA'] == loc]

    # Create income deciles for each dataset
    for df in [df_baseline, df_2020]:
        df['income_decile'] = pd.qcut(df['median_household_income'], q=10, labels=range(1, 11))
    
    # Initialize the figure
    fig, axs = plt.subplots(1, 3, figsize=(18, 6), sharey=False)
    fig.suptitle(f"Comparison Across Income Deciles {loc} (Baseline vs. 2020)")
    
    # Loop through each column and create a grouped bar chart
    bar_width = 0.4
    for i, (col, label) in enumerate(zip(columns, labels)):
        # Calculate mean values for each dataset
        decile_means_baseline = df_baseline.groupby('income_decile', observed=True)[col].mean()
        decile_means_2020 = df_2020.groupby('income_decile', observed=True)[col].mean()
        
        # Create positions for the bars
        x = range(1, 11)  # Deciles

        if col != 'overall_mortality_rate_over_25':
            x_baseline = [pos - bar_width / 2 for pos in x]
            x_2020 = [pos + bar_width / 2 for pos in x]
            axs[i].bar(x_2020, decile_means_2020.values, width=bar_width, label='2020', color='orange')
            axs[i].bar(x_baseline, decile_means_baseline.values, width=bar_width, label='Baseline', color='skyblue')

        else:
            x_baseline = [pos for pos in x]
            axs[i].bar(x_baseline, decile_means_baseline.values, width=bar_width*2,color='skyblue')

        # Set title, labels, and ticks
        axs[i].set_title(label)
        axs[i].set_xlabel("Income Deciles")
        axs[i].set_ylabel(f"Mean {label}")
        axs[i].set_xticks(range(1, 11))
        axs[i].set_xticklabels([str(d) for d in range(1, 11)])
        axs[i].legend()

    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

# Call the function with the baseline and 2020 datasets
income_decile_comparison(pm_mortality_data_df_baseline, pm_mortality_data_df_2020, loc="")


# %%
for cbsa in list(pm_mortality_data_df_2020.groupby('CBSA')[['Total', 'Attributable_Mortality']].sum().sort_values('Attributable_Mortality',ascending = False)[0:10].index):
    income_decile_comparison(pm_mortality_data_df_baseline, pm_mortality_data_df_2020, loc=cbsa)

# %% [markdown]
# ##### For Poster

# %%
for cbsa in list(pm_mortality_data_df_2020.groupby('CBSA')[['Total', 'Attributable_Mortality']].sum().sort_values('Attributable_Mortality',ascending = False)[0:10].index):
    print(cbsa)

# %%
import pandas as pd
import matplotlib.pyplot as plt

def income_decile_comparison(df_baseline, df_2020, col="Attributable_Mortality_Rate", 
                             cities = ["New York-Newark-Jersey City, NY-NJ-PA", "Los Angeles-Long Beach-Anaheim, CA", "Chicago-Naperville-Elgin, IL-IN-WI"]):
    """
    Creates a figure with three panels comparing:
    - PM2.5 concentration
    - Attributable mortality rate
    - Overall mortality rate
    across income deciles, showing baseline and 2020 data side by side.
    
    Parameters:
    - df_baseline: DataFrame containing the baseline data.
    - df_2020: DataFrame containing the 2020 data.
    - loc: Optional location filter for the data.
    - columns: Columns to compare across income deciles.
    """
        
    # Filter data by location if specified
    df_baseline = df_baseline.copy()
    df_2020 = df_2020.copy()
    

    # Initialize the figure
    fig, axs = plt.subplots(1, 3, figsize=(18, 6), sharey=False)
    #fig.suptitle(f"Comparison Across Income Deciles {loc} (Baseline vs. 2020)")

    i = 0
    for loc in cities:
        place = loc.split("-")[0]
        if place == "Los Angeles": place = "LA"
        if place == "New York": place = "NYC"
        
        df_baseline_loc = df_baseline[df_baseline['CBSA'] == loc]
        df_2020_loc = df_2020[df_2020['CBSA'] == loc]

        if "Attributable" in col:
            df_baseline_loc[col] = df_baseline_loc[col] * 100000
            df_2020_loc[col] = df_2020_loc[col] * 100000
            axs[i].set_title(f"{place} Income vs. Mortality", fontsize = 18)
            if i == 0:
                axs[i].set_ylabel(f"Attributable Annual Deaths per 100,000", fontsize = 14)
        else:
            axs[i].set_title(f"{place} Income vs. PM2.5 Concentration", fontsize = 18)
            if i == 0:
                axs[i].set_ylabel("PM2.5 Concentration (µg/m³)", fontsize = 14)
            
        # Create income deciles for each dataset
        for df in [df_baseline_loc, df_2020_loc]:
            df['income_decile'] = pd.qcut(df['median_household_income'], q=10, labels=range(1, 11))
        
        # Loop through each column and create a grouped bar chart
        bar_width = 0.4

        # Calculate mean values for each dataset
        decile_means_baseline = df_baseline_loc.groupby('income_decile', observed=True)[col].mean()
        decile_means_2020 = df_2020_loc.groupby('income_decile', observed=True)[col].mean()
            
        # Create positions for the bars
        x = range(1, 11)  # Deciles

        x_baseline = [pos - bar_width / 2 for pos in x]
        x_2020 = [pos + bar_width / 2 for pos in x]
        axs[i].bar(x_2020, decile_means_2020.values, width=bar_width, label='NEI 2020')
        axs[i].bar(x_baseline, decile_means_baseline.values, width=bar_width, label='NEI 2014')
        axs[i].tick_params(axis='y', labelsize=12)

        # Set title, labels, and ticks
        axs[i].set_xlabel("Median Household Income Decile", fontsize = 14)


        # Get decile bounds
        decile_bounds = df_baseline_loc['median_household_income'].quantile([i/10 for i in range(11)]).values
        decile_labels = [f"${decile_bounds[i]:,.0f}–${decile_bounds[i+1]:,.0f}" for i in range(10)]
        
        # # Set tick positions and labels
        # axs[i].set_xticks(range(1, 11))
        # axs[i].set_xticklabels(decile_labels, fontsize=12, rotation=0, ha="right")
        
        axs[i].set_xticks(range(1, 11))
        axs[i].set_xticklabels([str(d) for d in range(1, 11)], fontsize=12)

    
        axs[i].legend(fontsize = 14)
        axs[i].set_xlim([0.6, 10.4])
        
        i+=1

    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.subplots_adjust(wspace=0.15)  # Adjust horizontal and vertical spacing
    plt.subplots_adjust(bottom=0.1)  # Adjust bottom padding to 0.75 inches
    plt.show()

# Call the function with the baseline and 2020 datasets
income_decile_comparison(pm_mortality_data_df_baseline, pm_mortality_data_df_2020, col="Attributable_Mortality_Rate", 
                             cities = ["New York-Newark-Jersey City, NY-NJ-PA", "Los Angeles-Long Beach-Anaheim, CA", "Chicago-Naperville-Elgin, IL-IN-WI"])

# %%
income_decile_comparison(pm_mortality_data_df_baseline, pm_mortality_data_df_2020, col="2016_pm2.5", 
                             cities = ["New York-Newark-Jersey City, NY-NJ-PA", "Los Angeles-Long Beach-Anaheim, CA", "Chicago-Naperville-Elgin, IL-IN-WI"])

# %%
##

# %% [markdown]
# ### Emissions to Census Tract

# %%
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

# %%
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

# %%
# 3. Create gridded pollutant data
import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np

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

# %%
grid_gdf['PM25_total_tons'] = rwc_2020_new_surg['PM25_total_tons'].to_numpy().ravel()
grid_gdf = grid_gdf.reset_index().rename(columns={'index': 'iD'})



# %%
intersection = gpd.overlay(grid_gdf, census_tract_gdf, how='intersection')


# %%
intersection['fraction'] = intersection.geometry.area/(cell_size ** 2) # area of intersection / divided by area of gridcell
intersection['PM25_total_tons_per_polygon'] = intersection['PM25_total_tons'] * intersection['fraction'] # emissions of gridcell * area of intersection/area of gridcell

# %%
summed_df = intersection.groupby('GISJOIN')['PM25_total_tons_per_polygon'].sum().reset_index()

# %%
summed_df.shape

# %%
census_tract_gdf['PM25_total_tons_per_area'] = census_tract_gdf['PM25_total_tons_per_polygon']/census_tract_gdf.geometry.area * 1_000_000 * (1.609**2)

# %%
census_tract_gdf['PM25_total_tons_per_capita'] = census_tract_gdf['PM25_total_tons_per_polygon']/census_tract_gdf['Total'] * 100_000

# %%
census_tract_gdf = census_tract_gdf.merge(summed_df, on='GISJOIN')
census_tract_gdf.to_file("2020_rwc_census_tract_pm25_emissions.shp")

# %%
census_tract_gdf['White_percentage'] = census_tract_gdf['White alone']/census_tract_gdf['Total'] * 100 
census_tract_gdf['Black_percentage'] = census_tract_gdf["Black or African American alone"]/census_tract_gdf['Total'] * 100 
census_tract_gdf['American Indian_percentage'] = census_tract_gdf["American Indian and Alaska Native alone"]/census_tract_gdf['Total'] * 100 
census_tract_gdf['Asian_percentage'] = census_tract_gdf["Asian alone"]/census_tract_gdf['Total'] * 100 
census_tract_gdf["Native Hawaiian or Pacific Islander_percentage"] = census_tract_gdf["Native Hawaiian and Other Pacific Islander alone"]/census_tract_gdf['Total'] * 100 
census_tract_gdf["Hispanic_percentage"] = census_tract_gdf['Hispanic or Latino']/census_tract_gdf['Total'] * 100 
census_tract_gdf['Non-White_Fraction'] = 100 - census_tract_gdf['White_percentage']


# %%
import geopandas as gpd

# Load the CBSA data
cbsa = gpd.read_file('cbsa/tl_2019_us_cbsa.shp', engine="pyogrio")

# Convert CBSA data to the same CRS as census_tract_gdf
cbsa = cbsa.to_crs(census_tract_gdf.crs)

# Perform a spatial join to add CBSA information to census_tract_gdf
census_tract_gdf = census_tract_gdf.sjoin(cbsa[['geometry', 'NAME']], how="left", predicate="intersects")

# Rename the CBSA column for clarity
census_tract_gdf.rename(columns={'NAME': 'CBSA'}, inplace=True)


# %%
import pandas as pd

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
        if demo == ""
        for target in target_vars:
            msa_correlations[f"{demo.split('_')[0]}_vs_{target}"] = correlations.loc[demo, target]

    # Append to list
    emissions_correlation_results.append(msa_correlations)

# Convert results into a DataFrame
emissions_correlation_df = pd.DataFrame(emissions_correlation_results)
emissions_correlation_df

# %%
import matplotlib.pyplot as plt
import seaborn as sns

# Select relevant columns and rename them to keep only the race part
race_columns = [
    "Non_vs_PM25_total_tons_per_capita", 
    "Black_vs_PM25_total_tons_per_capita", 
    "American Indian_vs_PM25_total_tons_per_capita", 
    "Asian_vs_PM25_total_tons_per_capita", 
    "Native Hawaiian or Pacific Islander_vs_PM25_total_tons_per_capita", 
    "Hispanic_vs_PM25_total_tons_per_capita", 
    "White_vs_PM25_total_tons_per_capita"
]

# Create a new DataFrame with renamed columns
emissions_boxplot_df = emissions_correlation_df[race_columns].rename(
    columns=lambda x: x.split("_vs_")[0]  # Keep only the race part
)

# Create the boxplot
plt.figure(figsize=(10, 6))
sns.boxplot(data=emissions_boxplot_df)

# Add title and labels
plt.title('Boxplot of Correlations PM2.5 Emissions per Capita vs. Population Fraction')
plt.xticks(rotation=45)  # Rotate x-axis labels for readability
plt.ylabel("Pearson r Emissions vs. Population Fraction")
plt.xlabel("Race/Ethnicity")
plt.show()


# %%
emissions_correlation_df.round(2).to_csv("Emissions_race_correlations.csv")

# %% [markdown]
# ### Emissions by Species

# %%
rwc_2020_july =  xr.open_dataset('../SMOKE_sensitivity_analyses/rwc_2020_SMOKE_avg_201607.nc')


# %%
rwc_2020_july = add_pm25(rwc_2020_july)

seconds_in_a_day = 24 * 60 * 60  # Number of seconds in a day
days_in_month = 31  # Approximate days in a month (adjust if necessary)
seconds_in_a_month = seconds_in_a_day * days_in_month
grams_to_tons = 1 / 1_000_000  # Conversion factor from grams to tons
gps_to_tons = seconds_in_a_month * grams_to_tons
rwc_2020_july['PM25_total_tons'] = rwc_2020_july['PM25_total'] * gps_to_tons
pmgdf['PM25_total_tons_july'] = rwc_2020_july['PM25_total_tons'][0,0,:,:].to_numpy().flatten()
PWC_df_emissions = gpd.GeoDataFrame.merge(pmgdf, population_df, left_on=['row', 'col'], right_on=['ROWS', 'COLS'])


# %%
PWC_df_emissions['PM25_total_tons'].sum()

# %%
PWC_df_emissions['PM25_total_tons_july'].sum()/PWC_df_emissions['PM25_total_tons'].sum()

# %%
# july to jan comparison
jan_total = 0
for pm_component in ['PAL', 'PCA', 'PEC', 'PFE', 'PH2O', 'PK', 'PMG', 'PMN', 'PMOTHR', 'PNCOM', 'PNH4', 'PNO3', 'POC', 'PSI', 'PSO4', 'PTI']:
	try:
		jan_total += float(rwc_2020_new_surg[pm_component].sum())
		print(jan_total)
		#print(pm_component, round(float(rwc_2020_july[pm_component].sum()/rwc_2020_new_surg[pm_component].sum() * 100),4))
	except Exception as e:
		print(e)

# %%
jan_total = jan_total* gps_to_tons

# %%
jan_total

# %%
# july to jan comparison
dict = {}
for pm_component in ['PAL', 'PCA', 'PEC', 'PFE', 'PH2O', 'PK', 'PMG', 'PMN', 'PMOTHR', 'PNCOM', 'PNH4', 'PNO3', 'POC', 'PSI', 'PSO4', 'PTI']:
	try:
		july = float(rwc_2020_july[pm_component].sum() * gps_to_tons)
		jan = float(rwc_2020_new_surg[pm_component].sum() * gps_to_tons)
		dict[pm_component] = {'Jan': round(jan),
							  'July': round(july),
							  'Percent': round(100 * jan/jan_total, 1)}
	except Exception as e:
		print(e)


pd.DataFrame(dict).T

# %% [markdown]
# ### Emissions SMOKE Processed Plots

# %%
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

# %%
seconds_in_a_day = 24 * 60 * 60  # Number of seconds in a day
days_in_month = 31  # Approximate days in a month (adjust if necessary)
seconds_in_a_month = seconds_in_a_day * days_in_month
grams_to_tons = 1 / 1_000_000  # Conversion factor from grams to tons
gps_to_tons = seconds_in_a_month * grams_to_tons
rwc_2020_new_surg['PM25_total_tons'] = rwc_2020_new_surg['PM25_total'] * gps_to_tons

pmgdf['PM25_total_tons'] = rwc_2020_new_surg['PM25_total_tons'][0,0,:,:].to_numpy().flatten()
PWC_df_emissions = gpd.GeoDataFrame.merge(pmgdf, population_df, left_on=['row', 'col'], right_on=['ROWS', 'COLS'])


# %%
PWC_df_emissions

# %%
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

def plot_emissions(df, title, col="PM25", cmap="viridis", vs=False, counties=False, difference=False, state_color="white", ax=None, title_size=16):
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

	# # Add features to the map
	ax.add_feature(ocean, linewidth=0, zorder=0, facecolor="white")
	ax.add_feature(lakes, linewidth=0, zorder=0, facecolor="white")

	# Set visualization limits
	vmin = vs[0]
	vmax = vs[1]

	# Plot clipped data
	sc = df.plot(ax=ax, column=col, cmap=cmap, legend=False, vmin=vmin, vmax=vmax, zorder=-1, edgecolor="none", antialiased=False)

	# Create the ScalarMappable for the color bar
	norm = Normalize(vmin=vmin, vmax=vmax)
	sm = ScalarMappable(cmap=cmap, norm=norm)
	sm.set_array([])  # Required for colorbar

	# Add the title with adjusted position
	ax.set_title(title, fontsize=title_size, loc='center', pad=10)

	# Adjust layout to ensure that titles fit
	# print('done')
	ax.axis("off")
	return sc, sm

target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

fig, ax = plt.subplots(1, 1, figsize=(16, 14), subplot_kw={'projection': target_proj})

upper_bound = PWC_df_emissions['PM25_total_tons'].quantile(0.99)

sc, sm = plot_emissions(
	df=PWC_df_emissions,
	col="PM25_total_tons",
	title=f'CONUS January PM₂.₅ Emissions',
	cmap="cividis_r",
	vs=[0, upper_bound],
	counties=False,
	difference=False,
	ax=ax,
	title_size=24
)

 # Adjust colorbar
cbar = fig.colorbar(sm, ax=ax, orientation="vertical", shrink=0.7)
cbar.ax.set_ylabel("Janaury Total PM₂.₅ Emissions (metric tons)", rotation=90, fontsize=20,  labelpad=9)
cbar.ax.tick_params(labelsize=20)


# %%
rwc_2020_july =  xr.open_dataset('../SMOKE_sensitivity_analyses/rwc_2020_SMOKE_avg_201607.nc')


rwc_2020_july = add_pm25(rwc_2020_july)

seconds_in_a_day = 24 * 60 * 60  # Number of seconds in a day
days_in_month = 31  # Approximate days in a month (adjust if necessary)
seconds_in_a_month = seconds_in_a_day * days_in_month
grams_to_tons = 1 / 1_000_000  # Conversion factor from grams to tons
gps_to_tons = seconds_in_a_month * grams_to_tons
rwc_2020_july['PM25_total_tons_july'] = rwc_2020_july['PM25_total'] * gps_to_tons

pmgdf['PM25_total_tons_july'] = rwc_2020_july['PM25_total_tons_july'][0,0,:,:].to_numpy().flatten()
PWC_df_emissions = gpd.GeoDataFrame.merge(pmgdf, population_df, left_on=['row', 'col'], right_on=['ROWS', 'COLS'])


# %%
PWC_df_emissions['percent_july'] = PWC_df_emissions['PM25_total_tons_july']/PWC_df_emissions['PM25_total_tons'] * 100

# %%
target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

fig, ax = plt.subplots(1, 1, figsize=(16, 14), subplot_kw={'projection': target_proj})

upper_bound = PWC_df_emissions['PM25_total_tons_july'].quantile(0.99)

sc, sm = plot_emissions(
	df=PWC_df_emissions,
	col="PM25_total_tons_july",
	title=f'CONUS July PM₂.₅ Emissions',
	cmap="cividis_r",
	vs=[0, upper_bound],
	counties=False,
	difference=False,
	ax=ax,
	title_size=24
)

 # Adjust colorbar
cbar = fig.colorbar(sm, ax=ax, orientation="vertical", shrink=0.7)
cbar.ax.set_ylabel("July Total PM₂.₅ Emissions (metric tons)", rotation=90, fontsize=20,  labelpad=9)
cbar.ax.tick_params(labelsize=20)

# %%
target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

fig, ax = plt.subplots(1, 1, figsize=(16, 14), subplot_kw={'projection': target_proj})

upper_bound = 100 #PWC_df_emissions['percent_july'].quantile(0.99)

sc, sm = plot_emissions(
	df=PWC_df_emissions,
	col="percent_july",
	title=f'CONUS July % of Jan Emissions',
	cmap="magma_r",
	vs=[0, upper_bound],
	counties=False,
	difference=False,
	ax=ax,
	title_size=24
)

 # Adjust colorbar
cbar = fig.colorbar(sm, ax=ax, orientation="vertical", shrink=0.7)
cbar.ax.set_ylabel("CONUS July % of Jan PM₂.₅ Emissions", rotation=90, fontsize=20,  labelpad=9)
cbar.ax.tick_params(labelsize=20)

# %% [markdown]
# ### Emissions - Concentrations - Mortality Plots

# %%
import matplotlib.pyplot as plt
from pyproj import CRS, Transformer
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
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




# %%
import matplotlib.pyplot as plt
from pyproj import CRS, Transformer
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import box


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


# %%
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

pm_mortality_data_df_2020['Attributable_Mortality_Rate_100000'] = pm_mortality_data_df_2020['Attributable_Mortality_Rate'] * 100000

# %%
ny_df = pmgdf.loc[pmgdf['CBSA'] == cbsa_list[0]]
ny_df.loc[ny_df['PM25_RWC_2020'] >= 2].shape[0]/ny_df.shape[0] * 100

# %%
100 * ny_df.loc[ny_df['PM25_RWC_2020'] >= 2]['PM25_total_tons'].sum()/ny_df['PM25_total_tons'].sum()

# %%
import matplotlib.pyplot as plt
import numpy as np

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


# %%

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

for MSA in cbsa_list[0:1]:
    df = pmgdf

    upper_bound = 5
    if df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(0.98) >= 5:
        upper_bound = np.ceil(df.loc[df['CBSA'] == MSA]['PM25_RWC_2020'].quantile(0.98))

    emissions_upper_bound = 1
    if df.loc[df['CBSA'] == MSA]['PM25_total_tons'].quantile(0.98) >= 1:
        emissions_upper_bound = np.ceil(df.loc[df['CBSA'] == MSA]['PM25_total_tons'].quantile(0.98))
    
    # Create a Matplotlib figure with 1 row and 3 columns
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(28, 10), subplot_kw={'projection': target_proj}, dpi = 150)
    
    place = MSA.split('-')[0]
    
    # Panel A: Emissions map
    sc, sm = plot_MSA(
            df=df,
            col="PM25_total_tons",
            cbsa=cbsa,
            MSA=MSA,
            title=f'A {place} RWC Emissions',
            cmap="cividis_r",
            vs=[0, emissions_upper_bound],
            counties=False,
            difference=False,
            ax=axes[0],
            title_size=24
        )
    cbar = fig.colorbar(sm, ax=axes[0], orientation="vertical", shrink=0)
    cbar.ax.tick_params(labelsize=20)
    cbar.ax.set_ylabel("PM₂.₅ Tons Emitted", rotation=90, fontsize=20, labelpad=7)

    # Panel B: Concentration map
    sc, sm = plot_MSA(
            df=df,
            col="PM25_RWC_2020",
            cbsa=cbsa,
            MSA=MSA,
            title=f'B {place} RWC Concentrations',
            cmap="magma_r",
            vs=[0, upper_bound],
            counties=False,
            difference=False,
            ax=axes[1],
            title_size=24
        )
    cbar = fig.colorbar(sm, ax=axes[1], orientation="vertical", shrink=0)
    cbar.ax.tick_params(labelsize=20)
    cbar.ax.set_ylabel("PM₂.₅ Concentration (µg/m³)", rotation=90, fontsize=20, labelpad=7)

# %%
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

# %%
for MSA in cbsa_list[1:]:
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


# %%
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
df = pmgdf

# List of MSAs to plot
MSAs = [
    'Seattle-Tacoma-Bellevue, WA',
    'Minneapolis-St. Paul-Bloomington, MN-WI',
    'New York-Newark-Jersey City, NY-NJ-PA',
]

# Create a Matplotlib figure with 3 rows and 3 columns
target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(30, 30), subplot_kw={'projection': target_proj})

# Adjust spacing between subplots
fig.subplots_adjust(hspace=0.15, wspace=0.15)

# Iterate through MSAs and plot data
for i, MSA in enumerate(MSAs):
    place = MSA.split('-')[0]

    # Plot for 2014 PM2.5 data
    sc = plot_MSA(
        df=df,
        col="PM25_total_tons",
        cbsa=cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 from RWC (NEI 2014)',
        cmap="magma",
        #vs=[0, 5],
        counties=False,
        difference=False,
        ax=axes[i, 0],  # Use the first column for 2014 data
        title_size=18
    )

    # Plot for 2020 PM2.5 data
    sc = plot_MSA(
        df=df,
        col="PM25_RWC_2020",
        cbsa=cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 from RWC (NEI 2020)',
        cmap="magma",
        vs=[0, 5],
        counties=False,
        difference=False,
        ax=axes[i, 1],  # Use the second column for 2020 data
        title_size=18
    )

    # Plot for PM2.5 difference data
    sc = plot_MSA(
        df=df,
        col="PM25_diff",
        cbsa=cbsa,
        MSA=MSA,
        title=f'{place} PM2.5 Difference (2020 - 2014)',
        cmap="RdBu",
        vs=[-3, 3],
        counties=False,
        difference=True,
        ax=axes[i, 2],  # Use the third column for difference data
        title_size=18
    )


# Show the plot
plt.show()


# %% [markdown]
# #### Emissions Removed Species Checking

# %%
rwc_2020_jan_og = xr.open_dataset('../SMOKE_sensitivity_analyses/emis_mole_rwc_201601_original_2020_RWC.nc')

# %%
vars_to_delete = ['AACD', 'APIN', 'FACD', 'GLY', 'GLYD', 'ISPD', 'IVOC', 'MGLY', 'NMOG', 'PACD']
for var in vars_to_delete:


# %%
rwc_2020_jan_og['PACD']

# %%
print(rwc_2020_jan_og['PACD'][0,0,:,:].to_numpy().sum() * seconds_in_a_month /1_000_000) # metric tons)


# %%
print(rwc_2020_jan_og['NMOG'][0,0,:,:].to_numpy().sum() * seconds_in_a_month /1_000_000) # metric tons)


# %%

print(rwc_2020_jan_og['MGLY'][0,0,:,:].to_numpy().sum() * seconds_in_a_month ) # metric tons)
print(rwc_2020_jan_og['MGLY'][0,0,:,:].to_numpy().sum() * seconds_in_a_month * 72.06266 / 1_000_000) # metric tons

# %%
seconds_in_a_day = 24 * 60 * 60  # Number of seconds in a day
days_in_month = 31  # Approximate days in a month (adjust if necessary)
seconds_in_a_month = seconds_in_a_day * days_in_month

print("IVOC percent of VOC_INV", rwc_2020_jan_og['IVOC'][0,0,:,:].to_numpy().sum() * seconds_in_a_month/( rwc_2020_jan_og['VOC_INV'][0,0,:,:].to_numpy().sum() * seconds_in_a_month)) #* 58.04 / 1_000_000 # metric tons

# %%
seconds_in_a_dy = 24 * 60 * 60  # Number of seconds in a day
days_in_month = 31  # Approximate days in a month (adjust if necessary)
seconds_in_a_month = seconds_in_a_day * days_in_month
print(rwc_2020_jan_og['GLY'][0,0,:,:].to_numpy().sum() * seconds_in_a_month ) # metric tons)
print(rwc_2020_jan_og['GLY'][0,0,:,:].to_numpy().sum() * seconds_in_a_month * 58.04 / 1_000_000) # metric tons

# %% [markdown]
# #### Emissions Table with Spatial Correlation Emissions vs Concentratinos

# %%
PWC_df_emissions['population'] = PWC_df_emissions[6]

# %%
PWC_df_emissions['PM25_total_tons_per_person'] = PWC_df_emissions['PM25_total_tons']/PWC_df_emissions['population']

# %%
# MSA data
cbsa = gpd.read_file('cbsa/tl_2019_us_cbsa.shp', engine="pyogrio")

# Convert the GeoDataFrame to the target CRS
cbsa = cbsa.to_crs("+proj=lcc +lat_1=33 +lat_2=45 +lon_0=-97 +lat_0=40")

# %%
cbsa

# %%
cbsa.loc[cbsa['NAME']]

# %%
import pandas as pd
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
results_df

# %%
cbsa.loc[cbsa['NAME'] == MSA].geometry.area.sum()

# %%
results_df['correlation'].median()

# %%
!(pip3 install scipy)

# %% [markdown]
# 

# %% [markdown]
# ### Decile Plotting

# %%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def decile_plot(df, loc, col, version="baseline"):
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
    place = loc.split("-")[0]
    df = df.copy(deep=True)
    if loc:
        df = df[df['CBSA'] == loc]
    
    # Define race columns and number of deciles
    race_columns = ["White", "Black", "American Indian", "Asian", "Native Hawaiian or Pacific Islander",
                    "Other", "Two or More Races", "Hispanic"]
    num_iles = 10
    
    # Create the figure and subplot
    fig, ax = plt.subplots(figsize=(7, 8))
    #fig.suptitle(f"Race Composition within Deciles for {col} in {loc} ({version})", fontsize=14)
    
    # Create decile bins and apply them to the DataFrame
    bins = pd.qcut(df[col], q=num_iles, retbins=True) 
    df['decile'] = pd.qcut(df[col], q=num_iles, labels=False)
    
    # Get decile ranges from the bins
    decile_ranges = bins[1]* 100000
    tick_labels = [f"{value:.2f}" for value in decile_ranges]
    
    # Calculate race percentages within each decile
    for race in race_columns:
        df[race + '_percentage'] = (df[race] / df['Total']) * 100

    # Plot the horizontal decile data
    bottom = np.zeros(len(decile_ranges) - 1)
    for race in race_columns:
        race_percentages = [df[df['decile'] == decile][race + '_percentage'].mean() for decile in range(num_iles)]
        ax.barh(range(len(race_percentages)), race_percentages, left=bottom, height=0.9, label=race)
        bottom += np.array(race_percentages)
    
    measure = ""
    if col == "Attributable_Mortality_Rate":
        measure = "Mortality"
        label = "Attributable Mortality Rate (Deaths per 100,000)"
    if col == "2016_pm2.5":
        measure = "PM2.5"
        label = "PM2.5 Concentration (µg/m³)"

    ax.tick_params(axis='x', labelsize=12)
    ax.set_title(f"{place} {measure}", fontsize=18)
    ax.set_ylabel(label, labelpad=35, fontsize = 14)
    ax.set_xlabel('Percentage of Population', fontsize = 14)
    ax.set_yticks([])
    ax.set_ylim(-0.45, len(race_percentages) - 0.55)  # Adjust to fit the bars tightly
    ax.set_xlim(0,100)
    ax.set_yticklabels([])
    ax.set_xticks(np.linspace(0, 100, 11))

    # Adjust the first and last bars
    ax.text(0.75, -0.45, "0.00-", ha='right', va='center', fontsize=14, color="black")  # Adjust first tick label
    for idx, label in enumerate(tick_labels):
        if idx < len(tick_labels) - 2:
            ax.text(0.75, idx + 0.5, f"{label}-", ha='right', va='center', fontsize=12, color="black")
    
    # Adjust last bar (if needed)
    ax.text(0.75, len(tick_labels) - 1 - 0.55, f"{tick_labels[-2]}-", ha='right', va='center', fontsize=12, color="black")

    
    # Add legend
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3, fontsize=12)
    
    # Adjust layout
    #plt.tight_layout(rect=[0, 0, 1, 0.94])  # Adjust layout to fit title and legend
    plt.tight_layout(pad=0.0, rect=[0, 0, 0, 0])  # Adjust layout to remove padding
    plt.show()

# Example usage:
# Example DataFrame creation and function call would go here.


# %%
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



# %%
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
plt.legend(loc='upper center', bbox_to_anchor=(-1.5, -.18), ncol=4, fontsize=16, borderaxespad=0.)

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



# %%
# Create proxy artists for the legend

legend_patches = [patches.Patch(color=color, label=label) for label, color in race_colors.items()]

# Create a figure for just the legend
fig_legend = plt.figure(figsize=(12, 2))
ax = fig_legend.add_subplot(111)
ax.axis('off')

# Add legend with bounding box
legend = ax.legend(
    handles=legend_patches,
    loc='center',
    ncol=2,
    fontsize=18,
    frameon=True,
    fancybox=True,
    edgecolor='black'
)

# Optional: Manually set the bounding box line width and face color
legend.get_frame().set_linewidth(1.5)
legend.get_frame().set_facecolor('white')

plt.tight_layout()
plt.show()
fig_legend.savefig("main_decile_legend.png", format = 'png', dpi = 600, bbox_inches='tight')

# %%
fig.savefig("main_decile_with_legend.png", format="png", dpi=600, bbox_inches='tight')

# %% [markdown]
# #### Decile plots MSAs

# %%
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

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

# %%
plot_msa_deciles(cbsa_list[4:8], "E")

# %%
plot_msa_deciles(cbsa_list[8:12], "I")

# %%
plot_msa_deciles(cbsa_list[12:16], "N")

# %%
plot_msa_deciles(cbsa_list[16:], "R")

# %% [markdown]
# ### Relative Disparity Plot

# %% [markdown]
# #### 90% percentile by race

# %%
import pandas as pd
import numpy as np
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

df['Attributable_Mortality_Rate_100000'] = df['Attributable_Mortality_Rate'] * 100000

# Define race columns and number of deciles
race_columns = ["White", "Black", "American Indian", "Asian", "Native Hawaiian or Pacific Islander",
                "Other", "Two or More Races", "Hispanic"]

#The DataFrame 'df' has the necessary columns for the calculations
for race in race_columns:
    df[race + '_percentage'] = (df[race] / df['Total']) * 100

# Calculate CONUS averages
CONUS_PM_avg = df['2016_pm2.5'].mean()
CONUS_BMR_avg = df['overall_mortality_rate_over_25'].mean()
CONUS_PM_mortality_avg = df['Attributable_Mortality_Rate_100000'].mean()

# Dictionary to store the race-specific data
race_decile_averages = {}

# Loop through each race group to process
for race in race_columns:
    # Filter rows where the race percentage is in the top 10%
    top_decile_df = df[df[f'{race}_percentage'] >= df[f'{race}_percentage'].quantile(0.9)]
    
    # Calculate the averages for the three metrics in the top decile for this race
    race_decile_averages[race] = {
        "PM_90th": top_decile_df['2016_pm2.5'].mean(),
        "BMR_90th": top_decile_df['overall_mortality_rate_over_25'].mean(),
        "PM_mortality_90th": top_decile_df['Attributable_Mortality_Rate_100000'].mean()
    }


# Compare race-specific decile averages with CONUS averages
comparison = {}

for race, averages in race_decile_averages.items():
    comparison[race] = {
        "PM_90th_diff": (averages["PM_90th"] - CONUS_PM_avg) / CONUS_PM_avg * 100,
        "BMR_90th_diff": (averages["BMR_90th"] - CONUS_BMR_avg)/ CONUS_BMR_avg * 100,
        "PM_mortality_90th_diff": (averages["PM_mortality_90th"] - CONUS_PM_mortality_avg) / CONUS_PM_mortality_avg * 100
    }

# %%
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
pm_differences = [comparison[race]["PM_90th_diff"] for race in race_colors.keys()]
bmr_differences = [comparison[race]["BMR_90th_diff"] for race in race_colors.keys()]
pm_mortality_differences = [comparison[race]["PM_mortality_90th_diff"] for race in race_colors.keys()]

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
ax.set_title("Differences in RWC PM2.5, Baseline Mortality Rates, \nand RWC PM2.5 Attributable Mortality Rates \nfor Upper Decile of Racial/Ethnic Classifications", fontsize=22, pad=20)

# Add grid lines
ax.xaxis.grid(True, linestyle='solid', alpha=0.7)
for i in range(len(race_colors.keys())):
    ax.axhline(i, color='gray', linestyle='solid', linewidth=0.5, alpha = 0.7)

# Format x-axis labels as percentages
ax.set_xticks([-25, 0, 25])
ax.set_xticklabels([f"{x}%" for x in [-25, 0, 25]], fontsize=20)

# Legend
legend_labels = ["RWC PM2.5", "BMR", "RWC PM2.5 Att. Mort. Rate"]
legend_patches = [
    plt.Rectangle((0, 0), 1, 1, fc="gray"),
    plt.Rectangle((0, 0), 1, 1, edgecolor="black", facecolor="none"),
    plt.Line2D([0], [0], marker="o", color="black", linestyle="None", markersize=15)
]

plt.legend(legend_patches, legend_labels, loc='lower center', bbox_to_anchor=(0.25, -0.2), fontsize=20, ncol = 3)


# %% [markdown]
# #### All CONUS

# %%
import pandas as pd
import numpy as np
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
    

# %%
(df['overall_mortality_rate_over_25'] * df['Total']).sum() / df['Total'].sum()

# %%
bmr_PW

# %%

CONUS_PM_mortality_avg

# %%
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

# %%
print(CONUS_PM_avg, CONUS_BMR_avg, CONUS_PM_mortality_avg)

# %%
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
fig.savefig("main_relative_disparity.png", format="png", dpi=600, bbox_inches='tight')

# %% [markdown]
# #### Relative Disaparity by MSA

# %%
MSA_df

# %%
import pandas as pd

# Create an empty list to store results
correlation_results = []

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

# %%
pm_mortality_data_df_2020['RWC_PM25_T']

# %%
correlation_df[0:20]['Non_vs_2016_pm2.5'].mean()
#correlation_df['Non_vs_2016_pm2.5']

# %%
correlation_df[0:20]['Non_vs_Attributable_Mortality_Rate_100000'].mean()


# %%
correlation_df = correlation_df.round(2)
correlation_df.to_csv('race_correlations.csv')

# %%
correlation_results = []

for msa in cbsa_list:
	# Filter data for the current MSA
	MSA_df = pm_mortality_data_df_2020.loc[pm_mortality_data_df_2020["CBSA"] == msa]
	MSA_df['Non_White_Fraction'] = 100 - MSA_df['White_percentage']

	# Compute Pearson correlation coefficients
	correlations = MSA_df[["Non_White_Fraction", "2016_pm2.5", "Attributable_Mortality_Rate_100000"]].corr(method="pearson")

	# Extract relevant correlations
	rwc_corr = correlations.loc["Non_White_Fraction", "2016_pm2.5"]
	mortality_corr = correlations.loc["Non_White_Fraction", "Attributable_Mortality_Rate_100000"]

	# Append to list
	correlation_results.append({"MSA": msa, "RWC_PM25_Correlation": rwc_corr, "Mortality_Rate_Correlation": mortality_corr})

# Convert results into a DataFrame
correlation_df = pd.DataFrame(correlation_results)
correlation_df

# %%
def MSA_relative_disparity(df = pm_mortality_data_df_2020, MSA = "CONUS", ax = None):
	import pandas as pd
	import numpy as np
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



# %%
for msa in cbsa_list:
	MSA_df = pm_mortality_data_df_2020.loc[pm_mortality_data_df_2020['CBSA'] == msa]
	MSA_relative_disparity(df = MSA_df, MSA = msa)

# %%
# Create a 5-row, 4-column figure
fig, axes = plt.subplots(2, 2, figsize=(20, 17))  # Adjust size as needed
axes = axes.flatten()  # Flatten to iterate easily

# Loop through each CBSA and its corresponding subplot
for i, msa in enumerate(cbsa_list[0:4]):
    MSA_df = pm_mortality_data_df_2020.loc[pm_mortality_data_df_2020['CBSA'] == msa]
    
    # Call function but use the current subplot
    MSA_relative_disparity(df=MSA_df, MSA=msa, ax=axes[i])  

# Adjust layout for better spacing
plt.tight_layout()
plt.subplots_adjust(hspace=0.25)  # Adjust vertical spacing
plt.show()


# %%
# Create a 5-row, 4-column figure
fig, axes = plt.subplots(2, 2, figsize=(20, 17))  # Adjust size as needed
axes = axes.flatten()  # Flatten to iterate easily

# Loop through each CBSA and its corresponding subplot
for i, msa in enumerate(cbsa_list[4:8]):
    MSA_df = pm_mortality_data_df_2020.loc[pm_mortality_data_df_2020['CBSA'] == msa]
    
    # Call function but use the current subplot
    MSA_relative_disparity(df=MSA_df, MSA=msa, ax=axes[i])  

# Adjust layout for better spacing
plt.tight_layout()
plt.subplots_adjust(hspace=0.25)  # Adjust vertical spacing
plt.show()


# %%
# Create a 5-row, 4-column figure
fig, axes = plt.subplots(2, 2, figsize=(20, 17))  # Adjust size as needed
axes = axes.flatten()  # Flatten to iterate easily

# Loop through each CBSA and its corresponding subplot
for i, msa in enumerate(cbsa_list[8:12]):
    MSA_df = pm_mortality_data_df_2020.loc[pm_mortality_data_df_2020['CBSA'] == msa]
    
    # Call function but use the current subplot
    MSA_relative_disparity(df=MSA_df, MSA=msa, ax=axes[i])  

# Adjust layout for better spacing
plt.tight_layout()
plt.subplots_adjust(hspace=0.25)  # Adjust vertical spacing
plt.show()


# %%
# Create a 5-row, 4-column figure
fig, axes = plt.subplots(2, 2, figsize=(20, 17))  # Adjust size as needed
axes = axes.flatten()  # Flatten to iterate easily

# Loop through each CBSA and its corresponding subplot
for i, msa in enumerate(cbsa_list[12:16]):
    MSA_df = pm_mortality_data_df_2020.loc[pm_mortality_data_df_2020['CBSA'] == msa]
    
    # Call function but use the current subplot
    MSA_relative_disparity(df=MSA_df, MSA=msa, ax=axes[i])  

# Adjust layout for better spacing
plt.tight_layout()
plt.subplots_adjust(hspace=0.25)  # Adjust vertical spacing
plt.show()

# %%
# Create a 5-row, 4-column figure
fig, axes = plt.subplots(2, 2, figsize=(20, 17))  # Adjust size as needed
axes = axes.flatten()  # Flatten to iterate easily

# Loop through each CBSA and its corresponding subplot
for i, msa in enumerate(cbsa_list[16:]):
    MSA_df = pm_mortality_data_df_2020.loc[pm_mortality_data_df_2020['CBSA'] == msa]
    
    # Call function but use the current subplot
    MSA_relative_disparity(df=MSA_df, MSA=msa, ax=axes[i])  

# Adjust layout for better spacing
plt.tight_layout()
plt.subplots_adjust(hspace=0.25)  # Adjust vertical spacing
# Legend
legend_labels = ["RWC PM₂.₅", "BMR", "RWC PM₂.₅ Att. Mort. Rate"]
legend_patches = [
    plt.Rectangle((0, 0), 1, 1, fc="gray"),
    plt.Rectangle((0, 0), 1, 1, edgecolor="black", facecolor="none"),
    plt.Line2D([0], [0], marker="o", color="black", linestyle="None", markersize=20)
]

plt.legend(legend_patches, legend_labels, loc='lower center', bbox_to_anchor=(-0.85, -0.3), fontsize=25, ncol = 3)
plt.show()

# %%
# Create a 5-row, 4-column figure
fig, axes = plt.subplots(10, 2, figsize=(20, 100))  # Adjust size as needed
axes = axes.flatten()  # Flatten to iterate easily

# Loop through each CBSA and its corresponding subplot
for i, msa in enumerate(cbsa_list):
    MSA_df = pm_mortality_data_df_2020.loc[pm_mortality_data_df_2020['CBSA'] == msa]
    
    # Call function but use the current subplot
    MSA_relative_disparity(df=MSA_df, MSA=msa, ax=axes[i])  

# Adjust layout for better spacing
plt.tight_layout()
plt.subplots_adjust(hspace=0.25)  # Adjust vertical spacing

# Legend
legend_labels = ["RWC PM₂.₅", "BMR", "RWC PM₂.₅ Att. Mort. Rate"]
legend_patches = [
    plt.Rectangle((0, 0), 1, 1, fc="gray"),
    plt.Rectangle((0, 0), 1, 1, edgecolor="black", facecolor="none"),
    plt.Line2D([0], [0], marker="o", color="black", linestyle="None", markersize=20)
]

plt.legend(legend_patches, legend_labels, loc='lower center', bbox_to_anchor=(-0.85, -0.3), fontsize=25, ncol = 3)
plt.savefig("../full_MSA_disparity.svg", dpi=300, bbox_inches="tight")

#plt.show()


# %%
plt.savefig("../full_MSA_disparity.svg", dpi=300, bbox_inches="tight")


# %%
plt.savefig("../full_MSA_disparity.png", dpi=300, bbox_inches="tight")


# %% [markdown]
# #### 90th Percentile by Pollution

# %%
df['STATE'].value_counts()

# %%
pm_mortality_data_df_2020['STATE'].value_counts()

# %%
import pandas as pd
import numpy as np
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

upper_decile_threshold = pm_mortality_data_df_2020['2016_pm2.5'].quantile(0.90)

# Select rows where 2016pm2.5 is in the upper decile
df = pm_mortality_data_df_2020[pm_mortality_data_df_2020['2016_pm2.5'] >= upper_decile_threshold]

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
    

# %%
# Calculate CONUS averages
CONUS_PM_avg = pm_mortality_data_df_2020['2016_pm2.5'].mean()
CONUS_BMR_avg = pm_mortality_data_df_2020['overall_mortality_rate_over_25'].mean()
CONUS_PM_mortality_avg = pm_mortality_data_df_2020['Attributable_Mortality_Rate_100000'].mean()

# Compare race-specific decile averages with CONUS averages
comparison = {}

for race, averages in race_decile_averages.items():
    comparison[race] = {
        "PM_diff": (averages["PM"] - CONUS_PM_avg) / CONUS_PM_avg * 100,
        "BMR_diff": (averages["BMR"] - CONUS_BMR_avg)/ CONUS_BMR_avg * 100,
        "PM_mortality_diff": (averages["PM_mortality"] - CONUS_PM_mortality_avg) / CONUS_PM_mortality_avg * 100
    }

# %%
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
ax.set_title("Differences in RWC PM₂.₅, Baseline Mortality Rates, \nand RWC PM₂.₅ Attributable Mortality Rates\n for Upper Decile of RWC PM2.5 Concentrations", fontsize=22, pad=20)

# Add grid lines
ax.xaxis.grid(True, linestyle='solid', alpha=0.7)
for i in range(len(race_colors.keys())):
    ax.axhline(i, color='gray', linestyle='solid', linewidth=0.5, alpha = 0.7)

# Format x-axis labels as percentages
ax.set_xticks([ 0, 50, 100])
ax.set_xticklabels([f"{x}%" for x in [ 0, 50, 100]], fontsize=20)

# Legend
legend_labels = ["RWC PM₂.₅", "BMR", "RWC PM₂.₅ Att. Mort. Rate"]
legend_patches = [
    plt.Rectangle((0, 0), 1, 1, fc="gray"),
    plt.Rectangle((0, 0), 1, 1, edgecolor="black", facecolor="none"),
    plt.Line2D([0], [0], marker="o", color="black", linestyle="None", markersize=15)
]

plt.legend(legend_patches, legend_labels, loc='lower center', bbox_to_anchor=(0.25, -0.2), fontsize=20, ncol = 3)


# %% [markdown]
# ### Old Decile plots

# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define racial groups and colors
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

# Dummy data for the three metrics
num_groups = len(race_colors.keys())
np.random.seed(42)  # For reproducibility
NO2_differences = np.random.uniform(-30, 40, num_groups)
BMR_values = np.random.uniform(-30, 40, num_groups)
NO2_mortality_rates = np.random.uniform(-30, 40, num_groups)

fig, ax = plt.subplots(figsize=(8, 8))

y_pos = np.arange(num_groups)

# Plot NO2 differences as filled bars
for i, race in enumerate(race_colors.keys()):
    ax.barh(y_pos[i], NO2_differences[i], color=race_colors[race], edgecolor='black', height=0.8)

# Overlay BMR values as hollow bars
for i, race in enumerate(race_colors.keys()):
    ax.barh(y_pos[i], BMR_values[i], edgecolor=race_colors[race], facecolor='none', height=0.8, linewidth=2)

# Plot NO2 attributable mortality rates as black dots
ax.scatter(NO2_mortality_rates, y_pos, color='black', zorder=3)

# Formatting
ax.set_yticks(y_pos)
ax.set_yticklabels(race_colors.keys(), fontsize=16)
ax.axvline(0, color='black', linestyle='dashed', linewidth=1)
ax.set_xlabel("Percentage Difference Compared to CONUS", fontsize=16)
ax.set_title("Differences in RWC PM2.5, Baseline Mortality Rates, \nand RWC PM2.5 Attributable Mortality Rates", fontsize=20, pad=20)

# Add grid lines
ax.xaxis.grid(True, linestyle='dashed', alpha=0.7)

# Format x-axis labels as percentages
ax.set_xticks(np.arange(-30, 50, 10))
ax.set_xticklabels([f"{x}%" for x in np.arange(-30, 50, 10)], fontsize=14)

# Legend for race colors
legend_patches = [Patch(color=color, label=race) for race, color in race_colors.items()]
plt.legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=4, fontsize=12)

# Adjust layout
plt.tight_layout()
plt.show()


# %%
# Example usage for plotting on subplots:
fig, axs = plt.subplots(2, 4, figsize=(18, 12))
plt.tight_layout(pad=2.0, rect=[0, 0, 1.0, 1.0])  # Adjust layout to remove padding
plt.subplots_adjust(hspace=0.25, wspace=0.35)  # Adjust horizontal and vertical spacing
plt.subplots_adjust(bottom=0.1)  # Adjust bottom padding to 0.75 inches

# # Call decile_plot for each subplot
# decile_plot(pm_mortality_data_df_baseline, "Los Angeles-Long Beach-Anaheim, CA", "Attributable_Mortality_Rate", ax=axs[0, 0], first = True, year = "2014")
# decile_plot(pm_mortality_data_df_baseline, 'Minneapolis-St. Paul-Bloomington, MN-WI', "Attributable_Mortality_Rate", ax=axs[0, 1], year = "2014")
# decile_plot(pm_mortality_data_df_baseline, "New York-Newark-Jersey City, NY-NJ-PA", "Attributable_Mortality_Rate", ax=axs[0, 2], year = "2014")
# decile_plot(pm_mortality_data_df_baseline, "CONUS", "Attributable_Mortality_Rate", ax=axs[0, 3], year = "2014")


decile_plot(pm_mortality_data_df_2020, "Los Angeles-Long Beach-Anaheim, CA", "Attributable_Mortality_Rate", ax=axs[1, 0], first = True, year = "2020")
decile_plot(pm_mortality_data_df_2020, "Minneapolis-St. Paul-Bloomington, MN-WI", "Attributable_Mortality_Rate", ax=axs[1, 1], year = "2020")
decile_plot(pm_mortality_data_df_2020, "New York-Newark-Jersey City, NY-NJ-PA", "Attributable_Mortality_Rate", ax=axs[1, 2], year = "2020")
decile_plot(pm_mortality_data_df_2020, "CONUS", "Attributable_Mortality_Rate", ax=axs[1, 3], year = "2020")

plt.legend(loc='upper center', bbox_to_anchor=(-1.5, -.18), ncol=8, fontsize=12, borderaxespad=0.)
plt.show()

# %%
# Example usage for plotting on subplots:
fig, axs = plt.subplots(2, 4, figsize=(18, 12))
plt.tight_layout(pad=2.0, rect=[0, 0, 1.0, 1.0])  # Adjust layout to remove padding
plt.subplots_adjust(hspace=0.25, wspace=0.35)  # Adjust horizontal and vertical spacing
plt.subplots_adjust(bottom=0.1)  # Adjust bottom padding to 0.75 inches

# Call decile_plot for each subplot
decile_plot(pm_mortality_data_df_2020, "Los Angeles-Long Beach-Anaheim, CA", "Attributable_Mortality_Rate", ax=axs[0, 0], first = True, year = "2020")
decile_plot(pm_mortality_data_df_2020, 'Minneapolis-St. Paul-Bloomington, MN-WI', "Attributable_Mortality_Rate", ax=axs[0, 1], year = "2020")
decile_plot(pm_mortality_data_df_2020, "New York-Newark-Jersey City, NY-NJ-PA", "Attributable_Mortality_Rate", ax=axs[0, 2], year = "2020")
decile_plot(pm_mortality_data_df_2020, "CONUS", "Attributable_Mortality_Rate", ax=axs[0, 3], year = "2020")


decile_plot(pm_mortality_data_df_2020, "Los Angeles-Long Beach-Anaheim, CA", "2016_pm2.5", ax=axs[1, 0], first = True, year = "2020")
decile_plot(pm_mortality_data_df_2020, "Minneapolis-St. Paul-Bloomington, MN-WI", "2016_pm2.5", ax=axs[1, 1], year = "2020")
decile_plot(pm_mortality_data_df_2020, "New York-Newark-Jersey City, NY-NJ-PA", "2016_pm2.5", ax=axs[1, 2], year = "2020")
decile_plot(pm_mortality_data_df_2020, "CONUS", "2016_pm2.5", ax=axs[1, 3], year = "2020")

plt.legend(loc='upper center', bbox_to_anchor=(-1.5, -.18), ncol=8, fontsize=12, borderaxespad=0.)
plt.show()

# %%
cbsa_list

# %%
# Example usage for plotting on subplots:
fig, axs = plt.subplots(2, 4, figsize=(18, 12))
plt.tight_layout(pad=2.0, rect=[0, 0, 1.0, 1.0])  # Adjust layout to remove padding
plt.subplots_adjust(hspace=0.25, wspace=0.35)  # Adjust horizontal and vertical spacing
plt.subplots_adjust(bottom=0.1)  # Adjust bottom padding to 0.75 inches

# Call decile_plot for each subplot
decile_plot(pm_mortality_data_df_baseline, "Los Angeles-Long Beach-Anaheim, CA", "2016_pm2.5", ax=axs[0, 0], first = True, year = "2014")
decile_plot(pm_mortality_data_df_baseline, 'Minneapolis-St. Paul-Bloomington, MN-WI', "2016_pm2.5", ax=axs[0, 1], year = "2014")
decile_plot(pm_mortality_data_df_baseline, "New York-Newark-Jersey City, NY-NJ-PA", "2016_pm2.5", ax=axs[0, 2], year = "2014")
decile_plot(pm_mortality_data_df_baseline, "CONUS", "2016_pm2.5", ax=axs[0, 3], year = "2014")


decile_plot(pm_mortality_data_df_2020, "Los Angeles-Long Beach-Anaheim, CA", "2016_pm2.5", ax=axs[1, 0], first = True, year = "2020")
decile_plot(pm_mortality_data_df_2020, "Minneapolis-St. Paul-Bloomington, MN-WI", "2016_pm2.5", ax=axs[1, 1], year = "2020")
decile_plot(pm_mortality_data_df_2020, "New York-Newark-Jersey City, NY-NJ-PA", "2016_pm2.5", ax=axs[1, 2], year = "2020")
decile_plot(pm_mortality_data_df_2020, "CONUS", "2016_pm2.5", ax=axs[1, 3], year = "2020")

plt.legend(loc='upper center', bbox_to_anchor=(-1.5, -.18), ncol=8, fontsize=12, borderaxespad=0.)
plt.show()

# %%
# Example usage for plotting on subplots:
fig, axs = plt.subplots(2, 3, figsize=(18, 12))
plt.tight_layout(pad=2.0, rect=[0, 0, 0, 0])  # Adjust layout to remove padding
plt.subplots_adjust(hspace=0.25, wspace=0.35)  # Adjust horizontal and vertical spacing
plt.subplots_adjust(bottom=0.1)  # Adjust bottom padding to 0.75 inches

# Call decile_plot for each subplot
decile_plot(pm_mortality_data_df_baseline, "Los Angeles-Long Beach-Anaheim, CA", "2016_pm2.5", ax=axs[0, 0], first = True, year = "2014")
decile_plot(pm_mortality_data_df_baseline, 'Chicago-Naperville-Elgin, IL-IN-WI', "2016_pm2.5", ax=axs[0, 1], year = "2014")
decile_plot(pm_mortality_data_df_baseline, "New York-Newark-Jersey City, NY-NJ-PA", "2016_pm2.5", ax=axs[0, 2], year = "2014")
decile_plot(pm_mortality_data_df_2020, "Los Angeles-Long Beach-Anaheim, CA", "2016_pm2.5", ax=axs[1, 0], first = True, year = "2020")
decile_plot(pm_mortality_data_df_2020, "Chicago-Naperville-Elgin, IL-IN-WI", "2016_pm2.5", ax=axs[1, 1], year = "2020")
decile_plot(pm_mortality_data_df_2020, "New York-Newark-Jersey City, NY-NJ-PA", "2016_pm2.5", ax=axs[1, 2], year = "2020")
plt.legend(loc='upper center', bbox_to_anchor=(-1.0, -.18), ncol=8, fontsize=12, borderaxespad=0.)
plt.show()

# %%
# Example usage for plotting on subplots:
fig, axs = plt.subplots(2, 3, figsize=(18, 12))
plt.tight_layout(pad=2.0, rect=[0, 0, 0, 0])  # Adjust layout to remove padding
plt.subplots_adjust(hspace=0.25, wspace=0.35)  # Adjust horizontal and vertical spacing
plt.subplots_adjust(bottom=0.1)  # Adjust bottom padding to 0.75 inches

# Call decile_plot for each subplot
decile_plot(pm_mortality_data_df_baseline, "Los Angeles-Long Beach-Anaheim, CA", "2016_pm2.5", ax=axs[0, 0], first = True, year = "2014")
decile_plot(pm_mortality_data_df_baseline, 'Chicago-Naperville-Elgin, IL-IN-WI', "2016_pm2.5", ax=axs[0, 1], year = "2014")
decile_plot(pm_mortality_data_df_baseline, "New York-Newark-Jersey City, NY-NJ-PA", "2016_pm2.5", ax=axs[0, 2], year = "2014")
decile_plot(pm_mortality_data_df_2020, "Los Angeles-Long Beach-Anaheim, CA", "2016_pm2.5", ax=axs[1, 0], first = True, year = "2020")
decile_plot(pm_mortality_data_df_2020, "Chicago-Naperville-Elgin, IL-IN-WI", "2016_pm2.5", ax=axs[1, 1], year = "2020")
decile_plot(pm_mortality_data_df_2020, "New York-Newark-Jersey City, NY-NJ-PA", "2016_pm2.5", ax=axs[1, 2], year = "2020")
plt.legend(loc='upper center', bbox_to_anchor=(-1.0, -.18), ncol=8, fontsize=12, borderaxespad=0.)
plt.show()

# %%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def decile_plot_3_panels(df, loc, columns, version="baseline"):
    """
    Creates a figure with three side-by-side decile bar plots, each showing the
    racial composition within deciles of different data columns.
    
    Parameters:
    - df: DataFrame containing the data.
    - loc: Location filter for the 'CBSA Title' column. If empty, the whole DataFrame is used.
    - columns: List of three columns to create decile plots for.
    - version: Version label to include in the plot title.
    """
    
    # Copy and filter data if loc is specified
    df = df.copy(deep=True)
    if loc:
        df = df[df['CBSA'] == loc]
    
    # Define race columns and number of deciles
    race_columns = ["White", "Black", "American Indian", "Asian", "Native Hawaiian or Pacific Islander",
                    "Other", "Two or More Races", "Hispanic"]
    num_iles = 10
    
    fig, axs = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    fig.suptitle(f"Race Composition within Deciles for Different Indicators in {loc} ({version})")
    
    for i, col in enumerate(columns):
        # Prepare data for the decile plot
        df['decile'] = pd.qcut(df[col], q=num_iles, labels=False)
        
        # Calculate race percentages within each decile
        for race in race_columns:
            df[race + '_percentage'] = (df[race] / df['Total']) * 100

        # Plot the decile data in the corresponding subplot
        bottom = np.zeros(len(df['decile'].unique()))
        for race in race_columns:
            race_percentages = [df[df['decile'] == decile][race + '_percentage'].mean() for decile in range(num_iles)]
            axs[i].bar(range(1, num_iles + 1), race_percentages, bottom=bottom, label=race)
            bottom += np.array(race_percentages)
        
        axs[i].set_title(f"{col}")
        axs[i].set_xlabel(f'Deciles of {col}')
        axs[i].set_xticks(range(1, num_iles + 1))
        if i == 0:
            axs[i].set_ylabel('Percentage of Population')
    
    # Add legend to the first subplot only to avoid repetition
    axs[0].legend(loc="upper right", bbox_to_anchor=(1.4, 1))
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to fit title
    plt.show()

# Example usage:


# %%
decile_plot_3_panels(pm_mortality_data_df_baseline, "", ["Attributable_Mortality_Rate", "2016_pm2.5", "overall_mortality_rate_over_25"], "baseline CONUS pm2.5")
decile_plot_3_panels(pm_mortality_data_df_2020, "", ["Attributable_Mortality_Rate", "2016_pm2.5", "overall_mortality_rate_over_25"], "2020 CONUS pm2.5")

# %%

for cbsa in list(pm_mortality_data_df_2020.groupby('CBSA')[['Total', 'Attributable_Mortality']].sum().sort_values('Attributable_Mortality',ascending = False)[0:10].index):
    decile_plot_3_panels(pm_mortality_data_df_baseline, loc=cbsa, columns = ["Attributable_Mortality_Rate", "2016_pm2.5", "overall_mortality_rate_over_25"], version = "baseline CONUS pm2.5")
    decile_plot_3_panels(pm_mortality_data_df_2020, loc=cbsa,  columns = ["Attributable_Mortality_Rate", "2016_pm2.5", "overall_mortality_rate_over_25"],version =  "2020 CONUS pm2.5")

# %% [markdown]
# ### Spatial plotting

# %%
import matplotlib.pyplot as plt
from pyproj import CRS, Transformer
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import box


def plot_MSA_mortality(df, MSA, title, validation_df, col="PM25", cmap="viridis", vs=False, counties=False, difference=False, state_color="white", ax=None):
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
    ax.add_feature(ocean, linewidth=0.2, zorder=0)
    ax.add_feature(lakes, linewidth=0.2, zorder=0)
    ax.add_feature(states_provinces, edgecolor=state_color, linewidth=0.5, zorder=0)
    msa_bounds = df.loc[df["CBSA"] == MSA].total_bounds  # [minx, miny, maxx, maxy]
    
    buffer_size = 4000  # Example buffer size in the unit of the projection (e.g., meters)
    
    # Apply the buffer to expand the bounding box
    buffered_bounds = [
        msa_bounds[0] - buffer_size,  # minx
        msa_bounds[1] - buffer_size,  # miny
        msa_bounds[2] + buffer_size,  # maxx
        msa_bounds[3] + buffer_size   # maxy
    ]
    
    # Create a box polygon with the buffered bounds
    plot_area = box(*buffered_bounds)
    ax.set_extent([msa_bounds[0] -buffer_size , msa_bounds[2] + buffer_size, msa_bounds[1] - buffer_size, msa_bounds[3] + buffer_size], crs=target_proj)

    
    # Filter the main dataframe for data within the buffered plot area
    sub_df = df[df.geometry.intersects(plot_area)]
    
    if vs:
        vmin = vs[0]
        vmax = vs[1]
    elif not difference:
        vmin = np.quantile(sub_df[col], 0.025)
        vmax = np.quantile(sub_df[col], 0.975)
    else:
        vmin = -max( [abs(np.quantile(sub_df[col], 0.025)), abs(np.quantile(sub_df[col], 0.975))] )
        vmax =-vmin
        
    # # Plot MSA area data
    sc = sub_df.plot(ax=ax, column=col, cmap=cmap, legend=True, vmin=vmin, vmax=vmax, zorder=-1, legend_kwds={'label': col})
    
    # cbar = plt.colorbar(sc, ax=ax, orientation="vertical", fraction=0.02, pad=0.04)
    # cbar.set_label("Deaths per 10000")
    
    msa_geometry = df.loc[df["CBSA"] == MSA, "geometry"].unary_union  # Get combined geometry for the MSA
    
    # Simplify the geometry to reduce jagged edges
    tolerance = 4000  # Adjust this value as needed for smoother lines (in coordinate units)
    #simplified_msa_geometry = msa_geometry.simplify(tolerance, preserve_topology=True)
    
    # Plot the simplified geometry
   # ax.add_geometries([msa_geometry], crs=target_proj, edgecolor=state_color, facecolor='none', linewidth=1.5)
    
    # Add title and return the subset of the data
    ax.set_title(title)
    return sc

def plot_multiple_MSA(df, df2, MSA, figsize = (15,8)):
    # Create a Matplotlib figure with 1 row and 3 columns
    target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=figsize, subplot_kw={'projection': target_proj})

    place = MSA.split('-')[0]

    sc1 = plot_MSA_mortality(
        df=df,
        col="Attributable_Mortality_Rate_10000",
        MSA=MSA,
        title=f'{place} PM2.5 RWC baseline',
        validation_df=pd.DataFrame(),  # No validation points for difference
        cmap="viridis",
        counties=False,
        difference=False,
        ax=axes[0],
        state_color="black"# Pass the third axis
    )

    sc2 = plot_MSA_mortality(
        df=df2,
        col="Attributable_Mortality_Rate_10000",
        MSA=MSA,
        title=f'{place} PM2.5 RWC 2020',
        validation_df=pd.DataFrame(),  # No validation points for difference
        cmap="viridis",
        vs=False,
        counties=False,
        difference=False,
        ax=axes[1],
        state_color="black"# Pass the third axis
    )

    df['Attributable_Mortality_Rate_10000_difference'] = -df['Attributable_Mortality_Rate_10000'] + df2['Attributable_Mortality_Rate_10000']

    sc3 = plot_MSA_mortality(
        df=df,
        col="Attributable_Mortality_Rate_10000_difference",
        MSA=MSA,
        title=f'{place} PM2.5 RWC 2020 - baseline',
        validation_df=pd.DataFrame(),  # No validation points for difference
        cmap="RdBu",
        vs=False,
        counties=False,
        difference=True,
        ax=axes[2],
        state_color="black"# Pass the third axis
    )


    baseline_deaths = round(df.loc[df["CBSA"] == MSA]['Attributable_Mortality'].sum())
    baseline_population = round(df.loc[df["CBSA"] == MSA]['Total'].sum())
    print("Baseline Premature Mortality:", baseline_deaths)
    print("Baseline Total Population:",baseline_population)

    baseline_per_10000 = round((baseline_deaths/baseline_population * 10000),2)
    print("Baseline Per 10000:",baseline_per_10000)

    deaths_2020 = round(df2.loc[df2["CBSA"] == MSA]['Attributable_Mortality'].sum())
    population_2020 = round(df2.loc[df2["CBSA"] == MSA]['Total'].sum())
    print("2020     Premature Mortality:", deaths_2020)
    print("2020     Total Population:",population_2020)

    per_10000_2020 = round((deaths_2020/population_2020 * 10000),2)
    print("2020     Per 10000:",per_10000_2020)

    plt.show()

# %%
list(pm_mortality_data_df_2020.groupby('CBSA')[['Total', 'Attributable_Mortality']].sum().sort_values('Attributable_Mortality',ascending = False)[0:10].index)

# %%
pm_mortality_data_df_baseline['2016_pm2.5']

# %%
result = (
    pm_mortality_data_df_baseline
    .groupby('CBSA')
    .agg({
        'Total': 'sum',
        'Attributable_Mortality': 'sum',
        'Attributable_Mortality_Rate': 'mean' * 100000
    })
    .sort_values('Total', ascending=False)
    .head(20)  # Top 10 rows
)

result

# %%
def get_results_table(pm_mortality_data_df_baseline):
    pm_mortality_data_df_baseline['area'] = pm_mortality_data_df_baseline.geometry.area
    pm_mortality_data_df_baseline['Weighted_PM2.5'] = pm_mortality_data_df_baseline['2016_pm2.5'] * pm_mortality_data_df_baseline['area']
    pm_mortality_data_df_baseline['population_weighted_total'] = pm_mortality_data_df_baseline['2016_pm2.5'] * pm_mortality_data_df_baseline['Total']
    pm_mortality_data_df_baseline['Attributable_Mortality_Rate_per_100000'] = pm_mortality_data_df_baseline['Attributable_Mortality_Rate'] * 100000
    
    # Group by CBSA and calculate necessary metrics
    result = (
        pm_mortality_data_df_baseline
        .groupby('CBSA')
        .agg({
            'Total': 'sum',  # Total population
            'Attributable_Mortality': 'sum',  # Attributable mortality
            'Attributable_Mortality_lower': 'sum',  # Lower bound of attributable mortality
            'Attributable_Mortality_upper': 'sum',  # Upper bound of attributable mortality
            'Weighted_PM2.5': 'sum',  # Sum of weighted PM2.5
            'area': 'sum',  # Total area for weighting
            'population_weighted_total': 'sum'  # Population weighted total PM2.5
        })
    )
    
    # Calculate weighted average 2016pm2.5
    result['Average_RWC_PM2.5'] = result['Weighted_PM2.5'] / result['area']
    
    # Calculate mortality rates per 100,000 people
    result['Attributable_Mortality_per_100000'] = (result['Attributable_Mortality'] / result['Total']) * 100000
    result['Attributable_Mortality_lower_per_100000'] = (result['Attributable_Mortality_lower'] / result['Total']) * 100000
    result['Attributable_Mortality_upper_per_100000'] = (result['Attributable_Mortality_upper'] / result['Total']) * 100000
    
    # Sort by Attributable Mortality and select top 10 rows
    result = result.sort_values('Total', ascending=False)
    
    # Drop intermediate columns to clean up
    result = result.drop(columns=['Weighted_PM2.5', 'area', 'population_weighted_total'])
    
    return result

# %%
result[0:20]['Attributable_Mortality'].sum() / pm_mortality_data_df_baseline['Attributable_Mortality'].sum()

# %%
result[0:20]['Total'].sum() / pm_mortality_data_df_baseline['Total'].sum()

# %%
pm_mortality_data_df_2020['Attributable_Mortality_Rate_lower_per_100000'] = pm_mortality_data_df_2020['Attributable_Mortality_Rate_lower'] * 100000
pm_mortality_data_df_2020['Attributable_Mortality_Rate_upper_per_100000'] = pm_mortality_data_df_2020['Attributable_Mortality_Rate_upper'] * 100000

# %%
result = get_results_table(pm_mortality_data_df_2020)
# Display the result
result[['Attributable_Mortality_per_100000', 'Attributable_Mortality_lower_per_100000', 'Attributable_Mortality_upper_per_100000']].head(20)

# %%
result.head(20)


# %%
result[0:100]['Attributable_Mortality'].sum() / pm_mortality_data_df_2020['Attributable_Mortality'].sum()

# %%
result[0:100]['Total'].sum() / pm_mortality_data_df_2020['Total'].sum()

# %%
pm_mortality_data_df_2020.groupby('CBSA')[['Total', 'Attributable_Mortality']].sum().sort_values('Attributable_Mortality',ascending = False)[0:10]

# %%
pm_mortality_data_df_baseline['Attributable_Mortality_Rate_10000'] = pm_mortality_data_df_baseline['Attributable_Mortality_Rate'] * 10000
pm_mortality_data_df_2020['Attributable_Mortality_Rate_10000'] = pm_mortality_data_df_2020['Attributable_Mortality_Rate'] * 10000

# %%
plot_multiple_MSA(pm_mortality_data_df_baseline, pm_mortality_data_df_2020,  "Boston-Cambridge-Newton, MA-NH", figsize = (17, 6))

# %%
plot_multiple_MSA(pm_mortality_data_df_baseline, pm_mortality_data_df_2020,  "Detroit-Warren-Dearborn, MI", figsize = (17, 6))

# %%
plot_multiple_MSA(pm_mortality_data_df_baseline, pm_mortality_data_df_2020,  "Houston-The Woodlands-Sugar Land, TX", figsize = (17, 6))

# %%
plot_multiple_MSA(pm_mortality_data_df_baseline, pm_mortality_data_df_2020,  "Los Angeles-Long Beach-Anaheim, CA", figsize = (17, 6))

# %%
plot_multiple_MSA(pm_mortality_data_df_baseline, pm_mortality_data_df_2020,  "Philadelphia-Camden-Wilmington, PA-NJ-DE-MD", figsize = (17, 6))

# %%
plot_multiple_MSA(pm_mortality_data_df_baseline, pm_mortality_data_df_2020,  "Seattle-Tacoma-Bellevue, WA", figsize = (17, 6))
#plot_multiple_MSA(pm_mortality_data_df_2020, "Seattle-Tacoma-Bellevue, WA", figsize = (20,10))

# %%
plot_multiple_MSA(pm_mortality_data_df_baseline, pm_mortality_data_df_2020, 'New York-Newark-Jersey City, NY-NJ-PA', figsize = (17, 6))

# %%
plot_multiple_MSA(pm_mortality_data_df_baseline, pm_mortality_data_df_2020, 'Chicago-Naperville-Elgin, IL-IN-WI', figsize = (16, 7))

# %%
plot_multiple_MSA(pm_mortality_data_df_baseline, pm_mortality_data_df_2020,  'Minneapolis-St. Paul-Bloomington, MN-WI', figsize = (16, 7))


# %% [markdown]
# ## EJ Index

# %%


# %%
cbsa_counties['CBSA Title'].unique()

# %%
box_whisker(rwc_census_tract_pm25, location = 'San Jose-Sunnyvale-Santa Clara, CA')
decile_plot(rwc_census_tract_pm25, loc = 'San Jose-Sunnyvale-Santa Clara, CA')

# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def box_whisker(df, location = "", concentration_col = "concentrat", index_col = "SPL_EJI", num_iles = 10):

    if location:
        city_index  = df[df['CBSA Title'] == location].index
        city_tracts = df.loc[city_index]
        df = city_tracts
    
    # Create deciles for pollution concentrations
    num_iles = 10
    df['decile'] = pd.qcut(df[concentration_col], q=num_iles, labels=False)

    # create colors
    mean_ej_index = df.groupby('decile')[index_col].median()
    
    # Normalize the mean EJ_index values to get colors
    norm = plt.Normalize(mean_ej_index.min() - 0.2, mean_ej_index.max() + 0.2)
    cmap = plt.get_cmap('Reds')
    colors = [cmap(norm(value)) for value in mean_ej_index]
    
    # Create the box and whisker plot
    plt.figure(figsize=(12, 8))
    sns.boxplot(x='decile', y=index_col, data=df, palette = colors)
    plt.xlabel('Decile of PM2.5 Pollution Concentration')
    plt.ylabel('EJ Index')
    plt.title(f"Box and Whisker Plot of EJ Index for Each Decile of Pollution Concentrations in {location}")
    plt.show()
    
box_whisker(rwc_census_tract_pm25)

# %%
box_whisker(rwc_census_tract_pm25, location = "Seattle-Tacoma-Bellevue, WA")
decile_plot(rwc_census_tract_pm25, loc = "Seattle-Tacoma-Bellevue, WA")

# %%
box_whisker(rwc_census_tract_pm25, location = 'New York-Newark-Jersey City, NY-NJ')
decile_plot(rwc_census_tract_pm25, loc = 'New York-Newark-Jersey City, NY-NJ')

# %%
box_whisker(rwc_census_tract_pm25, location = 'Las Vegas-Henderson-North Las Vegas, NV')
decile_plot(rwc_census_tract_pm25, loc = 'Las Vegas-Henderson-North Las Vegas, NV')

# %%
box_whisker(rwc_census_tract_pm25, location = 'Minneapolis-St. Paul-Bloomington, MN-WI')
decile_plot(rwc_census_tract_pm25, loc = 'Minneapolis-St. Paul-Bloomington, MN-WI')

# %%
box_whisker(rwc_census_tract_pm25, location = 'Philadelphia-Camden-Wilmington, PA-NJ-DE-MD')
decile_plot(rwc_census_tract_pm25, loc = 'Philadelphia-Camden-Wilmington, PA-NJ-DE-MD')

# %%
box_whisker(rwc_census_tract_pm25, location = 'Chicago-Naperville-Elgin, IL-IN')
decile_plot(rwc_census_tract_pm25, loc = 'Chicago-Naperville-Elgin, IL-IN')

# %%
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#
num_iles = 10

# Calculate deciles of pollution concentrations
merged['EJI_decile'] = pd.qcut(merged['SPL_EJI'], q=num_iles, labels=False)

# Calculate the percentage of each race within each decile
EJI_deciles = [0,1,2,3,4,5,6,7,8,9]
    
#Create a bar plot with each bar representing a decile and colored by race percentages
plt.figure(figsize=(10, 6))
bottom = np.zeros(10)



for EJ_decile in EJI_deciles:
    EJ_percentages = []
    for decile in range(num_iles):
        pm25_decile_df = (merged.loc[merged['pm25_decile'] == decile])
        EJ_percentages.append(pm25_decile_df.loc[pm25_decile_df['EJI_decile'] == EJ_decile].shape[0] / pm25_decile_df.shape[0] * 100) 
    
    plt.bar(range(num_iles), EJ_percentages, bottom=bottom, label=EJ_decile)
    bottom += np.array(EJ_percentages)
    
plt.xlabel('Deciles of Pollution Concentrations')
plt.ylabel('Percentage of Population in each EJ index')
plt.title('EJ Index Composition within Deciles of Pollution Concentrations')
plt.legend()
plt.show()

# %% [markdown]
# ## 2020 total to census tract

# %%
### Converting rwc pollutant varaible to census tract

#0. Import RWC data
#import necessary libraries
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
cell_size = 4000
variables = ['PM25_TOT']
baseline = xr.open_dataset("avg_201601.nc")
rwc_2020 = xr.open_dataset("avg_201601_2020_RWC.nc")
no_rwc = xr.open_dataset("no_rwc_avg_201601.nc")

# 1. Import census tract and ethnicity data csv
import pandas as pd
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


# 3. Create gridded pollutant data
import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np

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


# add RWC variables to the grid shapefile of CONUS grid
for var in variables:
    RWC_var = (rwc_2020[var])[0,0,:,:].to_numpy().ravel()
    grid_gdf[var] = RWC_var
    
grid_gdf = grid_gdf.reset_index().rename(columns={'index': 'iD'})


# 5. Create intersections shapefile
merged_gdf = gpd.read_file("../SMOKE_sensitivity_analyses/merged_gdf/merged_gdf.shp")
merged_gdf = merged_gdf.to_crs(grid_gdf.crs) # convert to the same SMOKE CRS projection 
merged_shapes = merged_gdf[['GISJOIN','geometry']]
intersection = gpd.overlay(grid_gdf, merged_shapes, how='intersection')
print(intersection.columns)
#intersection.to_file("intersection_2020.shp")

from shapely.geometry import Polygon

# Define a function to calculate the area of a polygon
def calculate_polygon_area(polygon):
    return polygon.area

for var in variables:
    intersection['fraction'] = intersection.geometry.area/(cell_size ** 2) # area of intersection / divided by area of gridcell
    intersection['emissions_per_polygon'] = intersection[var] * intersection['fraction'] # concentration * area of gridcell
    summed_df = intersection.groupby('GISJOIN')['emissions_per_polygon'].sum().reset_index()
    merged_gdf = merged_gdf.merge(summed_df, on='GISJOIN')
    merged_gdf[var] = merged_gdf['emissions_per_polygon']/merged_gdf.geometry.area * 1000000
    census_tract_gdf = census_tract_gdf.merge(summed_df, on='GISJOIN')

census_tract_gdf.to_file("2020_tot_census_tract_pm25.shp")



