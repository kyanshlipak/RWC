
"""
This script performs data analysis and visualization on Residential Wood Combustion (RWC) 
emission data and PM2.5 concentrations. The analysis includes:
- Loading and processing emission data from different years (2016 and 2020).
- Plotting county-level RWC emissions and PM2.5 concentrations using geospatial data.
- Aggregating emissions data by appliance categories (e.g., woodstoves, hydronic heaters, firelogs, etc.).
- Comparing emissions across different months (January vs. July) for RWC.
- Visualizing the differences in PM2.5 emissions between baseline, RWC, and non-RWC scenarios.
- Converting emissions from grams to tons and handling the data in both pandas DataFrames and xarray datasets.
- Performing statistical analysis and generating output for understanding the contributions of RWC to total PM2.5 emissions.

Libraries Used:
- pandas: Data manipulation and processing.
- geopandas: Geospatial analysis and plotting.
- xarray: Handling multi-dimensional data arrays.
- matplotlib, cartopy: Visualization of geospatial data and plots.
- numpy: Array manipulation and mathematical operations.
- pyproj: Handling coordinate projections.
- shapely: Geometry manipulation for geospatial data.

Data Sources:
- RWC emission data files (2016 and 2020).
- SMOKE processed emissions data for January and July 2016.
- U.S. county shapefile for mapping.
"""


#import necessary libraries
import xarray as xr
from netCDF4 import Dataset, MFDataset, num2date
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import Proj
import geopandas as gpd

# ------------------------------------------------------------------
# Plot and anaysis of NEI RWC inputs (county level resolution)
# from the 2020 version released in 2023.


my_dir = '../SMOKE_sensitivity_analyses/'
rwc_2016 = pd.read_csv(my_dir +'2016ff.csv',skiprows=29)
rwc_2020 = pd.read_csv(my_dir +'rwc_2020.csv',skiprows=18)

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

dataframe, gdf_counties, merged_gdf = plot_counties_matplotlib(rwc_2020)

outdoor_wb = rwc_2020.loc[(rwc_2020['scc'] ==2104008700) & (rwc_2020['emission'] == "PM25-PRI")]['ann_value'].sum()
sum = rwc_2020.loc[(rwc_2020['emission'] == "PM25-PRI")]['ann_value'].sum()
print(f"Percent total emissions from outdoor woodburning {outdoor_wb/sum * 100}")


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

for (appliance_sccs, appliance_string) in appliance_categories[:]:
    agg_data = aggregate(rwc_2020, appliance_sccs)
    agg_data['region_cd'] = agg_data['region_cd'].astype(str)
    _, _, _ = plot_counties_matplotlib(agg_data, title = f"Annual Total {appliance_string} Primary PM₂.₅ Emission")
    

print("Total pm2.5 primary emissions", rwc_2020.loc[rwc_2020['emission'] == "PM25-PRI"]['ann_value'].sum())


emis_state = rwc_2020.loc[rwc_2020['emission'] == "PM25-PRI"].groupby('state')['ann_value'].sum().reset_index().set_index('state').sort_values(by='ann_value', ascending = False)
emis_state['ann_value'] = (emis_state['ann_value'].astype(int))
print("\n\n EMISSIONS by state \n", emis_state)




# ------------------------------------------------------------------
# anaysis of SMOKE processed CONUS January 2016 emissions from RWC

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

rwc_2020_july =  xr.open_dataset('../SMOKE_sensitivity_analyses/rwc_2020_SMOKE_avg_201607.nc')

rwc_2020_july = add_pm25(rwc_2020_july)
seconds_in_a_day = 24 * 60 * 60  # Number of seconds in a day
days_in_month = 31  # Approximate days in a month (adjust if necessary)
seconds_in_a_month = seconds_in_a_day * days_in_month
grams_to_tons = 1 / 1_000_000  # Conversion factor from grams to tons
gps_to_tons = seconds_in_a_month * grams_to_tons
rwc_2020_july['PM25_total_tons'] = rwc_2020_july['PM25_total'] * gps_to_tons


# july to jan comparison
jan_total = 0
for pm_component in ['PAL', 'PCA', 'PEC', 'PFE', 'PH2O', 'PK', 'PMG', 'PMN', 'PMOTHR', 'PNCOM', 'PNH4', 'PNO3', 'POC', 'PSI', 'PSO4', 'PTI']:
	try:
		jan_total += float(rwc_2020_new_surg[pm_component].sum())
		print(pm_component, jan_total)
		#print(pm_component, round(float(rwc_2020_july[pm_component].sum()/rwc_2020_new_surg[pm_component].sum() * 100),4))
	except Exception as e:
		print(e)
          
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


print("\njan vs july emissions of PM primary\n", pd.DataFrame(dict).T)



# ------------------------------------------------------------------
# Plotting of SMOKE processed CONUS January 2016 emissions from RWC
#
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


pmgdf['PM25_total_tons'] = rwc_2020_new_surg['PM25_total_tons'][0,0,:,:].to_numpy().flatten()
PWC_df_emissions = gpd.GeoDataFrame.merge(pmgdf, population_df, left_on=['row', 'col'], right_on=['ROWS', 'COLS'])
pmgdf['PM25_total_tons_july'] = rwc_2020_july['PM25_total_tons'][0,0,:,:].to_numpy().flatten()
PWC_df_emissions = gpd.GeoDataFrame.merge(pmgdf, population_df, left_on=['row', 'col'], right_on=['ROWS', 'COLS'])


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


target_proj = ccrs.LambertConformal(central_longitude=-97, central_latitude=40, standard_parallels=(33, 45))

fig, ax = plt.subplots(1, 1, figsize=(16, 14), subplot_kw={'projection': target_proj})

upper_bound = 100 #PWC_df_emissions['percent_july'].quantile(0.99)
PWC_df_emissions['percent_july'] = PWC_df_emissions['PM25_total_tons_july']/PWC_df_emissions['PM25_total_tons'] * 100
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
