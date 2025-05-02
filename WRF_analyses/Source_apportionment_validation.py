"""
===============================================================================
This script processes gridded PM2.5 concentration data from NetCDF files
to evaluate the contribution of Residential Wood Combustion (RWC) to PM2.5
pollution across the continental U.S and specificlaly compare to source apportionment
studied sites.

The script:
- Loads PM2.5 datasets with and without RWC emissions.
- Constructs a 4 km Lambert Conformal Conic grid and associates it with PM2.5 values.
- Computes RWC-related PM2.5 concentrations and percentage contributions.
- Defines a function to extract RWC PM2.5 values for specific geographic locations
  and their 3×3 neighboring grid cells.
- Applies this function to a list of urban and rural monitoring sites for comparative analysis.

Requirements:
- xarray
- numpy
- pandas
- geopandas
- shapely
- pyproj
- cartopy

Output:
- Console prints of RWC PM2.5 values and percent contributions at key sites.
===============================================================================
"""


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

def check_rwc_with_neighbords(location, lat, lon): 
	row, col = latlon_to_grid(lat, lon)

	# Select the 3×3 neighborhood
	neighbors = pmgdf.loc[
		(pmgdf['row'].between(row - 1, row + 1)) & 
		(pmgdf['col'].between(col - 1, col + 1))
	]

	# Compute the average PM2.5
	print(location, 
        round(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col), 'PM25_RWC_2020'].iloc[0],1), \
		round(pmgdf.loc[(pmgdf['row'] == row) & (pmgdf['col'] == col), 'PM25_RWC_percent_2020'].iloc[0],1) )

	print(round(neighbors['PM25_RWC_2020'].min(), 1), 
		round(neighbors['PM25_RWC_2020'].max(), 1), 
		round(neighbors['PM25_RWC_percent_2020'].min(), 1), 
		round(neighbors['PM25_RWC_percent_2020'].max(), 1))



#Phillips neighborhood
lat, lon = 44.9543561,-93.2691361  # Replace with actual coordinates
check_rwc_with_neighbords("Phillips", lat, lon)

## NYC 
lat, lon = 40.73708, -73.82158
check_rwc_with_neighbords("NYC", lat, lon)

lat, lon  = 33.8128, -112.2381
check_rwc_with_neighbords("Phoenix", lat, lon)


lat, lon  = 33.7750752, -84.417395
check_rwc_with_neighbords("jefferson st. atlanta", lat, lon)

lat, lon = 33.5516709,-86.8399714
# north birmingham
check_rwc_with_neighbords("north birmingham", lat, lon)


lat, lon = 32.9649163,-87.1426895
check_rwc_with_neighbords("centreville", lat, lon)

lat, lon = 48.3000793,-124.636648
check_rwc_with_neighbords("ccheaka peak ", lat, lon)


lat, lon = 47.5857665,-122.3331929
check_rwc_with_neighbords("Seattle ", lat, lon)


# brownseville
lat, lon = 26.06973, -97.16222
check_rwc_with_neighbords("brownseville", lat, lon)

#Los angeles
lat, lon = 34.067186, -118.227172
check_rwc_with_neighbords("Los Angeles", lat, lon)


# roubidoux
lat, lon = 34.002296, -117.417796
check_rwc_with_neighbords("roubidoux", lat, lon)


# Detroit
lat, lon = 42.228611, -83.20833
check_rwc_with_neighbords("Detroit", lat, lon)


# Chicago
lat, lon = 41.7514, -87.713488
check_rwc_with_neighbords("Chicago", lat, lon)


# Chicago
lat, lon = 41.7514, -87.713488
check_rwc_with_neighbords("Chicago", lat, lon)


### NYC STUDY
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
	check_rwc_with_neighbords(key, lat, lon)


### SOUTHEAST STUDY
dictionary = {
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

for key, value in dictionary.items():
	lat, lon = value[0], value[1]
	check_rwc_with_neighbords(key, lat, lon)


# NORTHWEST STUDY
dictionary = {
    #"Fairbanks": [64.8407, -147.7225],
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

for key, value in dictionary.items():
	lat, lon = value[0], value[1]
	check_rwc_with_neighbords(key, lat, lon)
