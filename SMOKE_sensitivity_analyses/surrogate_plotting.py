
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Polygon
import geopandas as gpd
from pyproj import Proj
import numpy as np
import pandas as pd


#counties shapefile
import cartopy.io.shapereader as shpreader
dirr = r""
reader = shpreader.Reader(dirr + 'counties/countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())


detached_housing_df = pd.read_csv(dirr + 'USA_135_NOFILL.txt', delimiter='\t', header=None, skiprows=16)
unit_housing_df = pd.read_csv(dirr + 'USA_136_NOFILL.txt', delimiter='\t', header=None, skiprows=16)
NLCD_low_intensity_df = pd.read_csv(dirr + 'USA_300_NOFILL.txt', delimiter='\t', header=None, skiprows=2)

def pandas_to_ndarray(df, num_rows, num_cols, data = 6, row_col = 2, col_col = 3):

    # Create an empty numpy array
    numpy_array = np.zeros((num_rows, num_cols), dtype=df[data].dtype)

    # Populate the numpy array with values from the DataFrame
    for index, row in df.iterrows():
        numpy_array[row[col_col]-1, row[row_col]-1] = row[data]

    return numpy_array


NLCD_low_intensity_ndarray_surg = pandas_to_ndarray(NLCD_low_intensity_df, num_rows = 1008, num_cols = 1332, data = 4)
detached_housing_ndarray_surg = pandas_to_ndarray(detached_housing_df, num_rows = 1332, num_cols = 1548, data = 4)
unit_housing_ndarray_surg = pandas_to_ndarray(unit_housing_df, num_rows = 1332, num_cols = 1548, data = 4)

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

lon_2016, lat_2016 = map_lat_lon(
    proj_params = {'proj': 'lcc',
                   'lat_1': 33,
                   'lat_2': 45,
                   'lon_0': -97,
                   'lat_0': 40},
    xorig = -2736000, yorig= -2088000,  # Coordinates of the origin
    ncols = 1332, nrows = 1008  # Number of columns and rows
    )

lon_2020, lat_2020 = map_lat_lon(
    proj_params = {'proj': 'lcc',
                   'lat_1': 33,
                   'lat_2': 45,
                   'lon_0': -97,
                   'lat_0': 40},
    xorig = -2952000, yorig = -2772000,  # Coordinates of the origin
    ncols = 1548, nrows = 1332  # Number of columns and rows
    )


counties_gdf = gpd.read_file(dirr + 'counties/tl_2024_us_county.shp')


def plot_gdf_in_county(data, bounds, lat_2016, lon_2016, counties_gdf, county_name, statefp, vs, label, title):
    """
    Plot a GeoDataFrame clipped to a specific county.

    Parameters:
    - data: 2D numpy array of grid data
    - bounds: List of row/column indices for slicing the lat/lon data
    - lat_2016, lon_2016: Latitude and Longitude arrays
    - counties_gdf: GeoDataFrame of county geometries
    - county_name: Name of the county to filter by
    - statefp: State FIPS code of the county
    - vs: [vmin, vmax] for color scale
    - label: Label for the colorbar
    - title: Title for the plot
    """
    lats = lat_2016[ bounds[0][0]:bounds[0][1], bounds[1][0]:bounds[1][1] ]
    lons = lon_2016[ bounds[0][0]:bounds[0][1], bounds[1][0]:bounds[1][1] ]
    
    # Initialize an empty list for polygons and data
    polygons = []
    values = []
    
    # Iterate over the grid cells
    for i in range(lats.shape[0] - 1):  # Iterate rows
        for j in range(lons.shape[1] - 1):  # Iterate columns
            # Get the corner coordinates of the grid cell
            cell_corners = [
                (lons[i, j], lats[i, j]),         # Top-left
                (lons[i, j + 1], lats[i, j + 1]), # Top-right
                (lons[i + 1, j + 1], lats[i + 1, j + 1]), # Bottom-right
                (lons[i + 1, j], lats[i + 1, j]), # Bottom-left
            ]
            # Create a polygon from the corners
            polygons.append(Polygon(cell_corners))
            # Append the data value for the cell
            values.append(data[i, j])
    
    # Create a GeoDataFrame
    gdf = gpd.GeoDataFrame(
        {'data': values},
        geometry=polygons,
        crs="EPSG:4326"  # WGS84 coordinate system
    )

    # Ensure the CRS of counties_gdf matches that of the grid (EPSG:4326)
    counties_gdf = counties_gdf.to_crs("EPSG:4326")

    # Extract the county geometry
    county_geom = counties_gdf.loc[
        (counties_gdf['NAME'] == county_name) & (counties_gdf['STATEFP'] == statefp),
        'geometry'
    ].values[0]

    # Filter and clip the GeoDataFrame
    gdf_filtered = gdf[gdf.intersects(county_geom)]
    gdf_clipped = gpd.clip(gdf_filtered, county_geom)

    # Plot
    plt.figure(figsize=(8, 6))
    ax = plt.subplot(111, projection=ccrs.PlateCarree())
    gdf_clipped.plot(
        column='data', cmap='viridis', legend=False, ax=ax, transform=ccrs.PlateCarree(),
        vmin=vs[0], vmax=vs[1]  # Set min and max values for the colormap
    )
    # Add county boundary
    #ax.add_geometries(
    #    [county_geom], crs=ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=0.7
    #)

    # Add map features
    ax.coastlines(resolution='10m')
    #ax.add_feature(cfeature.BORDERS, linestyle='--', alpha=0.5)
    #ax.add_feature(cfeature.STATES, edgecolor='gray', linewidth=0.5)
    lakes = cfeature.NaturalEarthFeature('physical', 'lakes', scale='10m', edgecolor='b', facecolor=cfeature.COLORS['water'])
    ax.add_feature(lakes, linewidth = 0)

    # Set color scale and labels
    plt.title(title, fontsize = 16)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')

    import matplotlib.colors as mcolors
    norm = mcolors.Normalize(vmin=vs[0], vmax=vs[1])

    cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='viridis', norm=norm), ax=ax, orientation='vertical')
    cbar.ax.tick_params(labelsize=14)  # Set the tick label size
    cbar.set_label(label, fontsize=14)  # Add a label with a specified font size
    plt.show()
    return ax


coords = [[570,600],[865,886]]
data1 = unit_housing_ndarray_surg[coords[0][0]+171: coords[0][1]+171, coords[1][0]+54: coords[1][1]+54] * 100
data2 = NLCD_low_intensity_ndarray_surg[coords[0][0]:coords[0][1], coords[1][0]:coords[1][1]]
data = data1 - data2
data3 = detached_housing_ndarray_surg[coords[0][0]+171: coords[0][1]+171, coords[1][0]+54: coords[1][1]+54] * 100


ax  = plot_gdf_in_county(data = data1, 
                   bounds = coords, 
                   lat_2016 = lat_2016, 
                   lon_2016 = lon_2016, 
                   counties_gdf = counties_gdf, 
                   county_name = "Cook", 
                   statefp = "17", 
                    vs=[0, 3],
                    label="Surrogate Weighting (%)",
                    title="A Cook County Single/Dual Unit Housing Surrogate")


ax  = plot_gdf_in_county(data = data3, 
                   bounds = coords, 
                   lat_2016 = lat_2016, 
                   lon_2016 = lon_2016, 
                   counties_gdf = counties_gdf, 
                   county_name = "Cook", 
                   statefp = "17", 
                    vs=[0, 3],
                    label="Surrogate Weighting (%)",
                    title="B Cook County Detached Housing Surrogate")

