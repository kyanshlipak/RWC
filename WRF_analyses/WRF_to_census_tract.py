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
    merged_gdf[var] = merged_gdf['emissions_per_polygon']/merged_gdf.geometry.area * 1_000_000
    census_tract_gdf = census_tract_gdf.merge(summed_df, on='GISJOIN')

census_tract_gdf.to_file("2020_tot_census_tract_pm25.shp")
