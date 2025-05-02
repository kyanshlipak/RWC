"""
Grid Boundary Visualization Script
----------------------------------

This script generates and visualizes the boundary of the CONUS 4 km domain using a Lambert Conformal Conic (LCC) projection.
It performs the following steps:

1. Defines the LCC projection and grid geometry.
2. Computes the projected coordinates of the grid boundary.
3. Transforms those coordinates to latitude/longitude.
4. Visualizes the boundary on a map using Cartopy and Matplotlib.
5. Optionally overlays a background raster image (e.g., Natural Earth terrain).

Dependencies:
- numpy
- pyproj
- matplotlib
- cartopy
- rasterio

Make sure "NE2_50M_SR.tif" is present in the working directory for the raster background to display properly.
"""

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
