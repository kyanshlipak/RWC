#import necessary libraries
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from wrf import getvar, interplevel, to_np, latlon_coords
import cartopy.crs as crs
import cartopy.feature as cfeature
import xarray as xr

# Variables of interest
runname = 'output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852'
dirToWRF = f'/projects/b1045/wrf-cmaq/output/CONUS4K/{runname}/'
filenames_d02 = sorted(glob.glob(os.path.join(dirToWRF, "wrfout_*")))
variables_of_interest = ['T2', 'Q2', 'U10', 'V10']

datasets =[]

# Load and combine all datasets along a new dimension (e.g., 'day')
for file in filenames_d02:
    dataset = xr.open_dataset(file)
    selected_data = dataset[variables_of_interest]
    datasets.append(selected_data)

combined_dataset = xr.concat(datasets, dim='Time')

average_dataset = combined_dataset.mean(dim="Time", keep_attrs=True)

# Manually copy over non-Time coordinates
for coord in combined_dataset.coords:
    if coord != "Time":  # Exclude the dimension being reduced
        average_dataset.coords[coord] = combined_dataset.coords[coord]

average_dataset = average_dataset.drop_dims("Time")

#average_dataset.to_netcdf(output_file)
average_dataset = average_dataset.assign_coords({
    "XLONG": combined_dataset["XLONG"].isel(Time=0),
    "XLAT": combined_dataset["XLAT"].isel(Time=0)
})

# Compute the average across the 'time' dimension
output_file = "baseline_wrf_avg.nc"

# Save the averaged dataset to a new NetCDF file
average_dataset.to_netcdf(output_file)

print(f"Averaged file saved as {output_file}")
