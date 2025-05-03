# Print total tons of RWC emissions for removed species
# showing relative insignificance


#import necessary libraries
import xarray as xr
from netCDF4 import Dataset, MFDataset, num2date
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import Proj
import geopandas as gpd

rwc_2020_jan_og = xr.open_dataset('../SMOKE_sensitivity_analyses/emis_mole_rwc_201601_original_2020_RWC.nc')

seconds_in_a_month = 60 * 60 * 24 * 31

vars_to_delete = ['AACD', 'APIN', 'FACD', 'GLY', 'GLYD', 'ISPD', 'IVOC', 'MGLY', 'NMOG', 'PACD']
for var in vars_to_delete:
	print(rwc_2020_jan_og[var][0,0,:,:].to_numpy().sum() * seconds_in_a_month /1_000_000) # metric tons)
