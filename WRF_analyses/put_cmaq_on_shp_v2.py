
# Updated 2022-06-07: Make box was projecting cells slightly off center. Fixed by using MCIP files to create shapes. 
#---------------------------------------------------------------------
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import netCDF4
from geopandas.tools import sjoin
from shapely.geometry import Point, shape, Polygon
import geopandas as gpd
#---------------------------------------------------------------------
print(gpd.__version__)
#exit()

# User defined functions
#---------------------------------------------------------------------
def make_box(X0,Y0,X1,Y1): # OLD 
    # Create the lat lon values given corner indices
    # For creating outer edge values of the array
    #       L
    #       lu----ru
    # H     |   c  |
    #       |      |
    #       ll----rl
    #
    L = np.abs(X0 - X1)
    H = np.abs(Y0 - Y1)
    lo_lu, la_lu = X0-L/2,Y0+H/2
    lo_ru, la_ru = X0+L/2,Y0+H/2
    lo_rl, la_rl = X0+L/2,Y0-H/2
    lo_ll, la_ll = X0-L/2,Y0-H/2
    return lo_lu, la_lu, lo_ru, la_ru,lo_rl, la_rl,lo_ll, la_ll

#---------------------------------------------------------------------
def make_shp(datas,datas_headers,lon3,lat3,tract_shapes):
    # 
    # datas = list of data in nvar x 1 array
    # datas_headers = the name of the variable that will go into shp file
    #
    # Begin:
    # create dictionary using air quality concs for future geodataframe
    d = {}
    for i in range(len(datas)):
        d[datas_headers[i]]=datas[i]
    #
    # pull locations and put into shapely point so we can do geographical transformations
    pointline = [Point(lon3.ravel()[i],lat3.ravel()[i]) for i in range(len(lon3.ravel()))]
    # add in geometry into dictionary
    d['geometry']=pointline
    #
    # make geodataframe using dictionary and transform lon/lats projections into the same geometry as the tracts
    gdf = gpd.GeoDataFrame(d, crs="EPSG:4326")
    gdf = gdf.to_crs(tract_shapes.crs)
    # pull out these lat lon points that have the tract crs
    x,y = gdf['geometry'].x,gdf['geometry'].y
    #
    # Start to make little polygons to surround the centroids
    # Create outer limits
    #       L
    #       lu----ru
    # H     |   c  |
    #       |      |
    #       ll----rl
    #
    # set up empty shape to hold the future polygons
    shape=np.zeros(lon.shape).tolist()
    # go through the indices in the lat lon
    for i in range(len(lon)):
            for j in range(len(lon[0])):
                X0,Y0 = lat[i][j],lon[i][j]
                # get corners
                lo_lu, la_lu = latd[i][j],lond[i][j]
                lo_ru, la_ru = latd[i+1][j],lond[i+1][j]
                lo_rl, la_rl = latd[i+1][j+1],lond[i+1][j+1]
                lo_ll, la_ll = latd[i][j+1],lond[i][j+1]
                # create points of corners
                pointList = [Point(lo_lu, la_lu),Point(lo_ru, la_ru), Point(lo_rl, la_rl), Point(lo_ll, la_ll),Point(lo_lu, la_lu)]
                # create a polygon from the corners
                poly = Polygon([[p.y, p.x] for p in pointList])
                # put the shape of the pixel into the shape list
                shape[i][j] = poly 
    # put list of shapes into final gdf, needs to be 1D
    shape = np.array(shape).ravel()
    d['geometry'] = shape
    d['lat'],d['lon'] = lat.ravel(),lon.ravel()
    d['lat_m'],d['lon_m'] = x,y
    #
    # make geodataframe
    gdf = gpd.GeoDataFrame(d, crs="EPSG:4326")
    gdf = gdf.to_crs(tract_shapes.crs)
    gdf.to_file('tmp.shp')
    gdf['geom_saved']=gdf['geometry']
    print('create polygon shapes')
    #
    # Uncomment/edit this area if things are too big 
    # Like if there are too many polygons -- more efficient to go through states
    #tractjoinlist = []
    #ntracks = []
    #combien tracts with the geodataframe -- there are many repeating variables
    #for tractstate in ['17','18','55','26']: # selecting only midwest states, can change to all
    #       tractsubset = tract_shapes[tract_shapes['STATEFP']==tractstate]
    #       for tractcounty in tractsubset.COUNTYFP.unique():
    #           tractjoin = tractsubset[tractsubset['COUNTYFP']==tractcounty]
    #           tractjoin = gpd.sjoin(tractjoin,gdf,how="inner",op='intersects')
    #           if len(tractjoin)!=0: 
    #               tractjoinlist.append(tractjoin)
    #               ntracks.append(tractjoin.NAME.unique().shape)
    #
    #tractjoinfinal = gpd.GeoDataFrame(tractjoinlist[0])
    #for i in range(len(tractjoinlist)):
    #       tractjoinfinal = tractjoinfinal.append(tractjoinlist[i])
    #
    # need to combine AQ by distances from centroids to remove repeating variables
    v=[];geo = [];gisjoin = [];geoid=[]
    #
    tractjoinfinal =  gpd.sjoin(tract_shapes,gdf,how="inner",op='intersects')
    grouped_tracts = tractjoinfinal.groupby('GISJOIN')
    #
    print('joined tracts to polygon shapes')
    # go through each join code and do computations to get final weight AQ
    for group in grouped_tracts.groups:
            # pull the section of the tract with repeating values
            tmp = grouped_tracts.get_group(group)
            # save geometry data from the census shape
            geo.append(tmp.iloc[0]['geometry'])
            # save the gidjoin coide
            gisjoin.append(tmp.iloc[0]['GISJOIN'])
            geoid.append(tmp.iloc[0]['GEOID'])
            # Now process data for weighting
            # get the area of intersection of cells inside of the shape
            area = tmp.geometry.intersection(tmp.geom_saved).area
            # get distance between tract and cmaq centroids
            # now weight each air quality constiteunt to get average depending on distances
            for dh in datas_headers:
                v.append(np.sum(1/np.array(area)*tmp[dh])/np.sum(1/np.array(area)))
    #
    v = np.array(v)
    v2 = v.reshape(int(v.shape[0]/len(datas_headers)),len(datas_headers))
    v2 = gpd.GeoDataFrame(v2)
    v2.columns = datas_headers
    v2['GISJOIN'] = gisjoin
    v2['geometry'] = geo
    v2.crs = tractjoinfinal.crs
    #
    print('outputting to file')
    v2.to_file(fdirout+fnameout)
    print('complete')
#---------------------------------------------------------------------



#---------------------------------------------------------------------
# Start User Input

# pull in lat lon data -- doesn't matter dates
d03 = Dataset('/projects/b1045/wrf-cmaq/output/Chicago_LADCO/mcip/PXLSM/ChicagoLADCO_d03/GRIDCRO2D_Chicago_LADCO_2018-08-16.nc')
lat3,lon3 = np.array(d03['LAT'][0][0]),np.array(d03['LON'][0][0])
lat,lon = lat3,lon3
print('lat',lat)
print('lon',lon)

# corners from MCIP file -- doesn't matter dates
f = '/projects/b1045/wrf-cmaq/output/Chicago_LADCO/mcip/PXLSM/ChicagoLADCO_d03/GRIDDOT2D_Chicago_LADCO_2018-08-16.nc'
lond,latd = np.array(Dataset(f)['LOND'][0][0]),np.array(Dataset(f)['LATD'][0][0])
print('latd',latd)
print('lond',lond)
#exit()

# get shape files -- this is the target shape you want the AQ data to be averaged on
tract_shapes = gpd.GeoDataFrame.from_file('/projects/b1045/montgomery/nhgis0006_shape/US_tract_2018.shp')
#tract_shapes = tract_shapes.to_crs("EPSG:4326") #kills

'''
# pull in the data that you want to go into the final shape
d = '/projects/b1045/montgomery/ForMaxime/nc'
no2_base= Dataset(d+'/no2_avg_base_annual.nc')['NO2'][0][0]
mo3_base = Dataset(d+'/dailymaxozone_annual_base.nc')['O3'][0][0]
pm_base = Dataset(d+'/pm25_avg_base_annual.nc')['PM25_TOT'][0][0]
o3_base = Dataset(d+'/o3_avg_base_annual.nc')['O3'][0][0]
'''

# pull in the data that you want to go into the final shape
#STEP1 Load BASELINE 
d_1 = '/projects/b1045/scamilleri'
#d_3 = '/projects/b1045/wrf-cmaq/output/Chicago_LADCO/output_BASE_FINAL_1.33km_sf_rrtmg_5_8_1_v3852/postprocess' #SUMMER
#d_3 = '/projects/b1045/wrf-cmaq/output/Chicago_LADCO/output_BASE_FINAL_fall_1.33km_sf_rrtmg_5_8_1_v3852/postprocess' #FALL
#d_3 = '/projects/b1045/wrf-cmaq/output/Chicago_LADCO/output_BASE_FINAL_wint_1.33km_sf_rrtmg_5_8_1_v3852/postprocess' #WINT
#d_3 = '/projects/b1045/wrf-cmaq/output/Chicago_LADCO/output_BASE_FINAL_spring_1.33km_sf_rrtmg_5_8_1_v3852/postprocess' #SPRING
#d_3 = '/projects/b1045/vlang/WH_shapefiles/' #Annual WHbase/75WH EDF Victoria
#d_4 = '/projects/b1045/vlang/ACT_shapefiles/' #Annual ACT EDF Victoria

no2_base= Dataset(d_1+'/Baseline_all.nc')['NO2'] #all Time 744 lay 0
mo3_base = Dataset(d_1+'/Baseline_Daily_MDA8O3.nc')['O3'] #orig 31,1,288,315
pm_base = Dataset(d_1+'/Baseline_all.nc')['PM25_TOT'] #all TIME 744
o3_base = Dataset(d_1+'/Baseline_all.nc')['O3'] #all TIME 744
pmEC_base = Dataset(d_1+'/Baseline_all.nc')['PM25_EC']
#no2_base= Dataset(d_1+'/Base_SumWint.nc')['NO2'] #all Time 744 lay 0
#mo3_base = Dataset(d_1+'/Base_SumWint_MDA8O3_all.nc')['O3'] #orig 31,1,288,315
#pm_base = Dataset(d_1+'/Base_SumWint.nc')['PM25_TOT'] #all TIME 744
#o3_base = Dataset(d_1+'/Base_SumWint.nc')['O3'] #all TIME 744
#pmEC_base = Dataset(d_1+'/Base_SumWint.nc')['PM25_EC']
#print('no2_b',no2_base)
#SUMMER
#no2_base= Dataset(d_3+'/all.nc')['NO2'] #all Time 744 lay 0
#mo3_base = Dataset(d_3+'/dailymaxozone_201808.nc')['O3'] #orig 31,1,288,315
#pm_base = Dataset(d_3+'/all.nc')['PM25_TOT'] #all TIME 744
#o3_base = Dataset(d_3+'/all.nc')['O3'] #all TIME 744
#pmEC_base = Dataset(d_3+'/all.nc')['PM25_EC']
#FALL
#no2_base= Dataset(d_3+'/all.nc')['NO2'] #all Time 744 lay 0
#mo3_base = Dataset(d_3+'/dailymaxozone_201810.nc')['O3'] #orig 31,1,288,315
#pm_base = Dataset(d_3+'/all.nc')['PM25_TOT'] #all TIME 744
#o3_base = Dataset(d_3+'/all.nc')['O3'] #all TIME 744
#pmEC_base = Dataset(d_3+'/all.nc')['PM25_EC']
#WINT
'''no2_base= Dataset(d_3+'/all.nc')['NO2'] #all Time 744 lay 0
mo3_base = Dataset(d_3+'/dailymaxozone_201901.nc')['O3'] #orig 31,1,288,315
pm_base = Dataset(d_3+'/all.nc')['PM25_TOT'] #all TIME 744
o3_base = Dataset(d_3+'/all.nc')['O3'] #all TIME 744
pmEC_base = Dataset(d_3+'/all.nc')['PM25_EC']'''
#SPRING
#no2_base= Dataset(d_3+'/all.nc')['NO2'] #all Time 744 lay 0
#mo3_base = Dataset(d_3+'/dailymaxozone_201904.nc')['O3'] #orig 31,1,288,315
#pm_base = Dataset(d_3+'/all.nc')['PM25_TOT'] #all TIME 744
#o3_base = Dataset(d_3+'/all.nc')['O3'] #all TIME 744
#pmEC_base = Dataset(d_3+'/all.nc')['PM25_EC']
#######WAREHOUSE
'''no2_base= Dataset(d_3+'/whbase_ncfiles/AnnualAvg_WHbase_FromAllTimeSteps.nc')['NO2'] #ANNUAL mean #lay 0
mo3_base = Dataset(d_3+'/whbase_ncfiles/AnnualAvgMDA8O3_WHbase_FromAllTimeSteps.nc')['O3'] # ANNUAL mean #orig 31,1,288,315
pm_base = Dataset(d_3+'/whbase_ncfiles/AnnualAvg_WHbase_FromAllTimeSteps.nc')['PM25_TOT'] # ANNUAL mean #all TIME 744
o3_base = Dataset(d_3+'/whbase_ncfiles/AnnualAvg_WHbase_FromAllTimeSteps.nc')['O3'] #ANNUAL MEA
pmEC_base = Dataset(d_3+'/whbase_ncfiles/AnnualAvg_WHbase_FromAllTimeSteps.nc')['PM25_EC']'''

#d_1 = '/projects/b1045/wrf-cmaq/input/emis/Chicago_LADCO/ChicagoLADCO_d03_spring'
#base_no2= Dataset(d_1+'/AvgBase_emis_mole_all_Seasons.nc')['NO2']
#base_no= Dataset(d_1+'/AvgBase_emis_mole_all_Seasons.nc')['NO']
#base_nox = base_no2[0,0,:,:] + base_no[0,0,:,:]


#STEP2 Load eHDV
#d_2 = '/projects/b1045/wrf-cmaq/output/Chicago_LADCO/output_Maxime_eHDV_1.33km_sf_rrtmg_5_8_1_v3852/postprocess/'
#d_2 = '/projects/b1045/scamilleri/'
#d_4 = '/projects/b1045/wrf-cmaq/output/Chicago_LADCO/output_Maxime_allTransport_1.33km_sf_rrtmg_5_8_1_v3852/postprocess' #SUMMER
#d_4 = '/projects/b1045/wrf-cmaq/output/Chicago_LADCO/output_Maxime_allTransport_fall_1.33km_sf_rrtmg_5_8_1_v3852/postprocess' #Fall
#d_4 = '/projects/b1045/wrf-cmaq/output/Chicago_LADCO/output_Maxime_allTransport_wint_1.33km_sf_rrtmg_5_8_1_v3852/postprocess' #Wint
#d_4 = '/projects/b1045/wrf-cmaq/output/Chicago_LADCO/output_Maxime_allTransport_spring_1.33km_sf_rrtmg_5_8_1_v3852/postprocess' #Spring


#no2_eHDV= Dataset(d_2+'/eHDV_all.nc')['NO2'] #all TIME 2952
#mo3_eHDV = Dataset(d_2+'/eHDV_Daily_MDA8O3.nc')['O3'] #TSTEP is 123
#pm_eHDV = Dataset(d_2+'/eHDV_all.nc')['PM25_TOT']
#o3_eHDV = Dataset(d_2+'/eHDV_all.nc')['O3']
#pmEC_eHDV = Dataset(d_2+'/eHDV_all.nc')['PM25_EC']
#no2_eHDV= Dataset(d_2+'/eHDVnoEGU_SumWint.nc')['NO2'] #all TIME 2952
#mo3_eHDV = Dataset(d_2+'/eHDVnoEGU_SumWint_MDA8O3_all.nc')['O3'] #TSTEP is 123
#pm_eHDV = Dataset(d_2+'/eHDVnoEGU_SumWint.nc')['PM25_TOT']
#o3_eHDV = Dataset(d_2+'/eHDVnoEGU_SumWint.nc')['O3']
#pmEC_eHDV = Dataset(d_2+'/eHDVnoEGU_SumWint.nc')['PM25_EC']
#no2_eHDV= Dataset(d_2+'/eHDV_SumWint.nc')['NO2'] #all TIME 2952
#mo3_eHDV = Dataset(d_2+'/eHDV_SumWint_MDA8O3_all.nc')['O3'] #TSTEP is 123
#pm_eHDV = Dataset(d_2+'/eHDV_SumWint.nc')['PM25_TOT']
#o3_eHDV = Dataset(d_2+'/eHDV_SumWint.nc')['O3']
#pmEC_eHDV = Dataset(d_2+'/eHDV_SumWint.nc')['PM25_EC']
#no2_eHDV= Dataset(d_2+'/eHDVnoEGU_all.nc')['NO2'] #all TIME 2952
#mo3_eHDV = Dataset(d_2+'/eHDVnoEGU_all_MDA8O3.nc')['O3'] #TSTEP is 123
#pm_eHDV = Dataset(d_2+'/eHDVnoEGU_all.nc')['PM25_TOT']
#o3_eHDV = Dataset(d_2+'/eHDVnoEGU_all.nc')['O3']
#pmEC_eHDV = Dataset(d_2+'/eHDVnoEGU_all.nc')['PM25_EC']
'''no2_eHDV= Dataset(d_2+'/allTrans_allTime_NewWint.nc')['NO2'] #all TIME 2952
mo3_eHDV = Dataset(d_2+'/allTrans_allTime_MDA8O3_NewWint.nc')['O3'] #TSTEP is 123
pm_eHDV = Dataset(d_2+'/allTrans_allTime_NewWint.nc')['PM25_TOT']
o3_eHDV = Dataset(d_2+'/allTrans_allTime_NewWint.nc')['O3']
pmEC_eHDV = Dataset(d_2+'/allTrans_allTime_NewWint.nc')['PM25_EC']'''
#no2_eHDV= Dataset(d_2+'/allTransnoEGU_all.nc')['NO2'] #all TIME 2952
#mo3_eHDV = Dataset(d_2+'/eATnoEGU_all_MDA8O3.nc')['O3'] #TSTEP is 123
#pm_eHDV = Dataset(d_2+'/allTransnoEGU_all.nc')['PM25_TOT']
#o3_eHDV = Dataset(d_2+'/allTransnoEGU_all.nc')['O3']
#pmEC_eHDV = Dataset(d_2+'/allTransnoEGU_all.nc')['PM25_EC']
'''no2_eHDV= Dataset(d_2+'/eHDV100p_noEGUs_all.nc')['NO2'] #all TIME 2952
mo3_eHDV = Dataset(d_2+'/eHDV100p_noEGUs_allMDA8O3.nc')['O3'] #TSTEP is 123
pm_eHDV = Dataset(d_2+'/eHDV100p_noEGUs_all.nc')['PM25_TOT']
o3_eHDV = Dataset(d_2+'/eHDV100p_noEGUs_all.nc')['O3']
pmEC_eHDV = Dataset(d_2+'/eHDV100p_noEGUs_all.nc')['PM25_EC']'''
#SUMMER
'''no2_eHDV= Dataset(d_4+'/all_201808.nc')['NO2'] #all TIME 2952
mo3_eHDV = Dataset(d_4+'/dailymaxozone_201808.nc')['O3'] #TSTEP is 123
pm_eHDV = Dataset(d_4+'/all_201808.nc')['PM25_TOT']
o3_eHDV = Dataset(d_4+'/all_201808.nc')['O3']
pmEC_eHDV = Dataset(d_4+'/all_201808.nc')['PM25_EC']'''
#FALL
'''no2_eHDV= Dataset(d_4+'/all_201810.nc')['NO2'] #all TIME 2952
mo3_eHDV = Dataset(d_4+'/dailymaxozone_201810.nc')['O3'] #TSTEP is 123
pm_eHDV = Dataset(d_4+'/all_201810.nc')['PM25_TOT']
o3_eHDV = Dataset(d_4+'/all_201810.nc')['O3']
pmEC_eHDV = Dataset(d_4+'/all_201810.nc')['PM25_EC']'''
#Wint
'''no2_eHDV= Dataset(d_4+'/all_201901.nc')['NO2'] #all TIME 2952
mo3_eHDV = Dataset(d_4+'/dailymaxozone_201901.nc')['O3'] #TSTEP is 123
pm_eHDV = Dataset(d_4+'/all_201901.nc')['PM25_TOT']
o3_eHDV = Dataset(d_4+'/all_201901.nc')['O3']
pmEC_eHDV = Dataset(d_4+'/all_201901.nc')['PM25_EC']'''
#SPRING
'''no2_eHDV= Dataset(d_4+'/all_201904.nc')['NO2'] #all TIME 2952
mo3_eHDV = Dataset(d_4+'/dailymaxozone_201904.nc')['O3'] #TSTEP is 123
pm_eHDV = Dataset(d_4+'/all_201904.nc')['PM25_TOT']
o3_eHDV = Dataset(d_4+'/all_201904.nc')['O3']
pmEC_eHDV = Dataset(d_4+'/all_201904.nc')['PM25_EC']'''
#######WAREHOUSE
'''no2_75WH= Dataset(d_3+'/75wh_ncfiles/AnnualAvg_75WH_FromAllTimeSteps.nc')['NO2'] #ANNUAL mean #lay 0
mo3_75WH = Dataset(d_3+'/75wh_ncfiles/AnnualAvgMDA8O3_75WH_FromAllTimeSteps.nc')['O3'] # ANNUAL mean #orig 31,1,288,315
pm_75WH = Dataset(d_3+'/75wh_ncfiles/AnnualAvg_75WH_FromAllTimeSteps.nc')['PM25_TOT'] # ANNUAL mean #all TIME 744
o3_75WH = Dataset(d_3+'/75wh_ncfiles/AnnualAvg_75WH_FromAllTimeSteps.nc')['O3'] #ANNUAL MEA
pmEC_75WH = Dataset(d_3+'/75wh_ncfiles/AnnualAvg_75WH_FromAllTimeSteps.nc')['PM25_EC']'''
########eLDV
no2_eLDV= Dataset(d_1+'/eLDV_allTime.nc')['NO2'] #all TIME 2952
mo3_eLDV = Dataset(d_1+'/eLDV_MDA8O3_allTime.nc')['O3'] #TSTEP is 123
pm_eLDV = Dataset(d_1+'/eLDV_allTime.nc')['PM25_TOT']
o3_eLDV = Dataset(d_1+'/eLDV_allTime.nc')['O3']
pmEC_eLDV = Dataset(d_1+'/eLDV_allTime.nc')['PM25_EC']

#print(no2_eHDV.shape)
#print(no2_eHDV[:][0].shape)
#print(no2_75WH.shape) #(1, 1, 288, 315)
#print(mo3_75WH.shape) #(1, 1, 288, 315)
#print(no2_base.shape) #(1, 1, 288, 315)
#print(mo3_base.shape) #(1, 1, 288, 315)
#exit()

#######ACT
'''no2_ACT= Dataset(d_4+'/ACT_ncfiles/AnnualAvg_ACT_FromAllTimeSteps.nc')['NO2'] #ANNUAL mean #lay 0
mo3_ACT = Dataset(d_4+'/ACT_ncfiles/AnnualAvgMDA8O3_ACT_FromAllTimeSteps.nc')['O3'] # ANNUAL mean #orig 31,1,288,315
pm_ACT = Dataset(d_4+'/ACT_ncfiles/AnnualAvg_ACT_FromAllTimeSteps.nc')['PM25_TOT'] # ANNUAL mean #all TIME 744
o3_ACT = Dataset(d_4+'/ACT_ncfiles/AnnualAvg_ACT_FromAllTimeSteps.nc')['O3'] #ANNUAL MEA
pmEC_ACT = Dataset(d_4+'/ACT_ncfiles/AnnualAvg_ACT_FromAllTimeSteps.nc')['PM25_EC']'''

#STEP3 TIME AVG
no2_base_avg=np.average(no2_base,axis=0) #avg time axis
pm_base_avg=np.average(pm_base,axis=0)
o3_base_avg=np.average(o3_base,axis=0)
mo3_base_avg=np.average(mo3_base,axis=0)
pmEC_base_avg=np.average(pmEC_base,axis=0)

'''no2_eHDV_avg=np.average(no2_eHDV,axis=0) #avg time axis
pm_eHDV_avg=np.average(pm_eHDV,axis=0)
o3_eHDV_avg=np.average(o3_eHDV,axis=0)
mo3_eHDV_avg=np.average(mo3_eHDV,axis=0)
pmEC_eHDV_avg=np.average(pmEC_eHDV,axis=0)'''

no2_eLDV_avg=np.average(no2_eLDV,axis=0) #avg time axis
pm_eLDV_avg=np.average(pm_eLDV,axis=0)
o3_eLDV_avg=np.average(o3_eLDV,axis=0)
mo3_eLDV_avg=np.average(mo3_eLDV,axis=0)
pmEC_eLDV_avg=np.average(pmEC_eLDV,axis=0)

'''
print(no2_base_avg.shape)
print(pm_base_avg.shape)
print(o3_base_avg.shape)
print(mo3_base_avg.shape)'''
#exit()

#STEP4 DIFF
'''no2_Diff=no2_eHDV_avg[0]-no2_base_avg[0]
mo3_Diff=mo3_eHDV_avg[0]-mo3_base_avg[0]
pm_Diff=pm_eHDV_avg[0]-pm_base_avg[0]
o3_Diff=o3_eHDV_avg[0]-o3_base_avg[0]
pmEC_Diff=pmEC_eHDV_avg[0]-pmEC_base_avg[0]'''
####WAREHOUSE
'''no2_Diff=no2_75WH[0] - no2_base[0]
print(no2_Diff.shape) #(1, 288, 315)
#exit()
mo3_Diff=mo3_75WH[0] - mo3_base[0] #mo3_eHDV_avg[0]-mo3_base_avg[0]
pm_Diff=pm_75WH[0] - pm_base[0] #pm_eHDV_avg[0]-pm_base_avg[0]
o3_Diff=o3_75WH[0] - o3_base[0] #o3_eHDV_avg[0]-o3_base_avg[0]
pmEC_Diff=pmEC_75WH[0] - pmEC_base[0] #pmEC_eHDV_avg[0]-pmEC_base_avg[0]'''

####ACT-75%WH
'''no2_ACTminus75WH=no2_ACT[0] - no2_75WH[0]
print(no2_ACTminus75WH.shape) #(1, 288, 315)
#exit()
mo3_ACTminus75WH=mo3_ACT[0] - mo3_75WH[0] #mo3_eHDV_avg[0]-mo3_base_avg[0]
pm_ACTminus75WH=pm_ACT[0] - pm_75WH[0] #pm_eHDV_avg[0]-pm_base_avg[0]
o3_ACTminus75WH=o3_ACT[0] - o3_75WH[0] #o3_eHDV_avg[0]-o3_base_avg[0]
pmEC_ACTminus75WH=pmEC_ACT[0] - pmEC_75WH[0] #pmEC_eHDV_avg[0]-pmEC_base_avg[0]

print('no2_ACTminus75WH',no2_ACTminus75WH.shape,np.max(no2_ACTminus75WH),np.min(no2_ACTminus75WH))
print('pm_ACTminus75WH',pm_ACTminus75WH.shape,np.max(pm_ACTminus75WH),np.min(pm_ACTminus75WH))
print('o3_ACTminus75WH',o3_ACTminus75WH.shape,np.max(o3_ACTminus75WH),np.min(o3_ACTminus75WH))
print('mda8o3_ACTminus75WH',mo3_ACTminus75WH.shape,np.max(mo3_ACTminus75WH),np.min(mo3_ACTminus75WH))
print('pmEC_ACTminus75WH',pmEC_ACTminus75WH.shape,np.max(pmEC_ACTminus75WH),np.min(pmEC_ACTminus75WH))'''
#exit()

no2_Diff=no2_eLDV_avg[0]-no2_base_avg[0]
mo3_Diff=mo3_eLDV_avg[0]-mo3_base_avg[0]
pm_Diff=pm_eLDV_avg[0]-pm_base_avg[0]
o3_Diff=o3_eLDV_avg[0]-o3_base_avg[0]
pmEC_Diff=pmEC_eLDV_avg[0]-pmEC_base_avg[0]

'''no2_base_av = no2_base_avg[0]
pm_base_av = pm_base_avg[0]
o3_base_av = o3_base_avg[0]
mo3_base_av = mo3_base_avg[0]
pmEC_base_av = pmEC_base_avg[0]'''
'''no2_base_av = no2_base[0,0] #was (1,1,288,315)
#print(no2_base_av.shape) #(288, 315)
#exit()
pm_base_av = pm_base[0,0]
o3_base_av = o3_base[0,0]
mo3_base_av = mo3_base[0,0]
pmEC_base_av = pmEC_base[0,0]

no2_75WH_av = no2_75WH[0,0]
pm_75WH_av = pm_75WH[0,0]
o3_75WH_av = o3_75WH[0,0]
mo3_75WH_av = mo3_75WH[0,0]
pmEC_75WH_av = pmEC_75WH[0,0]

no2_ACT_av = no2_ACT[0,0]
pm_ACT_av = pm_ACT[0,0]
o3_ACT_av = o3_ACT[0,0]
mo3_ACT_av = mo3_ACT[0,0]
pmEC_ACT_av = pmEC_ACT[0,0]

no2_ACTminus75WH_av = no2_ACTminus75WH[0]
pm_ACTminus75WH_av = pm_ACTminus75WH[0]
o3_ACTminus75WH_av = o3_ACTminus75WH[0]
mo3_ACTminus75WH_av = mo3_ACTminus75WH[0]
pmEC_ACTminus75WH_av = pmEC_ACTminus75WH[0] '''

# put all the data you want to put into the shape into 1 list and the header that will go into the shape
#datas = [no2_base.ravel(),o3_base.ravel(),pm_base.ravel(),mo3_base.ravel()]
#datas_headers= ['no2','o3','pm25','mdao3']
datas = [no2_base_avg.ravel(),o3_base_avg.ravel(),pm_base_avg.ravel(),mo3_base_avg.ravel(),pmEC_base_avg.ravel(),no2_Diff.ravel(),mo3_Diff.ravel(),pm_Diff.ravel(),o3_Diff.ravel(),pmEC_Diff.ravel()]
#datas = [no2_base_av.ravel(),o3_base_av.ravel(),pm_base_av.ravel(),mo3_base_av.ravel(),pmEC_base_av.ravel(),no2_75WH_av.ravel(),o3_75WH_av.ravel(),pm_75WH_av.ravel(),mo3_75WH_av.ravel(),pmEC_75WH_av.ravel(),no2_diff_av.ravel(),mo3_diff_av.ravel(),pm_diff_av.ravel(),o3_diff_av.ravel(),pmEC_diff_av.ravel()]
#datas = [base_nox.ravel()]
#datas = [no2_75WH_av.ravel(),o3_75WH_av.ravel(),pm_75WH_av.ravel(),mo3_75WH_av.ravel(),pmEC_75WH_av.ravel(),no2_ACT_av.ravel(),o3_ACT_av.ravel(),pm_ACT_av.ravel(),mo3_ACT_av.ravel(),pmEC_ACT_av.ravel(),no2_ACTminus75WH_av.ravel(),mo3_ACTminus75WH_av.ravel(),pm_ACTminus75WH_av.ravel(),o3_ACTminus75WH_av.ravel(),pmEC_ACTminus75WH_av.ravel()]


print(datas)
datas_headers= ['no2_base','o3_base','pm25_base','mo3_base','pEC_base','no2_diff','mdao3_diff','pm25_diff','o3_diff','pEC_diff']
#datas_headers= ['no2_75WH','o3_75WH','pm_75WH','mo3_75WH','pmEC_75WH','no2_ACT','o3_ACT','pm_ACT','mo3_ACT','pmEC_ACT','no2_ACTminus75WH','mo3_ACTminus75WH','pm_ACTminus75WH','o3_ACTminus75WH','pmEC_ACTminus75WH']
#datas_headers= ['noxEmis_base']

#exit()

# The filename and directory
#fdirout = '/projects/b1045/montgomery/'
#fnameout= 'baseline.shp'
fdirout = '/projects/b1045/scamilleri/'
#fnameout= '75WH_ACT_ACTminus75WH_Annual_allPollEC_NEW.shp'
fnameout= 'eLDV_BaseDiff_Annual_allPoll_NEW.shp'

make_shp(datas,datas_headers,lon3,lat3,tract_shapes)

