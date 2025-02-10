# Check wrf-cmaq output
# Need EPA annual AQS data
# Will append the corresponding CMAQ value to the file

import pandas as pd
import numpy as np
from netCDF4 import Dataset
from datetime import timedelta, date,datetime;

#---------------------------------------------------------#
# Before running, pull in files from here: https://aqs.epa.gov/aqsweb/airdata/download_files.html
#---------------------------------------------------------#

# Configure these inputs

domain='d03_WHbase'
time='hourly'
year='2019'
month='01'
#epa_code=['42401','42602','44201','42101','88101']; var=['SO2','NO2','O3','CO','PM25_TOT']
stdate,enddt = 1,30
yr = 2019
mo = 1
name = '20190101'
start_dt = datetime(yr, mo, stdate); end_dt = datetime(yr,mo,enddt,23)
fstart='COMBINE_ACONC_'

#CMAQ output
#dir_cmaq='E:/CMAQ_test/'
#dir_cmaq ='/projects/b1045/wrf-cmaq/output/CONUS4K/output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852/postproces
#dir_cmaq ='/projects/b1045/wrf-cmaq/output/CONUS4K/output_CONUS4K_d02_Grace_0pctLB_4km_sf_rrtmg_10_10_1_v38
dir_cmaq ='/projects/b1045/wrf-cmaq/output/Chicago_LADCO/output_ChicagoLADCO_d03.2019-01_WHbase_1.33km_sf_rrtmg_5_8_1_v3852/postprocess/'

# Output of this script
dir_epa='/projects/b1045/vlang/CMAQ_EPA/'

#grid='/projects/b1045/wrf-cmaq/output/CONUS4K/output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852/lat_lon_CONUS4K
grid='latlon_ChicagoLADCO_d03.nc'

# EPA codes for the diff species to make the file names of AQS fies
epa_code=['42401','42602','44201','42101','88101','81102']; var=['SO2','NO2','O3','CO','PM25_TOT','PM10'] #numerical identi vars
epa_files = [dir_epa+'%s_%s_%s.csv'%(time,epa_code[i],year,) for i in range(len(epa_code))]


# What is latitude longitude called in the grid file and get corners
#la,lo='LAT','LON'
la,lo='lat','lon'
#cmaq_lon,cmaq_lat=np.array(Dataset(grid)[lo][0][0]),np.array(Dataset(grid)[la][0][0])
cmaq_lon,cmaq_lat=np.array(Dataset(grid)[lo]),np.array(Dataset(grid)[la])
llat,ulat,llon,ulon=cmaq_lat.min(), cmaq_lat.max(), cmaq_lon.min(), cmaq_lon.max()

def find_index(stn_lon, stn_lat, wrf_lon, wrf_lat):
# stn -- points
# wrf -- list
#for iz in range(1):
   xx=[];yy=[]
   for i in range(len(stn_lat)):
   #for i in range(1):
      abslat = np.abs(wrf_lat-stn_lat[i])
      abslon= np.abs(wrf_lon-stn_lon[i])
      c = np.maximum(abslon,abslat)
      latlon_idx = np.argmin(c)
      x, y = np.where(c == np.min(c))
      #add indices of nearest wrf point station
      xx.append(x)
      yy.append(y)
   #
   xx=[xx[i][0] for i in range(len(xx))];yy=[yy[i][0] for i in range(len(yy))]
   #return indices list
   return xx, yy

def make_cmaq_fnames(fstart,year,month,stdate,enddt):
        fnames = [fstart+year+month+str(i).zfill(2)+'.nc' for i in range(stdate,enddt)]
        return fnames

def crop_epa(file, start_dt, end_dt,llat,ulat,llon,ulon):
#for t in range(1):
        df = pd.read_csv(file)
        df = df[(df['Latitude'] >= llat) & (df['Latitude'] <= ulat) & (df.Longitude >= llon) & (df.Longitude <= ulon)].reset_index()
        df['Datetime GMT']=pd.to_datetime(df['Date GMT']+ ' ' + df['Time GMT'])
        df= df[(df['Datetime GMT'] >= pd.to_datetime(start_dt) ) & (df['Datetime GMT'] <= pd.to_datetime(end_dt))]
        lon,lat=df['Longitude'].unique(), df['Latitude'].unique()
        #return_files =[dir_epa+'%s_%s_%s.csv'%(time,epa_code[i],year,) for i in range(len(epa_code))]
        return lon,lat,df

# Start scripts
#---------------------------------------------------------#
#cmaq_lon,cmaq_lat=np.array(Dataset(grid)[lo][0][0]),np.array(Dataset(grid)[la][0][0])
#llat,ulat,llon,ulon=cmaq_lat.min(), cmaq_lat.max(), cmaq_lon.min(), cmaq_lon.max()

fnames = make_cmaq_fnames(fstart,year,month,stdate,enddt+1)
cmaq=[Dataset(dir_cmaq+fnames[i]) for i in range(len(fnames))]
#t_index = pd.DatetimeIndex(start=start_dt, end=end_dt, freq='1h')

for v in range(len(var)):
        print(epa_files[v], start_dt, end_dt,llat,ulat,llon,ulon)
        epa_lon,epa_lat,df = crop_epa(epa_files[v], start_dt, end_dt,llat,ulat,llon,ulon)
        cmaqv = [np.array(cmaq[d][var[v]]) for d in range(len(cmaq))]
        cmaqu10 = [np.array(cmaq[d]['U10']) for d in range(len(cmaq))]
        cmaqv10 = [np.array(cmaq[d]['V10']) for d in range(len(cmaq))]

        for i in range(len(epa_lon)):
        #       print(tmp)
                tmp = df[(df['Longitude'] == epa_lon[i]) & (df['Latitude'] == epa_lat[i])]
                tsn = tmp['State Name'].unique()[0]
                tcn = tmp['County Name'].unique()[0]
                um = tmp['Units of Measure'].unique()[0]
                x,y=find_index([epa_lon[i]],[epa_lat[i]],cmaq_lon,cmaq_lat);
                x,y= x[0],y[0]
                if (x != 0) | (y != 0) | (x != cmaq_lon.shape[0]) | (y != cmaq_lat.shape[1]): #remove edge s
                        tmp.index = tmp['Datetime GMT']
                        t_index = pd.date_range(start=start_dt, end=end_dt, freq='1h')
                        tmp = tmp.resample('1H').mean().reindex(t_index).fillna(np.nan)
                        tmp['x'],tmp['y']=x,y
                        tmp['State Name'] = tsn
                        tmp['County Name'] = tcn
                        tmp['Units of Measure'] = um
                        #tmp['CMAQ_'+var[v]] = [cmaq[d][var[v]][t][0][x][y] for d in range(len(cmaq)) for t

                        z = [cmaqv[d][:,0,x,y] for d in range(len(cmaq))]
                        zi = []
                        z_u = [cmaqu10[d][:,0,x,y] for d in range(len(cmaq))]
                        zi_u = []
                        z_v = [cmaqv10[d][:,0,x,y] for d in range(len(cmaq))]
                        zi_v = []

                        for q in range(len(z)):
                                zi = zi + list(z[q])
                                zi_u = zi_u + list(z_u[q])
                                zi_v = zi_v + list(z_v[q])

                        z = np.array(zi).ravel(); z_u = np.array(zi_u).ravel(); z_v = np.array(zi_v).ravel()
                        print(len(z))
                        if len(z) < len(tmp): z = list(z) + [np.nan]*(len(tmp)-len(z))
                        tmp['CMAQ'] = z
                        tmp['CMAQ_u10'] = z_u
                        tmp['CMAQ_v10'] = z_v

                        # remove dt
                        tmp['dt'] = tmp.index
                        tmp['level_0'] = tmp.index
                        tmp=tmp.reset_index(drop=True)
                        if i == 0:
                                print('here')
                                final_df = tmp
                        else:
                                final_df = final_df.append(tmp,ignore_index=True)
        # output epa with cmaq final_df.to_csv(dir_epa+var[v]+'_'+domain+'_'+year+'_withCMAQ.csv')
        final_df.loc[final_df['Sample Measurement'] <= 0] = np.nan # remove zeros
        final_df.loc[final_df['Sample Measurement'] < 0] = np.nan # remove negatives
        final_df.to_csv(dir_epa+var[v]+'_'+domain+'_'+str(year)+'_'+str(mo)+'_EPA_CMAQ_Combine.csv'); print('done '+var[v])
        #fout = #CO_d03_2018_10_EPA_CMAQ_Combine.csv

