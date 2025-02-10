#!/usr/bin/env python

# -----------------
# Step 3
# -----------------
# Stacy Montgomery, May 2019
# Find correlations between WRF and Station data, output figures

# ------------------------------------------------------------------------------------
# ~~~~~~ START USER INPUT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# variables of interest
minTemp = 242; maxTemp = 294;

#runname='wrf_pure_PXLSM'
#runname='BASE_PXLSM_v0'
#runnamed02='output_CONUS4K_d02_2020_RWC_4km_sf_rrtmg_10_10_1_v3852'
runnamed02='output_CONUS4K_d02_4km_sf_rrtmg_10_10_1_v3852'
runname = runnamed02
#dirout='/projects/b1045/vlang/CMAQ_LCD/'+runname+'/'
dirout='/home/ksz4578/Heat_Pump_Project/WRF_analyses/CMAQ_LCD/'+runnamed02+'/'

# Processed US data, from previous file
comp_dataset_name = dirout+'wrfcheck_withstations_'+runname+'_012016.csv'                     # name and directory to write out to
comp_dataset_extra = dirout+'completeddata_mini_extras2.csv'             # has dates 
station_names = dirout+'station_out_removedmissing.csv'     # names of stations, lat, lon of stations 
comp_dataset_name2= dirout+'wrfcheck_withstations_complete_rain.csv' 

# Location of WRF output
#dirToWRF='/projects/b1045/wrf-cmaq/output/Chicago_LADCO/'+ runname +'/'
#dirToWRF='/projects/b1045/wrf-cmaq/output/SAVEUR_201304/'+ runname +'/'
dirToWRF = '/projects/b1045/wrf-cmaq/output/CONUS4K/'+runnamed02+'/'

#where to output figures
picdir=dirout+'/figures/'

Chatty= True       # false if you want to remove print statements

if Chatty: print('Starting ....')

#name=['t2d01.csv', 't2d02.csv', 't2d03.csv', 'raind01.csv', 'raind02.csv', 'raind03.csv', 'rainncd01.csv', 'rainncd02.csv', 'rainncd03.csv']
name=['t2d02.csv', 'winddird02.csv', 'windsd02.csv', 'rhd02.csv']


# -------------------------------------------------------------------------------------
def getWRFfromIND(ncfile,indxy, filenames,varname):
    t2d01=[ncfile[z][varname][i] for z in range(len(filenames)) for i in range(24)]
    t2d01_xx= [[t2d01[t][indxy[l]] for t in range(24*len(filenames_d01))] for l in range(len(indxy))]
    del t2d01
    return t2d01_xx

#t2d01=[(filenames[z],varname,i) for z in range(len(filenames)) for i in range(24)] #check days n shit
#t2d01_xx=[[(t2d01[t], indxy[l]) for t in range(24*len(filenames_d01))] for l in range(len(indxy))]

def kelvinToF(K):
   F=(K - 273.15)*9/5 + 32
   return F 

def mmToInch(mm):
   return 0.039370078740157*mm

def find(lst, a):
    return [i for i, x in enumerate(lst) if x==a]

# get weekly temperature in d03
def hourlyTempOverweek(numdays, dates,day, stationData,indz,l,q,ind0x,name):
#for zz in range(1):
   st_real=[]; hourly_wrf=[[] for i in range(len(l))]; hourly_station =[]; hourly_r=[]; real_list=[]; wrf_list=[]
   indOfDate =find(dates,day)[0]
   start=24*indOfDate
   end=24*indOfDate+24*numdays
   for t in range(len(stationData)):
         tmps=[stationData[indz[start+u]][t] for u in range(len(indz[start:end]))]
         st_real.append(tmps)
   tmp_wrf=l[q][indz[start:end]]
   tmp_wrf_transpose=tmp_wrf.T
   st_realClip=[]
   if (q == 1) or (q== 2):
       le=len(st_real[0])
       st_realClip=[np.compress(stationData[ind0x],pd.DataFrame(st_real)[qq].tolist()) for qq in range(le)]
       st_realClip=pd.DataFrame(st_realClip)
       print('not d01')
   else:
       le=len(st_real[0])
       st_realClip=[pd.DataFrame(st_real)[qq].tolist() for qq in range(le)]
       st_realClip=pd.DataFrame(st_realClip)
   # Subplots are organized in a Rows x Cols Grid
   # Tot and Cols are known
   Tot = len(st_realClip.columns)
   Cols = 3
   # Compute Rows required
   Rows = Tot // Cols 
   Rows += Tot % Cols
   # Create a Position index
   Position = range(1,Tot + 1)
   fig = plt.figure(1,figsize=(8, 11))
   for i in range(Tot):
       print(i)
       ax = fig.add_subplot(Rows,Cols,Position[i])
       plt.plot(range(0,len(st_realClip[i])),st_realClip[i], label='real')
       plt.plot(range(0,len(st_realClip[i])), kelvinToF(np.array(tmp_wrf_transpose[i])), label='wrf')
       ax.set_title('%s'%(name[i],))
   #outof loop
   fig.subplots_adjust(hspace=0.5, wspace=0.5)
   ax.legend(bbox_to_anchor=(-.5, -.5), loc='lower left', borderaxespad=0.1)
   plt.savefig(picdir+'/%s_%s_%s_hourly_station_and_WRF.png'%(ind0x,day,'T2',))
 #  plt.show()
   plt.close()

#find correlation at hourly level -----------------------------------------------------
#q is for d01,d02,d03

def getHourlyStation(dates,day, stationData,indz,l,q,ind0x):
#for i in range(1):
   st_real=[]; hourly_wrf=[[] for i in range(len(l))]; hourly_station =[]; hourly_r=[]; real_list=[]; wrf_list=[]
   indOfDate =find(dates,day)[0]
   start=24*indOfDate
   end=24*indOfDate+24
   for t in range(len(stationData)):
         #tmps=[stationData[indz[start+u]][t] for u in range(len(indz[start:end]))]
         st_real.append([stationData[indz[start+u]][t] for u in range(len(indz[start:end]))])
   #
   wrf_indz=[int(indz[i]) for i in range(len(indz))]
   tmp_wrf=l[q][wrf_indz[start:end]]
   tmp_wrf_transpose=tmp_wrf.T
   st_realClip=[]
   if (q == 1) or (q== 2):
       le=len(st_real[0])
       st_realClip=[np.compress(stationData[ind0x],pd.DataFrame(st_real)[qq].tolist()) for qq in range(le)]
       st_realClip=pd.DataFrame(st_realClip)
       print('not d01')
   else:
       le=len(st_real[0])
       st_realClip=[pd.DataFrame(st_real)[qq].tolist() for qq in range(le)]
       st_realClip=pd.DataFrame(st_realClip)
   #
   for z in range(len(tmp_wrf)):
      #print(z)
      stn1=kelvinToF(np.array(tmp_wrf_transpose[z]))
      bad = ~np.logical_or(np.isnan(stn1), np.isnan(st_realClip[z]))
      wrf=np.compress(bad, stn1); real=np.compress(bad, st_realClip[z])
      try:
         r=linregress(wrf,real).rvalue
         hourly_r.append(r)
      except:
         hourly_r.append(np.nan)
      real_list.append(list(real))
      wrf_list.append(list(wrf))   
   #
   return hourly_r,real_list,wrf_list


#dynamic plotting-----------------------------------------
# Need to find station names lollllzzz

def dynamicDaily(hourly_r,real_list,wrf_list,name, domain,day,var):
   # Subplots are organized in a Rows x Cols Grid
   # Tot and Cols are known
   Tot = len(wrf_list)
   Cols = 3
   # Compute Rows required
   Rows = Tot // Cols 
   Rows += Tot % Cols
   # Create a Position index
   Position = range(1,Tot + 1)
   fig = plt.figure(1,figsize=(8, 11))
   for i in range(Tot):
       print(i)
       ax = fig.add_subplot(Rows,Cols,Position[i])
       plt.plot(range(0,len(wrf_list[i])),wrf_list[i],label='wrf')
       plt.plot(range(0,len(wrf_list[i])),real_list[i],label='real')
       ax.set_title('%s, R: %.2f'%(name[i],hourly_r[i],))
   ax.legend(bbox_to_anchor=(-1.0, -.5), loc='lower left', borderaxespad=0.)
   fig.subplots_adjust(hspace=0.5, wspace=0.5)
   plt.savefig(picdir+'/%s_%s_%s_hourly_station_and_WRF.png'%(domain,day,var,))
  # plt.show()
   plt.close()

#regular station plotting-----------------------------------------
#get rid of interpolation

def plotStationFig(ncfile,rho, x, y, grdx, grdy, dates, day, coef, domain, var):
   from scipy.interpolate import griddata
   #fig = plt.figure(num=1)
   #ax = plt.axes(projection= ccrs.PlateCarree())
   #cmap = plt.get_cmap('Spectral_r')
   #bad = ~np.logical_or(np.isnan(rho), np.isnan(x))
   #rhoz,x1,y1=np.compress(bad, rho),np.compress(bad, x),np.compress(bad, y)
   #--- do interpolation
   #zi = griddata((x1,y1),rhoz,(grdx,grdy),method='cubic')
   #plt.pcolormesh(grdx,grdy,zi, transform= ccrs.PlateCarree(), cmap= cmap)
   #plt.scatter(x1,y1,c=rhoz, transform= ccrs.PlateCarree(), cmap= cmap, edgecolors= 'k')
   #ax.coastlines()
   #state_boundaries = cfeature.NaturalEarthFeature(category='cultural',
   #                name='admin_1_states_provinces_lines',
   #                scale='110m', facecolor='none')
   #ax.add_feature(cfeature.OCEAN, zorder=10, edgecolor='k')
   #ax.add_feature(cfeature.LAKES, zorder=10, edgecolor='k')
   #ax.add_feature(state_boundaries, edgecolor='black')
   #ax.add_feature(cfeature.BORDERS, linewidth=1, edgecolor='black')
   #plt.title('%s Real Temp on %s, R=%.2f'%(domain,day,coef,))
   #plt.colorbar()
   #plt.show()
   #plt.savefig('./%s_%s_%s_wrf_station.png'%(domain,day,var,)
   #plt.close(fig)
   #from scipy.interpolate import griddata
   #---------------
   #WRF plot   #----------------
   fig = plt.figure()
   sns.set()
   bad = ~np.logical_or(np.isnan(rho), np.isnan(x))
   rhoz,x1,y1=np.compress(bad, rho),np.compress(bad, x),np.compress(bad, y)
   #fig = plt.figure(num=2)
   ind=find(dates,day)[0]
   if var == 'T2':
       print('True')
       averageday=[ncfile[ind][var][i] for i in range(24)]
       grddata = kelvinToF(sum(averageday)/24)
   elif var =='Precip':
       averageday=[ncfile[ind]['RAINC'][i]+ncfile[ind]['RAINNC'][i] for i in range(24)]
       grddata = mmToInch(np.array(sum(averageday)/24))
   else:
       averageday=[ncfile[ind]['RAINC'][i]+ncfile[ind]['RAINNC'][i] for i in range(24)]
       grddata = sum(averageday)/24
   ax = plt.axes(projection= ccrs.PlateCarree())
   cmap = plt.get_cmap('Spectral_r')
   plt.pcolormesh(grdx,grdy, grddata, transform= ccrs.PlateCarree(), cmap= cmap, vmin=min(rhoz), vmax=max(rhoz))
   plt.scatter(x1,y1,c=rhoz, transform= ccrs.PlateCarree(), cmap= cmap, edgecolors= 'k', vmin=min(rhoz), vmax=max(rhoz),s=5)
   ax.coastlines()
   state_boundaries = cfeature.NaturalEarthFeature(category='cultural',
                   name='admin_1_states_provinces_lines',
                   scale='110m', facecolor='none')
   ax.add_feature(state_boundaries, edgecolor='black')
   ax.add_feature(cfeature.BORDERS, linewidth=1, edgecolor='black')
   ax.add_feature(cfeature.LAKES, linewidth=1, edgecolor='black', facecolor='none')
   plt.title('%s avg. daily %s on %s, R=%.2f'%(domain,var,day,coef,))
   plt.colorbar()
   plt.savefig(picdir+'/%s_%s_%s_wrf_station.png'%(domain,day,var,))
   plt.show()
   plt.close()



def correlationplot(hourly_r, x,y,grdx, grdy, dates,day,coef,domain,var):
   from scipy.interpolate import griddata
   fig = plt.figure(1,figsize=(8, 11))
   ax = plt.axes(projection= ccrs.PlateCarree())
   ind=find(dates,day)[0]
   if var=='T2':
       averageday=[ncfile[ind][var][i] for i in range(24)]
   if var =='Precip':
       averageday=[ncfile[ind]['RAINC'][i]+ncfile[ind]['RAINNC'][i] for i in range(24)]
   grddata = sum(averageday)/24
   cmap = plt.get_cmap('Spectral_r')
   cmap2= plt.get_cmap('seismic')
   pa=plt.pcolormesh(grdx,grdy, grddata, transform= ccrs.PlateCarree(), cmap= cmap)
   cba = plt.colorbar(pa,shrink=0.5)
   pb=plt.scatter(x,y,c= hourly_r, transform= ccrs.PlateCarree(),cmap=cmap2,s=5)
   cbb = plt.colorbar(pb,shrink=0.5)
   ax.coastlines()
   state_boundaries = cfeature.NaturalEarthFeature(category='cultural',
                   name='admin_1_states_provinces_lines',
                   scale='110m', facecolor='none')
   ax.add_feature(cfeature.OCEAN, zorder=10, edgecolor='k')
   ax.add_feature(state_boundaries, edgecolor='black')
   ax.add_feature(cfeature.BORDERS, linewidth=1, edgecolor='black')
   plt.title('Hourly R on %s, Daily R=%.2f'%(day,coef,))
   if var=='T2':
       cba.set_label('Temperature (F)')
   elif var=='Precip':
       cba.set_label('Inches')
   else:
       cba.set_label('units needed')
   cbb.set_label('Correlation')
   plt.savefig(picdir+'/%s_%s_correlation.png'%(domain,day,))
#   plt.show()
#   plt.close(fig)
   plt.close()

#remove missing files
def rm_missing(filenames_d01):
   testrm=[]
   for i in filenames_d01:
      try:
         test=Dataset(dirToWRF+i)
      except FileNotFoundError:
         print('File not found for ' + i)
         testrm.append(i)
#
   for i in testrm:
      filenames_d01.remove(i)
   #return
   return filenames_d01

# -------------------------------------------------------------------------------------
# ~~~~~~ IMPORT PACKAGES ~~~~~~~~~~~~
#Station
import glob, os
import pandas as pd, numpy as np, matplotlib.pyplot as plt, cartopy.crs as crs, cartopy.feature as cpf
from netCDF4 import Dataset; from matplotlib.cm import get_cmap
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords)
import time; from timezonefinder import TimezoneFinder; from pytz import timezone
import pytz; from datetime import datetime,date, timedelta; import dateutil.parser as dparser

tf = TimezoneFinder(in_memory=True)

import cartopy.crs as ccrs;import scipy.interpolate
import cartopy.feature as cfeature
import seaborn as sns
#------------------------------------------------------------------------------
# ~~~~~~ START MAIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#-------------------------------------------------------------------
# compute correlation ----------------------------------------------
#-------------------------------------------------------------------


# $1 Get WRF file names
filenames_d01=[] 
os.chdir(dirToWRF)
for file in glob.glob("wrfout_d01_*"):
    filenames_d01.append(file)

filenames_d01.sort() #files are now sorted by date and time

# $1 Get WRF file names
filenames_d02=[] 
os.chdir(dirToWRF)
for file in glob.glob("wrfout_d02_*"):
    filenames_d02.append(file)

filenames_d02.sort() #files are now sorted by date and time

# $1 Get WRF file names
filenames_d03=[] 
os.chdir(dirToWRF)
for file in glob.glob("wrfout_d03_*"):
    filenames_d03.append(file)

filenames_d03.sort() #files are now sorted by date and time

# remove missing files
filenames_d01=rm_missing(filenames_d01)
filenames_d02=rm_missing(filenames_d02)
filenames_d03=rm_missing(filenames_d03)

# pull out viable dates for analysis
dates=[filenames_d01[z].split("wrfout_d01_")[1].split("_00:00:00")[0] for z in range(len(filenames_d01))]
dates = dates[9:-1]

ncfiled01 = [Dataset(filenames_d01[i]) for i in range(len(filenames_d01))]
ncfiled02 = [Dataset(filenames_d02[i]) for i in range(len(filenames_d02))]
ncfiled03 = [Dataset(filenames_d03[i]) for i in range(len(filenames_d03))]
 
ncfiled02 = ncfiled01 
ncfiled03 = ncfiled01 

# Editted this to only look at d01 
wrf_latd01, wrf_lond01 = latlon_coords(getvar(Dataset(filenames_d01[1]),"RAINNC"))
wrf_latd02, wrf_lond02 = latlon_coords(getvar(Dataset(filenames_d01[1]),"RAINNC"))
wrf_latd03, wrf_lond03 = latlon_coords(getvar(Dataset(filenames_d01[1]),"RAINNC"))

# ------------------------------- Daily Correlations ---------------------------------------

l=[pd.read_csv(dirout+name[i],header=0,index_col=0) for i in range(len(name))] # each l is a diff variable 
                                                        # from wrf file, made from file2

# Get columns headers to be integers
for i in range(len(l)):
   l[i].columns=l[i].columns.astype(int)



# Adjust from accumulation to hourly
for z in range(0,3):
   l_hourly=np.array([l[z][i]-l[z][(i-1)] if i>0 else l[z][i] for i in range(len(list(l[z])))])
   l[z]=pd.DataFrame(l_hourly).T
   # set first hour as no rain --- would probably be wrong but
   l[z][0]=0 

# This is an artifact from the fact that when i first wrote the code I couldn't make the 
# header not be string, so I made strings, but now I can't use strings
indz=[i for i in range(24*len(dates))]
indz_str=[str(i) for i in range(24*len(dates))]

stationData=pd.read_csv(comp_dataset_name,header=0,index_col=0)
stationDataRain= pd.read_csv(comp_dataset_name2,header=0,index_col=0)

#stationData.columns=stationData.columns.astype(int)
#stationDataRain.columns=stationDataRain.columns.astype(int)

#set up vaxxriables for rain
daily_wrf=[[] for i in range(3)] # 3 is num domains
daily_st =[]; daily_std02=[]; daily_std03=[]

#daily_wrf_rain=[[] for i in range(3)] # num domains is 3
daily_st_rain =[]; daily_std02_rain=[]; daily_std03_rain=[]

#calculate daily averages for WRF
for q in range(3): # 3 because there are 3 domains in our list of L
   for z in range(len(dates)):
      start=24*z
      end=24*z+24
      # pull out 24 hours of temp
      #print(l, q, [indz[start:end]])
      tmp=l[q][indz[start:end]]
      # add rainnc + rainn over 24 hours
      #tmpR=l[q+3][indz[start:end]]+l[q+6][indz[start:end]]
      # get average, assume there are no missing values
      daily_wrf[q].append(list(tmp.sum(axis=1)/24))
      #daily_wrf_rain[q].append(list(tmpR.sum(axis=1)/24))

# calculate daily average for stations
for z in range(len(dates)):
   start=24*z
   end=24*z+24
   station=[]; stationR=[]
   for t in range(len(stationData)):
         try:
            tmps=[stationData[indz_str[start+u]][t] for u in range(len(indz[start:end]))]
            tmpsR=[stationDataRain[indz_str[start+u]][t] for u in range(len(indz[start:end]))]
         except Exception as e:
            if end >= 740:
               end = 740
               tmps=[stationData[indz_str[start+u]][t] for u in range(len(indz[start:end]))]
               tmpsR=[stationDataRain[indz_str[start+u]][t] for u in range(len(indz[start:end]))]
         try:
            daily=np.nanmean(tmps)
         except:
            daily=np.nan
         try:
            dailyR=np.nanmean(tmpsR)
         except:
            dailyR=np.nan
         station.append(daily)
         stationR.append(dailyR)
   daily_st.append(station)
   daily_st_rain.append(stationR)
         

from scipy.stats.stats import pearsonr
from scipy.stats import linregress

#get daily correlations
rd01_daily=[];rd02_daily=[];rd03_daily=[];


for i in range(len(daily_st)):
   daily_std02= np.compress(stationData['in_d02'],daily_st[i])
   x= daily_std02
   y= daily_wrf[1][i]
   bad = ~np.logical_or(np.isnan(x), np.isnan(y))
   rd02_daily.append(linregress(np.compress(bad, x),np.compress(bad, y)).rvalue)


worst_d02=[(dates[rd02_daily.index(sorted(rd02_daily)[i])],rd02_daily[rd02_daily.index(sorted(rd02_daily)[i])]) for i in range(3)]
best_d02=[(dates[rd02_daily.index(sorted(rd02_daily)[i])],rd02_daily[rd02_daily.index(sorted(rd02_daily)[i])]) for i in range(-3,0)]

#----------------------------------------------------------------
names = pd.read_csv(station_names)
names=names['stn']
names=[names[i].split('.csv')[0] for i in range(len(names))]
name_d02=np.compress(stationData['in_d02'],np.array(names))

lat_std02= np.compress(stationData['in_d02'], stationData['lat'])
lon_std02= np.compress(stationData['in_d02'], stationData['lon'])

#----------------------------------------------------------------

# Get correlations and plots
#rho = values for grid; x = stationlat, y= staitony; grdxgridy make sense, coef

var='T2'
for t in range(len(worst_d02)):
    ind0x='in_d02'
    x= lon_std02; y= lat_std02
    worst=worst_d02; best=best_d02
    day_worst,coef_worst= worst[t][0], worst[t][1]; 
    indOfDate_worst =find(dates, day_worst)[0]
    day_best,coef_best= best[t][0], best[t][1]; 
    indOfDate_best =find(dates, day_best)[0]
    grdx= wrf_lond02; grdy= wrf_latd02
    domain='d02'
    ncfile= ncfiled02
    rho_best= np.compress(stationData[ind0x], daily_st[indOfDate_best])
    rho_worst= np.compress(stationData[ind0x], daily_st[indOfDate_worst])
    hourly_r_best,real_list_best,wrf_list_best =getHourlyStation(dates,day_best, stationData,indz_str,l,q,ind0x)
    hourly_r_worst,real_list_worst,wrf_list_worst =getHourlyStation(dates,day_worst, stationData,indz_str,l,q,ind0x)
    plotStationFig(ncfile,rho_best, x, y, grdx, grdy, dates, day_best, coef_best, domain,var)
    plotStationFig(ncfile,rho_worst, x, y, grdx, grdy, dates, day_worst, coef_worst, domain,var)
    correlationplot(hourly_r_best, x,y,grdx, grdy, dates, day_best, coef_best,domain,var)
    correlationplot(hourly_r_worst, x,y,grdx, grdy, dates, day_worst, coef_worst,domain,var)
    # if q ==2:
    #   dynamicDaily(hourly_r_best,real_list_best,wrf_list_best, name_d03, domain, day_best,var)
    #   dynamicDaily(hourly_r_worst,real_list_worst,wrf_list_worst, name_d03, domain, day_worst,var)
#          hourlyTempOverweek(7, dates, day_best, stationData,indz,l,q,'in_d03', name_d03)

#plot hourly temp and wrf hourly temp
#numdays=7
#numdays, dates,day, stationData,indz,l,q,ind0x,name = 7,dates,day_best, stationData,indz,l,q,ind0x,name_d03

# for q in range(3):
#    for t in range(len(worst_d01)):
#       if q==2:
#          ind0x='in_d03'
#          x= lon_std03; y= lat_std03
#          worst=worst_d03; best=best_d03
#          day_worst,coef_worst= worst[t][0], worst[t][1]; 
#          indOfDate_worst =find(dates, day_worst)[0]
#          day_best,coef_best= best[t][0], best[t][1]; 
#          indOfDate_best =find(dates, day_best)[0]
#          grdx= wrf_lond03; grdy= wrf_latd03
#          domain='d03'
#          ncfile= ncfiled03
#          rho_best= np.compress(stationData[ind0x], daily_st[indOfDate_best])
#          rho_worst= np.compress(stationData[ind0x], daily_st[indOfDate_worst])
#          hourlyTempOverweek(7, dates, day_worst, stationData,indz,l,q,'in_d03', name_d03)


print("Done with step 3")


