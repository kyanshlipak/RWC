# Model evauation
import pandas as pd
from scipy.stats import pearsonr
import numpy as np
from sklearn.metrics import mean_squared_error
import scipy.stats as st
import matplotlib.pyplot as plt

# functions
def stats_normalized(data,prediction):
    x,y=data[~np.isnan(data)],prediction[~np.isnan(data)] # get rid of NaNs
    mb = np.mean(y - x) #predicted - observed
    rmse= np.sqrt(mean_squared_error(y,x))
    mu_d,mu_p = np.mean(x),np.mean(y)
    nmb = np.sum(y-x)/np.sum(x)*100
    nme = np.sum(np.abs(y-x))/np.sum(x)*100
    r,p = st.pearsonr(x,y)
    return mb,rmse,mu_d,mu_p,nmb,nme,r,p

# CMAQ RUN things
domain='d03_WHbase'  #"CMAP"
time='hourly'
dir_epa='/projects/b1045/vlang/CMAQ_EPA/'
epa_code=['42401','42602','44201','42101','88101','81102']; var=['SO2','NO2','O3','CO','PM25_TOT','PM10'] #numerical identifiers and corresponding vars
#epa_code=  #numerical identifiers and corresponding vars
years= ['2019']  #['2018','2018','2019','2019']
months=['1']         #['8','10','1','04']

epa_files = []
names=[]
for i in range(len(years)):
        year,month = years[i],months[i]
        epa_files.append([dir_epa+'%s_%s_%s_%s_EPA_CMAQ_Combine.csv'%(var[loop],domain,year,month) for loop in range(len(var))])
        names.append(['%s_%s_%s_%s'%(var[loop],domain,year,month) for loop in range(len(var))])


epa_files=[]
names=[]
for i in range(len(var)):
    for m in range(len(months)):
        month,year = months[m],years[m]
        epa_files.append([dir_epa+'%s_%s_%s_%s_EPA_CMAQ_Combine.csv'%(var[i],domain,year,month)])
        names.append('%s_%s_%s_%s'%(var[i],domain,year,month))

epa_files = np.array(epa_files).ravel(); names= np.array(names).ravel()


# Create hourly correlation
corrs = []
for i in range(len(epa_files)):
        f = pd.read_csv(epa_files[i])
        if f['Units of Measure'][10] == 'Parts per million': x,y = f['Sample Measurement']*1000,f['CMAQ']
        else: x,y = f['Sample Measurement'],f['CMAQ']
        f['Hourly Bias'] = (y-x)
        mb,rmse,mu_d,mu_p,nmb,nme,r,p = stats_normalized(x,y)
        corrs.append([mb,rmse,mu_d,mu_p,nmb,nme,r,p])
        f.to_csv(dir_epa+var[i]+'_hourlybias_'+domain+'_'+str(year)+'_EPA_CMAQ_Combine.csv');

corrs = pd.DataFrame(corrs)
corrs.columns = ['mb','rmse','mu_d','mu_p','nmb','nme','r','p']
corrs['Simulation'] = names
corrs.to_csv(dir_epa+'hourly_stn_coefficients_%s.csv'%(domain))


# create daily correlation
daily_bias = pd.DataFrame() # Make empty dataframe to retain daily bias for spatial analysis
daily_corrs = []
for i in range(len(epa_files)):
        f = pd.read_csv(epa_files[i])
        f['SiteID'] = f['State Code'].astype(str) + '-' + f['County Code'].astype(str) + '-' + f['Site Num'].astype(str)
        f['level_0']=pd.to_datetime(f['level_0'])
        if f['Units of Measure'][10] == 'Parts per million': f['Sample Measurement'],f['CMAQ'] = f['Sample Measurement']*1000,f['CMAQ']
        f.index = f['level_0']
        #favg = f.groupby('Site Num').resample('D').mean()
        favg = f.groupby('SiteID').resample('D').mean()
        favg=favg.reset_index(level=0, drop=True).reset_index()
        daily_bias = favg
        x,y = favg['Sample Measurement'],favg['CMAQ']
        daily_bias['daily_bias']= y-x
        mb,rmse,mu_d,mu_p,nmb,nme,r,p = stats_normalized(x,y)
        daily_corrs.append([mb,rmse,mu_d,mu_p,nmb,nme,r,p])
        daily_bias.to_csv(dir_epa+var[i]+'_dailybias_'+domain+'_EPA_CMAQ_Combine.csv');

daily_corrs = pd.DataFrame(daily_corrs)
daily_corrs.columns = ['mb','rmse','mu_d','mu_p','nmb','nme','r','p']
daily_corrs['Simulation'] = names
daily_corrs.to_csv(dir_epa+'daily_stn_coefficients_%s.csv'%(domain))

# create daily correlation
daily_max_corrs = []
for i in range(len(epa_files)):
        f = pd.read_csv(epa_files[i])
        f['SiteID'] = f['State Code'].astype(str) + '-' + f['County Code'].astype(str) + '-' + f['Site Num'].astype(str)
        f['level_0']=pd.to_datetime(f['level_0'])
        if f['Units of Measure'][10] == 'Parts per million': f['Sample Measurement'],f['CMAQ'] = f['Sample Measurement']*1000,f['CMAQ']
        f.index = f['level_0']
        #favg = f.groupby('Site Num').resample('D')
        favg = f.groupby('SiteID').resample('D')
        x,y = favg['Sample Measurement'].max(),favg['CMAQ'].max()
        x,y = x.reset_index(),y.reset_index()
        x,y = x['Sample Measurement'],y['CMAQ']
        mb,rmse,mu_d,mu_p,nmb,nme,r,p = stats_normalized(x,y)
        daily_max_corrs.append([mb,rmse,mu_d,mu_p,nmb,nme,r,p])

daily_max_corrs = pd.DataFrame(daily_max_corrs)
daily_max_corrs.columns = ['mb','rmse','mu_d','mu_p','nmb','nme','r','p']
daily_max_corrs['Simulation'] = names
daily_max_corrs.to_csv(dir_epa+'daily_max_stn_coefficients_%s.csv'%(domain))


plt.scatter(corrs.Simulation, corrs.r,label='hourly')
plt.scatter(daily_corrs.Simulation,daily_corrs.r,label='daily')
plt.scatter(daily_max_corrs.Simulation,daily_max_corrs.r,label='daily max')
#plt.errorbar(corrs.Simulation,corrs.r,label='hourly',yerr=corrs.nme/100,marker='o')
#plt.errorbar(daily_corrs.Simulation,daily_corrs.r,label='daily',yerr=daily_corrs.nme/100,marker='o')
#plt.errorbar(daily_max_corrs.Simulation,daily_max_corrs.r,yerr=daily_max_corrs.nme/100,label='daily max',marker='o')
plt.xticks(rotation = 90)
plt.grid()
plt.legend()
plt.ylabel('pearson r')
plt.tight_layout()
plt.savefig('/projects/b1045/vlang/CMAQ_EPA/2019_01_1.3m_75wh_modelstats.png')
plt.show()
