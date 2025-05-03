# This script validates January 2016 as a typical winter month 
# based on EPA AQS data from 2010-2024.

# Analyses and figures are generated that show Janaury 2016
# looks relatively normal


# %% [markdown]
# ## Import AQS Data

# %%
import os
import pandas as pd

# Directory containing the CSV files
data_dir = "AQS_pm2.5"

# Initialize an empty list to store DataFrames
dataframes = []

# Loop through all files in the directory
for file_name in os.listdir(data_dir):
    if file_name.endswith(".csv"):
        file_path = os.path.join(data_dir, file_name)
        # Read each CSV file into a DataFrame
        df = pd.read_csv(file_path)
        dataframes.append(df)

# Combine all DataFrames into a single DataFrame
combined_data = pd.concat(dataframes, ignore_index=True)

# %%
# Load the CSV file into a DataFrame
file_path = 'AQS_data_2016/PM25_TOT_d02_WHbase_2016_1_EPA_CMAQ_Combine.csv'
df_cmaq_baseline = pd.read_csv(file_path)

file_path = 'AQS_data_2016/PM25_TOT_d02_WHbase_2016_1_EPA_CMAQ_Combine_2020_RWC.csv'
df_cmaq_2020 = pd.read_csv(file_path)

# %% [markdown]
# ## Filter for site numbers consistent across the years

# %%
### create dt columns 

# Convert 'Date Local' to datetime
combined_data['Date Local'] = pd.to_datetime(combined_data['Date Local'])

# Extract year and month for filtering
combined_data['Year'] = combined_data['Date Local'].dt.year
combined_data['Month'] = combined_data['Date Local'].dt.month

df_cmaq_2020.loc[:,'Year'] = (pd.to_datetime(df_cmaq_2020['dt'])).dt.year
df_cmaq_baseline.loc[:,'Year'] = (pd.to_datetime(df_cmaq_baseline['dt'])).dt.year

# %%
### Create a unique site identifier

combined_data['Site ID'] = (
    combined_data['State Code'].astype(str) + "-" +
    combined_data['County Code'].astype(str) + "-" +
    combined_data['Site Num'].astype(str)
)

df_cmaq_baseline['Site ID'] = (
    df_cmaq_baseline['State Code'].dropna().astype(int).astype(str) + "-" +
    df_cmaq_baseline['County Code'].dropna().astype(int).astype(str) + "-" +
    df_cmaq_baseline['Site Num'].dropna().astype(int).astype(str)
)

df_cmaq_2020['Site ID'] = (
    df_cmaq_2020['State Code'].dropna().astype(int).astype(str) + "-" +
    df_cmaq_2020['County Code'].dropna().astype(int).astype(str) + "-" +
    df_cmaq_2020['Site Num'].dropna().astype(int).astype(str)
)



len(combined_data['Site ID'].unique())

# %%
#### FILTERING FOR COMMON SITES

# Get the count of unique years for each site
site_year_counts = combined_data.groupby('Site ID')['Year'].nunique()

# Find sites present in all years
common_sites = site_year_counts[site_year_counts == combined_data['Year'].nunique()].index


filtered_data = combined_data[combined_data['Site ID'].isin(common_sites)]
filtered_cmaq_baseline = df_cmaq_baseline[df_cmaq_baseline['Site ID'].isin(common_sites)]
filtered_cmaq_2020 = df_cmaq_2020[df_cmaq_2020['Site ID'].isin(common_sites)]

data_lost =round( 100 * (combined_data.shape[0] - filtered_data.shape[0])/ combined_data.shape[0] , 2)
print(f"Data lost from filtering {data_lost}%")

# %%
filtered_data.loc[filtered_data['Sample Duration'] == "1 HOUR"][['Sample Duration','Local Site Name', 'Date Local']]

# %%
### DAY NIGHT FILTERING

filtered_cmaq_2020['dt'] = pd.to_datetime(filtered_cmaq_2020['dt'])
filtered_cmaq_2020['hour'] = filtered_cmaq_2020['dt'].dt.hour
filtered_cmaq_baseline['dt'] = pd.to_datetime(filtered_cmaq_baseline['dt'])
filtered_cmaq_baseline['hour'] = filtered_cmaq_baseline['dt'].dt.hour
day_bound = 6
night_bound = 17

df_day_2020 = filtered_cmaq_2020.loc[(filtered_cmaq_2020['hour'] >= day_bound) & (filtered_cmaq_2020['hour'] <= night_bound)]
df_night_2020 = filtered_cmaq_2020.loc[(filtered_cmaq_2020['hour'] < day_bound) | (filtered_cmaq_2020['hour'] > night_bound)]
df_cmaq_baseline_day = filtered_cmaq_baseline.loc[(filtered_cmaq_baseline['hour'] >= day_bound) & (filtered_cmaq_baseline['hour'] <= night_bound)]
df_cmaq_baseline_night= filtered_cmaq_baseline.loc[(filtered_cmaq_baseline['hour'] < day_bound) | (filtered_cmaq_baseline['hour'] > night_bound)]

# filtered_baseline_day =  df_cmaq_baseline_day#[df_cmaq_baseline_day['Site ID'].isin(common_sites)]
# filtered_baseline_night =  df_cmaq_baseline_night[df_cmaq_baseline_night['Site ID'].isin(common_sites)]
# filtered_2020_day =  df_day_2020[df_day_2020['Site ID'].isin(common_sites)]
# filtered_2020_night =  df_night_2020[df_night_2020['Site ID'].isin(common_sites)]



# %%
df_cmaq_baseline_day['CMAQ'].mean()

# %%
df_cmaq_baseline_night['CMAQ'].mean()

# %%
filtered_cmaq_baseline['CMAQ'].mean()

# %%
filtered_data.columns

# %%
filtered_cmaq_2020

# %% [markdown]
# ## Show is a typical ish year

# %%
general_stats = filtered_data.groupby('Year')['Arithmetic Mean'].agg(['mean', 'median', 'std'])
cmaq_baseline_stats = filtered_cmaq_baseline.groupby('Year')['CMAQ'].agg(['mean', 'median', 'std'])
cmaq_2020_stats =   filtered_cmaq_2020.groupby('Year')['CMAQ'].agg(['mean', 'median', 'std'])

cmaq_baseline_stats_day = df_cmaq_baseline_day.groupby('Year')['CMAQ'].agg(['mean', 'median', 'std'])
cmaq_baseline_stats_night = df_cmaq_baseline_night.groupby('Year')['CMAQ'].agg(['mean', 'median', 'std'])
cmaq_2020_stats_day = df_day_2020.groupby('Year')['CMAQ'].agg(['mean', 'median', 'std'])
cmaq_2020_stats_night =  df_night_2020.groupby('Year')['CMAQ'].agg(['mean', 'median', 'std'])

# %%
# Group by year and calculate mean, median, and standard deviation for January PM2.5
january_data = filtered_data[filtered_data['Date Local'].dt.month == 1]
january_stats = january_data.groupby('Year')['Arithmetic Mean'].agg(['mean', 'median', 'std'])

# Statistics for January 2016
jan_2016_stats = january_stats.loc[2016]
print("January 2016 Statistics:")
print(jan_2016_stats)

# Display all January statistics
print("All January Statistics:")
print(january_stats)


# %%
import matplotlib.pyplot as plt

# Plot the mean PM2.5 levels for January
plt.figure(figsize=(10, 6))
plt.plot(january_stats.index, january_stats['mean'], marker='o', label='Mean January PM2.5')
plt.plot(general_stats.index, general_stats['mean'], marker='o', label='Mean Annual PM2.5')
plt.plot(cmaq_baseline_stats.index, cmaq_baseline_stats['mean'], marker = 'o', label='Baseline modeled PM2.5')
plt.plot(cmaq_baseline_stats_night.index, cmaq_baseline_stats_night['mean'], marker = 'o', label='Baseline Night modeled PM2.5')
plt.plot(cmaq_baseline_stats_day.index, cmaq_baseline_stats_day['mean'], marker = 'o', label='Baseline Day modeled PM2.5')

# plt.plot(january_stats.index, january_stats['median'], marker='o', label='Median January PM2.5')
# plt.plot(general_stats.index, general_stats['median'], marker='o', label='Median Annual PM2.5')

# plt.plot(january_stats.index, january_stats['mean'], marker='o', label='Mean January PM2.5')
# plt.plot(general_stats.index, general_stats['mean'], marker='o', label='Mean Annual PM2.5')

plt.axvline(2016, color='red', linestyle='--', label='January 2016')
plt.title('January Mean PM2.5 Levels by Year')
plt.xlabel('Year')
plt.ylabel('Mean PM2.5')
plt.legend()
plt.grid(True)
plt.show()


# %%
import matplotlib.pyplot as plt

# Plot the mean PM2.5 levels for January
plt.figure(figsize=(10, 6))
plt.plot(january_stats.index, january_stats['mean'], marker='o', label='Mean January PM2.5')
plt.plot(general_stats.index, general_stats['mean'], marker='s', label='Mean Annual PM2.5')
# plt.plot(cmaq_2020_stats.index, cmaq_2020_stats['mean'], marker = 'o', label='2020 modeled PM2.5')
# plt.plot(cmaq_2020_stats_night.index, cmaq_2020_stats_night['mean'], marker = 'o', label='2020 Night modeled PM2.5')

# plt.plot(cmaq_2020_stats_day.index, cmaq_2020_stats_day['mean'], marker = 'o', label='2020 Day modeled PM2.5')

# plt.plot(january_stats.index, january_stats['median'], marker='o', label='Median January PM2.5')
# plt.plot(general_stats.index, general_stats['median'], marker='o', label='Median Annual PM2.5')

# plt.plot(january_stats.index, january_stats['mean'], marker='o', label='Mean January PM2.5')
# plt.plot(general_stats.index, general_stats['mean'], marker='o', label='Mean Annual PM2.5')

#plt.axvline(2016, color='red', linestyle='--', label='January 2016')
plt.title('January vs. Annual Average PM₂.₅ Levels by Year',fontsize = 16)
plt.xlabel('Year', fontsize = 14)
plt.ylabel('Average PM₂.₅ Concentrations (µg/m³)', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.legend(fontsize = 14)
plt.grid(True)
plt.show()


# %%
df_cmaq_baseline_day['Year'] = "CMAQ_Baseline"
filtered_cmaq_2020['Year'] = "WRF-CMAQ"

# %%
import matplotlib.pyplot as plt
import seaborn as sns

# Create the box plot
january_data_clean = january_data.dropna(subset=['Arithmetic Mean', 'Year'])

plt.figure(figsize=(12, 6))
sns.boxplot(x='Year', y='Arithmetic Mean', data=january_data_clean)

#sns.boxplot(x='Year', y='CMAQ', data=df_cmaq_baseline_day)
sns.boxplot(x='Year', y='CMAQ', data=filtered_cmaq_2020)


#plt.axvline(2016, color='red', linestyle='--', label='January 2016')
plt.title('Boxplots of January PM₂.₅ Levels by Year', fontsize = 16)
plt.ylabel('Average January PM₂.₅ Concentrations (µg/m³)', fontsize = 14)
plt.xlabel('Year', fontsize = 14)
plt.ylim([0,30])
plt.xticks(rotation=45, fontsize = 14)  # Rotate x-axis labels for readability
plt.yticks(fontsize = 14) 
plt.show()


# %%
# Filter data for Winter months (December, January, February)
winter_data = filtered_data[filtered_data['Date Local'].dt.month.isin([12, 1, 2])]

# Add a column for the month name (optional, for clearer plotting)
winter_data['Month'] = winter_data['Date Local'].dt.month_name()

# %%
import matplotlib.pyplot as plt
import seaborn as sns
filtered_cmaq_2020['Year'] = "WRF-CMAQ"
# Create the box plot comparing PM2.5 levels in December, January, and February
plt.figure(figsize=(12, 6))
sns.boxplot(x='Month', y='Arithmetic Mean', data=winter_data, order=['December', 'January', 'February'])
sns.boxplot(x='Year', y='CMAQ', data=filtered_cmaq_2020)
#sns.boxplot(x='Year', y='CMAQ', data=df_day_2020)

#plt.axvline('January', color='red', linestyle='--', label='January')
plt.title('Comparison of Winter PM₂.₅ Levels (2010-2024)', fontsize = 16)
plt.ylabel('Average PM₂.₅ Concentration (µg/m³)', fontsize= 14)
plt.xlabel('Month', fontsize = 16)
plt.ylim([0,25])
plt.xticks(fontsize = 14) 
plt.yticks(fontsize = 14) 

plt.legend()
plt.show()

# %%
import matplotlib.pyplot as plt
import pandas as pd

# Calculate the mean PM2.5 by year and month
monthly_stats = winter_data.groupby(['Year', 'Month'])['Arithmetic Mean'].mean().unstack()

# Plot the mean PM2.5 levels for January, February, and December
plt.figure(figsize=(10, 6))
plt.plot(monthly_stats.index, monthly_stats['January'], marker='o', label='January')
plt.plot(monthly_stats.index, monthly_stats['February'], marker='s', label='February')
plt.plot(monthly_stats.index, monthly_stats['December'], marker='v', label='December')
# plt.plot(cmaq_baseline_stats_night.index, cmaq_baseline_stats_night['mean'], marker = 'o', label='Baseline Night modeled PM2.5')
# plt.plot(cmaq_baseline_stats_day.index, cmaq_baseline_stats_day['mean'], marker = 'o', label='Baseline Day modeled PM2.5')

# Highlight January 2016 with a vertical line
#plt.axvline(2016, color='red', linestyle='--', label='January 2016')

plt.title('Comparing PM₂.₅ Levels for January, February, and December by Year', fontsize = 16)
plt.xlabel('Year', fontsize = 14)
plt.ylabel('Average PM₂.₅ Concentration (µg/m³)', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.legend(fontsize = 14 )
plt.grid(True)
plt.show()


# %%
import matplotlib.pyplot as plt
import pandas as pd

# Calculate the mean PM2.5 by year and month
monthly_stats = winter_data.groupby(['Year', 'Month'])['Arithmetic Mean'].mean().unstack()

# Plot the mean PM2.5 levels for January, February, and December
plt.figure(figsize=(10, 6))
plt.plot(monthly_stats.index, monthly_stats['January'], marker='o', label='January')
plt.plot(monthly_stats.index, monthly_stats['February'], marker='s', label='February')
plt.plot(monthly_stats.index, monthly_stats['December'], marker='v', label='December')
plt.plot(cmaq_baseline_stats.index, cmaq_baseline_stats['mean'], marker = '*', markersize = 10, label='WRF-CMAQ')
# plt.plot(cmaq_baseline_stats_night.index, cmaq_baseline_stats_night['mean'], marker = 'o', label='Baseline Night modeled PM2.5')
# plt.plot(cmaq_baseline_stats_day.index, cmaq_baseline_stats_day['mean'], marker = 'o', label='Baseline Day modeled PM2.5')

# Highlight January 2016 with a vertical line
#plt.axvline(2016, color='red', linestyle='--', label='January 2016')

plt.title('Comparing PM₂.₅ Levels for January, February, and December by Year', fontsize = 16)
plt.xlabel('Year', fontsize = 14)
plt.ylabel('Average PM₂.₅ Concentration (µg/m³)', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.legend(fontsize = 14 )
plt.grid(True)
plt.show()


# %%


# %%



