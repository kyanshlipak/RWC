import pandas as pd
import geopandas as gpd

# EJ Index csv I was using for some analyses, probably not relevant for others
EJ_index = pd.read_csv('EJ_index.csv')


def import_rwc_shapefile(loc = 'census_tract_data/2016_rwc_census_tract_pm25.shp'):
    """"
	Reads shapefile into geopandas dataframe
	minor formatting to merge EJ_index df and rename/change datatypes
	"2016_rwc_census_tract_pm25.shp" is a shapefile including pm2.5, ACS population characteristics,
	and geographic data
	"""
    rwc_census_tract_baseline = gpd.read_file(loc)
    rwc_census_tract_baseline = rwc_census_tract_baseline.rename(columns={
               "White alon":    "White",
               "Black or A":    "Black",
               "American I":    "American Indian",
               "Asian alon":    "Asian",
               "Native Haw":    "Native Hawaiian or Pacific Islander",
               'Some Other': 'Other', 
               'Two or Mor': 'Two or More Races',
               'Hispanic o': 'Hispanic'
                       })
    
    rwc_census_tract_baseline['GEOID_x'] = rwc_census_tract_baseline['GEOID_x'].astype(int)
    rwc_census_tract_baseline = rwc_census_tract_baseline.merge(EJ_index[['GEOID', 'SPL_EJI']], left_on='GEOID_x', right_on='GEOID', how='left')
    rwc_census_tract_baseline['COUNTYFP'] = rwc_census_tract_baseline['COUNTYFP'].astype(int)
    rwc_census_tract_baseline['STATEFP'] = rwc_census_tract_baseline['STATEFP'].astype(float)
    return rwc_census_tract_baseline

# read the rwc census tract level pollution and population data
rwc_census_tract_2020 = import_rwc_shapefile(loc = 'census_tract_data/2020_rwc_census_tract_pm25.shp')
#rwc_census_tract_baseline = import_rwc_shapefile(loc = 'census_tract_data/2016_rwc_census_tract_pm25.shp')

# Combined statistical area shapefile data
cbsa = gpd.read_file('cbsa/tl_2019_us_cbsa.shp', engine="pyogrio")

# Convert the GeoDataFrame to the target CRS
cbsa = cbsa.to_crs(rwc_census_tract_2020.crs)
cbsa['CBSA'] = cbsa['NAME']

def add_CBSA(rwc_census_tract_baseline, cbsa):
    """
    Add the CBSA information to the rwc census tract geodataframe.
    CBSA info was obtained from the IPUMS NIGHS.
    Uses centroids, but in reality, census tract lines align with conuty --> CBSA
    """
    
    # Convert centroid to GeoDataFrame to use in spatial join
    rwc_census_tract_baseline = rwc_census_tract_baseline.copy(deep = True)
    rwc_census_tract_baseline['centroid'] = rwc_census_tract_baseline.geometry.centroid
    centroids_gdf = rwc_census_tract_baseline.set_geometry('centroid')
    
    
    # Perform spatial join based on centroid within cbsa geometry
    rwc_census_tract_baseline_merged = gpd.sjoin(centroids_gdf, cbsa[['geometry', 'CBSA']], predicate='within', how = 'left')
    
    # Merge the CBSA column back to the original rwc_census_tract_baseline dataframe
    rwc_census_tract_baseline['CBSA'] = rwc_census_tract_baseline_merged['CBSA']
    return rwc_census_tract_baseline

# add CBSA data
rwc_census_tract_2020 = add_CBSA(rwc_census_tract_2020, cbsa)
#rwc_census_tract_baseline = add_CBSA(rwc_census_tract_baseline, cbsa)

#______________________________________________________________________________________

"""
Okay so now we have a gdf with population characteristics, 
census tract shapes, pollution information, and CBSA info.
Next step is to read and clean/process the mortality data
"""

# read our mortality data
mortality_data = pd.read_csv("BenMAP_mortality.csv")



###
### Now we need to make a GISJOIN column to allow for merges with the big census tract gdf
###
mortality_data['Row'] = mortality_data['Row'].astype(str)
mortality_data['Column'] = mortality_data['Column'].astype(str)

# Step 1: Format 'Row' to be a 5-digit County FIPS code
mortality_data['County_FIPS'] = mortality_data['Row'].astype(str).str.zfill(5)

# Step 2: Format 'Column' to be a 6-digit Tract FIPS code
mortality_data['Tract_FIPS'] = mortality_data['Column'].astype(str).str.zfill(6)

# Step 3: Derive State FIPS Code (first 2 digits of the Tract FIPS code)
mortality_data['State_FIPS'] = mortality_data['County_FIPS'].str[:2]
mortality_data['County_FIPS'] = mortality_data['County_FIPS'].str[2:]

# Step 4: Create GISJOIN by concatenating 'G' + State FIPS + County FIPS + Tract FIPS + Suffix
mortality_data['GISJOIN'] = 'G' + mortality_data['State_FIPS'] + '0' +  mortality_data['County_FIPS'] + '0' + mortality_data['Tract_FIPS']


### Next we  need to deal with the age ranges
# 1: Make age range table
mortality_data["Age Range"] = mortality_data["Start Age"].astype(str) + "-" + mortality_data["End Age"].astype(str)

# 2: Pivot the table
pivoted_df = mortality_data.pivot_table(index=["GISJOIN", "Endpoint"],
                            columns="Age Range", values="Value").reset_index()

# 3: Aggregate age ranges of mortality data to match ACS age information
age_mapping = {
    "0-0": ["0-5"],
    "1-4": ["0-5"],
    "5-14": ["5-9", "10-14"],
    "15-24": ["15-17", "18-19", "20-21", "22-24"],
    "25-34": ["25-29", "30-34"],
    "35-44": ["35-39", "40-44"],
    "45-54": ["45-49", "50-54"],
    "55-64": ["55-59", "60-64"],
    "65-74": ["65-69", "70-74"],
    "75-84": ["75-79", "80-84"],
    "85-99": ["85+"]
}


# sum columns based on the mapping
mapped_df = pd.DataFrame()
mapped_df['GISJOIN'] = pivoted_df['GISJOIN']

for mapped_range, age_ranges in age_mapping.items():
    # Sum the columns corresponding to the age ranges in the mapping
    for col in age_ranges:
        mapped_df[col] = pivoted_df[mapped_range]

# Add non-age columns (GISJOIN, Endpoint) to the new DataFrame
mapped_df['Endpoint'] = pivoted_df['Endpoint']
mortality_data = mapped_df


#______________________________________________________________________________________

"""
Process age data from the ACS to match health dataframe
"""

age_data = pd.read_csv("nhgis_ages/nhgis0006_ds244_20195_tract.csv")

# Define age range groups for summing, with column ranges for males and females
age_ranges = {
    '0-5': ['ALT0M003', 'ALT0M027'],
    '5-9': ['ALT0M004', 'ALT0M028'],
    '10-14': ['ALT0M005', 'ALT0M029'],
    '15-17': ['ALT0M006', 'ALT0M030'],
    '18-19': ['ALT0M007', 'ALT0M031'],
    '20-21': ['ALT0M008', 'ALT0M009', 'ALT0M032', 'ALT0M033'],
    '22-24': ['ALT0M010', 'ALT0M034'],
    '25-29': ['ALT0M011', 'ALT0M035'],
    '30-34': ['ALT0M012', 'ALT0M036'],
    '35-39': ['ALT0M013', 'ALT0M037'],
    '40-44': ['ALT0M014', 'ALT0M038'],
    '45-49': ['ALT0M015', 'ALT0M039'],
    '50-54': ['ALT0M016', 'ALT0M040'],
    '55-59': ['ALT0M017', 'ALT0M041'],
    '60-64': ['ALT0M018', 'ALT0M019', 'ALT0M042', 'ALT0M043'],
    '65-69': ['ALT0M020', 'ALT0M021', 'ALT0M044', 'ALT0M045'],
    '70-74': ['ALT0M022', 'ALT0M046'],
    '75-79': ['ALT0M023', 'ALT0M047'],
    '80-84': ['ALT0M024', 'ALT0M048'],
    '85+': ['ALT0M025', 'ALT0M049']
}

# Calculate total population for each age range
for age_range, columns in age_ranges.items():
    age_data[age_range] = age_data[columns].sum(axis=1)

# Drop the original columns (the ones that go ALTOM0...)
age_data = age_data.drop(columns=[col for cols in age_ranges.values() for col in cols])

# Get rid of other columns we don't need (sorry there might be a better way to do this)
for col in [
    "YEAR", "STUSAB", "REGIONA", "DIVISIONA", "STATE", "STATEA", "COUNTY", 
    "COUNTYA", "COUSUBA", "PLACEA", "TRACTA", "BLKGRPA", "CONCITA", "AIANHHA", 
    "RES_ONLYA", "TRUSTA", "AIHHTLI", "AITSCEA", "ANRCA", "CBSAA", "CSAA", 
    "METDIVA", "NECTAA", "CNECTAA", "NECTADIVA", "UAA", "CDCURRA", "SLDUA", 
    "SLDLA", "ZCTA5A", "SUBMCDA", "SDELMA", "SDSECA", "SDUNIA", "PCI", 
    "PUMAA", "GEOID", "BTTRA", "BTBGA", "NAME_E", "ALT0E001", "ALT0E002", 
    "ALT0E003", "ALT0E004", "ALT0E005", "ALT0E006", "ALT0E007", "ALT0E008", 
    "ALT0E009", "ALT0E010", "ALT0E011", "ALT0E012", "ALT0E013", "ALT0E014", 
    "ALT0E015", "ALT0E016", "ALT0E017", "ALT0E018", "ALT0E019", "ALT0E020", 
    "ALT0E021", "ALT0E022", "ALT0E023", "ALT0E024", "ALT0E025", "ALT0E026", 
    "ALT0E027", "ALT0E028", "ALT0E029", "ALT0E030", "ALT0E031", "ALT0E032", 
    "ALT0E033", "ALT0E034", "ALT0E035", "ALT0E036", "ALT0E037", "ALT0E038", 
    "ALT0E039", "ALT0E040", "ALT0E041", "ALT0E042", "ALT0E043", "ALT0E044", 
    "ALT0E045", "ALT0E046", "ALT0E047", "ALT0E048", "ALT0E049", "NAME_M", 
    "ALT0M001", "ALT0M002", "ALT0M026"
]:
    del age_data[col]


#______________________________________________________________________________________

"""
Merge age and mortality data to calculate baseline all-cause mortality by census tract
"""

merged_data = pd.merge(age_data, mortality_data, on='GISJOIN', suffixes=('_population', '_mortality'))

###
###  Compute mortality * population for each age bin
###
age_bins = ['0-5', '5-9', '10-14', '15-17', '18-19', '20-21', '22-24',
            '25-29', '30-34', '35-39', '40-44', '45-49', '50-54', '55-59',
            '60-64', '65-69', '70-74', '75-79', '80-84']#, '85+']

# Initialize a list to store the results
mortality_rates = []

# Iterate over each GISJOIN census tract
for _, row in merged_data.iterrows():
    total_population = row[[col + '_population' for col in age_bins]].sum()  # Sum of population for all age bins. ig this is ineffiicent but im not changing it now :(
    total_mortality = sum(row[age_bin + '_mortality'] * row[age_bin + '_population'] for age_bin in age_bins)  # Sum of mortality * population for each tract
    overall_mortality_rate = total_mortality / total_population if total_population != 0 else 0  # Avoid division by zero (not sure if there are any zero pop, but leavin it)
    mortality_rates.append(overall_mortality_rate) # 

# Add the calculated mortality rates to the DataFrame
merged_data['overall_mortality_rate'] = mortality_rates


###
###  Compute mortality * population for each age bin > 25
###  for my specific concentration response function
###

# Age bins to consider for age > 25
age_bins_over_25 = ['25-29', '30-34', '35-39', '40-44', '45-49', '50-54', '55-59',
                    '60-64', '65-69', '70-74', '75-79', '80-84']#,'85+']

# Initialize a list to store the results
mortality_rates_over_25 = []

# Iterate over each GISJOIN census tract
for _, row in merged_data.iterrows():
    total_population_over_25 = row[[col + '_population' for col in age_bins_over_25]].sum()  # Sum of population for age bins > 25
    total_mortality_over_25 = sum(row[age_bin + '_mortality'] * row[age_bin + '_population'] for age_bin in age_bins_over_25)  # Sum of mortality * population for age bins > 25
    overall_mortality_rate_over_25 = total_mortality_over_25 / total_population_over_25 if total_population_over_25 != 0 else 0  # Avoid division by zero
    mortality_rates_over_25.append(overall_mortality_rate_over_25)

# Add the calculated mortality rates for age > 25 to the DataFrame
merged_data['overall_mortality_rate_over_25'] = mortality_rates_over_25


#______________________________________________________________________________________

"""
We have "merged_data" with all-cause mortality information by census tract
and we have gdf with population characteristics, 
census tract shapes, pollution information, and CBSA info.
Put it all together!
"""

# merge population & pollution information with mortality information
pm_mortality_data_df_2020 =  pd.merge(rwc_census_tract_2020, mortality_data, on = "GISJOIN")

# calculations!
def mortality_function(pm_mortality_data_df, relative_risk = 1.08):
	pm_mortality_data_df = pm_mortality_data_df.copy(deep = True) # full copy so changes don't happen to the inputted df

	# divide by four because only three months of the year
	# 1.08 is per 10µg/m^3 so we divide by 10
	# real realtive risk as opposed to the relative risk function(concentration)
	pm_mortality_data_df['real_relative_risk'] = relative_risk ** (pm_mortality_data_df['2016_pm2.5'] / 10 /4) 

	# attributable fraction from relative risk
	pm_mortality_data_df['AF'] = (pm_mortality_data_df['real_relative_risk'] - 1) / pm_mortality_data_df['real_relative_risk']

	# attributable mortality using (tract population) * (tract all cause mortality rate over 25) * (RWC PM attributable fraction)
	pm_mortality_data_df['Attributable_Mortality'] = pm_mortality_data_df['Total'] * pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']

	# rate (population not included)
	pm_mortality_data_df['Attributable_Mortality_Rate'] = pm_mortality_data_df['overall_mortality_rate_over_25'] * pm_mortality_data_df['AF']

	# printing total estimated attributable mortalities for domain (CONUS) for the given RR
	print("RR:", relative_risk, " - ", round(pm_mortality_data_df['Attributable_Mortality'].sum()))
	return pm_mortality_data_df, pm_mortality_data_df[['Attributable_Mortality', 'Attributable_Mortality_Rate']]



# Calculate mortalities for RR and 95% confidence intervals
pm_mortality_data_df_2020, _ = mortality_function(pm_mortality_data_df_2020, relative_risk = 1.08)
_, pm_mortality_data_df_2020[['Attributable_Mortality_lower', 'Attributable_Mortality_Rate_lower']] = mortality_function(pm_mortality_data_df_2020, relative_risk  = 1.06)
_, pm_mortality_data_df_2020[['Attributable_Mortality_upper', 'Attributable_Mortality_Rate_upper']] = mortality_function(pm_mortality_data_df_2020, relative_risk  = 1.09)
_, _ = mortality_function(pm_mortality_data_df_2020, relative_risk = 1.17) # just to try






