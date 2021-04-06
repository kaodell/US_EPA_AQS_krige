#process_sfc_obs_update.py
#	pull in csv files of PM2.5 data from the AQS and save as array of daily PM for each monitor
#written by: Kate O'Dell 
#version history: 1.0 initial 09/30/16 (version in krige code)
#		  2.0 version with CDPHE added 05/11/18
#		  3.0 updated version to account for double-counted sites 
#		  4.0 updated to python3
#		  4.1 updated monID and siteID to strings
#                 4.2 added check for duplicate match between 'None' and 'Included' event labels
###################################################################################################
# load modules
import datetime as dt
import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as pl
import cartopy as cartopy
###################################################################################################
# set year
year = 2017
# define minimum observation pct for each day
min_obs_pct = 25.0 #%
# aqs file location
aqs_fp = '/fischer-scratch/kaodell/surface_obs/01-15-21dwnload/'
# out file location and description
out_fp = '/home/kaodell/nasa_fires/Lipner_fullUSkrige/datafiles/krigedPM25_06-19_v2/'
out_desc = '_v2'

###################################################################################################
#Time period of interest, this has to be full year. Code won't adjust for partial years.
###################################################################################################
t0 = dt.datetime( year = year, month = 1, day = 1 )
tf = dt.datetime( year = year, month = 12, day = 31 )
NT = (tf - t0).days
alldays = np.empty( NT+1, dtype = object )
alldays[0] = t0
for i in range(1, NT+1):
	alldays[i] = alldays[i-1] + dt.timedelta(days = 1)

###################################################################################################
# Identify files to read in
###################################################################################################
# files can be downloaded here: https://aqs.epa.gov/aqsweb/airdata/download_files.html
# it is important to keep track of the file download date because they are updated twice a year

# load the FRM (88101) and nFRM (88502) data
fn1 = aqs_fp + 'daily_88101_' + str(year) +'.csv'
fn2 = aqs_fp + 'daily_88502_' + str(year) + '.csv'

###################################################################################################
# Load in data, add site ID and monitor ID
###################################################################################################
print('loading AQS files ',year)
PM_df_FRM_all = pd.read_csv(fn1)
PM_df_nFRM_all = pd.read_csv(fn2)

# add a column for Monitor ID (site explaining IDs: https://aqs.epa.gov/aqsweb/airdata/FileFormats.html#_monitors)
# parameter code (we use 88101 and 88502) is part of the ID
# 2005 & 2006 files have data from Canada with a state code of 'CC', let's change to 99
# file 1
if year in [2005,2006]:
	CanInds = np.where(PM_df_FRM_all['State Code'] == 'CC')[0]
	PM_df_FRM_all['State Code'][CanInds] = 99
	PM_df_FRM_all['State Code'] = PM_df_FRM_all['State Code'].astype('str')

# make monitor identifiers strings with leading zeros
PM_df_FRM_all['State Code'] = PM_df_FRM_all['State Code'].apply(lambda x: '{0:0>2}'.format(x))
PM_df_FRM_all['County Code'] = PM_df_FRM_all['County Code'].apply(lambda x: '{0:0>3}'.format(x))
PM_df_FRM_all['Site Num'] = PM_df_FRM_all['Site Num'].apply(lambda x: '{0:0>4}'.format(x))
PM_df_FRM_all['POC'] = PM_df_FRM_all['POC'].apply(lambda x: '{0:0>2}'.format(x))

monIDs_FRM1 = PM_df_FRM_all['State Code'] + '-' + PM_df_FRM_all['County Code'] + '-' + PM_df_FRM_all['Site Num'] +'-'+ '88101' + '-' + PM_df_FRM_all['POC'] 
monIDs_FRM = np.array(monIDs_FRM1,dtype=str)
PM_df_FRM_all['Monitor ID'] = monIDs_FRM

siteIDs_FRM1 = PM_df_FRM_all['State Code']+'-' + PM_df_FRM_all['County Code'] +'-'+ PM_df_FRM_all['Site Num']
siteIDs_FRM = np.array(siteIDs_FRM1,dtype=str)
PM_df_FRM_all['Site ID'] = siteIDs_FRM

# file 2
if year in [2005,2006]:
	CanInds = np.where(PM_df_nFRM_all['State Code'] == 'CC')[0]
	PM_df_nFRM_all['State Code'][CanInds] = 99
	PM_df_nFRM_all['State Code'] = PM_df_nFRM_all['State Code'].astype('str')

# make monitor identifiers strings with leading zeros
PM_df_nFRM_all['State Code'] = PM_df_nFRM_all['State Code'].apply(lambda x: '{0:0>2}'.format(x))
PM_df_nFRM_all['County Code'] = PM_df_nFRM_all['County Code'].apply(lambda x: '{0:0>3}'.format(x))
PM_df_nFRM_all['Site Num'] = PM_df_nFRM_all['Site Num'].apply(lambda x: '{0:0>4}'.format(x))
PM_df_nFRM_all['POC'] = PM_df_nFRM_all['POC'].apply(lambda x: '{0:0>2}'.format(x))

monIDs_nFRM1 = PM_df_nFRM_all['State Code'] + '-' + PM_df_nFRM_all['County Code'] + '-' + PM_df_nFRM_all['Site Num'] +'-'+ '88502' + '-' + PM_df_nFRM_all['POC'] 
monIDs_nFRM = np.array(monIDs_nFRM1,dtype=str)
PM_df_nFRM_all['Monitor ID'] = monIDs_nFRM

siteIDs_nFRM1 = PM_df_nFRM_all['State Code']+'-' + PM_df_nFRM_all['County Code'] +'-'+ PM_df_nFRM_all['Site Num']
siteIDs_nFRM = np.array(siteIDs_nFRM1,dtype=str)
PM_df_nFRM_all['Site ID'] = siteIDs_nFRM


# now combine arrays
PM_df_all = pd.concat([PM_df_FRM_all,PM_df_nFRM_all],ignore_index=True)
del PM_df_FRM_all, PM_df_nFRM_all
###################################################################################################
# remove observations with events excluded (event flag = Excluded) and low observation %
###################################################################################################
# remove days based on event flag
#	these are either 'none' (the ones we keep), 'included' a copy of the 'none' line (also keep incase
#	'none' isn't always copied) indicating the PM with the exceptional event included, 
#	and 'excluded' giving the PM with the exceptional event not included
#	(so these values are a bit lower typically and we don't want to average them in)
exl_events_inds = np.where(PM_df_all['Event Type'] == 'Excluded')[0]
n_obs_exl_events = len(exl_events_inds)
# drop the excluded events from file and re-lable indicies
PM_df1 = PM_df_all.drop(exl_events_inds)
PM_df1 = PM_df1.reset_index(drop=True)

# also remove days where %obs is less than specified amount
low_pct_obs_inds = np.where(PM_df1['Observation Percent'] < min_obs_pct)[0]
n_obs_low_pct = len(low_pct_obs_inds)
PM_df = PM_df1.drop(low_pct_obs_inds)
PM_df = PM_df.reset_index(drop=True)

print('low % obs removed: ',n_obs_low_pct, ', excluded events obs removed: ',n_obs_exl_events,', sum: ',n_obs_exl_events + n_obs_low_pct)
print('total obs before: ',len(PM_df_all),', after: ',len(PM_df),', difference: ', len(PM_df_all) - len(PM_df))
print('check all excluded events removed (should only say included and none here):',np.unique(PM_df['Event Type']))
###################################################################################################
# check that duplicate observations in each file actually are duplicates (within rounding errors)
###################################################################################################
# file 1
print('checking that duplicate obs match (e.g. - same monID and day) before avg')
mon_count = 0
diffs = []
for ID in np.unique(PM_df['Monitor ID']):
        mon_count +=1
        ID_inds = np.where(PM_df['Monitor ID']==ID)[0]
        dates = PM_df['Date Local'].iloc[ID_inds].values
        for date in np.unique(dates):
                date_inds = np.where(date==dates)[0]
                if len(date_inds) > 1:
                        dup_pm = []
                        #print(len(date_inds))
                        for dind in date_inds:
                                dup_pm.append(PM_df['Arithmetic Mean'].iloc[ID_inds[dind]])
                        if np.any(np.diff(dup_pm) > 0.1): # cut-off at 0.1 because this is the lowest sig figs in the AQS file
                                print('duplicate mismatch',dup_pm)
                        if np.any(np.diff(dup_pm) < -0.1):
                                print('duplicate mismatch',dup_pm)
        
###################################################################################################
# get table of date vs ID filled in with PM2.5 observations for each file
###################################################################################################
print('making sPM matrix')
# use pivot table, will average places where there is a duplicate ID/date combo
#	duplicates are due to:
#		1) flagged events 
#		2) different avg method calculated by EPA (typically same or very similar value, rounded differently)
#		see: https://aqs.epa.gov/aqsweb/airdata/FileFormats.html#_daily_summary_files
sPM = PM_df.pivot_table(index = 'Date Local', columns = 'Monitor ID', values = 'Arithmetic Mean', aggfunc = 'mean')
# get lat lons for the monitor IDs
slons = []
slats = []
sIDs = []
for ID in sPM.columns:
		PM_IDloc = np.where(PM_df['Monitor ID']==ID)[0][0]
		slons.append( PM_df.loc[PM_IDloc,['Longitude']][0])
		slats.append( PM_df.loc[PM_IDloc,['Latitude']][0])
		sIDs.append( PM_df.loc[PM_IDloc, ['Site ID']][0])
slons = np.array(slons)
slats = np.array(slats)
sIDs = np.array(sIDs,dtype='str')
monitor_IDs = sPM.columns
monitor_IDs = np.array(monitor_IDs,dtype='str')
sPM = np.array(sPM)
# check all IDs got a row
print('# unique IDs: ',mon_count,', sPM shape: ', sPM.shape)

###################################################################################################
# save ouput
###################################################################################################
out_fn = out_fp + 'sfc_pm' + str(year)[2:] + out_desc +'.npz'
np.savez(out_fn, sPM = sPM, lons = slons, lats = slats, monIDs = monitor_IDs, siteIDs = sIDs)
print('file saved, making plots')

###################################################################################################
# plot sites in December to check for file completeness (if there are there sites missing for an
# entire state .. red flag that the data is incomplete) states have 6 months to report data to the EPA but can take longer for nFRM data
###################################################################################################
mean_dec_obs = np.nanmean(sPM[-32:,:],axis=0)
tinds = np.where(np.isnan(mean_dec_obs)==False)
def plot_background(ax):
    ax.set_extent([235., 290., 23., 50.])
    ax.add_feature(cartopy.feature.COASTLINE.with_scale('50m'), linewidth=0.5)
    ax.add_feature(cartopy.feature.STATES, linewidth=0.5)
    ax.add_feature(cartopy.feature.BORDERS, linewidth=0.5)

fig, ax = pl.subplots(nrows=1, ncols=1, figsize=(20, 13), constrained_layout=True,
                      subplot_kw={'projection': cartopy.crs.PlateCarree()})
plot_background(ax)
c=ax.scatter(slons[tinds], slats[tinds], c=mean_dec_obs[tinds], transform=cartopy.crs.PlateCarree())
ax.set_title('Mean December PM2.5 ' + str(year), fontsize=16)
cb1 = fig.colorbar(c, ax=ax, orientation='horizontal', shrink=0.74, pad=0)
fig.show()
