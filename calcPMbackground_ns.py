# calcPMbackground_ns.py
#	A python script to calculate the background PM concentration using a seasonal average of days with no 'smoke' according
#	to HMS. 
# written by: Katelyn O'Dell
# 08/08/17
# version 2: 12/16/20, updated to python 3 and changed file i/o to load the npz files instead of netCDF
#			this version comments out the site background check, but it should be uncommented 
#			for running new years. Plotting methods need to be updated to cartopy. also
#                       changed seasons to JFM, AMJ, JAS, OND (from DJF, MAM, JJA, SON). renamed calcPMbackground_ns.py
###################################################################################################
#load important modules
import numpy as np
import pylab as pl
import netCDF4
import csv
import datetime as dt
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mplt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.feature import NaturalEarthFeature
###################################################################################################
# user inputs
###################################################################################################
# set year
year = 2009

# select background averaging method
avg_method = 'mean' 
#avg_method = 'median' # used in publications

# specify file locations
# location of folder with sPM npz files and kPM npz files output from process_sfc_obs_update.py and kriqe_aqs_US.py, respectively
in_fp = '/home/kaodell/nasa_fires/Lipner_fullUSkrige/datafiles/krigedPM25_06-19_v2/'
# where to put output file
out_fp = '/home/kaodell/nasa_fires/Lipner_fullUSkrige/datafiles/krigedPM25_06-19_v2/'
# where to put output figures
fig_fp = '/home/kaodell/nasa_fires/Lipner_fullUSkrige/datafiles/krigedPM25_06-19_v2/'

# Hazard Mapping System (HMS) smoke plume polygons. Acess here: https://www.ospo.noaa.gov/Products/land/hms.html#data
# shapefiles regridded to the krigging grid, curtosey of Bonne Ford
HMS_fp = '/pierce-scratch/bford/satellite_data/HMS/4kate/'

# choose version of kPM and sPM data to load
kPM_file_desc = '_v2_statfix'
sPM_file_desc = '_v2'

# set description of out files and figures
out_desc = '_v2_statfix'+avg_method

###################################################################################################
# user-defined functions
###################################################################################################
def haversine(lon0,lon1,lat0,lat1):
    # from Will Lassman
    r = 6371000. #m
    lon0 = lon0*np.pi/180

    lon1 = lon1*np.pi/180

    lat0 = lat0*np.pi/180

    lat1 = lat1*np.pi/180

    return 2*r*np.arcsin(np.sqrt(np.sin((lat1 - lat0)/2.)**2 +\
		 np.cos(lat0)*np.cos(lat1)*np.sin((lon1 - lon0)/2.)**2))

# for plotting
def mk_map(ax):
    ax.patch.set_visible(False)
    shapename = 'admin_1_states_provinces_lakes_shp'
    states_shp = shpreader.natural_earth(resolution='110m',category='cultural', name=shapename)
    facecolor = (0, 0, 0, 0)
    edgecolor = 'black'
    for state in shpreader.Reader(states_shp).geometries():
        ax.add_geometries([state], ccrs.PlateCarree(),
                      facecolor=facecolor, edgecolor=edgecolor)
###################################################################################################
# make array of days
###################################################################################################
# make date time this year
t0 = dt.datetime( year = year, month = 1, day = 1 )
tf = dt.datetime( year = year, month = 12, day = 31 )
NT = (tf - t0).days + 1
alldays = np.empty( NT, dtype = object )
date = ['na']*( NT)
months = ['na']*( NT)
alldays[0] = t0
date[0] = str(alldays[0])[:10]
months[0] = 1
for i in range(1, NT):
	alldays[i] = alldays[i-1] + dt.timedelta(days = 1)
	date[i] = str(alldays[i])[:10]
	months[i] = alldays[i].month

###################################################################################################
# Load Data
###################################################################################################
# kPM file out from krige_aqs_US.py
fn1 = in_fp + 'kpm' +str(year)+ kPM_file_desc + '.npz'
f1 = np.load(fn1)
kPM = f1['kPM']
glats = f1['glats']
glons = f1['glons']
f1.close()

# sPM output from process_sfc_obs_update.py
fn2 = in_fp + 'sfc_pm' + str(year)[2:] + sPM_file_desc + '.npz'
f2 = np.load(fn2)
sPM = f2['sPM']
slon = f2['lons']
slat = f2['lats']
f2.close()

# HMS smoke plume polygons re-gridded to the kriging grid curtosey of Bonne Ford
file = HMS_fp + 'hms_smoke_' + str(year) + '.nc'
nc_fid = Dataset(file, 'r')
#print nc_fid.variables.keys()
HMS_smoke = nc_fid.variables['HMS_Smoke'][:]
nc_fid.close()

###################################################################################################
# Calculate Background PM (seasonal median or mean)
###################################################################################################
# we're going to do a seasonal background average: [JFM AMJ JAS OND]
seasons = [[1,2,3],[4,5,6],[7,8,9],[10,11,12]]
months = np.array(months)

# make arrays to fill
Nlon,Nlat = glons.shape
background_PM = np.zeros([4,Nlon,Nlat])
background_PM[:] = -999.9
season_smoky_days = np.zeros([4,Nlon,Nlat])
season_smoky_days[:] = -999.9
numdays_nosmoke = np.zeros([4,Nlon,Nlat])
numdays_nosmoke[:] = -999.9

for k in range(0,4):
    dayinds_season = []
    for n in range(len(seasons[k])):
        dayinds = np.where(months == seasons[k][n])
        dayinds_season = np.hstack([dayinds_season, dayinds[0]])
    dayinds_season = np.array(dayinds_season,dtype=int)
    for i in range(Nlon):
        for j in range(Nlat):
            no_smoke_PM_inds = np.where(HMS_smoke[dayinds_season,i,j] == 0)
            no_smoke_PM = kPM[dayinds_season[no_smoke_PM_inds[0]],i,j]
            if avg_method == 'median':
                #print('taking seasonal median')
                background_PM[k,i,j] = np.median(no_smoke_PM)
            elif avg_method == 'mean':
                #print('taking seasonal mean')
                background_PM[k,i,j] = np.mean(no_smoke_PM)
            if len(no_smoke_PM) < 20:
                # check that there are sufficient days for the background estimate
                # still assign background here, but check figures to make sure it looks reasonable
                print('few nosmoke days at inds:' + str(k) + ' ' + str(i) + ' '+ str(j))
                print('estimated background: ',background_PM[k,i,j])
                print('# nosmoke days: ',len(no_smoke_PM))
            #also find monthly sum of smoky days
            season_smoky_days[k,i,j] = len(dayinds_season) - len(no_smoke_PM)
            numdays_nosmoke[k,i,j] = len(no_smoke_PM)

# find sum of smoky days in each grid cell
# simply sum along time axis of HMS_smoke
smoky_days = np.sum(HMS_smoke,axis=0)

print('background max: ',background_PM.max())
print('background min: ',background_PM.min())
# note, sometimes minimum is negative. I make any < 0, 0 in my analysis. It doesn't make a difference because there are few negative points.
print('background median: ',np.median(background_PM))
print('background mean: ',background_PM.mean())

###################################################################################################
# Save Background Data
###################################################################################################
fn = out_fp + 'backgroundPM' + str(year) + out_desc + '.npz'
np.savez(fn, background_PM=background_PM, glons=glons,glats=glats,season_smoky_days=season_smoky_days, smoky_days=smoky_days,
         numdays_nosmoke=numdays_nosmoke)
print('data saved, making figures')

###################################################################################################
# Make Figures
###################################################################################################

# calculate background at sPM sites to overlay in a figure
background_PM_sites = np.zeros([4,len(slon)])
background_PM_sites[:] = -999.9
season_smoky_days_sites = np.zeros([4,len(slon)])
season_smoky_days_sites[:] = -999.9
numdays_nosmoke_sites = np.zeros([4,len(slon)])
numdays_nosmoke_sites[:] = -999.9
nobs_nosmoke_sites = np.zeros([4,len(slon)])
nobs_nosmoke_sites[:] = -999.9

for k in range(0,4):
    dayinds_season = []
    for n in range(3):
        dayinds = np.where(months == seasons[k][n])
        dayinds_season = np.hstack([dayinds_season, dayinds[0]])
        dayinds_season = np.array(dayinds_season,dtype=int)
        for i in range(len(slon)):
            HMSdists = haversine(glons, slon[i], glats, slat[i])
            HMSind = np.where(HMSdists == HMSdists.min())
            no_smoke_PM_inds = np.where(HMS_smoke[dayinds_season,HMSind[0],HMSind[1]] == 0)
            no_smoke_PM = sPM[dayinds_season[no_smoke_PM_inds[0]],i]
            ntruevals = len(np.where(np.isnan(no_smoke_PM)==False)[0])
            if avg_method == 'median':
                background_PM_sites[k,i] = np.nanmedian(no_smoke_PM)
            elif avg_method == 'mean':
                background_PM_sites[k,i] = np.nanmean(no_smoke_PM)
            if ntruevals < 0.1*(len(no_smoke_PM)):
                background_PM_sites[k,i] = np.nan
            #also find monthly sum of smoky days
            season_smoky_days_sites[k,i] = len(dayinds_season) - len(no_smoke_PM)
            numdays_nosmoke_sites[k,i] = len(no_smoke_PM)
            nobs_nosmoke_sites[k,i] = ntruevals

season_abr = ['JFM','AMJ','JAS','OND']
#1) plot seasonal background PM
fig, axarr = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection':ccrs.PlateCarree()})
plt.tight_layout()
axlist = axarr.flatten()
for i in range(4):
    ax = axlist[i]
    mk_map(ax)
    cs = ax.pcolormesh(glons,glats,background_PM[i,:,:],cmap = 'OrRd',shading='nearest')
    cax,kw = mplt.colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.65)
    cbar = fig.colorbar(cs,cax=cax,**kw)
    ax.set_title(season_abr[i])
fig.suptitle('background PM2.5 ' +str(year))
pl.savefig(fig_fp +'bkPM_'+ str(year) + out_desc+'.png')

#2) plot seasonal background PM with overlayed site background
fig, axarr = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection':ccrs.PlateCarree()})
plt.tight_layout()
axlist = axarr.flatten()
for i in range(4):
    ax = axlist[i]
    mk_map(ax)
    cs = ax.pcolormesh(glons,glats,background_PM[i,:,:],cmap = 'OrRd',shading='nearest',vmin=0,vmax=25)
    tinds = np.where(np.isfinite(background_PM_sites[i,:]))[0]
    ninds = np.where(np.isnan(background_PM_sites[i,:]))[0]
    ax.scatter(slon[tinds],slat[tinds],c=background_PM_sites[i,tinds],s=10,vmin=0,vmax=25,cmap='OrRd')
    #ax.scatter(slon[ninds],slat[ninds],c='k',s=10)# make nan sites black
    cax,kw = mplt.colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.65)
    cbar = fig.colorbar(cs,cax=cax,**kw)
    ax.set_title(season_abr[i])
fig.suptitle('background PM2.5 ' +str(year))
pl.savefig(fig_fp +'bkPMwsites_'+ str(year) + out_desc+'.png')
pl.show()

#3) plot season smoky days
fig, axarr = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection':ccrs.PlateCarree()})
plt.tight_layout()
axlist = axarr.flatten()
for i in range(4):
    ax = axlist[i]
    mk_map(ax)
    cs = ax.pcolormesh(glons,glats,season_smoky_days[i,:,:],cmap = 'Greys',shading='nearest')
    cax,kw = mplt.colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.65)
    cbar = fig.colorbar(cs,cax=cax,**kw)
    ax.set_title(season_abr[i])
fig.suptitle('# of smoke days ' +str(year))
pl.savefig(fig_fp +'season_smokey_days_'+ str(year) + out_desc+'.png')

#4) check number of days feeding nosmoke estimates
fig, axarr = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection':ccrs.PlateCarree()})
plt.tight_layout()
axlist = axarr.flatten()
for i in range(4):
    ax = axlist[i]
    mk_map(ax)
    cs = ax.pcolormesh(glons,glats,numdays_nosmoke[i,:,:],shading='nearest')
    cax,kw = mplt.colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.65)
    cbar = fig.colorbar(cs,cax=cax,**kw)
    ax.set_title(season_abr[i])
fig.suptitle('# of no-smoke days ' +str(year))
pl.savefig(fig_fp +'num_nosmoke_days_'+ str(year) + out_desc+'.png')
pl.show()

