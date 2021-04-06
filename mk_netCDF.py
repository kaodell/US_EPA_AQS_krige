#mk_netCDF.py
#	process_ouput from krige_aqs_US.py and calcPMbackground_ns.py, make figures,check stats, and save as netCDF file
#written by: Katelyn O'Dell 08/08/17
# version 2.0: update to python 3, and renamed mk_netCDF, updated figures for plotting
###################################################################################################
#load important modules
import numpy as np
import netCDF4
import datetime as dt
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
year = 2018

# file paths for input data
# output from krige_aqs_US.py
kPM_in_fp = '/home/kaodell/nasa_fires/Lipner_fullUSkrige/datafiles/krigedPM25_06-19_v2/'
# output from calcPMbackground_ns.py
bkPM_in_fp = '/home/kaodell/nasa_fires/Lipner_fullUSkrige/datafiles/krigedPM25_06-19_v2/'
# Hazard Mapping System (HMS) smoke plume polygons. Acess here: https://www.ospo.noaa.gov/Products/land/hms.html#data
# shapefiles regridded to the krigging grid, curtosey of Bonne Ford
hms_fp = '/pierce-scratch/bford/satellite_data/HMS/4kate/'

# file paths for output data and figures
out_fp = '/home/kaodell/nasa_fires/Lipner_fullUSkrige/datafiles/krigedPM25_06-19_v2/'
fig_fp = '/home/kaodell/nasa_fires/Lipner_fullUSkrige/datafiles/krigedPM25_06-19_v2/'

# file descriptions
in_bk_desc = '_v2_statfixmean'
in_kPM_desc = '_v2_statfix'
out_desc = '_v2_statfix_meanbk'

###################################################################################################
# user defined functions
###################################################################################################
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
# make date time array
###################################################################################################
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
# Open datafiles
###################################################################################################
fn1 = kPM_in_fp +'kpm' + str(year) + in_kPM_desc +'.npz'
fn2 = bkPM_in_fp + 'backgroundPM' + str(year) + in_bk_desc +'.npz'

# kriged PM
f1 = np.load(fn1)
kPM = f1['kPM']
glats = f1['glats']
glons = f1['glons']
rsq = f1['rsq']
slope = f1['slope']
MB = f1['MB']
MAE = f1['MAE']
site_rsq = f1['ksRSQ']
site_slope = f1['ksslope']
site_MB = f1['ksMB']
site_MAE = f1['ksMAE']
site_nobs = f1['knobs']
slat = f1['ilat'] # lats of sites within kriging domain
slon = f1['ilon']
sPM = f1['sPMin']
f1.close()

# background estimates
f2 = np.load(fn2)
background_PM = f2['background_PM']
f2.close()

# HMS
file = hms_fp + 'hms_smoke_' + str(year) + '.nc'
nc_fid = netCDF4.Dataset(file, 'r')
HMS_smoke = nc_fid.variables['HMS_Smoke'][:]
nc_fid.close()

# grab and save edges of grid boxes
file = '/fischer-scratch/kaodell/WRF_output/summer16/wrfout_d01_2016-06-23_00:00:00'
nc_fid = netCDF4.Dataset(file, 'r')
we_lats  = nc_fid.variables['XLAT_U'][:]
we_lons = nc_fid.variables['XLONG_U'][:]
ns_lats = nc_fid.variables['XLAT_V'][:]
ns_lons = nc_fid.variables['XLONG_V'][:]
nc_fid.close()

###################################################################################################
# create background array
###################################################################################################
seasons = [[1,2,3],[4,5,6],[7,8,9],[10,11,12]]
months = np.array(months)
# create grid
kPMbackground = np.empty(kPM.shape)
kPMbackground[:] = -999.9
for k in range(0,4):
	dayinds_season = []
	for n in range(len(seasons[k])):		
		dayinds = np.where(months == seasons[k][n])
		dayinds_season = np.hstack([dayinds_season, dayinds[0]])
	dayinds_season = np.array(dayinds_season,dtype=int)
	kPMbackground[dayinds_season,:,:] = background_PM[k,:,:]

###################################################################################################
# Save data as a netCDF
###################################################################################################
outfn = out_fp + 'krigedPM25_' + str(year) + out_desc + '.nc'
nc_w_fid = netCDF4.Dataset(outfn, 'w', format='NETCDF4',clobber=True)
nc_w_fid.description = 'Krigged PM2.5 concentrations from EPA AQS surface monitors for Jan-Dec ' + str(year) + 'created: '+dt.datetime.now().strftime('%Y-%m-%d')

# Need to define dimensions that will be used in the file
nc_w_fid.createDimension('date',NT)
nc_w_fid.createDimension('lonx',189)
nc_w_fid.createDimension('lony',309)

nc_w_fid.createDimension('welonx',189)
nc_w_fid.createDimension('welony',310)

nc_w_fid.createDimension('nslonx',190)
nc_w_fid.createDimension('nslony',309)

nc_w_fid.createDimension('nsites',None)

# define variables for netcdf file
                          #      name  format  dimensions
doy_nc = nc_w_fid.createVariable('doy', 'i4', ('date',))

lon_nc = nc_w_fid.createVariable('lon', 'f8', ('lonx','lony'))

lat_nc = nc_w_fid.createVariable('lat', 'f8', ('lonx','lony'))

we_lon_nc = nc_w_fid.createVariable('we_lon', 'f8', ('welonx','welony'))

we_lat_nc = nc_w_fid.createVariable('we_lat', 'f8', ('welonx','welony'))

ns_lon_nc = nc_w_fid.createVariable('ns_lon', 'f8', ('nslonx','nslony'))

ns_lat_nc = nc_w_fid.createVariable('ns_lat', 'f8', ('nslonx','nslony'))

PM25_nc = nc_w_fid.createVariable('PM25', 'f8', ('date', 'lonx','lony'))

PM25_background_nc = nc_w_fid.createVariable('Background_PM25', 'f8', ('date', 'lonx', 'lony'))

HMS_smoke_nc = nc_w_fid.createVariable('HMS_Smoke', 'f8', ('date', 'lonx', 'lony'))


site_lon_nc = nc_w_fid.createVariable('testing_sites_longitudes', 'f8', ('nsites'))

site_lat_nc = nc_w_fid.createVariable('testing_sites_latitudes', 'f8', ('nsites'))

site_rsq_nc = nc_w_fid.createVariable('r_squared', 'f8', ('nsites'))

site_MB_nc = nc_w_fid.createVariable('mean_bias', 'f8', ('nsites'))

site_MAE_nc = nc_w_fid.createVariable('mean_absolute_error', 'f8', ('nsites'))

site_slope_nc = nc_w_fid.createVariable('slope', 'f8', ('nsites'))

nobs_nc = nc_w_fid.createVariable('nobs', 'f8', ('nsites'))

# add variable attributes 
doy_nc.setncatts({'units':'date','long_name':'day of year',\
               'var_desc':'Day of year for 24-hr avg PM2.5'})

lat_nc.setncatts({'units':'degrees','var_desc':'degrees latitude for data grid centers',\
               'long_name':'Latitude (centers) [degrees]'})
lon_nc.setncatts({'units':'degrees','var_desc':'degrees longitude for data grid centers',\
               'long_name':'Longitude (centers) [degrees]'})
we_lat_nc.setncatts({'units':'degrees','var_desc':'degrees latitude for data grid west-east borders',\
               'long_name':'Latitude (WE border) [degrees]'})
we_lon_nc.setncatts({'units':'degrees','var_desc':'degrees longitude for data grid west-east borders',\
               'long_name':'Longitude (WE border) [degrees]'})
ns_lat_nc.setncatts({'units':'degrees','var_desc':'degrees latitude for data grid north-south borders',\
               'long_name':'Latitude (NS border) [degrees]'})
ns_lon_nc.setncatts({'units':'degrees','var_desc':'degrees longitude for data grid north-south borders',\
               'long_name':'Longitude (NS border) [degrees]'})

PM25_nc.setncatts({'units':'ug/m3','var_desc':'24 hour average PM2.5 concentration',\
               'long_name':'PM2.5 concentration [ug/m3]'})
PM25_background_nc.setncatts({'units':'ug/m3','var_desc':'Seasonal PM25 Background (JFM, AMJ, JAS, OND)',\
               'long_name':'PM2.5 background concentration [ug/m3]'})
HMS_smoke_nc.setncatts({'units':'binary','var_desc':'Binary HMS Smoke: 1 = smoke, 0 = no smoke',\
                        'long_name':'HMS Smoke Plumes'})

site_lat_nc.setncatts({'units':'degrees','long_name':'testing sites latitudes',\
               'var_desc':'latitudes for EPA AQS sites used for kriging LOOCV'})
site_lon_nc.setncatts({'units':'degrees','long_name':'testing sites longitudes',\
               'var_desc':'longitudes for EPA AQS sites used for kriging LOOCV'})

site_rsq_nc.setncatts({'units':'unitless','long_name':'r-squared',\
               'var_desc':'r-squared for LOOCV against sfc monitors'})
site_MB_nc.setncatts({'units':'ug m-3','long_name':'mean bias',\
               'var_desc':'mean bias for LOOCV against sfc monitors'})
site_MAE_nc.setncatts({'units':'ug m-3','long_name':'mean absolute error',\
               'var_desc':'mean absolute error for LOOCV against sfc monitors'})
site_slope_nc.setncatts({'units':'unitless','long_name':'linear regression slope',\
               'var_desc':'slope for LOOCV against sfc monitors'})
nobs_nc.setncatts({'units':'unitless','long_name':'number of site observations',\
               'var_desc':'number of monitor observations available for LOOCV against sfc monitors'})

# assign data to the variables
doy_nc[:] = np.arange(0,kPM.shape[0])
lon_nc[:] = glons
lat_nc[:] = glats
we_lon_nc[:] = we_lons[0,:,:]
we_lat_nc[:] = we_lats[0,:,:]
ns_lon_nc[:] = ns_lons[0,:,:]
ns_lat_nc[:] = ns_lats[0,:,:]
PM25_nc[:] = kPM
PM25_background_nc[:] = kPMbackground
HMS_smoke_nc[:] = HMS_smoke
site_lat_nc[:] = slat
site_lon_nc[:] = slon
site_rsq_nc[:] = site_rsq
site_MB_nc[:] = site_MB
site_MAE_nc[:] = site_MAE
site_slope_nc[:] = site_slope
nobs_nc[:] = site_nobs
nc_w_fid.close()

print('data saved, making figures')
###################################################################################################
# Make figures
###################################################################################################

# 1) plot map of R2, MB, MAE, and slope at sites from LOOCV w title giving overall stats (add slope once we have it)
# only plot sites with over 60 obs
tinds = np.where(site_nobs > 60.0)[0] 

plot_data = [site_rsq[tinds], site_MB[tinds], site_MAE[tinds],site_slope[tinds]]
titles = ['site R2 \n overall: '+str(rsq)[:4],'site MB \n overall: '+str(MB)[:5],'site MAE \n overall: '+str(MAE)[:5],
          'site slope \n overall: '+str(slope)[:4]]
fig, axarr = plt.subplots(nrows=2,ncols=2,subplot_kw={'projection':ccrs.PlateCarree()})
plt.tight_layout()
axlist = axarr.flatten()
# rsq
ax = axlist[0]
mk_map(ax)
cs = ax.scatter(slon[tinds],slat[tinds],c=plot_data[0],s=10,cmap='jet',vmin=0,vmax=1)
cax,kw = mplt.colorbar.make_axes(ax,location='bottom',shrink=0.35)
cbar = fig.colorbar(cs,cax=cax,**kw)
ax.set_title(titles[0])
# mean bias
ax = axlist[1]
norm = mplt.colors.TwoSlopeNorm(vmin=-5, vcenter=0, vmax=5)
mk_map(ax)
cs = ax.scatter(slon[tinds],slat[tinds],c=plot_data[1],s=10,cmap='bwr',norm=norm)
cax,kw = mplt.colorbar.make_axes(ax,location='bottom',shrink=0.35,extend='both')
cbar = fig.colorbar(cs,cax=cax,**kw)
ax.set_title(titles[1])
# mean absolute error
ax = axlist[2]
mk_map(ax)
cs = ax.scatter(slon[tinds],slat[tinds],c=plot_data[2],s=10,vmin=0,vmax=5)
cax,kw = mplt.colorbar.make_axes(ax,location='bottom',shrink=0.35,extend='max')
cbar = fig.colorbar(cs,cax=cax,**kw)
ax.set_title(titles[2])
# slope
ax = axlist[3]
mk_map(ax)
cs = ax.scatter(slon[tinds],slat[tinds],c=plot_data[3],s=10,vmin=0.5,vmax=1.5)
cax,kw = mplt.colorbar.make_axes(ax,location='bottom',shrink=0.35,extend='both')
cbar = fig.colorbar(cs,cax=cax,**kw)
ax.set_title(titles[3])
fig.suptitle('kriging stats ' +str(year))
#plt.tight_layout()
plt.savefig(fig_fp +'kriging_stats_'+ str(year) + out_desc+'.png')

# 2) kPM annual average, annual smoke, annual nosmoke
# some quick calculations (note not removing negative backgrounds or total PM here, which we typically do in the full smoke calc)

kPM_anavg = np.nanmean(kPM,axis=0)
nosmoke = np.where(HMS_smoke==0,kPM,kPMbackground)
nosmoke_anavg = np.nanmean(nosmoke,axis=0)
smoke = np.where(HMS_smoke==1,kPM-kPMbackground,0)
smoke_anavg = np.nanmean(smoke,axis=0)

plot_data = [kPM_anavg,nosmoke_anavg,smoke_anavg]
titles = ['total PM','nosmoke PM', 'smoke PM']
fig, axarr = plt.subplots(nrows=3,ncols=1,subplot_kw={'projection':ccrs.PlateCarree()})
plt.tight_layout()
axlist = axarr.flatten()
for i in range(3):
    ax = axlist[i]
    mk_map(ax)
    cs = ax.pcolormesh(glons,glats,plot_data[i],shading='nearest')
    cax,kw = mplt.colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.65)
    cbar = fig.colorbar(cs,cax=cax,**kw)
    ax.set_title(titles[i])
fig.suptitle('annual mean PM2.5 ' +str(year))
plt.savefig(fig_fp +'kriging_smoke_nosmoke_annavg_'+ str(year) + out_desc+'.png')
fig.show()

# 3-4) kriging one day of each month with sites to check it looks OK
# split into two figures, 6 panels each
# Jan - June
dayinds_plot = []
titles = ['01-01','02-01','03-01', '04-01', '05-01', '06-01']
date = np.array(date)
for ds in titles:
    date_str = str(year) +'-'+ds
    dind = np.where(date==date_str)[0][0]
    dayinds_plot.append(dind)
fig, axarr = plt.subplots(nrows=3,ncols=2,subplot_kw={'projection':ccrs.PlateCarree()})
plt.tight_layout()
axlist = axarr.flatten()
for i in range(6):
    ax = axlist[i]
    mk_map(ax)
    cs = ax.pcolormesh(glons,glats,kPM[dayinds_plot[i],:,:],shading='nearest',cmap='OrRd',
                       vmin=kPM[dayinds_plot[i],:,:].min(),vmax = kPM[dayinds_plot[i],:,:].max())
    ax.scatter(slon,slat,c=sPM[dayinds_plot[i],:],s=5,cmap='OrRd',
               vmin=kPM[dayinds_plot[i],:,:].min(),vmax = kPM[dayinds_plot[i],:,:].max())
    cax,kw = mplt.colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.65)
    cbar = fig.colorbar(cs,cax=cax,**kw)
    ax.set_title(titles[i])
fig.suptitle('daily PM2.5 ' +str(year))
plt.savefig(fig_fp +'kriging_daily_plot_examplesA_'+ str(year) + out_desc+'.png')

# July - December
dayinds_plot = []
titles = ['07-01','08-01', '09-01', '10-01', '11-01', '12-01']
for ds in titles:
    date_str = str(year) +'-'+ds
    dind = np.where(date==date_str)[0][0]
    dayinds_plot.append(dind)
fig, axarr = plt.subplots(nrows=3,ncols=2,subplot_kw={'projection':ccrs.PlateCarree()})
plt.tight_layout()
axlist = axarr.flatten()
for i in range(6):
    ax = axlist[i]
    mk_map(ax)
    cs = ax.pcolormesh(glons,glats,kPM[dayinds_plot[i],:,:],shading='nearest',cmap='OrRd',
                       vmin=kPM[dayinds_plot[i],:,:].min(),vmax = kPM[dayinds_plot[i],:,:].max())
    ax.scatter(slon,slat,c=sPM[dayinds_plot[i],:],s=5, cmap='OrRd',
               vmin=kPM[dayinds_plot[i],:,:].min(),vmax = kPM[dayinds_plot[i],:,:].max())
    cax,kw = mplt.colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.65)
    cbar = fig.colorbar(cs,cax=cax,**kw)
    ax.set_title(titles[i])
fig.suptitle('daily PM2.5 ' +str(year))
plt.savefig(fig_fp +'kriging_daily_plot_examplesB_'+ str(year) + out_desc+'.png')
plt.show()

# Distribution of total PM2.5 and background PM2.5

# remove over ocean values
state_grid = np.load('state_grid.npy') # assigns grid cells to states
US_inds = np.where(state_grid!='NA')
kPM_plot = kPM[:,US_inds[0],US_inds[1]].flatten()
bkPM_plot = kPMbackground[:,US_inds[0],US_inds[0]].flatten()

# plot
fig, axs = plt.subplots(2,1)
ax=axs.flatten()
if np.max(kPM_plot) < 100:
    if np.min(kPM_plot) > -5:
        ax[0].hist(kPM_plot, bins=[np.min(kPM_plot),0,5,10,20,30,40,50,np.max(kPM_plot)])
    else:
        ax[0].hist(kPM_plot, bins=[np.min(kPM_plot),-5,0,5,10,20,30,40,50,np.max(kPM_plot)])
else:
    if np.min(kPM_plot) > -5:
        ax[0].hist(kPM_plot, bins=[np.min(kPM_plot),0,5,10,20,30,40,50,100,np.max(kPM_plot)])
    else:
        ax[0].hist(kPM_plot, bins=[np.min(kPM_plot),-5,0,5,10,20,30,40,50,100,np.max(kPM_plot)])
ax[0].set_xlabel('daily-average PM2.5 [ug m-3]')
ax[0].set_ylabel('count')
ax[0].set_yscale('log')
ax[0].set_title('total PM2.5 ' +str(year) + ' histogram')
if np.min(bkPM_plot)< 0:
    if np.max(bkPM_plot)<12:
        ax[1].hist(bkPM_plot, bins=[np.min(bkPM_plot),0,4,6,8,10,np.max(bkPM_plot)])
    else:
        ax[1].hist(bkPM_plot, bins=[np.min(bkPM_plot),0,4,6,8,10,12,np.max(bkPM_plot)])
elif np.min(bkPM_plot)< 4:
    if np.max(bkPM_plot)<12:
        ax[1].hist(bkPM_plot, bins=[np.min(bkPM_plot),4,6,8,10,np.max(bkPM_plot)])
    else:
        ax[1].hist(bkPM_plot, bins=[np.min(bkPM_plot),4,6,8,10,12,np.max(bkPM_plot)])
else:
    if np.max(bkPM_plot)<12:
        ax[1].hist(bkPM_plot, bins=[np.min(bkPM_plot),6,8,10,np.max(bkPM_plot)])
    else:
        ax[1].hist(bkPM_plot, bins=[np.min(bkPM_plot),6,8,10,12,np.max(bkPM_plot)])
ax[1].set_xlabel('background PM2.5 [ug m-3]')
ax[1].set_ylabel('count')
ax[1].set_title('background PM2.5 ' +str(year) + ' histogram')
plt.tight_layout()
plt.savefig(fig_fp +'kriging_histrogram_'+ str(year) + out_desc+'.png')
plt.show()

# In addition to these plots, it is reccommended to look at the output netCDF file variables using ncview

