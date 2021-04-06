#Krige_aqs_US.py
#	Reading aqs data from given year and kriging to full US wrf-chem grid 
#written by: Kate O'Dell 
#version history: 1.0 initial 09/30/16 
#		 2.0 re-written for 2015 data 10/25/16
#		 3.0 re-written for 2013 Oregon (renamed from krige_aqs_2015 to krige_aqs_OR2013) 1/16/17
#		 4.0 re-written for 2015 West for comp. methods class project, remaned krige_aqs_West2015
#		 5.0 re-written for kriging over the entire US for any year, renamed krige_aqs_US 
#		 5.1 slightly edited (sections added, file save names changed) to use different sites for trends, 01/09/18
#		 6.0 reverted to version used for jim crooks data
#		 6.1 added section for LOOCV instead of k-fold for validation
#		 6.2 removed section to create sPM array since this is done in another code in a different way now
########################################################################
#load important modules
import numpy as np
import pylab as pl
from netCDF4 import Dataset
import csv
import datetime as dt
from pykrige.ok import OrdinaryKriging
from scipy.stats import linregress
from random import sample
import time
import sys
from random import shuffle
########################################################################
# user inputs
########################################################################
year = 2014
# to run multiple years at once
'''
if len(sys.argv) > 1:
	print(sys.argv[1],sys.argv[2])
	#sys.exit()
	year = int(sys.argv[1])
	sPM_fn = str(sys.argv[2])
else:
	print('no')
'''
#variogram parameters to use for kriging
vgp = [2.6,8.5,0.1]		#West 2015 params found with kfold method (tested June - September)

# sPM file
sPM_fn = '/home/kaodell//nasa_fires/Lipner_fullUSkrige/datafiles/krigedPM25_06-15_update/sfc_pm' + str(year)[2:] + '_oldsPM_test_github_updateIDs.npz'
# out fp
out_fp = '/home/kaodell//nasa_fires/Lipner_fullUSkrige/datafiles/krigedPM25_06-15_update/'
out_desc = '_code_update_check_oldsPM_oldstats'
########################################################################
#user-defined functions
########################################################################
def haversine(lon0,lon1,lat0,lat1):
    
    r = 6371000.
#m                                                                                                                                                                                                                                                 

    lon0 = lon0*np.pi/180

    lon1 = lon1*np.pi/180

    lat0 = lat0*np.pi/180

    lat1 = lat1*np.pi/180



    return 2*r*np.arcsin(np.sqrt(np.sin((lat1 - lat0)/2.)**2 +\
		 np.cos(lat0)*np.cos(lat1)*np.sin((lon1 - lon0)/2.)**2))
########################################################################
#Time period of interest
########################################################################
t0 = dt.datetime( year = year, month = 1, day = 1 )
tf = dt.datetime( year = year, month = 12, day = 31 )
NT = (tf - t0).days
alldays = np.empty( NT+1, dtype = object )
alldays[0] = t0
for i in range(1, NT+1):
	alldays[i] = alldays[i-1] + dt.timedelta(days = 1)

############################################################################################################
#Load in data and process
############################################################################################################
fn = sPM_fn
f = np.load(fn,'r')
sPM = f['sPM']  
lons = f['lons']
lats = f['lats']
monIDs = f['monIDs']
siteIDs = f['siteIDs']
f.close()

file = '/home/kaodell/nasa_fires/processed_datafiles/US_multiyear_krige/wrf_season_wes_usgs_d01.nc'
nc_fid = Dataset(file, 'r')
glats = nc_fid.variables['XLAT'][:] #grid center lat/lon
glons = nc_fid.variables['XLONG'][:]
nc_fid.close()	

NT, NS = sPM.shape
NX, NY = glons.shape

############################################################################################################
#Part 4: Trim and sort data
############################################################################################################
#sPMin is sPM for sites inside glon and glat grid  
#sPMout sPM for sites close to glon and glat (but not inside)

NT, NS = sPM.shape

dmax_in = 15000. #m - max distance from domain to be considered inside the domain and use for validation stats
dmax_out = 150000. #m - max distance from domain to use in kriging but not in stats 

iIND = []
oIND = []
other_IND = []
for siteIND in range(0,NS):
	dists = haversine(lons[siteIND],glons, lats[siteIND], glats)
	dmin = dists.min()
	if dmin <= dmax_in:
		iIND.append(siteIND)
	elif dmin <= dmax_out:
		oIND.append(siteIND)
		#print dmin

sPMin = sPM[:,iIND]
sPMout = sPM[:,oIND]
ilat = lats[iIND]
ilon = lons[iIND]
olat = lats[oIND]
olon = lons[oIND]

NIS = len(iIND)
NOS = len(oIND)
print('number of in sites: ', NIS)
print('number of out sites: ', NOS)

############################################################################################################
#Part 5: LOOCV Cross Validation
############################################################################################################
#preallocate matricies
dsx = np.array([-999])
dsy = np.array([-999])

ksRSQ = np.empty(NIS)
ksMAE = np.empty(NIS)
ksMB = np.empty(NIS)
knobs = np.empty(NIS)

ksRSQ[:] = np.nan
ksMAE[:] = np.nan
ksMB[:] = np.nan
knobs[:] = np.nan

kPM_site = np.zeros([NT,NIS])
INVALIDS = []				
yout_all = np.zeros([NT,NIS])
for si in range(0,NIS):
#for si in [7]: # for testing code
	#remove monitor being tested
	sPMinTEMP = np.delete(sPMin, si, axis = 1)
	ilatTEMP = np.delete(ilat, si)
	ilonTEMP = np.delete(ilon, si)

	#save data from monitor being tested
	dsx = np.hstack([dsx, sPMin[:, si]])
	clat = ilat[si]
	clon = ilon[si]

	#recombine monitors
	sPMtemp = np.hstack([sPMinTEMP, sPMout])
	latTEMP = np.hstack([ilatTEMP, olat])
	lonTEMP = np.hstack([ilonTEMP, olon])
	
	#remove monitors at same site
	dups = np.where( latTEMP == ilat[si])[0]
	sPMtemp2 = np.delete(sPMtemp, dups, axis = 1)
	latTEMP2 = np.delete(latTEMP, dups)
	lonTEMP2 = np.delete(lonTEMP,dups)

	#loop through time
	yout = np.empty(NT)	
	for ti in range(NT):		
		sPMtemp3 = sPMtemp2[ti,:]
		
		#find and remove all nans				
		trueINDS = np.where(np.isnan(sPMtemp3) == False )[0]
		sPMtemp4 = sPMtemp3[trueINDS]
		latTEMP3 = latTEMP2[trueINDS]
		lonTEMP3 = lonTEMP2[trueINDS]

		#krige to removed monitor location
		#print(lonTEMP3.shape, latTEMP3.shape, sPMtemp4.shape)
		ok = OrdinaryKriging( lonTEMP3, latTEMP3, sPMtemp4, variogram_model='spherical', \
			variogram_parameters = vgp, verbose=False, enable_plotting=False)
		z, ss = ok.execute('points', clon, clat)
		yout[ti] = z
	
	dsy = np.hstack([dsy, yout])

	# calculate kriging stats at the site of the removed monitor (ks)
	xout = sPMin[:,si] 
	kPM_site[:,si] = yout
	nanINDS = np.where( np.isnan(xout) == True )[0]		
	xout = np.delete(xout, nanINDS)
	yout = np.delete(yout, nanINDS)
	ksRSQ[si] = np.corrcoef( xout, yout )[1,0]**2.
	ksMB[si] = np.mean(yout - xout)
	ksMAE[si] = np.mean(np.abs(yout - xout))
	knobs[si] = len(xout) - len(nanINDS)

# calculate kriging stats for full LOOCV (all monitors vs. LOO krige)
#remove invalids (-999s and nans)	
dsx = dsx[1:]
dsy = dsy[1:]

trueINDS = np.where(np.isnan(dsx) == False)[0]
dsx = dsx[trueINDS]
dsy = dsy[trueINDS]

slope, intercept, rvalue, pvalue, stderr = linregress( dsx,dsy )
rsq = rvalue**2.

MB = np.mean(dsy - dsx)
MAE = np.mean(np.abs(dsy - dsx))

############################################################################################################
#Part 6: Krige to grid
############################################################################################################
lonK = np.hstack( [ilon, olon] )
latK = np.hstack( [ilat, olat] )

#preallocate for kriging output
kPM = np.empty( [NT, NX, NY] )
kPM[:,:,:] = np.nan
NED = 0 # number of days with < 30% of sites available
EDind = [] # array for inds of these days, won't do anything with them here but will save
print('total kriging sites: ', NIS + NOS)
for ti in range(NT):
#for ti in range(1): # for testing
	datatemp = np.hstack([sPMin[ti,:], sPMout[ti,:]])
	INDS = np.where( np.isnan( datatemp ) == False )
	if len(INDS[0]) <= 0.3 * (NIS + NOS):
		print('less than 30% of sites available on DOY ',ti)
		#print NED
		NED += 1
		EDind.append(ti)
	datatemp = datatemp[INDS]
	templon = lonK[INDS]
	templat = latK[INDS]
	
	ok = OrdinaryKriging( templon, templat, datatemp, variogram_model='spherical', \
		variogram_parameters = vgp, verbose=False, enable_plotting=False)
	z, ss = ok.execute( 'points', glons.reshape((1,-1)), glats.reshape((1,-1)) )
	kPM[ti, :, :] = z.reshape( (NX, NY) )

print('number of days with <30% of sites reporting: ',NED)

f = open(out_fp + 'kpm' + str(year) + out_desc+'.npz', 'wb')
np.savez(f, alldays=alldays, kPM=kPM,monIDs=monIDs,siteIDs=siteIDs,glats=glats,glons=glons,sPMin=sPMin,\
	sPMout=sPMout,ilon=ilon,ilat=ilat,olat=olat,olon=olon,slope=slope,kPM_site=kPM_site,\
	intercept=intercept,rvalue=rvalue,pvalue=pvalue,stderr=stderr,ksRSQ=ksRSQ,ksMB=ksMB,ksMAE=ksMAE,knobs=knobs,rsq=rsq,MB=MB,MAE=MAE)
f.close()

