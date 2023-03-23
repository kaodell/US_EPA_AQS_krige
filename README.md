# US_EPA_AQS_krige readme
written by Katelyn O'Dell 08.20.21

This repository contains python code used to create observation-based smoke PM2.5 data available at https://doi.org/10.25675/10217/230602. 
Package versions used with this code to create the kriged smoke PM2.5 dataset are availabe in the readme located in the data repository.

Scripts in this folder (run in this order to create the kriged PM2.5 data):

1. process_sfc_obs_update.py
written by Katelyn O'Dell
This script reads in the pre-generated data files from the EPA AQS Data Mart (https://aqs.epa.gov/aqsweb/airdata/download_files.html) for PM2.5 FRM 
(parameter code 88101) and PM2.5 non FRM (parameter code 88502) monitors. It is important to keep track of the download date of these files as they
are updated twice a year. 
The script creates a date x monitor array of PM2.5 concentrations for a single year. Monitors are defined using a monitor ID as suggested by the 
EPA AQS system (for details see: https://aqs.epa.gov/aqsweb/airdata/FileFormats.html#_monitors). Output data includes monitor lat, monitor lon, monitor ID,
site ID (again, see link above for details), and the date x monitor PM2.5 array. The code also creates a figure of monitor PM2.5 concentrations 
in December as a first check for EPA AQS Data Mart file completeness. If an entire state is missing monitors or there are a significanly lower number of
monitors than previous years, this idicates there is likely a significant amount of data not yet incorporated into the pre-generated files. 

2. krige_aqs_US.py
modified from kriging and LOOCV python code written by William Lassman
This script reads in the npz output file from process_sfc_obs_update.py and kriges the surface PM2.5 concentrations to a 15 x 15 km WRF-chem grid 
across the full contiguous US. Grid details are available in the data repository. For details on this kriging process and how the kriging 
parameters were determined see Lassman et al. (2017) and O'Dell et al. (2019) (or if you really want all the gory details, 
refer to my masters thesis - https://hdl.handle.net/10217/193183 ). This code outputs a date x grid_X x grid_Y array of spatially 
interpolated surface PM2.5 concentrations and leave-one-out cross validation (LOOCV) statistics. PM2.5 concentrations on the grid outside 
the contiguous US should be discarded as they are extrapolated. With the LOOCV this code takes several days to run for one year. 
Without the LOOCV it takes less than a day. 

3. calcPMbackground_ns.py
written by Katelyn O'Dell
This script reads in the output from process_sfc_obs_update.py and krige_aqs_US.py scripts and gridded NOAA HMS smoke from Bonne Ford. 
The script outputs an estimated non-smoke background PM2.5 concentration on the 15 x 15 km kriging grid. Background can be set as the 
mean or the median of the non-smoke impacted days. The code also creates several figures comparing the estimated kriged non-smoke background 
concentrations to estimated non-smoke background concentrations at the surface sites. For details about this method of estimating non-smoke PM2.5
please see O'Dell et al. (2019) and O'Dell et al. (2021). 

4. mk_netCDF.py
written by Katelyn O'Dell
This script reads in the output from krige_aqs_US.py and calcPMbackground_ns.py and the gridded NOAA HMS smoke from Bonne Ford. It outputs a netCDF
file containing the krigged PM2.5, gridded HMS, background PM2.5, and LOOCV statistics. These files were uploaded to the Mountain Scholar data repository
( https://doi.org/10.25675/10217/230602 ). This code also outputs figures of the data to check for obvious errors. 

# References
Lassman, W.; Ford, B.; Gan, R. W.; Pfister, G.; Magzamen, S.; Fischer, E. V.; Pierce, J. R. 
Spatial and Temporal Estimates of Population Exposure to Wildfire Smoke during the Washington 
State 2012 Wildfire Season Using Blended Model, Satellite, and In-Situ Data. GeoHealth 2017, 2017GH000049. 
https://doi.org/10.1002/2017GH000049.

O’Dell, K.; Ford, B.; Fischer, E. V.; Pierce, J. R. Contribution of Wildland-Fire Smoke to US PM2.5 
and Its Influence on Recent Trends. Environ. Sci. Technol. 2019, 53 (4), 1797–1804. https://doi.org/10.1021/acs.est.8b05430.

O'Dell, K., K. Bilsback, B. Ford, S. E. Martenies, S. Magzamen, E. V. Fischer and J. R. Pierce (in press), 
Estimated Mortality and Morbidity Attributable to Smoke Plumes in the US: Not Just a Western US Problem, GeoHealth, 
DOI: 10.1029/2021GH000457

Note for future code users

03.32.23 - I have since realized I have a few instances in the process_sfc_obs_update.py code of setting pandas values with chained indexing (which is ill-advised as it may not often work as exptected). I re-looked over this code and the instances in which I had those issues I had confirmed the code was working as expected or the output of the code indicates it worked as expected. In addition, this only impacted years 2005 and 2006. However, if anyone uses this code in the future for those years, those istances should be resolved. See here for details on this warning in pandas: https://pandas.pydata.org/docs/user_guide/indexing.html#indexing-view-versus-copy.

