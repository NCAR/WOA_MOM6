######################################################################

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

import gsw  # TEOS-10 equation of state

import fill  # land fill

###########################################################################
# specs for output files

path_out = '/glade/work/bryan/Observations/WOA18/'
file_out_unfilled = 'test_sfc_script_unfilled.nc'
file_out_filled = 'test_sfc_script_filled.nc'

dep_max = 50.0  # maximum depth to include in output

###########################################################################
# Specs for input data from NODC

# URL for NODC THREDSS server
ncei_thredds_url = 'https://data.nodc.noaa.gov/thredds/dodsC/ncei/woa/'

# Setup which WOA data set to download
# 1 deg. WOA
fext = '01'
res = '1.00'

# 0.25 deg WOA
#fext = '04'
#res = '0.25'

decade = 'decav'
date_mon = range(1,13)

###########################################################################
# Get the data from NODC

# Get monthly temperature & salinity
base_url_t=ncei_thredds_url + '/temperature/' + decade + '/' + res + '/'
file_root_t='woa18_' + decade + '_t'
files_t=['{0:s}{1:s}{2:02d}_{3:s}.nc'.format(base_url_t,file_root_t,date_str,fext) for date_str in date_mon]
print('Monthly temeperature files:',files_t)
dst_mon = xr.open_mfdataset(files_t,decode_times=False,data_vars='minimal',combine='by_coords')
dst_mon = dst_mon.sel(depth=slice(0,dep_max)).compute()

base_url_s=ncei_thredds_url + '/salinity/' + decade + '/' + res + '/'
file_root_s='woa18_' + decade + '_s'
files_s=['{0:s}{1:s}{2:02d}_{3:s}.nc'.format(base_url_s,file_root_s,date_str,fext) for date_str in date_mon]
print('Monthly salinity files:',files_s)
dss_mon = xr.open_mfdataset(files_s,decode_times=False,data_vars='minimal',combine='by_coords')
dss_mon = dss_mon.sel(depth=slice(0,dep_max)).compute()


# Merge seasonal data sets
# force download to overcome limitation of interpolation on chunked dims
ds_upper = dst_mon.merge(dss_mon)

###########################################################################
# Compute potential temperature
z1d = -ds_upper['depth']
lat1d = ds_upper['lat']
lon1d = ds_upper['lon']

z3d,lat3d,lon3d = xr.broadcast(z1d,lat1d,lon1d)

p3d = gsw.p_from_z(z3d,lat3d)
SR = gsw.SR_from_SP(ds_upper['s_an'].data)
PT0 = gsw.conversions.pt0_from_t(SR,ds_upper['t_an'],p3d)


ds_upper['theta0'] = xr.DataArray(PT0,dims=('time','depth','lat','lon'))                         
for var in ('time', 'depth','lat','lon') :
    ds_upper['theta0'].assign_coords({var:ds_upper[var]})
    
ds_upper['theta0'].attrs = ds_upper['t_an'].attrs
ds_upper['theta0'].attrs['long_name']='Potential temperature from objectively analyzed in situ temperature and salinty'
ds_upper['theta0'].attrs['standard_name']='sea_water_potential_temperature'
ds_upper['theta0'].encoding = ds_upper['t_an'].encoding

ds_upper['time'].attrs = dst_mon['time'].attrs
ds_upper['depth'].attrs = dst_mon['depth'].attrs
ds_upper['lat'].attrs = dst_mon['lat'].attrs
ds_upper['lon'].attrs = dst_mon['lon'].attrs

###########################################################################
# Output unfilled arrays

ds_upper.to_netcdf(path_out+file_out_unfilled)

###########################################################################
# Do land fill and output filled arrays

mask_all_points = xr.ones_like(ds_upper['t_an']).astype('bool')
ds_fill = xr.Dataset()

for var in ('t_an', 's_an', 'theta0') :
    print('Starting ',var,' ...')

    ds_fill[var] = fill.lateral_fill(ds_upper[var],mask_all_points,ltripole=False,
                                     tol=5.0e-3,use_sor=True,rc=1.88,max_iter=1000)

ds_fill['time'].attrs = dst_mon['time'].attrs

ds_fill.to_netcdf(path_out+file_out_filled)
