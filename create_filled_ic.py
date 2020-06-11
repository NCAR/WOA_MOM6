######################################################################

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

import gsw

import fill

path_out = '/glade/work/bryan/Observations/WOA18/'
file_out_unfilled = 'test_myscript_unfilled.nc'
file_out_filled = 'test_myscript_filled.nc'

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
date_mon = 1      # This is the month of the output IC file
date_seas = range(13,17)   # all four seasons


###########################################################################
# Get the data from NODC

# Get monthly temperature & salinity
base_url_t=ncei_thredds_url + '/temperature/' + decade + '/' + res + '/'
file_root_t='woa18_' + decade + '_t'
files_t='{0:s}{1:s}{2:02d}_{3:s}.nc'.format(base_url_t,file_root_t,date_mon,fext) 
print('Monthly temperature files:',files_t)
dst_mon = xr.open_dataset(files_t,decode_times=False) #,data_vars='minimal',combine='by_coords')

base_url_s=ncei_thredds_url + '/salinity/' + decade + '/' + res + '/'
file_root_s='woa18_' + decade + '_s'
files_s='{0:s}{1:s}{2:02d}_{3:s}.nc'.format(base_url_s,file_root_s,date_mon,fext)
print('Monthly salinity files:',files_s)
dss_mon = xr.open_dataset(files_s,decode_times=False) #,data_vars='minimal',combine='by_coords')

# Merge monthly input data sets
ds_mon = dst_mon.merge(dss_mon).squeeze()

# Get seasonal temperature & salinity
files_t=['{0:s}{1:s}{2:02d}_{3:s}.nc'.format(base_url_t,file_root_t,date_str,fext) for date_str in date_seas]
print('Seasonal temeperature files:',files_t)
dst_seas = xr.open_mfdataset(files_t,decode_times=False,data_vars='minimal',combine='by_coords')

files_s=['{0:s}{1:s}{2:02d}_{3:s}.nc'.format(base_url_s,file_root_s,date_str,fext) for date_str in date_seas]
print('Seasonal salinity files:',files_s)
dss_seas = xr.open_mfdataset(files_s,decode_times=False,data_vars='minimal',combine='by_coords')

# Merge seasonal data sets
# force download to overcome limitation of interpolation on chunked dims
ds_seas = dst_seas.merge(dss_seas).compute() 

###########################################################################
# Extend Monthly data to deepest levels by interpolating from seasonal

nseas = np.shape(ds_seas['time'])[0]
nmon = np.shape(ds_mon['time'])
print('length seasonal file=',nseas)
print('length monthly file=',nmon)

time_target = ds_mon['time'].values
print('target time of ic file=',time_target)

nlon = np.shape(ds_mon['lon'])[0]
nlat = np.shape(ds_mon['lat'])[0]
print(' nlon=',nlon,' nlat=',nlat)
ndep_upper = np.shape(ds_mon['depth'])[0]
ndep_full = np.shape(ds_seas['depth'])[0]
print('# levels seasonal file=',ndep_full)
print('# levels monthly file=',ndep_upper)

# Check that the right depths are being extracted and concatinated
#tmp_dep = xr.concat([ds_mon['depth'],ds_seas['depth'][ndep_upper:]],dim='depth')
#print('# of levels in extended file =',np.shape(tmp_dep)[0])
#print('concatinated full depth')
#for k in range(0,ndep_full) :
#    print('  k=',k,' z=',tmp_dep.values[k])

# Pull out seasons -1 and N for periodic b.c on time interpolation.
ds_seas_bef = ds_seas.isel(time=nseas-1)
ds_seas_aft = ds_seas.isel(time=0)

#  Adjust times to be monotonic
ds_seas_bef['time'] = ds_seas_bef['time'] - 12.
ds_seas_aft['time'] = ds_seas_aft['time'] + 12.

print('time before=',ds_seas_bef['time'].values, ' time after =',ds_seas_aft['time'].values)

###########################################################################
# Create a data set with monthly data in upper ocean, 
#  interpolated seasonal data in deeper ocean

ds_deep = xr.Dataset()
for var in ('t_an', 's_an') :
    print('Starting ',var,' ...')
    tmp_extend = xr.concat([ds_seas_bef[var][ndep_upper:,:,:],
                            ds_seas[var][:,ndep_upper:,:,:],
                            ds_seas_aft[var][ndep_upper:,:,:]],
                           dim='time').compute()
    print(' ... extended seasonal data in time. shape=',np.shape(tmp_extend))
    print(' ... extended seasonal data times=',tmp_extend['time'].values)

    print(' ... interpolating to time ',time_target)
    tmp_interp = tmp_extend.interp(time=time_target)
    print('... finished inteprolation to months. shape=',np.shape(tmp_interp))
    print( ' ... shape of monthly upper ocean array=',np.shape(ds_mon[var]))

    ds_deep[var] = xr.concat([ds_mon[var],tmp_interp],dim='depth')
    print('... finished inserting deep values into monthly arrays. shape = ',np.shape(ds_deep[var]))


###########################################################################
# Compute pressure, reference salinity and potential temperature using TEOS-10

z1d = -ds_deep['depth']
lat1d = ds_deep['lat']
lon1d = ds_deep['lon']

z3d,lat3d,lon3d = xr.broadcast(z1d,lat1d,lon1d)

p3d = gsw.p_from_z(z3d,lat3d)
SR = gsw.SR_from_SP(ds_deep['s_an'].data)
PT0 = gsw.conversions.pt0_from_t(SR,ds_deep['t_an'],p3d)


ds_deep['theta0'] = xr.DataArray(PT0,dims=('depth','lat','lon'))                         
for var in ('depth','lat','lon') :
    ds_deep['theta0'].assign_coords({var:ds_seas[var]})
    
ds_deep['theta0'].attrs = ds_seas['t_an'].attrs
ds_deep['theta0'].attrs['long_name']='Potential temperature from objectively analyzed in situ temperature and salinty'
ds_deep['theta0'].attrs['standard_name']='sea_water_potential_temperature'
ds_deep['theta0'].encoding = ds_seas['t_an'].encoding


ds_deep.to_netcdf(path_out+file_out_unfilled)


###########################################################################
# Fill all land with values interpolated from nearest ocean

mask_all_points=xr.DataArray(np.ones((ndep_full,nlat,nlon)).astype('bool'),
                            dims=('depth','lat','lon'),
                            coords=(ds_seas['depth'],ds_seas['lat'],ds_seas['lon']))

ds_fill = xr.Dataset()

for var in ('t_an', 's_an', 'theta0') :
    print('Starting ',var,' ...')

    ds_fill[var] = fill.lateral_fill(ds_deep[var],mask_all_points,ltripole=False,
                                     tol=5.0e-3,use_sor=True,rc=1.88,max_iter=1000)

ds_fill['depth'].attrs = ds_seas['depth'].attrs
ds_fill['lat'].attrs = ds_seas['lat'].attrs
ds_fill['lon'].attrs = ds_seas['lon'].attrs


ds_fill.to_netcdf(path_out+file_out_filled)
