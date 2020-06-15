######################################################################

import xarray as xr
import numpy as np
from datetime import datetime
import gsw
import fill
import os, subprocess

try: import argparse
except: raise Exception('This version of python is not new enough. python 2.7 or newer is required.')

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Creates filled full depth (T,S,PT) for one month (January) using WOA dataset.
      WOA dataset is retrived using  NODC THREDSS server
      ''')
  parser.add_argument('-path_out', type=str, default='./',
      help='''Path where to output netCDF files''')

  parser.add_argument('-file_out_unfilled', type=str, default='test_myscript_unfilled.nc',
      help='''File name for the unfilled output''')

  parser.add_argument('-file_out_filled', type=str, default='test_myscript_filled.nc',
      help='''File name for the filled output''')

  parser.add_argument('-resolution', type=str, default='01',
      help='''Setup which WOA data set to download. Valid entries are: 01 (1 deg) or 04 (0.25 deg)''')

  parser.add_argument('-author', type=str, default='Frank Bryan (bryan@ucar.edu)',
      help='''Name and email of person creating the dataset''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)
  return
#=======================================================
# This is where all the action happends, i.e., functions
# functions for each diagnostic are called.
#=======================================================
def driver(args):

  path_out = args.path_out
  file_out_unfilled = args.file_out_unfilled
  file_out_filled =  args.file_out_filled

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

  # Global attrs
  ds_fill.attrs['title'] = 'T and S from WOA filled over continents'
  ds_fill.attrs['WOA_resolution'] = args.resolution + ', 01 (1 deg), 04 (0.25 deg)'
  ds_fill.attrs['author'] = args.author
  ds_fill.attrs['date'] = datetime.now().isoformat()
  ds_fill.attrs['created_using'] = os.path.basename(__file__) + ' which can be found at https://github.com/NCAR/WOA_MOM6'
  ds_fill.attrs['git_hash'] = str(subprocess.check_output(["git", "describe","--always"]).strip())
  # save
  ds_fill.to_netcdf(path_out+file_out_filled)
  return

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
