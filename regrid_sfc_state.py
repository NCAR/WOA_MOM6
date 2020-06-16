#!/usr/bin/env python

######################################################################

import xarray as xr
import numpy as np
from datetime import datetime
import esmlab
import esmlab_regrid
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
      Regrid surface salinity and potential temperature from land filled WOA to the ocean model grid.
      This scrip relies on the WOA file created using create_filled_sfc.py.
      ''')
  parser.add_argument('-infile', type=str, default='./file_out_filled.nc',
      help='''Full path to the WOA file created using create_filled_sfc.py''')

  parser.add_argument('-path_out', type=str, default='./',
      help='''Path where to output netCDF files''')

  parser.add_argument('-src_grid_name', type=str, default='WOA_04',
      help='''Source grid name. Grids currently supported: WOA_04 and WOA_01.''')

  parser.add_argument('-dst_grid_name', type=str, default='tx0.66v1',
      help='''Destination grid name. Grids currently supported: tx0.66v1''')

  parser.add_argument('-author', type=str, default='Gustavo Marques (gmarques@ucar.edu)',
      help='''Name and email of person creating the dataset''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)
  return

#=======================================================
# This is where all the action happends, i.e., functions
# functions for each diagnostic are called.
#=======================================================

def driver(args):

  gridpath = './script_files/' #args.gridpath #'/glade/work/gmarques/cesm/mom6_input_files/tx0.66v1/salinity_restoring'
  if not os.path.isdir(gridpath):
    print('Creating a directory to place SCRIP files: {} ... \n'.format(gridpath))
    os.system('mkdir '+gridpath)

  esmlab.config.set({'regrid.gridfile-directory': gridpath})
  # src and dst grids
  if args.src_grid_name == 'WOA_01':
    src_grid_name = 'WOA_01_SCRIP'
    os.system('ln -s  /glade/work/gmarques/cesm/datasets/WOA18/WOA_01_SCRIP.nc '+gridpath)
  elif args.src_grid_name == 'WOA_04':
    src_grid_name = 'WOA_04_SCRIP'
    os.system('ln -s  /glade/work/gmarques/cesm/datasets/WOA18/WOA_04_SCRIP.nc '+gridpath)
  else:
    raise ValueError('The source grid name provided, {}, is not supported. Please use WOA_01 or WOA_04.'.format(args.src_grid_name))

  if args.dst_grid_name == 'tx0.66v1':
    dst_grid_name = 'tx0.66v1_SCRIP_190314'
    os.system('ln -s  /glade/work/altuntas/mom.input/tx0.66v1/gen_grid_190314/tx0.66v1_SCRIP_190314.nc '+gridpath)
    # prototype for the restoring file for the tx0.66v1 grid
    ds_out = xr.open_dataset('/glade/p/cesmdata/cseg/inputdata/ocn/mom/tx0.66v1/salt_restore_tx0.66v1_180828.nc',
             decode_times=False)
    ds_out['theta0'] = xr.Variable(dims=('TIME','LAT','LON'), data = np.zeros(ds_out.salt.shape))
  else:
    raise ValueError('The destination grid name provided, {}, is not supported. Only tx0.66v1 is supported at this point. '.format(args.dst_grid_name))

  # generate weights
  R_bilinear = esmlab_regrid.regridder(name_grid_src=src_grid_name, name_grid_dst=dst_grid_name,
                            method='bilinear', overwrite_existing=True)

  ###########################################################################
  # WOA salinity file with land fill, created using create_filled_sfc.py
  woa = xr.open_dataset(args.infile, decode_times=False)

  # average between two-layers (depth = 0 and depth = 10, depth indices 0 and 2)
  woa_s_an_surface_ave = woa.s_an.isel(depth=[0,2]).mean('depth')
  woa_theta0_surface_ave = woa.theta0.isel(depth=[0,2]).mean('depth')

  # regrid and compare against original
  for m in range(len(woa.time)):
    ds_out.salt[m,:] = R_bilinear(woa_s_an_surface_ave[m,:]).rename({'lat':'LAT', 'lon':'LON'})
    ds_out.theta0[m,:] = R_bilinear(woa_theta0_surface_ave[m,:]).rename({'lat':'LAT', 'lon':'LON'})

  ###########################################################################
  # Global attrs
  ds_out.attrs['title'] = 'surface salinity and potential temperature from WOA filled over continents'
  ds_out.attrs['src_file'] = args.infile
  ds_out.attrs['src_grid_name'] = args.src_grid_name
  ds_out.attrs['dst_grid_name'] = args.dst_grid_name
  ds_out.attrs['author'] = args.author
  ds_out.attrs['date'] = datetime.now().isoformat()
  ds_out.attrs['created_using'] = os.path.basename(__file__) + ' which can be found at https://github.com/NCAR/WOA_MOM6'
  ds_out.attrs['git_hash'] = str(subprocess.check_output(["git", "describe","--always"]).strip())
  # save
  fname = 'state_restore_{}_{}{}{}.nc'.format(args.dst_grid_name, datetime.now().isoformat()[0:4],datetime.now().isoformat()[5:7],
           datetime.now().isoformat()[8:10])
  ds_out.to_netcdf(args.path_out+fname)
  print('Done!')
  return

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
