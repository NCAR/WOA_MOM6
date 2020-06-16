# WOA_MOM6

About
=====
This repository contains Python scripts for generating initial conditions and restoring files for MOM6 using data from the World Ocean Atlas

Requirements
=====
Almost all dependencies can be installed using Conda. Please see the [list of requirements](requirements.txt).
Only esmlab-regrid needs to be installed manually. The following release has been used:

https://github.com/NCAR/esmlab-regrid/releases/tag/v2019.5.1

To work properly, esmlab-regrid requires esmf=8.0.0 and xesmf=0.1.1.

Usage
=====
### create_filled_ic.py

  ```
  create_filled_ic.py [-h] [-path_out PATH_OUT]
                             [-file_out_unfilled FILE_OUT_UNFILLED]
                             [-file_out_filled FILE_OUT_FILLED]
                             [-resolution RESOLUTION] [-author AUTHOR]

  Creates filled full depth (T,S,PT) for one month (January) using WOA dataset.
  WOA dataset is retrived using NODC THREDSS server

  optional arguments:
    -h, --help            show this help message and exit
    -path_out PATH_OUT    Path where to output netCDF files
    -file_out_unfilled FILE_OUT_UNFILLED
                          File name for the unfilled output
    -file_out_filled FILE_OUT_FILLED
                          File name for the filled output
    -resolution RESOLUTION
                          Setup which WOA data set to download. Valid entries
                          are: 01 (1 deg) or 04 (0.25 deg)
    -author AUTHOR        Name and email of person creating the dataset
```

### create_filled_sfc.py

   ```
   create_filled_sfc.py [-h] [-path_out PATH_OUT]
                            [-file_out_unfilled FILE_OUT_UNFILLED]
                            [-file_out_filled FILE_OUT_FILLED]
                            [-resolution RESOLUTION] [-author AUTHOR]

  Creates filled (T, S, PT) for upper 50m for 12 months using WOA dataset. WOA
  dataset is retrived using NODC THREDSS server

  optional arguments:
    -h, --help            show this help message and exit
    -path_out PATH_OUT    Path where to output netCDF files
    -file_out_unfilled FILE_OUT_UNFILLED
                          File name for the unfilled output
    -file_out_filled FILE_OUT_FILLED
                          File name for the filled output
    -resolution RESOLUTION
                          Setup which WOA data set to download. Valid entries
                          are: 01 (1 deg) or 04 (0.25 deg)
    -author AUTHOR        Name and email of person creating the dataset
```

### fill.py

  used in by the scipts above (modified version from pop_tools)

### regrid_sfc_state.py

  ```
  regrid_sfc_state.py [-h] [-infile INFILE] [-path_out PATH_OUT]
                             [-src_grid_name SRC_GRID_NAME]
                             [-dst_grid_name DST_GRID_NAME] [-author AUTHOR]

  Regrid surface salinity and potential temperature from land filled WOA to the
  ocean model grid. This scrip relies on the WOA file created using
  create_filled_sfc.py.

  optional arguments:
    -h, --help            show this help message and exit
    -infile INFILE        Full path to the WOA file created using
                          create_filled_sfc.py
    -path_out PATH_OUT    Path where to output netCDF files
    -src_grid_name SRC_GRID_NAME
                          Source grid name. Grids currently supported: WOA_04
                          and WOA_01.
    -dst_grid_name DST_GRID_NAME
                          Destination grid name. Grids currently supported:
                          tx0.66v1
    -author AUTHOR        Name and email of person creating the dataset
 ```

