Structure for https://tree.nathanfriend.com
```
Med12
  nemo-src
    arch
    ...
    src
    test
    tools
  static-data
    basin_mask.nc
    ...
    domain_cfg.nc
    ...
    weights_ERA5-MED7km_bicub.nc
    weights_ERA5-MED7km_bilin.nc
  dynamic-data
    era5/
    cmecs/
  era5-preprocessing
    grib/
      *.grib (raw data)
    run/
      *.run (download scripts)
    tools/
      *.F90 (process raw data into NetCDF)
    fetch-and-process.py
  cmecs-preprocessing
```
