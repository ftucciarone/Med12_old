Structure for https://tree.nathanfriend.com
```
Med12
  data-dynamic
    era5/
    cmecs/
  data-static
    basin_mask.nc
    ...
    domain_cfg.nc
    ...
    weights_ERA5-MED7km_bicub.nc
    weights_ERA5-MED7km_bilin.nc
  nemo-src
    arch
    ...
    src
    test
    tools
  preprocessing-era5
    grib/
      *.grib (raw data)
    run/
      *.run (download scripts)
    tools/
      *.F90 (process raw data into NetCDF)
    fetch-and-process.py
  preprocessing-cmecs
```
