# Folder setup
```
.
└── Med12/
    ├── nemo-src/
    │   ├── arch
    │   ├── ...
    │   ├── src
    │   ├── test
    │   └── tools
    ├── static-data/
    │   ├── basin_mask.nc
    │   ├── ...
    │   ├── domain_cfg.nc
    │   ├── ...
    │   ├── weights_ERA5-MED7km_bicub.nc
    │   └── weights_ERA5-MED7km_bilin.nc
    ├── dynamic-data/
    │   ├── era5/
    │   └── cmecs/
    ├── era5-preprocessing/
    │   ├── grib/
    │   │   └── *.grib (raw data)
    │   ├── run/
    │   │   └── *.run (download scripts)
    │   ├── tools/
    │   │   └── *.F90 (process raw data into NetCDF)
    │   └── fetch-and-process.py
    └── cmecs-preprocessing
```



# Med12 configuration setup

The configuration is described in Storto et. al. https://gmd.copernicus.org/articles/16/4811/2023/gmd-16-4811-2023.html

installation can be done with the following steps

```
mkdir NemoMed12
cd NemoMed12
svn co https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r4.0.7/
```
It is important to stick with version 4.0.7, version 4.0 does not work

## To create the configuration MED12
```
cd cfgs
cp -r AMM12 MED12
cd MED12
echo "bld::tool::fppkeys key_mpp_mpi key_asminc" > cpp_MED12.fcm
```
then remove `top` from `cpp_MED12.fcm`
copy `EXPREF` into 
compile


# Retrieving and processing forcing fields

### [Atmospheric forcing: ERA5 from Copernicus Climate Data Store](https://cds.climate.copernicus.eu/), [readme (instructions)](forcings/era5-atmos/readme-era5.md)
### [Lateral boundary conditions: CMEMS from Copernicus Marine Service](https://marine.copernicus.eu/), [readme (instructions)](forcings/cmems-latbnd/readme-cmems.md)

