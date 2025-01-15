# Folder setup
The code is structured with base folder `Med12` and inside it it is found
- nemo-src, containing the source code of nemo v4.0.7;
- data-static, containing all those input fields that are not required to change during execution;
- data-dynamic, containing conversely all those data that change from run to run (forcings, lateral boundary conditions and so on)
- preprocessing-[era5, cmecs], all the scripts necessary to retrieve data and process them to finally put them into data-dynamics.
This structure is represented in the following tree.
```
└── Med12/
    ├── data-dynamic/
    │   ├── era5/
    │   └── cmecs/
    ├── data-static/
    │   ├── basin_mask.nc
    │   ├── ...
    │   ├── domain_cfg.nc
    │   ├── ...
    │   ├── weights_ERA5-MED7km_bicub.nc
    │   └── weights_ERA5-MED7km_bilin.nc
    ├── nemo-src/
    │   ├── arch
    │   ├── ...
    │   ├── src
    │   ├── test
    │   └── tools
    ├── preprocessing-era5/
    │   ├── grib/
    │   │   └── *.grib (raw data)
    │   ├── run/
    │   │   └── *.run (download scripts)
    │   ├── tools/
    │   │   └── *.F90 (process raw data into NetCDF)
    │   └── fetch-and-process.py
    └── preprocessing-cmecs
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

