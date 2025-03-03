# Era5 overview
ERA5 is the fifth generation ECMWF reanalysis for the global climate and weather for the past 8 decades. Data is available from 1940 onwards. ERA5 replaces the ERA-Interim reanalysis.

Reanalysis combines model data with observations from across the world into a globally complete and consistent dataset using the laws of physics. This principle, called data assimilation, is based on the method used by numerical weather prediction centres, where every so many hours (12 hours at ECMWF) a previous forecast is combined with newly available observations in an optimal way to produce a new best estimate of the state of the atmosphere, called analysis, from which an updated, improved forecast is issued. Reanalysis works in the same way, but at reduced resolution to allow for the provision of a dataset spanning back several decades. Reanalysis does not have the constraint of issuing timely forecasts, so there is more time to collect observations, and when going further back in time, to allow for the ingestion of improved versions of the original observations, which all benefit the quality of the reanalysis product.

ERA5 provides hourly estimates for a large number of atmospheric, ocean-wave and land-surface quantities. An uncertainty estimate is sampled by an underlying 10-member ensemble at three-hourly intervals. Ensemble mean and spread have been pre-computed for convenience. Such uncertainty estimates are closely related to the information content of the available observing system which has evolved considerably over time. They also indicate flow-dependent sensitive areas. To facilitate many climate applications, monthly-mean averages have been pre-calculated too, though monthly means are not available for the ensemble mean and spread.

ERA5 is updated daily with a latency of about 5 days. In case that serious flaws are detected in this early release (called ERA5T), this data could be different from the final release 2 to 3 months later. In case that this occurs users are notified.

The data set presented here is a regridded subset of the full ERA5 data set on native resolution. It is online on spinning disk, which should ensure fast and easy access. It should satisfy the requirements for most common applications.

An overview of all ERA5 datasets can be found in [this article](https://confluence.ecmwf.int/display/CKB/The+family+of+ERA5+datasets). Information on access to ERA5 data on native resolution is provided in [these guidelines](https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5).

Data has been regridded to a regular lat-lon grid of 0.25 degrees for the reanalysis and 0.5 degrees for the uncertainty estimate (0.5 and 1 degree respectively for ocean waves). There are four main sub sets: hourly and monthly products, both on pressure levels (upper air fields) and single levels (atmospheric, ocean-wave and land surface quantities).

The present entry is "ERA5 hourly data on single levels from 1940 to present".

[ERA5 data documentation](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation)
Detailed information relating to the ERA5 data archive can be found in the web link above.


[The ERA5 global reanalysis](https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1002/qj.3803)
Journal article describing ERA5.


## Era5 preprocessing in python

We have three main directories (along with a `home` directory): 
- `work`: here we store temporary and intermediate files;
- `grib`: here we store the raw `.grib` files from the download;
- `arch`: here we store the processed files, those that will be fed to NEMO.
 
These directories are set in the `.json` file `directories.json` with the following syntax:
```json
{
    "home": "/home/ftucciar",
    "work": "/home/ftucciar/Med12/preprocessing-era5/work",
    "grib": "/home/ftucciar/Med12/preprocessing-era5/grib",
    "arch": "/home/ftucciar/Med12/preprocessing-era5/archive"
}
```
and are accessed in each processing file as
```python
import json
# Directories
dirs = json.load( open('directories.json') )
```


### Parameters dictionary for the processing
Refer to [ERA5 data documentation](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation) page to understand where to find the parameters for this dictionary
```python
param = {
         "long_name" : str,           # Variable name in CDS, e.g. '10m_u_component_of_wind'
         "var_name" : str,            # ShortName, e.g. '10u'
         "out_name" : str,            # Output name of the variable
         "chr_id" : str,              # paramID, e.g. '165'
         "nts" : int,                 # Number of time steps
         "nx" : int,                  # Dimensions in x
         "ny" : int                   # Dimensions in y
         "daymean": bool,             # Flag to average over one day
         "units_change": {
                     "change": bool,  # Flag to convert into another set of units
                     "ucf": float     # Units conversion factor
                     },
         "maskland": {
                     "mask": bool,    # Flag to mask the land
                     "exec": str,     # Path to executable (as in this case this is in fortran)
                     "maskfile": str  # Path to NetCDF file witht he mask 
                     }
         }
```
### 10m wind $u,v-$ components (`u10m` and `v10m`)
These parameter are the [eastward ( $u$ )](https://codes.ecmwf.int/grib/param-db/165) and [northward ($v$)](https://codes.ecmwf.int/grib/param-db/166) components of the 10m wind. They are the horizontal speed of air moving towards the east and the north, respectively, at a height of ten metres above the surface of the Earth, in metres per second.

```python
# -*- coding: utf-8 -*-
import json
from era5_process import era5_process

# Directories
dirs = json.load( open('directories.json') )

# Parameters of the input field
param = {
         # Parameters depending on the field processed
         # Refer to: https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation 
         #
         "long_name": "Variable_name_in_CDS", # Variable name in CDS, e.g. '10m_u_component_of_wind'
         "var_name": "ShortName",             # ShortName, e.g. '10u'
         "out_name": "output_name",           # Output name of the variable
         "chr_id": "paramID",                 # paramID, e.g. '165'
         # Common parameters
         "nts": 24,
         "nx": 1440,
         "ny": 721,
         "daymean": False,
         "units_change": {
                          "change": False,
                          "ucf": 1
                         },
         "maskland": {
                      "mask": True,
                      "exec": "/home/ftucciar/Med12/preprocessing-era5/tools/scr2/flandR.x",
                      "maskfile": "/home/ftucciar/Med12/preprocessing-era5/tools/lsm_ERA5_0.25.nc"
                     }
        }
# %%
# Parameters of the input field
param["long_name"] = "10m_u_component_of_wind"
param["var_name"] = "10u"
param["out_name"] = "u10m"
param["chr_id"] = "165"

era5_process(param, dirs, 2020, 1, 2020, 1, cleanup=False)

# Parameters of the input field
param["long_name"] = "10m_v_component_of_wind"
param["var_name"] = "10v"
param["out_name"] = "v10m"
param["chr_id"] = "166"
        
era5_process(param, dirs, 2020, 1, 2020, 1, cleanup=False)
```
### 2m Temperature (`t2m`)
### Total Precipitation (`precip`)
### Snowfall (`snow`)
### Surface Solar Radiation Downwards (`swrd`)
### Surface Thermal Radiation Downwards (`lwrd`)
