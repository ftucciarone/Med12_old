# Era5 overview
<img align="right" src="https://cds.climate.copernicus.eu/thumbnails/JaABbrZjt5iZG1MVCZFaW4gsA4E=/360x0/filters:format(webp)/object-store.os-api.cci2.ecmwf.int/cci2-prod-catalogue/resources/reanalysis-era5-single-levels/overview-detail_d2d128d66670fd68c467a008450654b2030747dba36dfa116bf69461747cc14a.png" width="400" height="400" />


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

## Downloading the forcing raw data
First, we need to download the necessary raw data. This can be done with a python API for the following fields

- [ ] 10 m $u-$ component of wind, [10m_u_component_of_wind](https://codes.ecmwf.int/grib/param-db/165)
- [ ] 10 m $v-$ component of wind, [10m_v_component_of_wind](https://codes.ecmwf.int/grib/param-db/166)
- [ ] 2 m temperature, [2m_temperature](https://codes.ecmwf.int/grib/param-db/167)
- [ ] 2 m dewpoint temperature, [2m_dewpoint_temperature](https://codes.ecmwf.int/grib/param-db/168)
- [ ] Mean sea level pressure, [mean_sea_level_pressure](https://codes.ecmwf.int/grib/param-db/151)
- [ ] Surface short-wave (solar) radiation downwards, [surface_solar_radiation_downwards](https://codes.ecmwf.int/grib/param-db/169)
- [ ] Surface long-wave (thermal) radiation downwards, [surface_thermal_radiation_downwards](https://codes.ecmwf.int/grib/param-db/175)
- [ ] Total precipitation, [total_precipitation](https://codes.ecmwf.int/grib/param-db/228)
- [ ] Snowfall, [snowfall](https://codes.ecmwf.int/grib/param-db/144)





## Process the raw data
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
These parameter are the [eastward](https://codes.ecmwf.int/grib/param-db/165) ($u$) and [northward](https://codes.ecmwf.int/grib/param-db/166) ($v$) components of the 10m wind. They are the horizontal speed of air moving towards the east and the north, respectively, at a height of ten metres above the surface of the Earth, in metres per second.

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
# Refer to: https://apps.ecmwf.int/codes/grib/param-db/165
param["long_name"] = "10m_u_component_of_wind"
param["var_name"] = "10u"
param["out_name"] = "u10m"
param["chr_id"] = "165"
era5_process(param, dirs, 2020, 1, 2020, 1, cleanup=False)

# Parameters of the input field
# Refer to: https://apps.ecmwf.int/codes/grib/param-db/166
param["long_name"] = "10m_v_component_of_wind"
param["var_name"] = "10v"
param["out_name"] = "v10m"
param["chr_id"] = "166"
era5_process(param, dirs, 2020, 1, 2020, 1, cleanup=False)
```
![wind](https://github.com/ftucciarone/Med12/blob/main/preprocessing-era5/python/figures/u10m_wind.gif?raw=true)

### 2m Temperature (`t2m`)
This parameter is the temperature of air at 2m above the surface of land, sea or in-land waters. This parameter has units of kelvin (K). Temperature measured in kelvin can be converted to degrees Celsius (°C) by subtracting 273.15.
```fortran
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
         "long_name": "Variable_name_in_CDS", # Variable name in CDS, e.g. '2m_temperature'
         "var_name": "ShortName",             # ShortName, e.g. r"\2t"
         "out_name": "output_name",           # Output name of the variable, e.g. 't2m'
         "chr_id": "paramID",                 # paramID, e.g. '167'
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
# Refer to: https://apps.ecmwf.int/codes/grib/param-db/167
param["long_name"] = "2m_temperature"
param["var_name"] = r"\2t"
param["out_name"] = "t2m"
param["chr_id"] = "167"
era5_process(param, dirs, 2020, 1, 2020, 1, cleanup=False)
```
### Total Precipitation (`precip`)
This parameter is the accumulated liquid and frozen water, comprising rain and snow, that falls to the Earth's surface. It is the sum of large-scale precipitation and convective precipitation. Large-scale precipitation is generated by the cloud scheme in the ECMWF Integrated Forecasting System (IFS). The cloud scheme represents the formation and dissipation of clouds and large-scale precipitation due to changes in atmospheric quantities (such as pressure, temperature and moisture) predicted directly by the IFS at spatial scales of the [grid box](https://confluence.ecmwf.int/display/CKB/Model+grid+box+and+time+step) or larger. Convective precipitation is generated by the convection scheme in the IFS, which represents convection at spatial scales smaller than the grid box. [See further information](https://confluence.ecmwf.int/display/CKB/Convective+and+large-scale+precipitation). This parameter does not include fog, dew or the precipitation that evaporates in the atmosphere before it lands at the surface of the Earth.

This parameter is the total amount of [water accumulated over a particular time period which depends on the data extracted](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Meanrates/fluxesandaccumulations). The units of this parameter are depth in metres of water equivalent. It is the depth the water would have if it were spread evenly over the grid box.
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
         "nts": 1,
         "nx": 1440,
         "ny": 721,
         "daymean": False,
         "units_change": {
                          "change": True,
                          "ucf": 3.6          
                         },
         "maskland": {
                      "mask": True,
                      "exec": "/home/ftucciar/Med12/preprocessing-era5/tools/scr2/flandR.x",
                      "maskfile": "/home/ftucciar/Med12/preprocessing-era5/tools/lsm_ERA5_0.25.nc"
                     }
        }
# %%
# Parameters of the input field
# Refer to: https://apps.ecmwf.int/codes/grib/param-db/228
param["long_name"] = "total_precipitation"
param["var_name"] = "tp"
param["out_name"] = "precip"
param["chr_id"] = "228"
era5_process(param, dirs, 2020, 1, 2020, 1, cleanup=False)
```
### Snowfall (`snow`)  
This parameter is the accumulated snow that falls to the Earth's surface. It is the sum of large-scale snowfall and convective snowfall. Large-scale snowfall is generated by the cloud scheme in the ECMWF Integrated Forecasting System (IFS). The cloud scheme represents the formation and dissipation of clouds and large-scale precipitation due to changes in atmospheric quantities (such as pressure, temperature and moisture) predicted directly by the IFS at spatial scales of the [grid box](https://confluence.ecmwf.int/display/CKB/Model+grid+box+and+time+step) or larger. Convective snowfall is generated by the convection scheme in the IFS, which represents convection at spatial scales smaller than the grid box. See further information.

This parameter is the total amount of [water accumulated over a particular time period which depends on the data extracted](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Meanrates/fluxesandaccumulations) The units of this parameter are depth in metres of water equivalent. It is the depth the water would have if it were spread evenly over the grid box.
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
         "nts": 1,
         "nx": 1440,
         "ny": 721,
         "daymean": False,
         "units_change": {
                          "change": True,
                          "ucf": 3.6
                         },
         "maskland": {
                      "mask": True,
                      "exec": "/home/ftucciar/Med12/preprocessing-era5/tools/scr2/flandR.x",
                      "maskfile": "/home/ftucciar/Med12/preprocessing-era5/tools/lsm_ERA5_0.25.nc"
                     }
        }
# %%
# Parameters of the input field
# Refer to: https://apps.ecmwf.int/codes/grib/param-db/144
param["long_name"] = "snowfall"
param["var_name"] = "sf"
param["out_name"] = "snow"
param["chr_id"] = "144"
era5_process(param, dirs, 2020, 1, 2020, 1, cleanup=False)
```
### Surface Solar Radiation Downwards (`swrd`)
This parameter is the amount of solar radiation (also known as shortwave radiation) that reaches a horizontal plane at the surface of the Earth. This parameter comprises both direct and diffuse solar radiation.

Radiation from the Sun (solar, or shortwave, radiation) is partly reflected back to space by clouds and particles in the atmosphere (aerosols) and some of it is absorbed. The rest is incident on the Earth's surface (represented by this parameter). See [further documentation](https://www.ecmwf.int/sites/default/files/elibrary/2015/18490-radiation-quantities-ecmwf-model-and-mars.pdf).

To a reasonably good approximation, this parameter is the model equivalent of what would be measured by a pyranometer (an instrument used for measuring solar radiation) at the surface. However, care should be taken when comparing model parameters with observations, because observations are often local to a particular point in space and time, rather than representing averages over a [model grid box](https://confluence.ecmwf.int/display/CKB/Model+grid+box+and+time+step).

This parameter [is accumulated over a particular time period which depends on the data extracted](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Meanrates/fluxesandaccumulations). The units are joules per square metre (J m-2). To convert to watts per square metre (W m-2), the accumulated values should be divided by the accumulation period expressed in seconds.
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
         "nts": 1,
         "nx": 1440,
         "ny": 721,
         "daymean": False,
         "units_change": {
                          "change": True,
                          "ucf": 3600
                         },
         "maskland": {
                      "mask": True,
                      "exec": "/home/ftucciar/Med12/preprocessing-era5/tools/scr2/flandR.x",
                      "maskfile": "/home/ftucciar/Med12/preprocessing-era5/tools/lsm_ERA5_0.25.nc"
                     }
        }
# %%
# Parameters of the input field
# Refer to: https://apps.ecmwf.int/codes/grib/param-db/169
param["long_name"] = "surface_solar_radiation_downwards"
param["var_name"] = "ssrd"
param["out_name"] = "swrd"
param["chr_id"] = "169"
era5_process(param, dirs, 2020, 1, 2020, 1, cleanup=False)
```
### Surface Thermal Radiation Downwards (`lwrd`)
This parameter is the amount of thermal (also known as longwave or terrestrial) radiation emitted by the atmosphere and clouds that reaches a horizontal plane at the surface of the Earth.

The surface of the Earth emits thermal radiation, some of which is absorbed by the atmosphere and clouds. The atmosphere and clouds likewise emit thermal radiation in all directions, some of which reaches the surface (represented by this parameter). See [further documentation](https://www.ecmwf.int/sites/default/files/elibrary/2015/18490-radiation-quantities-ecmwf-model-and-mars.pdf).

This parameter [is accumulated over a particular time period which depends on the data extracted](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Meanrates/fluxesandaccumulations). The units are joules per square metre (J m-2). To convert to watts per square metre (W m-2), the accumulated values should be divided by the accumulation period expressed in seconds.
```python
[comment]: <> ![Alt Text](https://media.giphy.com/media/vFKqnCdLPNOKc/giphy.gif)
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
         "nts": 1,
         "nx": 1440,
         "ny": 721,
         "daymean": False,
         "units_change": {
                          "change": True,
                          "ucf": 3600
                         },
         "maskland": {
                      "mask": True,
                      "exec": "/home/ftucciar/Med12/preprocessing-era5/tools/scr2/flandR.x",
                      "maskfile": "/home/ftucciar/Med12/preprocessing-era5/tools/lsm_ERA5_0.25.nc"
                     }
        }
# %%
# Parameters of the input field
# Refer to: https://apps.ecmwf.int/codes/grib/param-db/175
param["long_name"] = "surface_thermal_radiation_downwards"
param["var_name"] = "strd"
param["out_name"] = "lwrd"
param["chr_id"] = "175"
era5_process(param, dirs, 2020, 1, 2020, 1, cleanup=False)
```
