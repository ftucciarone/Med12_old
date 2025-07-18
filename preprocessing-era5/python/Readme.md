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
- `arch`: here we store the processed files, those that will be fed to NEMO;
- `figs`: here we store postprocessed images made with `makegif.py`.
 
These directories are set in the `.json` file `directories.json` with the following syntax:

https://github.com/ftucciarone/Med12_old/blob/0c60baaae4cfdeec9b4ef27d86b35d11ade75f03/preprocessing-era5/python/directories.json#L1-L7

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

https://github.com/ftucciarone/Med12_old/blob/6733a2eacd2244d05870f1bbe636c14d8648cd32/preprocessing-era5/python/runproc_10m_components_of_wind.py#L1-L50

<p align="center">
  <img src="https://github.com/ftucciarone/Med12_old/blob/main/preprocessing-era5/python/figures/u10m_wind.gif?raw=true" alt="animated" />
  <img src="https://github.com/ftucciarone/Med12_old/blob/main/preprocessing-era5/python/figures/v10m_wind.gif?raw=true" alt="animated" />
</p>


### 2m Temperature (`t2m`)
This parameter is the temperature of air at 2m above the surface of land, sea or in-land waters. This parameter has units of kelvin (K). Temperature measured in kelvin can be converted to degrees Celsius (°C) by subtracting 273.15.

https://github.com/ftucciarone/Med12_old/blob/6733a2eacd2244d05870f1bbe636c14d8648cd32/preprocessing-era5/python/runproc_2m_temperature.py#L1-L41

<p align="center">
  <img src="https://github.com/ftucciarone/Med12_old/blob/main/preprocessing-era5/python/figures/t2m.gif?raw=true" alt="animated" />
</p>

### Total Precipitation (`precip`)
This parameter is the accumulated liquid and frozen water, comprising rain and snow, that falls to the Earth's surface. It is the sum of large-scale precipitation and convective precipitation. Large-scale precipitation is generated by the cloud scheme in the ECMWF Integrated Forecasting System (IFS). The cloud scheme represents the formation and dissipation of clouds and large-scale precipitation due to changes in atmospheric quantities (such as pressure, temperature and moisture) predicted directly by the IFS at spatial scales of the [grid box](https://confluence.ecmwf.int/display/CKB/Model+grid+box+and+time+step) or larger. Convective precipitation is generated by the convection scheme in the IFS, which represents convection at spatial scales smaller than the grid box. [See further information](https://confluence.ecmwf.int/display/CKB/Convective+and+large-scale+precipitation). This parameter does not include fog, dew or the precipitation that evaporates in the atmosphere before it lands at the surface of the Earth.

This parameter is the total amount of [water accumulated over a particular time period which depends on the data extracted](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Meanrates/fluxesandaccumulations). The units of this parameter are depth in metres of water equivalent. It is the depth the water would have if it were spread evenly over the grid box.


https://github.com/ftucciarone/Med12_old/blob/6733a2eacd2244d05870f1bbe636c14d8648cd32/preprocessing-era5/python/runproc_total_precipitation.py#L1-L41

<p align="center">
  <img src="https://github.com/ftucciarone/Med12_old/blob/main/preprocessing-era5/python/figures/total_precip.gif?raw=true" alt="animated" />
</p>

### Snowfall (`snow`)  
This parameter is the accumulated snow that falls to the Earth's surface. It is the sum of large-scale snowfall and convective snowfall. Large-scale snowfall is generated by the cloud scheme in the ECMWF Integrated Forecasting System (IFS). The cloud scheme represents the formation and dissipation of clouds and large-scale precipitation due to changes in atmospheric quantities (such as pressure, temperature and moisture) predicted directly by the IFS at spatial scales of the [grid box](https://confluence.ecmwf.int/display/CKB/Model+grid+box+and+time+step) or larger. Convective snowfall is generated by the convection scheme in the IFS, which represents convection at spatial scales smaller than the grid box. See further information.

This parameter is the total amount of [water accumulated over a particular time period which depends on the data extracted](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Meanrates/fluxesandaccumulations) The units of this parameter are depth in metres of water equivalent. It is the depth the water would have if it were spread evenly over the grid box.

https://github.com/ftucciarone/Med12_old/blob/6733a2eacd2244d05870f1bbe636c14d8648cd32/preprocessing-era5/python/runproc_snowfall.py#L1-L41

<p align="center">
  <img src="https://github.com/ftucciarone/Med12_old/blob/main/preprocessing-era5/python/figures/snowfll.gif?raw=true" alt="animated" />
</p>

### Surface Solar Radiation Downwards (`swrd`)
This parameter is the amount of solar radiation (also known as shortwave radiation) that reaches a horizontal plane at the surface of the Earth. This parameter comprises both direct and diffuse solar radiation.

Radiation from the Sun (solar, or shortwave, radiation) is partly reflected back to space by clouds and particles in the atmosphere (aerosols) and some of it is absorbed. The rest is incident on the Earth's surface (represented by this parameter). See [further documentation](https://www.ecmwf.int/sites/default/files/elibrary/2015/18490-radiation-quantities-ecmwf-model-and-mars.pdf).

To a reasonably good approximation, this parameter is the model equivalent of what would be measured by a pyranometer (an instrument used for measuring solar radiation) at the surface. However, care should be taken when comparing model parameters with observations, because observations are often local to a particular point in space and time, rather than representing averages over a [model grid box](https://confluence.ecmwf.int/display/CKB/Model+grid+box+and+time+step).

This parameter [is accumulated over a particular time period which depends on the data extracted](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Meanrates/fluxesandaccumulations). The units are joules per square metre (J m-2). To convert to watts per square metre (W m-2), the accumulated values should be divided by the accumulation period expressed in seconds.

https://github.com/ftucciarone/Med12_old/blob/6733a2eacd2244d05870f1bbe636c14d8648cd32/preprocessing-era5/python/runproc_surface_solar_radiation_downwards.py#L1-L41

<p align="center">
  <img src="https://github.com/ftucciarone/Med12_old/blob/main/preprocessing-era5/python/figures/swrd.gif?raw=true" alt="animated" />
</p>

### Surface Thermal Radiation Downwards (`lwrd`)
This parameter is the amount of thermal (also known as longwave or terrestrial) radiation emitted by the atmosphere and clouds that reaches a horizontal plane at the surface of the Earth.

The surface of the Earth emits thermal radiation, some of which is absorbed by the atmosphere and clouds. The atmosphere and clouds likewise emit thermal radiation in all directions, some of which reaches the surface (represented by this parameter). See [further documentation](https://www.ecmwf.int/sites/default/files/elibrary/2015/18490-radiation-quantities-ecmwf-model-and-mars.pdf).

This parameter [is accumulated over a particular time period which depends on the data extracted](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Meanrates/fluxesandaccumulations). The units are joules per square metre (J m-2). To convert to watts per square metre (W m-2), the accumulated values should be divided by the accumulation period expressed in seconds.

https://github.com/ftucciarone/Med12_old/blob/6733a2eacd2244d05870f1bbe636c14d8648cd32/preprocessing-era5/python/runproc_surface_thermal_radiation_downwards.py#L1-L41

<p align="center">
  <img src="https://github.com/ftucciarone/Med12_old/blob/main/preprocessing-era5/python/figures/lwrd.gif?raw=true" alt="animated" />
</p>


#
https://github.com/ftucciarone/Med12_old/blob/c6f489633bd606f16a767a00ea3333fc11039eb3/preprocessing-era5/python/era5_process.py#L1-L147
