# Era5 preprocessing in python

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


## Parameters dictionary for the processing
```python
param = {
         "long_name" : str,           # Long name describing the field
         "var_name" : str,            # Variable ID
         "out_name" : str,            # Output name of the variable
         "chr_id" : str,              # Character version of numerical ID
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
## 10m wind $u,v-$ components
10m u component of wind
m s-1
Magnitude of the eastward component of the two-dimensional horizontal air velocity near the surface.
10m v component of wind
m s-1
Magnitude of the northward component of the two-dimensional horizontal air velocity near the surface.

```python
# -*- coding: utf-8 -*-
import json
from era5_process import era5_process

# Directories
dirs = json.load( open('directories.json') )

# Parameters of the input field
param = {
         # Parameters depending on the field processed
         "long_name": "variable_long_name",
         "var_name": "var_nam_in_grib_file",
         "out_name": "output_name",
         "chr_id": "variable_ID_in_grib_file",
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
