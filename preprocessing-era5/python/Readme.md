# Era5 preprocessing in python

We have three main directories (along with a `home` directory): 
- `work`: here we store temporary and intermediate files;
- `grib`: here we store the raw `.grib` files from the download;
- `arch`: here we store the processed files, those that will be fed to NEMO.
 
Each processing script will have the following declaration:
```python
# Directories
dirs = {
        "home": "/home/ftucciar",                                   # Home directory
        "work": "/home/ftucciar/Med12/preprocessing-era5/work",     # Work directory
        "grib": "/home/ftucciar/Med12/preprocessing-era5/grib",     # Grib directory
        "arch": "/home/ftucciar/Med12/preprocessing-era5/archive"   # Archive (processed)
    }
```
> [!IMPORTANT]  
> This could be set as a global variable somewhere

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
