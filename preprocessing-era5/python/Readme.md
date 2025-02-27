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

