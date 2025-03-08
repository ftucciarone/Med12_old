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
