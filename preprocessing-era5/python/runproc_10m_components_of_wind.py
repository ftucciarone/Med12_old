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
        
era5_process(param, dirs, 2020, 1, 2020, 12, cleanup=False)










