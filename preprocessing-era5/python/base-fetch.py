# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

import numpy as np
from datetime import date, timedelta
from calendar import monthrange

import os, shutil

import cdo

from nco import Nco
from nco.custom import Rename

param = {
         "long_name": "10m_u_component_of_wind",
         "short_name": "u10m",
         "num_id":  165,
         "chr_id": "165",
         "var_id": "var165",
         "ouv": "u10m",
         "nts": 24,
         "cfs": None,
         "nx": 1440,
         "ny": 721
    }

years_s = 2020 # Starting year
years_e = 2020 # Ending year

month_s = 1  # Starting month
month_e = 1  # Ending month, up to its last day (so 1-1 gives you January,
             #   1-12 gives you the whole year)

home="/home/ftucciar"                                   # Home directory
work="/home/ftucciar/Med12/preprocessing-era5/work"     # Work directory
grib="/home/ftucciar/Med12/preprocessing-era5/grib"     # Grib directory
arch="/home/ftucciar/Med12/preprocessing-era5/archive"  # Archive (processed)


date_s = date(years_s, month_s, 1)
date_e = date(years_e, month_e, monthrange(years_e, month_e)[-1]) + timedelta(days=1)
delta = date_e - date_s
print(delta.days)

cdo = cdo.Cdo(tempdir=work)
nco = Nco()

# process all the the years
for year in range(years_s, years_e+1):
    
    # Control print
    print(year)
    
    # Set in and out files
    fgrib = grib + "/" + param["long_name"] + "_" + str(year) + ".grib"
    fnCDF = work + "/era5_" + str(year) + "_" + param["chr_id"] + ".nc"
    
    # Convert grib file into netCDF file using CDO (Climate Data Operators)
    cdo.copy(input=fgrib, output=fnCDF, options="-f nc")
    
    first_day = 0
    for month in range(month_s, month_e+1):
        
        # Create output file (definitive one)
        out_file = "/" + "ERA5_" + param["ouv"] + "_y" + str(year) + "m" + str(month).zfill(2) + ".nc"
    
        # Control print
        print(arch + out_file)
        
        # Slabbing one month
        dim = "time"
        istart = "," + str(first_day+1)
        istop = "," + str(first_day+param["nts"]*monthrange(years_e, month_e)[-1])
        nco.ncks(input=fnCDF, output=work+out_file, options=["-h -F -d "+dim+istart+istop])
    
        # %% Rename variable
        nco.ncrename(input=work+out_file, options=["-h", Rename("d", {"."+param["var_id"]: param["ouv"]}) ])
    
        # Move variable to archive
        shutil.move(work+out_file, arch+out_file)
