# -*- coding: utf-8 -*-



def era5_process(param, dirs, s_year, s_month, e_year, e_month, cleanup=False):
    """ 
        Main processing function for era5 atmospheric data
            
        Parameters of the input field
        param = {
                 "long_name" : str,                # Long name describing the field
                 "var_name" : str,                 # Variable ID
                 "out_name" : str,                 # Output name of the variable
                 "chr_id" : str,                   # Character version of numerical ID
                 "nts" : int,                      # Number of time steps per day
                 "nx" : int,                       # Dimensions in x
                 "ny" : int                        # Dimensions in y
                 "daymean": bool,                  # Flag to average over one day
                 "units_change": {
                                  "change": bool,  # Flag to convert into another set of units
                                  "ucf": float     # Units conversion factor
                                  },
                 "maskland": {
                              "mask": bool,        # Flag to mask the land
                              "exec": str,         # Path to executable (as in this case this is in fortran)
                              "maskfile": str      # Path to NetCDF file witht he mask 
                              }
                 }        
        
        Directories 
        dirs = {
                "home" : str,          # Home directory
                "work" : str,          # Work directory
                "grib" : str,          # Grib directory
                "arch" : str           # Archive (processed)
            }     
       
        s_year : int,                  # First year to process
        e_year : int,                  # Last year to process
       
        month_s : int,                 # First month to process
        month_e : int,                 # Last (full) month to process
       
        
    """
    import os
    import shutil
    import subprocess
    from calendar import monthrange
    from datetime import date, timedelta

    # Climate Data Operators (CDO) wrappers
    import cdo

    # NetCDF Operator (NCO) python wrappers 
    from nco import Nco
    from nco.custom import Rename
    
    # Directories
    home = dirs["home"]  # Home directory
    work = dirs["work"]  # Work directory
    grib = dirs["grib"]  # Grib directory
    arch = dirs["arch"]  # Archive (processed)
    
    # Initialize wrappers for CDO and NCO
    cdo = cdo.Cdo(tempdir=work)
    nco = Nco()
    
    # Compute the number of days of the interval of months
    date_s = date(s_year, s_month, 1)
    date_e = date(e_year, e_month, monthrange(e_year, e_month)[-1]) + timedelta(days=1)
    delta = date_e - date_s
    # print(delta.days)

    # Cycle over the years
    for year in range(s_year, e_year+1):
        
        # Control print
        print(year)
        
        # Set in and out files
        fgrib = grib + "/" + param["long_name"] + "_" + str(year) + ".grib"
        fwork = work + "/" + param["long_name"] + "_" + str(year) + "_fc.nc"
        fnCDF = work + "/era5_" + str(year) + "_" + param["chr_id"] + ".nc"
        
        # Convert grib file into netCDF file using CDO (Climate Data Operators)
        cdo.copy(input=fgrib, output=fwork, options="-f nc")
        
        # Compute daily mean
        if "daymean" in param.keys() and param["daymean"] == True:
            cdo.daymean(input=fwork, output=fnCDF, options="")
        else:
            fnCDF=fwork
        
        # % Rescale the data 
        if "units_change" in param.keys() and param["units_change"]["change"] == True:
            operation = param["var_name"] + " = float(" + param["var_name"] + "/" + str(param["units_change"]["ucf"]) + ");"
            print(operation)
            nco.ncap2(input=fnCDF, output=fnCDF, options=["-h -s '"+operation + "' -O "])
        
        
        # Cleanup work file if a) we want and b) the file exists
        if cleanup and os.path.exists(fwork): os.remove(fwork)
    
        # Cycle over the months
        first_day = 0
        for month in range(s_month, e_month+1):
            
            # Create output file (definitive one)
            out_file = "/" + "ERA5_" + param["out_name"] + "_y" + str(year) + "m" + str(month).zfill(2) + ".nc"
        
            # Control print
            print(arch + out_file)
            
            # Slabbing one month
            dim = "time"
            istart = "," + str(first_day+1)
            istop = "," + str(first_day+param["nts"]*monthrange(e_year, e_month)[-1])
            nco.ncks(input=fnCDF, output=work+out_file, options=["-h -F -d "+dim+istart+istop])
        
            # % Rename variable
            # nco.ncrename(input=work+out_file, options=["-h", Rename("v", {""+param["var_name"]: param["out_name"]}) ])
            subprocess.run(["/usr/bin/ncrename -h -v " + param["var_name"] +","+param["out_name"] + " " + work+out_file], shell=True)
        
        
            # Mask land with fortran code
            if "maskland" in param.keys() and param["maskland"]["mask"] == True:
                fortexe = param["maskland"]["exec"]
                maskfile = param["maskland"]["maskfile"]
                era5_maskland(fortexe, work+out_file, param["out_name"], param["nx"], param["ny"], param["nts"]*monthrange(e_year, e_month)[-1], maskfile)
        
            # Move variable to archive
            shutil.move(work+out_file, arch+out_file)    

# %%
def era5_spHumidity(mslp, dewt, dirs, s_year, s_month, e_year, e_month, cleanup=False):
    """
    """
    
# %%
def era5_maskland(fortexe, out, varname, nx, ny, nts, maskfile):
    """
    """  
    import subprocess
    excomm = fortexe + " " + out + " " + varname + " " + str(nx) + " " + str(ny) + " " + str(nts) + " lsm"
    subprocess.run(["ln -sf " + maskfile + " lsm.nc"], shell=True) 
    subprocess.run([excomm], shell=True)
