#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np
import math
from netCDF4 import Dataset
import json

# Directories
dirs = json.load( open('directories.json') )

# Initialise MPI and get rank and nprocs
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
# Print out some info. Always print with a single process, otherwise the printing
# is messy as there is no order between the processes.
if rank ==0:
    print("mpi version is ", MPI.Get_version())
    print("Employed cores:", comm.Get_size(), "/", comm.Get_attr(MPI.UNIVERSE_SIZE))
    
# some proxy data
ndata = 31

slicelen = int(math.ceil(ndata/nprocs))
# some infos
if rank == 0:
    print("Every processor will process", slicelen, " elements")
# Define start and stop
start = rank * slicelen
stop = min(( rank + 1 ) * slicelen, ndata)
# Here I am printing with all processes. As you will see, it is messy.


# %%
filepath = "/home/ftucciar/Med12/preprocessing-era5/archive/ERA5_snow_y2020m01.nc"
varname = "snow"
ncfile = Dataset(filepath, "r", format="NETCDF4")

fig, ax = plt.subplots(figsize=(7, 4), dpi=300)
for s in range(start, stop):
    print(s)
    field = np.squeeze(ncfile.variables[varname][s, :, :].data)
    ax.imshow(field[:,:], cmap='RdBu_r', vmin=0, vmax=0.0002)
    plt.title(r"Snowfall")
    plt.axis('off')
    filename = dirs["figs"] + "/"+ varname + str(s).zfill(5) + ".png"
    #plt.show()
    plt.savefig(filename, format='png' ,
               bbox_inches='tight', 
               transparent=True,
               pad_inches=0)
