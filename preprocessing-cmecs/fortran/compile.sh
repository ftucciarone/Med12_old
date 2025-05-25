#!/bin/bash

set -x

NETCDF_HOME=/usr/local/apps/netcdf4/4.7.4/INTEL/170
NETCDF_HOME=/home/Andrea.Storto/libs/netcdf-c-4.7.4-f-4.5.3
NETCDF_HOME=/opt/intel/netcdf-c-4.7.4-f-4.5.3-impi

ifort -O2 -o bdy_T.exe create_bdy_T.F90 -L${NETCDF_HOME}/lib -lnetcdff -lnetcdf -I${NETCDF_HOME}/include
ifort -O2 -o bdy_U.exe create_bdy_U.F90 -L${NETCDF_HOME}/lib -lnetcdff -lnetcdf -I${NETCDF_HOME}/include
ifort -O2 -o bdy_V.exe create_bdy_V.F90 -L${NETCDF_HOME}/lib -lnetcdff -lnetcdf -I${NETCDF_HOME}/include

exit 0
