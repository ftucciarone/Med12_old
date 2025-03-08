#!/bin/ksh
set -x

COMP="ifort"
NCDF=/opt/intel/netcdf-c-4.7.4-f-4.5.3-impi
ADD=""

#
# Netcdf stuff: to know how to add NetCDF links just run the command "nc-config --all"
#               and then take the following entries:
#               --fflags
#               --flink
#
#  --fflags    -> -I/usr/local/include -I/usr/local/include
#  --flibs     -> -L/usr/local/Cellar/netcdf-fortran/4.6.1/lib -lnetcdff
#
#
nc_config__fflags="-I/usr/include "
nc_config__flibs="-L/usr/lib/x86_64-linux-gnu -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf -ldl -lm"

COMP="gfortran -Ofast -fbounds-check -Wall -Wno-uninitialized -ffree-line-length-512 "

str=""
for f in master_flandr.F90 read_field2dt_fillLand.F90 read_msk.F90 handlerr.F90 write_field2dt_fillLand.F90 fill_land.F90 ; do
        $COMP -c -I$NCDF/include $f
str="$str `echo $f | cut -d. -f1`.o"
done

$COMP $str -o flandR.x -L$NCDF/lib -lnetcdf
