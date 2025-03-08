#!/bin/ksh

NCDF=/usr/local/apps/netcdf/3.6.3/GNU/5.3.0
COMP="gfortran"

str=""
for f in master_flandr.F90 read_field2dt_fillLand.F90 read_msk.F90 handlerr.F90 write_field2dt_fillLand.F90 fill_land.F90 ; do
        $COMP -c -I$NCDF/include $f
str="$str `echo $f | cut -d. -f1`.o"
done

$COMP $str -o flandR.x -L$NCDF/lib -lnetcdf
