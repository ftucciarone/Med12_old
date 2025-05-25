#
# Netcdf stuff: to know how to add NetCDF links just run the command "nc-config --all"
#               and then take the following entries:
#               --fflags
#               --flink
#

#nc_config__which=/usr/bin/nc-config
nc_config__which=$( which nc-config ) || { echo "which nc-config failed: check syntax or specify directly the path."; exit ; }
nc_config__cc=$( $nc_config__which --cc) || { echo "nc-config --cc failed:check syntax or specify the compiler."; exit ; }
nc_config__cflags=$( $nc_config__which --cflags ) || { echo "nc-config --cflags failed: check syntax."; exit ; }


#nf_config__which=/usr/bin/nf-config
nf_config__which=$( which nf-config ) || { echo "which nf-config failed: check syntax or specify directly the path."; exit ; }
nf_config__fc=$( $nf_config__which --fc) || { echo "nf-config --fc failed:check syntax or specify the compiler."; exit ; }
nf_config__prefix=$( $nf_config__which --prefix) || { echo "nf-config --prefix failed:check syntax or specify the compiler."; exit ; }
nf_config__includedir=$( $nf_config__which --includedir) || { echo "nf-config --includedir failed:check syntax or specify the compiler."; exit ; }
nf_config__fflags=$( $nf_config__which --fflags ) || { echo "nf-config --fflags failed: check syntax."; exit ; }
nf_config__flibs=$( $nf_config__which --flibs) || { echo "nf-config --flibs failed:check syntax."; exit ; }


echo
echo "netCDF C version"
echo "nc-config directory: " $nc_config__which
echo "nc-config --cflags:  " $nc_config__cflags
echo "nc-config --cc:      " $nc_config__cc
echo
echo "netCDF Fortran version"
echo "nf-config directory:   " $nf_config__which
echo "nf-config --prefix:    " $nf_config__prefix
echo "nf-config --includedir:" $nf_config__includedir
echo "nf-config --fflags:    " $nf_config__fflags
echo "nf-config --flibs:     " $nf_config__flibs
echo "nf-config --fc:        " $nf_config__fc
echo


netCDF_Libs="$nf_config__fflags -L$nf_config__prefix/lib -lnetcdff"


#
# Compiler options and flags: 
#              
#
gFort_flags=" -Ofast -fbounds-check  -Wno-uninitialized -ffree-line-length-512 "
iFort_flags=""

comp_flags=$gFort_flags

COMP="$nf_config__fc $comp_flags" 

for f in T U V; do
   (set -x
   $COMP -o bdy_$f.exe create_bdy_$f.F90 $netCDF_Libs 
   )
done

