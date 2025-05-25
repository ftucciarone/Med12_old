source netCDF.macro
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
   $COMP -o create_bdy/bdy_$f.exe create_bdy/create_bdy_$f.F90 $netCDF_Libs 
   )
done

