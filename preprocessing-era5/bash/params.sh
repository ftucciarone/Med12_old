#!/bin/sh

# Temporal parameters 
#
#   ys: start year
#   ye: end year
#   ms: start month
#   me: end month
#
#   note: to do a single month you set me=ms
#         not me=ms+1, as we always go to the 
#         end of the month (me)
ys=2021
ye=2021
ms=1
me=1

# Spatial parameters
#
#   nx: points in x (longitude)
#   ny: points in y (latitude)
nx=1440
ny=721

# Fill land
#
#  0:  no
#  1: yes
fland=0

# Directories set
base_dir=
grib_dir=
work_dir=
arch_dir=
tools_dir=

# 
landseamask=$tools_dir/lsm_ERA5_0.25.nc
merge_exe=$tools_dir/merge.x
mergerad_exe=$tools_dir/merge_rad.x
mergeprec_exe=$tools_dir/merge_prec.x
qs_exe=$tools_dir/td2qs.x
flexe=$tools_dir/flandR.x

# CDO and NC specification
i_know_them=1
if [ $i_know_them == 1 ]; then
	cdo=/usr/local/bin/cdo
   ncap2=/usr/local/bin/ncap2
   ncatted=/usr/local/bin/ncatted
   ncbo=/usr/local/bin/ncbo
   ncchecker=/usr/local/bin/ncchecker
   ncclimo=/usr/local/bin/ncclimo
   ncdiff=/usr/local/bin/ncdiff
   ncdu=/usr/local/bin/ncdu
   ncea=/usr/local/bin/ncea
   ncecat=/usr/local/bin/ncecat
   nces=/usr/local/bin/nces
   ncflint=/usr/local/bin/ncflint
   ncks=/usr/local/bin/ncks
   ncra=/usr/local/bin/ncra
   ncrcat=/usr/local/bin/ncrcat
   ncremap=/usr/local/bin/ncremap
   ncrename=/usr/local/bin/ncrename
   ncview=/usr/local/bin/ncview
   ncwa=/usr/local/bin/ncwa
   ncz2psx=/usr/local/bin/ncz2psx

   return
else
   cdo=$( which cdo )             || { echo "which cdo failed: check syntax, installation or specify the path."; exit ; }
   ncap2=$( which ncap2 )         || { echo "which ncap2 failed: check syntax, installation or specify the path."; exit ; }
   ncatted=$( which ncatted )     || { echo "which ncatted failed: check syntax, installation or specify the path."; exit ; }
   ncbo=$( which ncbo )           || { echo "which ncbo failed: check syntax, installation or specify the path."; exit ; }
   ncchecker=$( which ncchecker ) || { echo "which ncchecker failed: check syntax, installation or specify the path."; exit ; }
   ncclimo=$( which ncclimo )     || { echo "which ncclimo failed: check syntax, installation or specify the path."; exit ; }
   ncdiff=$( which ncdiff )       || { echo "which ncdiff failed: check syntax, installation or specify the path."; exit ; }
   ncdu=$( which ncdu )           || { echo "which ncdu failed: check syntax, installation or specify the path."; exit ; }
   ncea=$( which ncea )           || { echo "which ncea failed: check syntax, installation or specify the path."; exit ; }
   ncecat=$( which ncecat )       || { echo "which ncecat failed: check syntax, installation or specify the path."; exit ; }
   nces=$( which nces )           || { echo "which nces failed: check syntax, installation or specify the path."; exit ; }
   ncflint=$( which ncflint )     || { echo "which ncflint failed: check syntax, installation or specify the path."; exit ; }
   ncks=$( which ncks )           || { echo "which ncks failed: check syntax, installation or specify the path."; exit ; }
   ncra=$( which ncra )           || { echo "which ncra failed: check syntax, installation or specify the path."; exit ; }
   ncrcat=$( which ncrcat )       || { echo "which ncrcat failed: check syntax, installation or specify the path."; exit ; }
   ncremap=$( which ncremap )     || { echo "which ncremap failed: check syntax, installation or specify the path."; exit ; }
   ncrename=$( which ncrename )   || { echo "which ncrename failed: check syntax, installation or specify the path."; exit ; }
   ncview=$( which ncview )       || { echo "which ncview failed: check syntax, installation or specify the path."; exit ; }
   ncwa=$( which ncwa )           || { echo "which ncwa failed: check syntax, installation or specify the path."; exit ; }
   ncz2psx=$( which ncz2psx )     || { echo "which ncz2psx failed: check syntax, installation or specify the path."; exit ; }

   echo "cdo=$( which cdo )"
   echo "ncap2=$( which ncap2 ) "
   echo "ncatted=$( which ncatted ) "
   echo "ncbo=$( which ncbo ) "
   echo "ncchecker=$( which ncchecker ) "
   echo "ncclimo=$( which ncclimo ) "
   echo "ncdiff=$( which ncdiff ) "
   echo "ncdu=$( which ncdu ) "
   echo "ncea=$( which ncea ) "
   echo "ncecat=$( which ncecat ) "
   echo "nces=$( which nces ) "
   echo "ncflint=$( which ncflint ) " 
   echo "ncks=$( which ncks ) "
   echo "ncra=$( which ncra ) "
   echo "ncrcat=$( which ncrcat ) "
   echo "ncremap=$( which ncremap ) "
   echo "ncrename=$( which ncrename ) "
   echo "ncview=$( which ncview ) "
   echo "ncwa=$( which ncwa ) "
   echo "ncz2psx=$( which ncz2psx ) "
   
   return
fi

