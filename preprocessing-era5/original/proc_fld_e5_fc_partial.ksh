#!/bin/ksh

ys=2024
ye=2024
ms=1
me=11
if [ $me -lt 12 ]; then
  ntmp=$( datediff ${ye}`printf "%02d" $(( $me + 1 ))`01 ${ye}`printf "%02d" $ms`01 )
else
  ntmp=$( datediff $(( $ye + 1 ))0101 ${ye}`printf "%02d" $ms`01 )
fi

nx=1440
ny=721
wrkdir=/home/Andrea.Storto/era5/postproc/work
cdo=cdo
ncod=/usr/bin
ncpdq=$ncod/ncpdq
ncks=$ncod/ncks
ncbo=$ncod/ncbo
ncdiff=$ncod/ncdiff
ncap2=$ncod/ncap2
ncwa=$ncod/ncwa
ncecat=$ncod/ncecat
ncrcat=$ncod/ncrcat
ncrename=$ncod/ncrename
ei_tools=/home/Andrea.Storto/era5/postproc/tools
landseamask=$ei_tools/lsm_ERA5_0.25.nc
flexe=$ei_tools/flandR.x
fland=0

nparams=4

set -A params snow   swrd   lwrd precip 
set -A inv     144    169    175    228
set -A ivv  var144 var169 var175 var228
set -A ouv    snow   swrd   lwrd precip
set -A nts       1      1      1      1   
set -A cfs     3.6  3600.  3600.    3.6   
set -A nampar snowfall surface_solar_radiation_downwards surface_thermal_radiation_downwards total_precipitation

mkdir -p $wrkdir
cd $wrkdir || exit -2
mkdir -p ../archive
set -x

nt=365
set -A ndds 0 31 28  31  30  31  30  31  31  30  31  30  31
set -A csum 0 0  31  59  90 120 151 181 212 243 273 304 334
if [ $(( $yy % 4 )) -eq 0 ]; then
	nt=366
    set -A ndds 0 31 29  31  30  31  30  31  31  30  31  30  31
    set -A csum 0 0  31  60  91 121 152 182 213 244 274 305 335
fi


for yy in `seq $ys $ye`; do

   year=$yy
   d=1

   # Merge analysis fields
   iv=0
   while [ $iv -lt $nparams ]; do

		par=${inv[$iv]}
		coe=${cfs[$iv]}
      namp=${nampar[$iv]}

		forig=/home/Andrea.Storto/era5/${namp}_${yy}.grib
      ff=era5_${yy}_${par}.grib
      fn=era5_${yy}_${par}_fc.nc
      fno=era5_${yy}_${par}.nc
		ln -sf $forig $ff
      $cdo -f nc copy $ff $fn || exit -1
		$cdo daymean $fn $fno || exit 3
		vt="var$par"
		$ncap2 -h -s "$vt = float ($vt / $coe );" $fno -O $fno || exit -1
		$ncks -h -F -d time,1,${ntmp} $fno -O $fno || exit -2
		rm -f $ff $fn

      iv=$(( $iv + 1 ))

   done



   ln -sf $landseamask lsm.nc

   iv=0
   while [ $iv -lt $nparams ]; do
      fin=era5_${yy}_${inv[$iv]}.nc
      ns=${nts[$iv]}
	   ien=0
	   for m in `seq $ms $me`; do
		   mm=`printf "%02d" $m`
         nm=${ndds[$m]}
			
         fou=ERA5_${ouv[$iv]}_y${yy}m${mm}.nc
         ist=$(( $ien + 1 ))
         ien=$(( $ien + ( $nm * $ns ) ))

         $ncks -h -F -d time,${ist},${ien} $fin $fou || exit -1
         $ncrename -h -v ${ivv[$iv]},${ouv[$iv]} $fou || exit -2
         vv=${ouv[$iv]}
		   if [ $fland -eq 1 ]; then
            $flexe $fou $vv $nx $ny $(( $nm * $ns )) lsm || exit -1
		   fi
		   mv $fou ../archive/ || exit 2
	   done
	   rm -f $fin
      iv=$(( $iv + 1 ))
   done


   rm lsm.nc

done

exit 0
