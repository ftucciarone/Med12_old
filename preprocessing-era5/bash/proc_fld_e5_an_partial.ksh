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
merge_exe=$ei_tools/merge.x
mergerad_exe=$ei_tools/merge_rad.x
mergeprec_exe=$ei_tools/merge_prec.x
qs_exe=$ei_tools/td2qs.x
flexe=$ei_tools/flandR.x
fland=0

nparams=5
set -A params  u10m   v10m    t2m   q2m    slp
set -A inv      165    166    167    qs    134
set -A ivv   var165 var166 var167    qs var134
set -A ouv     u10m   v10m    t2m   q2m    slp
set -A nts       24     24     24    24     24

mkdir -p $wrkdir

set -x

for yy in `seq $ys $ye`; do

year=$yy
d=1

cd $wrkdir || exit -2
mkdir -p ../archive

# Merge analysis fields
set -A nampar mean_sea_level_pressure 10m_u_component_of_wind 10m_v_component_of_wind 2m_temperature 2m_dewpoint_temperature

kp=0
for par in 134 165 166 167 168; do

	    namp=${nampar[$kp]}

            for t in an ; do
                forig=/home/Andrea.Storto/era5/${namp}_${yy}.grib
                ff=era5_${yy}_${par}_${t}.grb
                fn=era5_${yy}_${par}.nc
		ln -sf $forig $ff
                $cdo -f nc copy $ff $fn || exit -1
		if [ $par -eq 134 ]; then
		 $ncrename -h -v var151,var134 $fn || exit -2
		fi
		rm -f $ff
            done

	    kp=$(( $kp + 1 ))

done

nt=365
set -A ndds 0 31 28  31  30  31  30  31  31  30  31  30  31
set -A csum 0 0  31  59  90 120 151 181 212 243 273 304 334
if [ $(( $yy % 4 )) -eq 0 ]; then
	nt=366
    set -A ndds 0 31 29  31  30  31  30  31  31  30  31  30  31
    set -A csum 0 0  31  60  91 121 152 182 213 244 274 305 335
fi

# Compute specific humidity
f1=era5_${yy}_168.nc
f2=era5_${yy}_134.nc
fo=era5_${yy}_qs.nc
$ncks -h -d time,0,0 $f1 -O $fo
$ncrename -h -v var168,qs $fo
$qs_exe $f1 var168 $f2 var134 $fo qs $nx $ny $yy f f $ntmp 24 || exit -3

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
		mv $fou ../archive/

	done

	rm -f $fin
        iv=$(( $iv + 1 ))
done

rm lsm.nc
rm -f [eE]*_${yy}*nc

done

exit 0
