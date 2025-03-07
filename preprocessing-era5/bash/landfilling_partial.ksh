#!/bin/ksh

set -x
 
export OMP_NUM_THREADS=12
export KMP_STACKSIZE="200M"
export MP_TASK_AFFINITY=-1
export KMP_AFFINITY=compact

ulimit -s unlimited

flexe=/home/Andrea.Storto/tools/fill_land/flandR_OMP.x
archdir=/home/Andrea.Storto/era5/postproc/archive

ys=2024
ye=$ys
ms=1
me=11

nparams=9
set -A params  u10m   v10m    t2m   q2m    slp  snow   swrd   lwrd precip
set -A nts       24     24     24    24     24     1      1      1      1

nx=1440
ny=721
landseamask=/home/Andrea.Storto/era5/postproc/tools/lsm_ERA5_0.25.nc
copy=0

cd /home/Andrea.Storto/era5/postproc
ln -sf $landseamask lsm.nc

for yyyy in `seq $ys $ye`; do
for m in `seq $ms $me`; do

set -A ndds 0 31 28  31  30  31  30  31  31  30  31  30  31
if [ $(( $yyyy % 4 )) -eq 0 ]; then
    set -A ndds 0 31 29  31  30  31  30  31  31  30  31  30  31
fi

mm=`printf "%02d" $m`
nm=${ndds[$m]}

iv=0
while [ $iv -lt $nparams ]; do

     par=${params[$iv]}
     ns=${nts[$iv]}
     ff=ERA5_${par}_y${yyyy}m${mm}.nc
     echo "$ff (from $ys to $ye)" > doing_now
     time $flexe $archdir/$ff $par $nx $ny $(( $nm * $ns )) lsm || exit -2

     iv=$(( $iv + 1 ))

done

done
done

exit 0
