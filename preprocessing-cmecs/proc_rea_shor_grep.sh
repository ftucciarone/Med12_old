#!/bin/sh
#PBS -e /scratch/ms/it/ugy/NAT/input/
#PBS -o /scratch/ms/it/ugy/NAT/input/

set -x

#--- Environm. config.
ulimit -s unlimited
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/intel/netcdf4_hdf5/lib
# export LD_LIBRARY_PATH=/usr/local/apps/netcdf4/4.7.4/INTEL/170/lib:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$home/libs/netcdf-c-4.7.4-f-4.5.3/lib

#--- Executables
home=
exe=$home/preproc/merc2iclbc_med7km_oras/m2r_nn.x
exe_rot=$home/preproc/merc2iclbc_med7km_oras/sosie_new/bin/mycorr_vect.x

#--- Variables
set -A vars thetao so uo vo zos
set -A type 3  3  3  3   2
nvs=${#vars[*]}

#--- Input/Output dirs
wrkd=$home/preproc/merc2iclbc_med7km_oras/work_dir
outd=$wrkd/output
nams=`ls $home/preproc/merc2iclbc_med7km_oras/namelist*`

#--- Ancillary files
meshm=$home/models/inputs/MED7km/mesh_mask.nc
ftdir=ftp://my.cmems-du.eu/Core/GLOBAL_REANALYSIS_PHY_001_031/global-reanalysis-phy-001-031-grepv2-daily
ft2dir=my.cmems-du.eu/Core/GLOBAL_REANALYSIS_PHY_001_031/global-reanalysis-phy-001-031-grepv2-daily

#--- Starting/end dates
this_year=2017

#
# For each year to process
#
for this_year in {1993..2004}; do
   startm=1
   if [ $this_year -eq 2007 ]; then
      startm=10
   fi

   #
   # For each month
   #
   for imm in `seq $startm 12` ; do
      #
      # Create a variabe ym=YYYYMM
      #
      ym=${this_year}`printf "%02d" $imm`
      #
      # Create two variables yy=YYYY and mm=MM
      #
      yy=`echo $ym | cut -c 1-4`
      mm=`echo $ym | cut -c 5-6`
      #
      # Compute the number of days in the month 
      #
      if [ $(( $yy % 4 )) -eq 0 ]; then
         set -A mdays 0 31 29 31 30 31 30 31 31 30 31 30 31
      else
         set -A mdays 0 31 28 31 30 31 30 31 31 30 31 30 31
      fi
      #
      # Create two variables str=YYYYMM01 and end=YYYYMM31 (or whatever is the last day)
      # 
      str=${ym}01
      end_day=${mdays[$mm]}
      end=${ym}${end_day}

      #--- Go!
      # Create and move to output directory
      mkdir -p $outd
      cd $wrkd || exit -1
      # Link the mesh_mask file in the directory
      ln -sf $meshm .
      # Copy in the namelists
      cp $nams . || exit -1

      ymd=$str
      #
      # For all the days in the month, in format YYYMMDD
      #
      while [ $ymd -le $end ]; do

         yy=`echo $ymd | cut -c 1-4`
         mm=`echo $ymd | cut -c 5-6`
         dd=`echo $ymd | cut -c 7-8`
         #
         # get the files: but how.
         #
         wget --user=******* --password=******** \
         	-r -l1 --no-parent -A "grepv2_daily_${ymd}*" \
         	$ftdir/$yy/$mm/ || exit 3
         mv $ft2dir/$yy/$mm/grepv2_daily_${ymd}* . || exit -4

         fi=grepv2_daily_${ymd}.nc
         funrot=$outd/grepv2_MED7km_unrotated_$ymd.nc

         ln -sf $fi input.nc || exit -3
         fo=$outd/grepv2_ORAS5_MED7km_$ymd.nc
   
         for iv in `seq 0 $(( $nvs - 1 ))`; do
      
      	   time $exe namelist_${vars[$iv]} ${type[$iv]} || { echo "$exe failed, exiting!" ; exit 2 ; }
   
         	if [ $iv -eq 0 ]; then
         		mv output.nc $fo || { echo "No output, exiting!" ; exit 3 ; }
   	      else
   		      ncks -h output.nc -A $fo || { echo "No output, exiting!" ; exit 3 ; }
         	fi
   
         done
   
         rm -f input.nc output.nc

         # Rotate u/v

         cp $fo $funrot

         # Uses the software SOSIE https://github.com/brodeau/sosie
         #
         #  -m   <mesh_mask_file>  => Specify which mesh_mask file to use
         #  -G   <T/U,V>           => Specify grid points where to save the (un)rotated vector
         #  -t   <time_name>       => name of time variable in <x.nc> and <y.nc>'
         #  -i   <x.nc> <y.nc>     => unrotate vector fields given in these 2 files to the same grid
         #  -v   <nameU> <nameV>   => names for x and y components in intput files
         #
         #
         # I guess -o is an extra to define the file name
         #
         #
         time $exe_rot -m $meshm -G U -t time -i $funrot $funrot -v uo_oras vo_oras -o $fo $fo || exit -5

         rm -f $funrot

   #---

##   fo=$outd/grepv2_CGLORS_MED7km_$ymd.nc
##   ln -sf $fi input.nc || exit -3
##
##   for iv in `seq 0 $(( $nvs - 1 ))`; do
##
##        time $exe namelist2_${vars[$iv]} ${type[$iv]} || { echo "$exe failed, exiting!" ; exit 2 ; }
##
##        if [ $iv -eq 0 ]; then
##                mv output.nc $fo || { echo "No output, exiting!" ; exit 3 ; }
##        else
##                ncks -h output.nc -A $fo || { echo "No output, exiting!" ; exit 3 ; }
##        fi
##
##   done
##   rm -f input.nc output.nc 
##
##   # Rotate u/v
##
##   cp $fo $funrot
##   time $exe_rot -m $meshm -G U -t time -i $funrot $funrot -v uo_cglo vo_cglo -o $fo $fo || exit -5
##
##   rm -f $funrot
##
##   #---
##
##   fo=$outd/grepv2_FOAM_MED7km_$ymd.nc
##   ln -sf $fi input.nc || exit -3
##
##   for iv in `seq 0 $(( $nvs - 1 ))`; do
##
##        time $exe namelist3_${vars[$iv]} ${type[$iv]} || { echo "$exe failed, exiting!" ; exit 2 ; }
##
##        if [ $iv -eq 0 ]; then
##                mv output.nc $fo || { echo "No output, exiting!" ; exit 3 ; }
##        else
##                ncks -h output.nc -A $fo || { echo "No output, exiting!" ; exit 3 ; }
##        fi
##
##   done
##   rm -f input.nc output.nc
##
##   # Rotate u/v
##
##   cp $fo $funrot
##   time $exe_rot -m $meshm -G U -t time -i $funrot $funrot -v uo_foam vo_foam -o $fo $fo || exit -5
##
##   rm -f $funrot

   # --- 
         rm -f $fi

         ymd=`dateincr $ymd +1`

      done

# ksh $home/preproc/merc2iclbc_med7km_oras/create_bdy.ksh $str $end 2 || exit -2

   done

done

exit 0
