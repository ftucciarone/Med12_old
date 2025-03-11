#!/bin/bash
source /home/ftucciar/Stockage12T/Med12/preprocessing-era5/bash/params.sh

# 
long_name_msl=mean_sea_level_pressure
var_name_msl=msl
chr_id_msl=134
# 
long_name_dpt=2m_dewpoint_temperature
var_name_dpt=2d
chr_id_dpt=168

out_name=qs
nts=24
# Preset for skipdownload, will be overridden by commandline parameters
skipdownload=0
# Read optional parameters
# Careful that s does not require an entry parameter and 
# this is done by removing the colon
shift $((OPTIND + 4))
while getopts f:s flag
do
    case "${flag}" in
        f) filland=${OPTARG:-0};;
        s) skipdownload=1;;
    esac
done
dayavg=${dayavg:-0}
uchng=${uchng:-1.}
filland=${filland:-0}


echo " Mean sea level pressure " 
echo " "
echo "    Long name: $long_name_msl";
echo "Variable name: $var_name_msl";
echo "  Variable ID: $chr_id_msl";
echo " " 
echo " 2 metres dewpoint temperature " 
echo " "
echo "    Long name: $long_name_dpt";
echo "Variable name: $var_name_dpt";
echo "  Variable ID: $chr_id_dpt";
echo " " 
echo "    Fill Land: $filland (if 1, execute fortran code )";
echo " " 


# build the list of years and compute the numbers of years to process
lyears=($(seq $ys $(( $ye-1 > $ys ? $ye : $ys ))))
nyears=$(( $(( $ye + 1 )) - $ys ))

var_label_msl="var$chr_id_msl"
var_label_dpt="var$chr_id_dpt"


for y in `seq 1 1 $nyears`; do

   year=${lyears[$(($y-1))]}
   # Trick to handle multiple years and 
   # get the correct months
   if [ $y == 1 ]; then
      s=$ms
   else
      s=1
   fi
   if [ $y == $nyears ]; then
      e=$me
   else
      e=12
   fi

   # Compute the number of days in this month with calendar operations 
   ndays=$( $datediff ${year}-`printf "%02d-" $s`01  ${year}-`printf "%02d-" $(( $e + 1 ))`01 )

   # Here I am creating the list of month to be substituted
   # into the python dictionary of parameters, so it will be
   # something like "01", "02", "03", ...
   lm=""
   for m in `seq $s $(( $e - 1 ))`; do
      mm=`printf "%02d" $m`
      lm="${lm}\"${mm}\","
   done
   mm=`printf "%02d" $e`
   lm="${lm}\"${mm}\""

   # Now prepare download script based on python API   
   out=${grib_dir}/${long_name_msl}_${year}.grib
   scr=${run_dir}/${long_name_msl}_${year}.fetch
   sed -e "s|_VARIABLE_|$long_name_msl|g" \
       -e "s|_YEAR_|${year}|g" \
       -e "s|_LIST_MONTHS_|${lm}|g" \
       -e "s|_OUTPUT_|$out|g" \
       run_sedf_ep > $scr
   chmod +x $scr
   echo "Download script prepared: " $scr

   # Launch the download of the files  
   scr=${run_dir}/${long_name_msl}_${year}.fetch
   if [ $skipdownload -eq 0 ]; then
      $scr || exit -1
   fi

   # Now prepare download script based on python API   
   out=${grib_dir}/${long_name_dpt}_${year}.grib
   scr=${run_dir}/${long_name_dpt}_${year}.fetch
   sed -e "s|_VARIABLE_|$long_name_dpt|g" \
       -e "s|_YEAR_|${year}|g" \
       -e "s|_LIST_MONTHS_|${lm}|g" \
       -e "s|_OUTPUT_|$out|g" \
       run_sedf_ep > $scr
   chmod +x $scr
   echo "Download script prepared: " $scr

   # Launch the download of the files  
   scr=${run_dir}/${long_name_dpt}_${year}.fetch
   if [ $skipdownload -eq 0 ]; then
      $scr || exit -1
   fi


   # Conversion of the .grib file into netcdf
   in=${grib_dir}/${long_name}_${year}.grib
   out=${work_dir}/era5_${year}_${var_label}.nc

   echo $cdo -f nc copy ${grib_dir}/${long_name_msl}_${year}.grib ${work_dir}/era5_${year}_${var_label_msl}.nc
   $cdo -f nc copy ${grib_dir}/${long_name_msl}_${year}.grib ${work_dir}/era5_${year}_${var_label_msl}.nc

   echo $cdo -f nc copy ${grib_dir}/${long_name_dpt}_${year}.grib ${work_dir}/era5_${year}_${var_label_dpt}.nc
   $cdo -f nc copy ${grib_dir}/${long_name_dpt}_${year}.grib ${work_dir}/era5_${year}_${var_label_dpt}.nc


   # Compute specific humidity
   f1=${work_dir}/era5_${year}_${var_label_dpt}.nc
   f2=${work_dir}/era5_${year}_${var_label_msl}.nc
   fo=${work_dir}/era5_${year}_qs.nc

   echo $ncks -h -d time,0,0 $f1 -O $fo
   $ncks -h -d time,0,0 $f1 -O $fo

   echo $ncrename -h -v $var_name_dpt,$out_name $fo
   $ncrename -h -v $var_name_dpt,$out_name $fo

   echo $qs_exe $f1 $var_name_dpt $f2 $var_name_msl $fo $out_name $nx $ny $year f f $ndays 24 
   $qs_exe $f1 $var_name_dpt $f2 $var_name_msl $fo $out_name $nx $ny $year f f $ndays 24 || exit -3



   # make a link here to the landsea mask because it is statically coded in 
   # the fortran code
   ln -sf $landseamask lsm.nc
    
   ien=0
   # for each month from start to end
   for m in `seq $s $e`; do
       # Format the current month for both datediff and label
       mm=`printf "%02d" $m`
       # Compute the number of days in this month with calendar operations 
       nm=$( $datediff ${year}-$mm-01  ${year}-`printf "%02d-" $(( $m + 1 ))`01 )
       # Set up the provisionary output file (mind the CAPS)		
       fou=${arch_dir}/ERA5_${out_name}_y${year}m${mm}.nc
       ist=$(( $ien + 1 ))
       ien=$(( $ien + ( $nm * $nts ) ))

       echo $ncks -h -F -d time,${ist},${ien} $fo $fou
       ncks -h -F -d time,${ist},${ien} $fo $fou || exit -1
       
       echo $ncrename -h -v $out_name,q2m $fou 
       $ncrename -h -v $out_name,q2m $fou || exit -2

       if [ $filland -eq 1 ]; then
           echo " fill land"
           echo $flexe $fou q2m $nx $ny $(( $nm * $nts )) lsm
           $flexe $fou q2m $nx $ny $(( $nm * $nts )) lsm || exit -1
       fi
   done

done

exit 0

