#!/bin/bash
source /home/ftucciar/Stockage12T/Med12/preprocessing-era5/bash/params.sh


# Argument validation check
if [ "$#" -lt 5 ]; then
    echo " "
    echo "Usage: $0 <long_name> <var_name> <out_name> <chr_id> <nts> <daymean>            "
    echo " "
    echo "Parameters depending on the field processed to be passed                        "
    echo " " 
    echo "  Refer to: https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation "
    echo " "
    echo " long_name: Variable name in CDS, e.g. '10m_u_component_of_wind'                "
    echo "  var_name: ShortName, e.g. 'u10m'                                              "
    echo "  out_name: Output name of the variable, e.g. 'u10m'                            "
    echo "    chr_id: paramID, e.g. '167'                                                 "
    echo "       nts: number of time snapshots downloaded                                 "
    echo "   daymean: optional, set to 1 to average over one day                          "
    echo " "
    echo "Usage: $0 <long_name> <var_name> <out_name> <chr_id> <nts> <daymean>            "
    echo " "
    echo "  Examples:"
    echo " " 
    echo " 10m_u_component_of_wind 10u u10m 165 24 " 
    echo " 10m_v_component_of_wind 10v v10m 166 24 " 

    exit 1
fi

# Read positional parameters defining the field
long_name=$1
var_name=$2
out_name=$3
chr_id=$4
nts=$5

# Preset for skipdownload, will be overridden by commandline parameters
skipdownload=0
# Read optional parameters
# Careful that s does not require an entry parameter and 
# this is done by removing the colon
shift $((OPTIND + 4))
while getopts d:u:f:s flag
do
    case "${flag}" in
        d) dayavg=${OPTARG:-0};;
        u) uchng=${OPTARG:-1};;
        f) filland=${OPTARG:-0};;
        s) skipdownload=1;;
    esac
done
dayavg=${dayavg:-0}
uchng=${uchng:-1.}
filland=${filland:-0}


echo " " 
echo "    Long name: $long_name";
echo "Variable name: $var_name";
echo "  Output name: $out_name";
echo "  Variable ID: $chr_id";
echo "    Timesteps: $nts (downloaded with python script)";
echo " " 
echo "Daily average: $dayavg (if 1, True)";
echo " Change units: $uchng (if neq 1. performs var/coef, do not set to 0.)";
echo "    Fill Land: $filland (if 1, execute fortran code )";
echo " " 

# build the list of years and compute the numbers of years to process
lyears=($(seq $ys $(( $ye-1 > $ys ? $ye : $ys ))))
nyears=$(( $(( $ye + 1 )) - $ys ))

var_label="var$chr_id"

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
   out=${grib_dir}/${long_name}_${year}.grib
   scr=${run_dir}/${long_name}_${year}.fetch
   sed -e "s|_VARIABLE_|$long_name|g" \
       -e "s|_YEAR_|${year}|g" \
       -e "s|_LIST_MONTHS_|${lm}|g" \
       -e "s|_OUTPUT_|$out|g" \
       run_sedf_ep > $scr
   chmod +x $scr
   echo "Download script prepared: " $scr

   # Launch the download of the files  
   scr=${run_dir}/${long_name}_${year}.fetch
   if [ $skipdownload -eq 0 ]; then
      $scr || exit -1
   fi

   # Conversion of the .grib file into netcdf
   in=${grib_dir}/${long_name}_${year}.grib
   out=${work_dir}/era5_${year}_${var_label}.nc
   echo $cdo -f nc copy $in $out
   $cdo -f nc copy $in $out || exit -1

   # If necessary, average over each day
   if [ $dayavg -ne 0 ]; then
      echo " Daymean"
      echo $cdo daymean $out $out || exit 3
      nts=1
   fi

   # If necessary, apply a scale factor and change the units
   if [ "$(echo "$uchng != 1" | bc)" = 1 ] ; then
      echo " Scale and change units"
      vt=$var_name
      echo $ncap2 -h -s "$vt = float ($vt / $uchng );" $out -O $out
      $ncap2 -h -s "$vt = float ($vt / $uchng );" $out -O $out || exit -1

      echo $ncks -h -F -d time,1,${ndays} $out -O $out
      $ncks -h -F -d time,1,${ndays} $out -O $out || exit -2

   fi

   # Slice the netcdf into months and execute fortran 
   fin=${work_dir}/era5_${year}_${var_label}.nc

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

       echo $ncks -h -F -d time,${ist},${ien} $fin $fou
       ncks -h -F -d time,${ist},${ien} $fin $fou || exit -1
       
       echo $ncrename -h -v $var_name,$out_name $fou 
       $ncrename -h -v $var_name,$out_name $fou || exit -2

       if [ $filland -eq 1 ]; then
           echo " fill land"
           echo $flexe $fou $out_name $nx $ny $(( $nm * $nts )) lsm
           $flexe $fou $out_name $nx $ny $(( $nm * $nts )) lsm || exit -1
       fi
   done

done

exit 0

