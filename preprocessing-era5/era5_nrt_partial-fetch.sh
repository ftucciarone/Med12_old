#!/bin/sh

ys=2024
ms=1
me=11

ls_vars="10m_u_component_of_wind 10m_v_component_of_wind 2m_dewpoint_temperature 2m_temperature mean_sea_level_pressure surface_solar_radiation_downwards surface_thermal_radiation_downwards total_precipitation snowfall"

lm=""
for m in `seq $ms $(( $me - 1 ))`; do
   mm=`printf "%02d" $m`
   lm="${lm}\"${mm}\","
done
mm=`printf "%02d" $me`
lm="${lm}\"${mm}\""

for y in $ys; do
   for var in $ls_vars; do

      out=grib/${var}_${y}.grib
      scr=run/${var}_${y}.run
      sed -e "s!_VARIABLE_!$var!g" \
          -e "s!_YEAR_!$y!g" \
          -e "s!_LIST_MONTHS_!${lm}!g" \
          -e "s!_OUTPUT_!$out!g" \
      base-fetch.py > $scr

      chmod +x $scr
      ./$scr || exit -1
   done
done

exit 0
