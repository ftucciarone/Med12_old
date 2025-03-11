#!/bin/bash
source /home/ftucciar/Stockage12T/Med12/preprocessing-era5/bash/params.sh




./fetch-and-process.sh             10m_u_component_of_wind  10u   u10m 165 24 -d 0 -u  1.   -f 1 -s
./fetch-and-process.sh             10m_v_component_of_wind  10v   v10m 166 24 -d 0 -u  1.   -f 1 -s

./fetch-and-process.sh                      2m_temperature   2t    t2m 167 24 -d 0 -u  1.   -f 1 -s
 
./fetch-and-process.sh                 total_precipitation   tp precip 228  1 -d 1 -u  3.6  -f 1 -s
./fetch-and-process.sh                            snowfall   sf   snow 144  1 -d 1 -u  3.6  -f 1 -s
./fetch-and-process.sh   surface_solar_radiation_downwards ssrd   swrd 169  1 -d 1 -u 3600. -f 1 -s
./fetch-and-process.sh surface_thermal_radiation_downwards strd   lwrd 175  1 -d 1 -u 3600. -f 1 -s



