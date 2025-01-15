# Atmospheric forcing: ERA5 data

ERA5 (atmospheric forcing) is downloaded from Copernicus Climate Data Store https://cds.climate.copernicus.eu/

The atmospheric data will be processed with the following workflow (checked box means tested)
- [x] Install Climate Data Store API 
- [ ] Download using fetch* script
- [ ] process analysis using the proc_fld* script
- [ ] process forecast using the proc_fld* script
- [ ] landfilling with landfilling*script








## Climate Data Store API, [API](https://cds.climate.copernicus.eu/how-to-api)
First, [register](https://accounts.ecmwf.int/auth/realms/ecmwf/login-actions/registration?client_id=cds&tab_id=w2gu_uF6V0o) to ECMWF and go to the [API web site](https://cds.climate.copernicus.eu/how-to-api) to [log in](https://accounts.ecmwf.int/auth/realms/ecmwf/login-actions/authenticate?client_id=cds&tab_id=w2gu_uF6V0o). Once logged in, the box that previously contained
```
url: https://cds.climate.copernicus.eu/api
key: <PERSONAL-ACCESS-TOKEN>
```
is now completed with a personal alphanumeric token. Copy the code (the whole box) displayed to the file `$HOME/.cdsapirc` (in your Unix/Linux environment). You can now install the CDS API client via the package management system pip, by running on Unix/Linux the command below.
```
pip install "cdsapi>=0.7.2"
```
Once the CDS API client is installed, it can be used to request data from the datasets listed in the CDS, ADS, ECDS and CEMS Early Warning DS catalogues. The basic template to download the data is the following:
```python
import cdsapi
client = cdsapi.Client()
# definition
dataset = "<DATASET-SHORT-NAME>"
request = {<SELECTION-REQUEST>}
target = "<TARGET-FILE>"
# Request
client.retrieve(dataset, request, target)
```
An example of downloading procedure is the following.
```python
import cdsapi
client = cdsapi.Client()
dataset = 'reanalysis-era5-pressure-levels'
request = {
    'product_type': ['reanalysis'],
    'variable': ['geopotential'],
    'year': ['2024'],
    'month': ['03'],
    'day': ['01'],
    'time': ['13:00'],
    'pressure_level': ['1000'],
    'data_format': 'grib',
}
target = 'download.grib'
client.retrieve(dataset, request, target)
```
 Please note that there are several terms and conditions to agree to when downloading the data, so please run the example in order to agree to them. The first request will end in an error such as below
```console
requests.exceptions.HTTPError: 403 Client Error: Forbidden for url: https://cds.climate.copernicus.eu/api/retrieve/v1/processes/reanalysis-era5-pressure-levels/execution
required licences not accepted
Not all the required licences have been accepted; please visit https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels?tab=download#manage-licences to accept the required licence(s).
```
visit the listed website and accept the license(s), then retry.

### Basic python sript to download the data
```python
#!/usr/bin/python
import cdsapi
client = cdsapi.Client()
dataset = "reanalysis-era5-single-levels"
request = {"product_type": ["reanalysis"],
           "variable": ["_VARIABLE_"],
           "year": ["_YEAR_"],
           "month": [_LIST_MONTHS_],
           "day": ["01", "02", "03", "04", "05", "06", "07",
                   "08", "09", "10", "11", "12", "13", "14",
                   "15", "16", "17", "18", "19", "20", "21",
                   "22", "23", "24", "25", "26", "27", "28",
                   "29", "30", "31"
           ],
           "time": ["00:00", "01:00", "02:00", "03:00", "04:00", "05:00",
                    "06:00", "07:00", "08:00", "09:00", "10:00", "11:00",
                    "12:00", "13:00", "14:00", "15:00", "16:00", "17:00",
                    "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"
           ],
           "data_format": "grib",
           "download_format": "unarchived"
}
target = "_OUTPUT_"
client.retrieve(dataset, request, target)
```
bash file to retrieve the data
```shell
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

		out=${var}_${y}.grib
		scr=${var}_${y}.run
		sed -e "s!_VARIABLE_!$var!g" \
		    -e "s!_YEAR_!$y!g" \
		    -e "s!_LIST_MONTHS_!${lm}!g" \
		    -e "s!_OUTPUT_!$out!g" \
		run_sedf_ep > $scr

		chmod +x $scr
		./$scr || exit -1

	done
done

exit 0
```
