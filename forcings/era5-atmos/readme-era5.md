# Atmospheric forcing: ERA5 data

ERA5 (atmospheric forcing) is downloaded from Copernicus Climate Data Store https://cds.climate.copernicus.eu/


### Climate Data Store API, [API](https://cds.climate.copernicus.eu/how-to-api)
First, [register](https://accounts.ecmwf.int/auth/realms/ecmwf/login-actions/registration?client_id=cds&tab_id=w2gu_uF6V0o) to ECMWF and then [log in](https://accounts.ecmwf.int/auth/realms/ecmwf/login-actions/authenticate?client_id=cds&tab_id=w2gu_uF6V0o). Once logged in, copy the code displayed below to the file `$HOME/.cdsapirc`
(in your Unix/Linux environment)
```
url: https://cds.climate.copernicus.eu/api
key: <PERSONAL-ACCESS-TOKEN>
```
You can now install the CDS API client via the package management system pip, by running on Unix/Linux the command below.
```
pip install "cdsapi>=0.7.2"
```
Once the CDS API client is installed, it can be used to request data from the datasets listed in the CDS, ADS, ECDS and CEMS Early Warning DS catalogues.
```Python
import cdsapi
client = cdsapi.Client()
# definition
dataset = "<DATASET-SHORT-NAME>"
request = {<SELECTION-REQUEST>}
target = "<TARGET-FILE>"
# Request
client.retrieve(dataset, request, target)
```
