# Med12 configuration setup

The configuration is described in Storto et. al. https://gmd.copernicus.org/articles/16/4811/2023/gmd-16-4811-2023.html

installation can be done with the following steps

```
mkdir NemoMed12
cd NemoMed12
svn co https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/r4.0/r4.0.7/
```
It is important to stick with version 4.0.7, version 4.0 does not work

## To create the configuration MED12
```
cd cfgs
cp -r AMM12 MED12
cd MED12
mv cpp_AMM12.fcm cpp_MED12.fcm
```
then remove `top` from `cpp_MED12.fcm`
copy `EXPREF` into 
compile
