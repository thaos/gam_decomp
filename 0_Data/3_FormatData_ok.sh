#bin/sh
list_inputs=("r1i1p1" "r2i1p1" "r3i1p1" "r4i1p1")

for input in ${list_inputs[*]}
do
#input="hircp1"
cdo timmean -yearmean -seldate,1961-01-01T00:00:00,1990-12-31T:23:59:59  -selmon,6,7,8 ${input}_ok.nc ${input}_clim.nc
done
cdo ensmean r[1-9]i1p1_clim.nc rxi1p1_clim.nc

for input in ${list_inputs[*]}
do
cdo sub -yearmean -selmon,6,7,8 ${input}_ok.nc ${input}_clim.nc ${input}_smean.nc
cdo fldmean ${input}_smean.nc ${input}_gbl_sts.nc
cdo sellonlatbox,-10,40,30,50, ${input}_smean.nc ${input}_eur_smean.nc
cdo fldmean -sellonlatbox,-10,40,30,50, ${input}_smean.nc ${input}_eur_sts.nc
done


