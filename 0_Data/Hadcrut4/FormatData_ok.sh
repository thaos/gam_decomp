#bin/sh
input="HadCRUT4_median"
cdo -yearmean -selmon,6,7,8 ${input}.nc ${input}_smean.nc
cdo fldmean ${input}_smean.nc ${input}_gbl_sts.nc
cdo sellonlatbox,-10,40,30,50, ${input}_smean.nc ${input}_eur_smean.nc
cdo fldmean -sellonlatbox,-10,40,30,50, ${input}_smean.nc ${input}_eur_sts.nc



