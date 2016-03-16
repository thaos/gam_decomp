curf <- getwd()
l_models <- list.files("/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/")

for (f in l_models){
  print(f)
  setwd(f)
  list_inputs <- list.files("./", pattern="r[1-9]i1p1.nc")
  list_inputs <- unlist(lapply(list_inputs, function(x) unlist(strsplit(x, split=".nc"))[1]))
  for (input in list_inputs)
    system(paste("cdo timmean -yearmean -seldate,1961-01-01T00:00:00,1990-12-31T:23:59:59  -selmon,6,7,8 ", input,"_ok.nc ", input,"_clim.nc", sep=""))
  system("cdo ensmean r[1-9]i1p1_clim.nc rxi1p1_clim.nc")
  for (input in list_inputs){
    system(paste("cdo sub -yearmean -selmon,6,7,8 ", input,"_ok.nc ", input,"_clim.nc ", input,"_smean.nc", sep=""))
    system(paste("cdo fldmean ", input,"_smean.nc ", input,"_gbl_sts.nc", sep=""))
    system(paste("cdo sellonlatbox,-10,40,30,50, ", input,"_smean.nc ", input,"_eur_smean.nc", sep=""))
    system(paste("cdo fldmean -sellonlatbox,-10,40,30,50, ",input,"_smean.nc ", input,"_eur_sts.nc", sep=""))
  }
  setwd(curf)
}



