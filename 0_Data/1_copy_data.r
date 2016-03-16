copy_data <- function(model, mfolder, type="hist", res_folder=paste("~/Data/gam_decomp/0_Data/", model, sep="")){
  res_folder <- paste(res_folder, type, sep="/")
  dir.create(res_folder, recursive=TRUE)
  ffolder <- switch(type,
                    hist = paste(mfolder,"/historical/mon/atmos/", sep=""),
                    rcp = paste(mfolder,"/rcp85/mon/atmos/", sep=""),
                    nat = paste(mfolder,"/historicalNat/mon/atmos/", sep="")
                    )
  l_folders <- list.files(ffolder,"r[1-9]i1p1", full.names=TRUE)
  l_files <- unlist(lapply(l_folders, list.files, pattern="tas_Amon_*", full.names=TRUE))
  l_files_small <- unlist(lapply(l_folders, list.files, pattern="tas_Amon_*", full.names=FALSE))
  lapply(l_files, file.copy, to=res_folder)
  unlist(lapply(l_files_small, function(x) paste(res_folder, x, sep="/")))
}
copy_data("CNRM","/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/CNRM-CERFACS/CNRM-CM5")

copy_all <- function(model, mfolder, res_folder=paste("~/Data/gam_decomp/0_Data/", model, sep="")){
 types <- c("hist", "nat", "rcp")
 ans <- lapply(types, copy_data, model=model, mfolder=mfolder, res_folder=res_folder)
 names(ans) <- types
 ans
}
cfiles <- copy_all("CNRM","/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/CNRM-CERFACS/CNRM-CM5")

create_ts <- function(input){
  input <- unlist(strsplit(input, ".nc"))[1]
  system(paste("cdo -yearmean -selmon,6,7,8 ", input,".nc ", input, "_smean.nc", sep=""))
  system(paste("cdo fldmean ", input,"_smean.nc ", input,"_gbl_sts.nc", sep=""))
  system(paste("cdo fldmean -sellonlatbox,-10,40,30,50, ", input,"_smean.nc ", input,"_eur_sts.nc", sep=""))
  c(paste(input, "_eur_sts.nc", sep=""), paste(input, "_gbl_sts.nc", sep=""))
}

transfer_data <- function(model, l_files, type="hist", res_folder=paste("~/Data/gam_decomp/1_PrepareData/2_FormatData/1_Data/", model, sep="")){
  res_folder <- paste(res_folder, type, sep="/")
  dir.create(res_folder, recursive=TRUE)
  lapply(l_files, file.copy, to=res_folder)
}

get_and_transfer <- function(model, mfolder, type){
  cfiles <- copy_data(model, mfolder, type)
  print(cfiles)
  cfiles <- unlist(lapply(cfiles, create_ts))
  transfer_data(model, cfiles, type)
}
get_and_transfer("CNRM","/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/CNRM-CERFACS/CNRM-CM5", "hist")

get_and_transfer_alltypes <- function(model, mfolder){
 types <- c("hist", "nat", "rcp")
 lapply(types, get_and_transfer, model=model, mfolder=mfolder)
}
get_and_transfer_alltypes("CNRM","/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/CNRM-CERFACS/CNRM-CM5")
l_models[5]="CNRM"
lapply(c(3,5,16,21), function(i) get_and_transfer_alltypes(l_models[i], l_folders[i]))
lapply(seq_along(l_models), function(i) get_and_transfer_alltypes(l_models[i], l_folders[i]))

# copy_data("GFDL", 
#           "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/GFDL/GFDL-CM3/historical/mon/atmos/",
#           "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/GFDL/GFDL-CM3/rcp85/mon/atmos/"
# )
# 
# copy_data("GFDL", 
#           "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/GFDL/GFDL-CM3/historical/mon/atmos/",
#           "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/GFDL/GFDL-CM3/rcp85/mon/atmos/"
# )
# 
l_models <- list.files("/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/")
l_folders <- c("/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/BCC/bcc-csm1-1/",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/BNU/BNU-ESM",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/CCCMA/CanESM2",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/CMCC/CMCC-CESM",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/CNRM-CERFACS/CNRM-CM5",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/CSIRO-QCCCE/CSIRO-Mk3-6-0",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/CSIRO-BOM/ACCESS1-0",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/CSIRO-QCCCE/CSIRO-Mk3-6-0",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/FIO/FIO-ESM",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/GFDL/GFDL-CM3",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/GISS/GISS-E2-R",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/IAP/FGOALS-g2",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/ICHEC/EC-EARTH",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/CMCC/CMCC-CESM",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/INM/inmcm4",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/IPSL/IPSL-CM5A-LR",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/MIROC/MIROC5",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/MOHC/HadGEM2-CC",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/MPIM/MPI-ESM-LR",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/MRI/MRI-CGCM3",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/NCAR/CCSM4/",
               "/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/NCC/NorESM1-M")
