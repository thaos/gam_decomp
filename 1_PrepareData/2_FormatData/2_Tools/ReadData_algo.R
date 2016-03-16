library(ncdf)

get_tas <- function(nc.file, var="tas"){
		nc <- open.ncdf(nc.file)
		tas <- as.numeric(get.var.ncdf(nc, var))
    origin <- strsplit(nc$dim$time$units, " ")[[1]][3]
    years  <-  Vectorize(function(x)substr(x, 1,4))(as.Date(get.var.ncdf(nc, "time"),  origin=origin))
    close.ncdf(nc)
    data.frame("year"=years, "tas"=tas)
}
get_tas("~/Data/gam_decomp/1_PrepareData/2_FormatData/1_Data/CNRM/hist/tas_Amon_CNRM-CM5_historical_r1i1p1_185001-200512_eur_sts.nc")
get_tas("~/Data/gam_decomp/1_PrepareData/2_FormatData/1_Data/CNRM/rcp/tas_Amon_CNRM-CM5_rcp85_r1i1p1_200601-210012_eur_sts.nc")
get_tas("~/Data/gam_decomp/1_PrepareData/2_FormatData/1_Data/IPSL/rcp/tas_Amon_IPSL-CM5A-LR_rcp85_r1i1p1_200601-230012_eur_sts.nc")

get_run <- function(nc.file){
  splitted <- strsplit(nc.file, "/")[[1]]
  splitted <- strsplit(splitted[length(splitted)],  "_")[[1]][5]
  run <- as.numeric(substr(splitted, 2, 2))
  run
}
# get_run("~/Data/gam_decomp/1_PrepareData/2_FormatData/1_Data/CNRM/rcp/tas_Amon_CNRM-CM5_rcp85_r1i1p1_200601-210012_eur_sts.nc")

get_tas_and_run <- function(nc.file, var="tas"){
  tas <- get_tas(nc.file, var=var)
  cbind(tas, "run"=get_run(nc.file))
}
# get_tas_and_run("~/Data/gam_decomp/1_PrepareData/2_FormatData/1_Data/CNRM/rcp/tas_Amon_CNRM-CM5_rcp85_r1i1p1_200601-210012_eur_sts.nc")
get_tas_and_run("~/Data/gam_decomp/1_PrepareData/2_FormatData/1_Data/IPSL/rcp/tas_Amon_IPSL-CM5A-LR_rcp85_r1i1p1_200601-230012_eur_sts.nc")

  
get_pattern_tas <- function(model="CNRM", type="hist", pattern="_eur_sts.nc"){
  folder <- list.files("../1_Data/", pattern=model, full.names=TRUE) 
  folder <- list.files(folder, pattern=type, full.names=TRUE) 
  folder <- list.files(folder, pattern=pattern, full.names=TRUE) 
  print(folder)
  if(length(folder) == 0) return()
  df <- lapply(folder, get_tas_and_run)
  df <- do.call(rbind, df)
  cbind(df, type=type)
}
# get_pattern_tas()

get_mod_tas <-  function(model="CNRM", type="hist"){
  eur_tas <- get_pattern_tas(model=model, type=type, pattern="_eur_sts.nc")
  if(is.null(eur_tas)) return()
  names(eur_tas)[2] <- "eur_tas"
  gbl_tas <- get_pattern_tas(model=model, type=type, pattern="_gbl_sts.nc")
  names(gbl_tas)[2] <- "gbl_tas"
  tas <- merge(eur_tas, gbl_tas, by=c("year", "run", "type"))
  tas
}
# get_mod_tas()

get_tas_alltypes <- function(model="CNRM"){
  types=c("hist", "nat", "rcp")
  df <- lapply(types, get_mod_tas, model=model)
  df <- do.call(rbind, df)
  df
}
# get_tas_alltypes()

create_rds <- function(model="CNRM", file=paste("tas_", tolower(model), ".rds", sep=""), ...){
  tas <- get_tas_alltypes(model=model)
  tas <- merge(tas, nat, by=c("year"))
  tas <- merge(tas, ant, by=c("year"))
  file <- paste("../4_Results/", file, sep="")
  saveRDS(tas, file=file)
  file
}
# create_rds()

create_rds_obs <- function(file="tas_obs.rds"){
  model="OBS"
  type="hist"
  folder <- list.files("../1_Data/", pattern=model, full.names=TRUE) 
  folder <- list.files(folder, pattern=type, full.names=TRUE) 
  eur_f <- list.files(folder, pattern="_eur_sts.nc", full.names=TRUE) 
  gbl_f <- list.files(folder, pattern="_gbl_sts.nc", full.names=TRUE) 
  if(length(eur_f) > 1){
    eur_tas <- do.call(rbind, lapply(eur_f, get_tas, var="temperature_anomaly"))
    gbl_tas <- do.call(rbind, lapply(gbl_f, get_tas, var="temperature_anomaly"))
  } else{
    eur_tas <-get_tas(eur_f, var="temperature_anomaly")
    gbl_tas <-get_tas(gbl_f, var="temperature_anomaly")
  }
  names(eur_tas)[2] <- "eur_tas"
  names(gbl_tas)[2] <- "gbl_tas"
  tas <- merge(eur_tas, gbl_tas, by="year")
  tas <- cbind(tas, run=1, type=type)
  tas <- merge(tas, nat, by=c("year"))
  tas <- merge(tas, ant, by=c("year"))
  file <- paste("../4_Results/", file, sep="")
  saveRDS(tas, file=file)
  file
}
# create_rds_obs()

