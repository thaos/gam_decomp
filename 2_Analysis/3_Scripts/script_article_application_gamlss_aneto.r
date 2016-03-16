#!/usr/bin/env Rscript
.libPaths("/cnrm/vdr/USERS/thaos/R/x86_64-mageia-linux-gnu-library/3.0")

# load functions
packages <- c("scales",
              "gtable",
              "grid",
              "extRemes",
              "quantreg",
              "ncdf",
              "gamlss",
              "FARg")
lapply(packages, library, character.only=TRUE)

setwd("~/Data/gam_decomp/2_Analysis/3_Scripts/")
make_shortcut <- function(folder){
function(file) paste(folder, file, sep="")
}
res_path <- make_shortcut("../4_Results/")
data_path <- make_shortcut("../1_Data/")
tools_path <- make_shortcut("../2_Tools/")
source(tools_path("gamlss_allnat.r"), chdir=TRUE)
source(tools_path("gam_rmnat.r"), chdir=TRUE)

load_data <- function(model){
  tas <- readRDS(data_path(paste("tas_", model, ".rds", sep="")))
  tas_subset <- subset(tas, year >= 1961 & year <= 1990 & type != "nat")
  eur_bias <- tapply(as.numeric(tas_subset$eur_tas), tas_subset$run, mean)
  gbl_bias <- tapply(as.numeric(tas_subset$gbl_tas), tas_subset$run, mean)
  for(r in unique(tas_subset$run)){
    tas[tas$run==r, "eur_tas"] <- as.numeric(tas[ tas$run==r, "eur_tas"]) - eur_bias[r] 
    tas[tas$run==r, "gbl_tas"] <- as.numeric(tas[ tas$run==r, "gbl_tas"]) - gbl_bias[r] 
  }
  tas$hnat <- tas$type=="nat"
  tas$ant[tas$hnat] <- 0
  tas
}

compute_far <- function(model, xp=1.6, stat_model=gauss_fit, R=10, ...){
  mdata <- load_data(model)
  bsample <- boot_samples(mdata, R=R)
  gam_an <- boot_gam_allnat(bsample, "eur_tas", "ant", "nat", "year", "hnat")
  far <- boot_far_allnat_onperiod(gam_an, xp=xp)
  far <- add_param(far, operation=p1/p0, name="RR")
  far["RR",,] %<>% imput_aurel(.) 
  far %<>%  add_param(., operation=al_trans(RR), name="alRR") %>% 
  add_param(., operation=el_trans(RR), name="elRR")
  ic_far <- get_ic_onperiod(far, method_name=model, ci_p=0.9)  
  ic_far$param <- as.character(ic_far$param)
  list("data"=mdata, "samples"=bsample, "gam"=gam_an, "far"=far, "ic_far"=ic_far)
} 

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

model <- args[1]
far_model <- compute_far(model)

saveRDS(far_ipsl, file="ipsl.rds")
far_ncar <- readRDS("ncar.rds")

