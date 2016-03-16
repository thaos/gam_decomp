setwd("~/Data/gam_decomp/1_PrepareData/2_FormatData/3_Scripts/")
source("../2_Tools/ReadData_algo.R")

res_folder <- "../4_Results/"
nat <- read.table("../../1_EBM/4_Results/nat_effect.txt",
                  header=TRUE)
ant <- read.table("../../1_EBM/4_Results/ant_effect.txt",
                  header=TRUE)
names(nat)[2] <- "nat"
names(ant)[2] <- "ant"
#l_models <- list.files("../1_Data")[-c(1, 10)]
#l_models <- unlist(lapply(strsplit(l_models, "_"), "[", 2))
# l_models <- c("CNRM", "IPSL", "MOHC", "NCAR")

tas_cnrm <- create_rds(model="CNRM", sim="HIST")
tas_cnrm <- readRDS(tas_cnrm)
tas_nat_cnrm <- merge(tas_cnrm, nat, by=c("year"))
tas_nat_cnrm <- merge(tas_nat_cnrm, ant, by=c("year"))
saveRDS(tas_nat_cnrm, file=paste(res_folder, "tas_histnat_cnrm.rds", sep=""))
tas_rcp <- create_rds(model="CNRM/NAT")
tas_rcp <- readRDS(tas_rcp)
tas_nat_rcp <- merge(tas_rcp, nat, by=c("year"))
tas_nat_rcp <- merge(tas_nat_rcp, ant, by=c("year"))
saveRDS(tas_nat_rcp, file=paste(res_folder, "tas_nat_cnrm.rds", sep=""))
# 
# tas_ipsl <- create_rds(model="IPSL")
# tas_ipsl <- readRDS(tas_ipsl)
# tas_nat_ipsl <- merge(tas_ipsl, nat, by=c("year"))
# tas_nat_ipsl <- merge(tas_nat_ipsl, ant, by=c("year"))
# saveRDS(tas_nat_ipsl, file=paste(res_folder, "tas_nat_ipsl.rds", sep=""))

tas_obs <- create_rds(model="OBS", var="temperature_anomaly")
tas_obs <- readRDS(tas_obs)
tas_nat_obs <- merge(tas_obs, nat, by=c("year"))
tas_nat_obs <- merge(tas_nat_obs, ant, by=c("year"))
tas_nat_obs$run <- 1
saveRDS(tas_nat_obs, file=paste(res_folder, "tas_nat_obs.rds", sep=""))

 l_models <- c("CNRM", "IPSL", "CCCMA", "NCAR")
 lapply(l_models, create_rds)

 create_rds_obs()

 mapply(format_rds, model=l_models, output=l_ouputs_nat, sim="NAT")
