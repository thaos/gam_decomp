setwd("~/Data/FAR_paper/1_PrepareData/1_EBM/3_Scripts/")
source("../2_Tools/held_model_AS_algo.r")
FF_data <- read.table("../1_Data/Time_Prof_all_RF.txt", header=TRUE)
rcp85_all <- data.frame("year"= c(2000, 2005, 2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080,  2090, 2100), 
                        "ant"=c(1.723, 1.906, 2.154, 2.665 ,3.276, 3.993, 4.762, 5.539, 6.299, 7.020, 7.742, 8.388))
rcp85_all_int <- as.data.frame(approx(x=rcp85_all$year, y=rcp85_all$ant, xout=seq(2005, 2100, 1)))
names(rcp85_all_int) <- c("year", "ant")

volcanic_forcing <- FF_data$Volcanic
volcanic_forcing <- c(volcanic_forcing, numeric(length(2012:2100)))
volcanic_effect <- held_model(volcanic_forcing)
volcanic_df <- data.frame("year"=1750:2100, "volcanic_effect"=volcanic_effect[-1,"T"])
write.table(volcanic_df, file="../4_Results/volcanic_effect.txt", row.names=FALSE)

nat_forcing <- apply(subset(FF_data, year<=2005, select=c("Volcanic", "Solar")), 1, sum)
nat_forcing <- c(nat_forcing, numeric(length(2006:2100)))
nat_effect <- held_model(nat_forcing)
nat_df <- data.frame("year"=1750:2100, "nat_effect"=nat_effect[-1,"T"])
write.table(nat_df, file="../4_Results/nat_effect.txt", row.names=FALSE)

ant_forcing <- apply(subset(FF_data, year<=2005, select=!(names(FF_data) %in% c("year", "Volcanic", "Solar"))), 1, sum)
ant_forcing <- c(ant_forcing, unlist(subset(rcp85_all_int, year>=2006, select="ant")))
ant_effect <- held_model(ant_forcing)
ant_df <- data.frame("year"=1750:2100, "ant_effect"=ant_effect[-1,"T"])
write.table(ant_df, file="../4_Results/ant_effect.txt", row.names=FALSE)

aer_forcing <- apply(subset(FF_data, select=c("aerosolERF")), 1, sum)
aer_forcing <- c(aer_forcing, numeric(length(2012:2100)))
aer_effect <- held_model(aer_forcing)

ghg_forcing <- apply(subset(FF_data, select=c("CO2", "OtherWMGHG")), 1, sum)
ghg_forcing <- c(ghg_forcing, numeric(length(2012:2100)))
ghg_effect <- held_model(ghg_forcing)

all_forcing <- apply(FF_data[, -1], 1, sum)
all_forcing <- c(all_forcing, numeric(length(2012:2100)))
all_effect <- held_model(all_forcing)
all_df <- data.frame("year"=1750:2100, "all_effect"=all_effect[-1,"T"])

nat_aer_ghg_df <- data.frame("year"=1750:2100, "nat"=nat_effect[-1,"T"], "aer"=aer_effect[-1,"T"], "ghg"=ghg_effect[-1, "T"]) 
write.table(nat_aer_ghg_df, file="../4_Results/nat_aer_ghg_df.txt", row.names=FALSE)

write.table(data.frame(year=ant_df$year, nat=nat_forcing, ant=ant_forcing), file="../4_Results/antnat_forcings.txt")
