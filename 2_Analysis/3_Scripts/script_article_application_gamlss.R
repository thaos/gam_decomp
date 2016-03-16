# load functions
packages <- c("scales",
              "gtable",
              "grid",
              "extRemes",
              "quantreg",
              "ncdf",
              "gamlss")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages,
                           rownames(installed.packages())),
                   repos="http://cran.us.r-project.org")
}
lapply(packages, library, character.only=TRUE)
library(devtools)
devtools::load_all("~/Data/FARpackage/FARg")

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


# l_models <- tolower(list.files("../../1_PrepareData/2_FormatData/1_Data")[-c(1, 5, 3, 10)])
# l_models <- unlist(lapply(strsplit(l_models, "_"), "[", 2))
l_models <- c("cnrm", "obs")
l_models <- c("cnrm", "ipsl", "cccma", "ncar", "obs")
l_models <- c("cnrm", "ipsl", "cccma", "ncar")
# l_models <- c(l_models, "obs")
l_models <- tolower(list.files("~/Data/gam_decomp/1_PrepareData/2_FormatData/1_Data")[-c(1, 8, 9, 14)])
tas_cnrm <- load_data("cnrm")
tas_obs <- load_data("obs")
tas_ipsl <- load_data("ipsl")
tas_cccma <- load_data("cccma")
l_tas <- lapply(l_models, load_data)

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
far_cnrm <- compute_far("cnrm")
far_ipsl <- compute_far("ipsl", R=50)
far_obs <- compute_far("obs", R=2)
saveRDS(far_cnrm, file="cnrm.rds")
saveRDS(far_ipsl, file="ipsl.rds")
far_ncar <- readRDS("ncar.rds")

l_far <- character(length(l_models))
for(i in tail(seq_along(l_models), n=15)){
  cat("------------------------------------\n", l_models[i])
  ans <- compute_far(l_models[i], R=100) 
  l_far[i] <- save_fast(ans, file=paste(l_models[i],".rds", sep="")) 
}
l_far_file <- save_fast(l_far)
l_far <- readRDS("l_far.rds")

df_tas_nat <- do.call(rbind, lapply(l_models, function(x)cbind(readRDS(paste(x, ".rds", sep=""))$data, method=x)))
ic_allnat_far <- do.call(rbind, lapply(l_models, function(x)readRDS(paste(x, ".rds", sep=""))$ic_far))
ic_allnat_far$param <- as.character(ic_allnat_far$param)
ic_allnat_far$param[ic_allnat_far$param == "gam_ant"] <- "x_t,ant" 
ic_allnat_far$param[ic_allnat_far$param == "gam_nat"] <- "x_t,nat"
ic_allnat_far$param[ic_allnat_far$param == "gam_all"] <- "x_t" 
p1_obs <- as.numeric(subset(ic_allnat_far, time==2003 & method=="obs" & param == "p1", select="Estim"))
ic_far_merged <- merge(ic_allnat_far, df_tas_nat, by.x=c("time", "method"), by.y=c("year", "method"))

compute_q <- function(model, p=p1_obs, stat_model=gauss_fit, ...){
  mdata <- load_data(model)
  mdata_pi <- change_anomalie_period(mdata, y="eur_tas")
  bsample <- boot_samples(mdata, R=10)
  gam_an <- boot_gam_allnat(bsample, "eur_tas", "ant", "nat", "year", reuse_ant_bias=TRUE)
  mfit <- fit_and_boot_allnat(gam_an$l_gam_an_boot, stat_model, ...)
  q <- boot_q_allnat_onperiod(mfit, gam_an$l_gam_an_origin, l_time=sort(unique(gam_an$l_gam_an_origin[[1]]$time)), p=p)
  ic_q <- get_ic_onperiod(q, method_name=model, ci_p=0.9)  
  ic_q$param <- as.character(ic_q$param)
  list("data"=mdata, "data_pi"=mdata_pi, "samples"=bsample, "gam"=gam_an, "fit"=mfit, "q"=q, "ic_q"=ic_q)
} 
q_cnrm <- compute_q("cnrm")
plot_boot_time(q_cnrm$ic_q, param="q")  
l_q <- character(length(l_models))
for(i in seq_along(l_models)){
  cat("------------------------------------\n", l_models[i])
  ans <- compute_q(l_models[i]) 
  l_q[i] <- save_fast(ans, file=paste(l_models[i],"q_.rds", sep="")) 
}
l_q_file <- save_fast(l_q)
l_q <- readRDS("l_q.rds")
ic_allnat_q <- do.call(rbind, lapply(l_models, function(x)readRDS(paste(x, "q_.rds", sep=""))$ic_q))
ic_allnat_q$param <- as.character(ic_allnat_q$param)
ic_allnat_q$param[ic_allnat_q$param == "gam_ant"] <- "x_t,ant" 
ic_allnat_q$param[ic_allnat_q$param == "gam_nat"] <- "x_t,nat"
ic_allnat_q$param[ic_allnat_q$param == "gam_all"] <- "x_t" 
plot_boot_time(ic_allnat_q, param="q")

p_gg <- ggplot(l_gam_an$l_gam_an_origin, aes(x=ant, y=y)) +
geom_point(alpha=0.3) +
geom_line(aes(y=gam_all), color="red", size=2)+
geom_line(aes(y=gam_ant), color="green", size=1)+
xlab("ALL forcing response") +
ylab("T_reg") +
ggtitle("GAM fit : T_reg=f(ALL)")+
facet_grid(clim_mod~.)
p_gg
dev.print(pdf, file="gam_ant_resp.pdf")

p_gg <- ggplot(df_tas_nat, aes(x=time, y=ant)) +
geom_point(alpha=0.3) +
geom_line(aes(y=ant), color="red", size=2)+
geom_line(aes(y=ant), color="red", size=2)+
ylab("ALL forcing resp") +
xlab("year") +
ggtitle("C02 forcings")
p_gg
dev.print(pdf, file="ant_resp.pdf")

# p0 <- ggplot(data=subset(ic_allnat_far, param %in% c("x_t", "x_t,nat", "x_t,ant") & method %in% c("cnrm", "obs")), 
p0 <- ggplot(data=subset(ic_allnat_far, param %in% c("x_t")), 
             aes(x=time))+
# ggtitle(expression(paste("Covariate Estimation and Decomposition: ", x[t], " = ", x["t, ant"], " + ",  x["t, nat"]))) +
ggtitle(expression(paste("Covariate Estimation : ", x[t], " = ", x["t, ant"], " + ",  x["t, nat"]))) +
# geom_point(data=subset(df_tas_nat, method %in% c("cnrm", "obs")), aes(x=year, y=eur_tas), size=0.9, alpha=0.7) +
# geom_point(data=df_tas_pi, aes(x=year, y=eur_tas), color="yellow",size=0.9, alpha=0.3) +
geom_point(data=df_tas_nat, aes(x=year, y=eur_tas), color="black",size=0.9, alpha=0.3) +
geom_line(aes(x=time, y=Estim, group=param, color=param), size=1) +
geom_ribbon(aes(x=time, ymin=IC_inf, ymax=IC_sup, group=param, color=param), alpha=0.4) + 
# facet_grid(param~method) +
facet_wrap(~method) +
coord_cartesian(xlim=c(1850,2100))+
theme(legend.position="none")+
# geom_hline(aes(yintercept=1.6), linetype=2)+
ylab("Estimate")
# pdf(file=res_path("gam_decomp_app__obs.pdf"))
pdf(file=res_path("gam_decomp_app__obs.pdf"), width=12, height=10)
plot(p0)
dev.off()

p00 <- ggplot(data=subset(ic_far_merged, param %in% c("x_t")), 
             aes(x=ant))+
# ggtitle(expression(paste("Covariate Estimation and Decomposition: ", x[t], " = ", x["t, ant"], " + ",  x["t, nat"]))) +
ggtitle(expression(paste("Covariate Estimation : ", x[t], " = ", x["t, ant"], " + ",  x["t, nat"]))) +
# geom_point(data=subset(df_tas_nat, method %in% c("cnrm", "obs")), aes(x=ant, y=eur_tas), size=0.9, alpha=0.7) +
# geom_point(data=df_tas_pi, aes(x=ant, y=eur_tas), color="yellow",size=0.9, alpha=0.3) +
geom_point(data=df_tas_nat, aes(x=ant, y=eur_tas), color="black",size=0.9, alpha=0.3) +
geom_line(aes(x=ant, y=Estim, group=param, color=param), size=1) +
geom_ribbon(aes(x=ant, ymin=IC_inf, ymax=IC_sup, group=param, color=param), alpha=0.4) + 
# facet_grid(param~method) +
facet_wrap(~method) +
coord_cartesian(xlim=c(0, 2))+
theme(legend.position="none")+
# geom_hline(aes(yintercept=1.6), linetype=2)+
ylab("Estimate")
# pdf(file=res_path("gam_decomp_app__obs.pdf"))
plot(p00)

col=c(gg_color_hue(length(l_models)-1), "black")
# ic_subset <- ic_allnat_far[ic_allnat_far$param %in% c("p0", "p1") & ic_allnat_far$method=="ipsl",] %T>% print()
ic_subset <- ic_allnat_far[ic_allnat_far$param %in% c("p0", "p1"), ] %T>% print()
# p0001 <- ggplot(subset(ic_subset, method=="cnrm"),  aes(x=time))+
p0001 <- ggplot(ic_subset,  aes(x=time))+
ggtitle(expression(paste(p[0], " and ", p[1], " from 1870 to 2070")))+ 
geom_point(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+
geom_line(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+ 
geom_ribbon(aes(x=time,  ymin=IC_inf,  ymax=IC_sup,  group=method,  fill=method),  alpha=0.2)+ 
coord_cartesian(xlim=c(1850,2100))+
# facet_grid(param~., scales="free_y") +
facet_grid(param~. ) +
ylab("Estimate")+
coord_trans(y="log10")+
# ylim(c(0, 0.000002))+
scale_color_manual(values=col)+
scale_fill_manual(values=col)+
theme(legend.position = "bottom")
pdf(file=res_path("p0p1_app_obs.pdf"), width=12, height=10)
# pdf(file="p0p1_app_poster.pdf", height=15)
plot(p0001)
dev.off()

sample(l_models, 10)
plot_far(ic_allnat_far, axis_trans="al", main=expression(paste("Relative Risk = ", p[1]/p[0], " from 1850 to 2100" )), xlim=c(1850,2100), col=c(gg_color_hue(length(l_models)-1), "black"))
plot_far(ic_allnat_far, axis_trans="al", main=expression(paste("Relative Risk = ", p[1]/p[0], " from 1850 to 2100" )), xlim=c(1850,2100), col=c(gg_color_hue(length(l_models))))
dev.print(pdf, file=res_path("far_ic_app_obs.pdf"), height=10, width=12)
plot_far(subset(ic_allnat_far, method %in% sample(l_models, 5)), axis_trans="al", main=expression(paste("Relative Risk = ", p[1]/p[0], " from 1850 to 2100" )), xlim=c(1850,2100), col=c(gg_color_hue(5)))
dev.print(pdf, file=res_path("far_ic_app_cnrm.pdf"), height=9, width=12)
plot_far(subset(ic_allnat_far, method=="gauss_fit"), axis_trans="al", main=expression(paste("Relative Risk = ", p[1]/p[0], " from 1850 to 2100" )), xlim=c(1850,2100), col=c(gg_color_hue(2)))
dev.print(pdf, file="far_ic_app_poster.pdf", width=12, height=10)

