packages <- c("devtools", 
              "magrittr", 
              "mgcv",
              "ggplot2",
              "parallel")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())),repos="http://cran.us.r-project.org")
}
lapply(packages, library, character.only=TRUE)
load_all("~/Data/FARpackage/FARg")

gam_allnat <- function(y, ant, nat, time, hnat, knots=NULL, fixed=FALSE ){
  data <- data.frame("y"=y, "ant"=ant, "nat"=nat, "hnat"=as.factor(hnat), "time"=time)
  if(!is.null(knots))
    gam_fit <- gam(y ~ s(ant, k=knots, fx=fixed) + nat, data=data, control=list(keepData=TRUE))
  else
    gam_fit <- gam(y ~ s(ant) + nat, data=data, control=list(keepData=TRUE))
  gam_fit
}


change_anomalie_period <- function(data, y="y", time="year", period=c(1850, 1879)){
  i_pre_ind <- which(data[, time] >= period[1] & data[,time ] <= period[2])
  bias <- mean(data[i_pre_ind, y])
  data[, y] %<>% `-`(., bias)
  attr(data, "bias") <- bias
  data
}

predict_gam_allnat <- function(gam_fit, newdata=NULL){
  gam_terms <- attr(gam_fit$terms, "term.labels")
  gam_pterms <- attr(gam_fit$pterms, "term.labels")
  fix_terms <- intersect(gam_terms, gam_pterms)
  smooth_terms <- setdiff(gam_terms, gam_pterms)
  if(is.null(newdata))
    newdata <- gam_fit$data
  gam_pred <- predict(gam_fit, newdata=newdata, type="terms")
  gam_int <- attr(gam_pred, "constant")
  gam_param <- apply(gam_pred[, fix_terms, drop=FALSE], 1, sum) 
  gam_smooth <- apply(gam_pred[, paste("s(", smooth_terms, ")", sep=""), drop=FALSE], 1, sum)
  cbind(newdata, "gam_all"=gam_int + gam_param + gam_smooth, "gam_nat"=gam_param, "gam_ant"=gam_int+gam_smooth)
}

predict_gam_allnat <- function(gam_fit, newdata=NULL){
  if(is.null(newdata))
    newdata <- gam_fit$data
  data_all <- subset(newdata, hnat==0)
  data_ant <- data_all
  data_ant$nat <- 0
  data_nat <- data_all
  data_nat$ant <- 0
  data_nat$hnat <- 1
  gam_all <- predict(gam_fit, newdata=data_all)
  gam_ant <- predict(gam_fit, newdata=data_ant)
  gam_nat <- predict(gam_fit, newdata=data_nat)	
  ans <- cbind(data_all, "gam_all"=gam_all, "gam_nat"=gam_nat, "gam_ant"=gam_ant)
  data_all <- subset(newdata, hnat==1)
  data_ant <- data_all
  data_ant$nat <- 0
  gam_all <- predict(gam_fit, newdata=data_all)
  gam_ant <- predict(gam_fit, newdata=data_ant)
  ans <- rbind(ans,
               cbind(data_all, "gam_all"=gam_all, "gam_nat"=gam_all, "gam_ant"=gam_ant)
               )
}


boot_samples <- function(data, R=250){
  i_samples <- sapply(1:R, function(x) sample.int(nrow(data), size=nrow(data), replace=TRUE))
  i_samples[, 1] <- 1:nrow(data)
  l_samples <- apply(i_samples, 2, function(x) data[x,])
}
# l_samples <- boot_samples(tas_nat)


boot_gam_allnat <- function(l_samples, y, ant, nat, time, hnat, ...){
  original_data <- l_samples[[1]]
  l_gam <- lapply(l_samples, function(data) gam_allnat(data[, y], data[, ant],data[, nat], data[, time], data[, hnat], fixed=FALSE, knots=NULL))
  l_gam_an  <- lapply(l_gam, predict_gam_allnat)	
  l_gam_an_origin  <- lapply(l_gam, predict_gam_allnat, newdata=l_gam[[1]]$data)
  # print_bias <- function(original_data) print(with(subset(original_data, time >= pre_ind[1] & time <= pre_ind[2]), mean(gam_all - gam_nat)))
  # lapply(l_gam_an_origin, print_bias)
  list("l_gam"=l_gam, "l_gam_an_boot"=l_gam_an, "l_gam_an_origin"=l_gam_an_origin)
}
# l_gam_an <- boot_gam_allnat(l_samples, "eur_tas", "year", "nat", "year")


get_p_allnat_onperiod <- function(object, xp, l_time, original_data, ant_nat_or_all="all"){
  if(ant_nat_or_all == "nat")
    original_data$gam_all <- original_data$gam_nat
  if(ant_nat_or_all == "ant")
    original_data$gam_all <- original_data$gam_ant
  get_p_allnat <- function(time){
    pnt <- set_pnt(time, xp, time_var=object$time_var, data=original_data)
    if(class(object) == "gpd_fit")
      get_p(object, pnt, under_threshold=TRUE)
    else
      get_p(object, pnt) 
  }
  print(system.time(
  ans <- sapply(l_time, get_p_allnat)
  ))
  colnames(ans) <- l_time
  ans
}
# l_p_nat <- mcmapply(get_p_allnat_onperiod, object=l_gpd_fit, original_data=l_gam_an_origin, MoreArgs=list("xp"=1.6, "l_time"=1850:2100, "nat_or_all"="nat"), SIMPLIFY=FALSE)
# l_p_all <- mcmapply(get_p_allnat_onperiod, object=l_gpd_fit, original_data=l_gam_an_origin, MoreArgs=list("xp"=1.6, "l_time"=1850:2100), SIMPLIFY=FALSE)

get_q_allnat_onperiod <- function(object, p, l_time, original_data, ant_nat_or_all="all"){
  if(ant_nat_or_all == "nat")
    original_data$gam_all <- original_data$gam_nat
  if(ant_nat_or_all == "ant")
    original_data$gam_all <- original_data$gam_ant
  get_q_allnat <- function(time){
    pnt <- set_pnt(time, p, time_var=object$time_var, data=original_data)
    get_q(object, pnt) 
  }
  print(system.time(
  ans <- sapply(l_time, get_q_allnat)
  ))
  colnames(ans) <- l_time
  ans
}

get_far_allnat_onperiod <- function(object, xp, l_time, original_data){
  p_all <- get_p_allnat_onperiod(object, xp, l_time, original_data)
  p_nat <- get_p_allnat_onperiod(object, xp, l_time, original_data, ant_nat_or_all="nat")
  far <- 1 - p_nat[1,]/p_all[1,]
  ans <- rbind(far, p_nat, p_all)
  rownames(ans) <- c("FAR", paste(rownames(p_nat),"0", sep=""), paste(rownames(p_all), "1", sep=""))
  ans
}

fit_and_boot_allnat <- function(l_gam_an, fit, ...){
  fit0 <- fit(y, l_gam_an[[1]], ~gam_all, ~gam_all, "time", ...)
  l_fit <- lapply(l_gam_an, function(data) fit(y, data, ~gam_all, ~gam_all, "time", init=fit0$par, ...))
}
# qthreshold <- 0.9
# l_gpd_fit <- fit_and_boot_allnat(l_gam_an$l_gam_an_boot, gpd_fit, qthreshold=qthreshold)

boot_far_allnat_onperiod  <- function(l_fit, l_gam_an_origin, l_time, xp){
  l_gam_an_origin <- lapply(l_gam_an_origin, function(x) x[x$hnat==0,])
  l_far <- mapply(get_far_allnat_onperiod, object=l_fit, original_data=l_gam_an_origin, MoreArgs=list("xp"=xp, "l_time"=l_time), SIMPLIFY=FALSE)
  l_gam_an_origin <- lapply(l_gam_an_origin, function(x) unique(x[x$time %in% l_time, c("time", "gam_all", "gam_nat", "gam_ant")]))
  rbinding <- function(x, y){
    print(str(x))
    y <- t(y[,c("gam_all", "gam_nat", "gam_ant")]) 
    print(str(y))
    rbind(x, y)
  }
  l_far2 <- mapply(rbinding, l_far, l_gam_an_origin, SIMPLIFY=FALSE)
  a_far <- array(unlist(l_far2), dim=c(dim(l_far2[[1]]), length(l_far)))
  dimnames(a_far)[1:2] <- dimnames(l_far2[[1]])
  dimnames(a_far)[[3]] <- 1:length(l_far)
  a_far
}
# debug(boot_far_allnat_onperiod)

boot_q_allnat_onperiod  <- function(l_fit, l_gam_an_origin, l_time, p, ant_nat_or_all="all"){
  l_q <- mapply(get_q_allnat_onperiod, object=l_fit, original_data=l_gam_an_origin, MoreArgs=list("p"=p, "l_time"=l_time, "ant_nat_or_all"=ant_nat_or_all), SIMPLIFY=FALSE)
  l_gam_an_origin <- lapply(l_gam_an_origin, function(x) unique(x[x$time %in% l_time, c("time", "gam_all", "gam_nat", "gam_ant")]))
  rbinding <- function(x, y){
    print(str(x))
    y <- t(y[,c("gam_all", "gam_nat", "gam_ant")]) 
    print(str(y))
    rbind(x, y)
  }
  l_q2 <- mapply(rbinding, l_q, l_gam_an_origin, SIMPLIFY=FALSE)
  a_q <- array(unlist(l_q2), dim=c(dim(l_q2[[1]]), length(l_q)))
  dimnames(a_q)[1:2] <- dimnames(l_q2[[1]])
  dimnames(a_q)[[3]] <- 1:length(l_q)
  a_q
}

imput_aurel_byyear <- function(RR){
  l_na <- which(is.na(RR))
  n_na <- length(l_na)
  print(n_na/length(RR))
  n_inf <- floor(n_na/2)
  n_sup <- n_na - n_inf
  if(n_inf > 0) {
    i_inf <- l_na[1:n_inf]
    RR[i_inf]  <-  0
  }
  if(n_sup > 0) {
    i_sup <- l_na[(n_inf+1) : n_na]
    RR[i_sup]  <-  Inf
  }
  RR
}

imput_aurel <- function(boot_res_RR){
  RR <- apply(boot_res_RR, 1, imput_aurel_byyear)
  t(RR)
}

# a_far <- boot_far_allnat_onperiod(l_gpd_fit, l_gam_an$l_gam_an_origin, 1850:2100, 1.6)
