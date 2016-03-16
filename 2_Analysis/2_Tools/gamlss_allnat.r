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

gam_allnat <- function(y, ant, nat, time, hnat, start.from=NULL){
  data <- data.frame("y"=y, "ant"=ant, "nat"=nat, "hnat"=hnat, "time"=time)
  if(is.null(start.from)){
    gam_fit <-gamlss(y~pb(ant)+nat ,sigma.fo=~pb(ant)+nat,family=NO, data=data, method=mixed(1,50))
  }else{
    gam_fit <- tryCatch({ 
      gamlss(y~pb(ant)+nat ,sigma.fo=~pb(ant)+nat,family=NO, data=data, method=mixed(1,50), start.from=start.from)
    },
    error = function(e){print(e); start.from}, 
    finally = print("-------"))

  }
  gam_fit$data <- data
  gam_fit
}

boot_samples <- function(data, R=250){
  i_samples <- sapply(1:R, function(x) sample.int(nrow(data), size=nrow(data), replace=TRUE))
  i_samples[, 1] <- 1:nrow(data)
  l_samples <- apply(i_samples, 2, function(x) data[x,])
}
# l_samples <- boot_samples(tas_nat)


boot_gam_allnat <- function(l_samples, y, ant, nat, time, hnat, ...){
  data_0 <- l_samples[[1]]
  start.from <- gam_allnat(data_0[, y], data_0[, ant],data_0[, nat], data_0[, time], data_0[, hnat])
  l_gam <- lapply(l_samples, function(data) gam_allnat(data[, y], data[, ant],data[, nat], data[, time], data[, hnat], start.from=start.from))
  l_gam
}
# l_gam_an <- boot_gam_allnat(l_samples, "eur_tas", "year", "nat", "year")


get_p_allnat_onperiod <- function(object, xp, data, original_data, ant_nat_or_all="all"){
  original_data <- subset(original_data, hnat==0, select=c("time", "ant", "nat"))
  original_data <- unique(original_data)
  original_data <- original_data[order(original_data$time), ]
  if(ant_nat_or_all == "nat")
    original_data$ant <- 0
  if(ant_nat_or_all == "ant")
    original_data$nat <- 0
  gamlss_p  <- predictAll(object, type="response", newdata=original_data, data=data)
  ans <- pNO(xp, gamlss_p$mu, gamlss_p$sigma, lower.tail=FALSE)
  ans <- t(cbind(ans, gamlss_p$mu, gamlss_p$sigma))
  rownames(ans) <- c("p", "mu", "sigma")
  colnames(ans) <- original_data$time
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

get_far_allnat_onperiod <- function(object, xp, original_data){
  p_all <- get_p_allnat_onperiod(object, xp, object$data, original_data)
  p_nat <- get_p_allnat_onperiod(object, xp, object$data, original_data, ant_nat_or_all="nat")
  far <- 1 - p_nat[1,]/p_all[1,]
  ans <- rbind(far, p_nat, p_all)
  rownames(ans) <- c("FAR", paste(rownames(p_nat),"0", sep=""), paste(rownames(p_all), "1", sep=""))
  ans
}

boot_far_allnat_onperiod  <- function(l_gam, xp){
  l_far <- lapply(l_gam, get_far_allnat_onperiod,  xp=xp, original_data=l_gam[[1]]$data)
  a_far <- array(unlist(l_far), dim=c(dim(l_far[[1]]), length(l_far)))
  dimnames(a_far)[1:2] <- dimnames(l_far[[1]])
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
