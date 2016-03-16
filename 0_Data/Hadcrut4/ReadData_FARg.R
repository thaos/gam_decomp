library(ncdf)
library(memoise)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(FARg)

# FAR STOTT
#               1990-1999   2003-2012
#--------------------------------------
#RT_ALL (years) 52(14-444)  5(2.7-11)
#RT_NAT (years) >10^4       >10^3
#--------------------------------------
#P_NAT/P_ALL    <10^-3      <10^-3

l.gbl=list.files(pattern="gbl_sts")
l.eur=list.files(pattern="eur_sts")

getData  <- function(l.nc, origin="1850-01-01"){
	require(date)
	getTAS <- function(nc.file){
		nc=open.ncdf(nc.file)
		get.var.ncdf(nc, "temperature_anomaly")
	}
	nc=open.ncdf(l.nc[1])
	years = Vectorize(function(x)substr(x, 1,4))(as.Date(get.var.ncdf(nc, "time"),  origin="1850-01-01"))
	years = as.numeric(years)
	close.ncdf(nc)
	tas=sapply(l.nc, getTAS)
	colnames(tas)=paste("run", 1:length(l.nc), sep="")
	rownames(tas)=years
	tas
}
d.gbl=getData(l.gbl)
d.gbl.v=as.numeric(d.gbl)
d.gbl.avg=matrix(rowMeans(d.gbl), nrow=nrow(d.gbl), ncol=ncol(d.gbl))
d.gbl.avg.v=as.numeric(d.gbl.avg)
d.eur=getData(l.eur)
d.eur.v=as.numeric(d.eur)
year=as.numeric(rep(rownames(d.eur), ncol(d.eur)))
par(mfrow=c(2, 2))
plot(x=as.numeric(rep(rownames(d.eur), ncol(d.eur))), y=d.eur.v)
plot(x=as.numeric(rep(rownames(d.eur), ncol(d.eur))), y=d.eur.v, xlim=c(1850, 2010))
plot(x=as.numeric(rep(rownames(d.gbl), ncol(d.gbl))), y=d.gbl.v)
plot(x=as.numeric(rep(rownames(d.gbl), ncol(d.gbl))), y=d.gbl.v, xlim=c(1850, 2010))
cor(d.gbl.avg.v, d.eur.v)
plot(d.gbl.avg.v, d.eur.v)
lm.fitted=lm(d.gbl.avg.v~d.eur.v)
par(mfrow=c(2, 2))
plot(lm.fitted)
hist(residuals(lm.fitted))
tas  <- as.data.frame(cbind(year,  d.eur.v,  d.gbl.avg.v))
names(tas) <- c("year", "eur.tas", "avg.gbl.tas") 
tas  <- as.data.frame(cbind(year,  d.eur.v,  d.gbl.v))
names(tas) <- c("year", "eur.tas", "gbl.tas") 
save(tas, file="tas.rdata")

