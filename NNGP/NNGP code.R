library(sp)
library(rgdal)
library(lattice)
library(data.table)
library(gstat)
library(tidycensus)
library(tidyverse)
library(sf)
library(RColorBrewer)
library(spNNGP)
library(geepack)
library(geoR)
library(fields)
library(rgdal)

data_merge<-readRDS('/Users/wenbozhang/Desktop/uw course/596/aclima project/aclima data/full_data_merge.rds')


# ther are na values in race
data_merge_comp<-data_merge
st_geometry(data_merge_comp) <- NULL

data_merge_comp<-data_merge[complete.cases(data_merge_comp[,1:22]),]
data_merge_comp[data_merge_comp$value<5,'value']=5/sqrt(2)
data_merge_comp$log10_value<-log10(data_merge_comp$value)
data_merge_comp$log_value<-log(data_merge_comp$value)

#data_merge_comp<-data_merge_comp[data_merge_comp$value<50,]
n<-nrow(data_merge_comp)

# transform coordinates
# utm 10
xy <- data.frame(st_coordinates(data_merge_comp$geometry_aclima))
coordinates(xy) <- c("X", "Y")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
res <- spTransform(xy, CRS("+proj=utm +zone=10 ellps=WGS84"))
coords<-coordinates(res)/1000
data_merge_comp$coords<-coords

plot(data_merge_comp$coords[,1],data_merge_comp$coords[,2], col=tim.colors()[cut(data_merge_comp$value,64)], pch=19,ylab='y(km)',xlab='x(km)')

set.seed(1)

sigma.sq <- 5
tau.sq <- 1
phi <- 3/0.5
starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)
tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)
priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))
cov.model <- "exponential"
# phi alpha
theta.alpha <- as.matrix(expand.grid(seq(0.1,1,length.out=15), seq(3/10,3/0.1,length.out=15)))
colnames(theta.alpha) <- c("alpha","phi")

start.time <- Sys.time()

model1 <- spConjNNGP(log_value~race_aa_per+race_asian_per+race_hisp_per+race_other_per+age_median+`edu_>high_per`,data=data_merge_comp,coords=coords,sigma.sq.IG=c(2, 40), n.neighbors=10,cov.model=cov.model,theta.alpha=theta.alpha,
                     k.fold = 2, score.rule = "crps",fit.rep=TRUE, n.samples=200,n.omp.threads=18,verbose=FALSE)

model2 <- spConjNNGP(log_value~race_aa_per+race_asian_per+race_hisp_per+race_other_per+age_median+`hhinc_>100000_per`,data=data_merge_comp,coords=coords,sigma.sq.IG=c(2, 40), n.neighbors=10,cov.model=cov.model,theta.alpha=theta.alpha,
                     k.fold = 2, score.rule = "crps",fit.rep=TRUE, n.samples=200,n.omp.threads=18,verbose=FALSE)

model3 <- spConjNNGP(log_value~race_aa_per+race_asian_per+race_hisp_per+race_other_per+age_median+unemp_per,data=data_merge_comp,coords=coords,sigma.sq.IG=c(2, 40), n.neighbors=10,cov.model=cov.model,theta.alpha=theta.alpha,
                     k.fold = 2, score.rule = "crps",fit.rep=TRUE, n.samples=200,n.omp.threads=18,verbose=FALSE)

end.time <- Sys.time()