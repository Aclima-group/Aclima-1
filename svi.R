

library(dplyr)
library(spdep)
library(geepack)
library(stringr)
library(tidyverse)

setwd("/Users/carolliu/Downloads")
options(tigris_use_cache = TRUE)
svi<-read.csv("California (2).csv",header = TRUE)
options(timeout = 2000)

head(svi)
svi_alameda <- subset(svi, svi$COUNTY =="Alameda")
svi_alameda <- svi_alameda[,c("LOCATION","RPL_THEME1")]
head(svi_alameda)
names(svi_alameda )<-c('tractID','SVI')
head(svi_alameda)
svi_alameda$tractID <- gsub(".*Tract (.+) Alameda.*", "\\1", svi_alameda$tractID)
svi_alameda$tractID <- str_replace_all(svi_alameda$tractID, "[^[:alnum:]]", "")
svi_alameda$tractID <-ifelse(nchar(svi_alameda$tractID)==4,
                             paste(svi_alameda$tractID,'00'),svi_alameda$tractID)
svi_alameda$tractID <- gsub(" ", "", svi_alameda$tractID , fixed = TRUE)


## assign census tract ID to air quality dataset 

full_data <- readRDS('data_merge_median.rds')
head(full_data)

full_data$tractID <- substr(full_data$geoid, 6,11)
full_data$value<-ifelse(full_data$value< (-10), NA,full_data$value)
full_data<-na.omit(full_data)
head(full_data)


air_svi_merge <-nest_join(full_data, svi_alameda, by='tractID')
head(air_svi_merge)


## substitue <LOD values
## new_value = values after substitution 
## logno2  = log transfromed values 
air_svi_merge$new_value <- ifelse(air_svi_merge$value < 5, 5/sqrt(2), air_svi_merge$value)
head(air_svi_merge)
# saveRDS(air_svi_merge, 'air_svi_merge.RDS')

air_svi_merge$logno2 <- log(air_svi_merge$new_value)


unlist_fn <- function(x){
  return (as.numeric(unlist(x[[1]])))
}


unlist_fn_1 <- Vectorize(unlist_fn)
new_data <- unlist_fn_1(air_svi_merge$svi_alameda)
air_svi_merge <- cbind(air_svi_merge,new_data)
colnames(air_svi_merge)[27]<-'svi_score'
head(air_svi_merge)


## OLS 

lm_model <- lm(logno2~ race_aa_per+race_asian_per+race_hisp_per+race_other_per+
                 age_median+svi_score, data = air_svi_merge)
summary(lm_model)
lm_model$coefficients

## gee
gee_model <- geeglm(logno2~ race_aa_per+race_asian_per+race_hisp_per+race_other_per+
                      age_median+svi_score, data = air_svi_unique,
                    id = tractID,
                    family = gaussian,corstr = 'independence') 

summary(gee_model)
cc <- coef(summary(gee_model))
citab <- with(as.data.frame(cc),
              cbind(lwr=Estimate-1.96*Std.err,
                    upr=Estimate+1.96*Std.err))
rownames(citab) <- rownames(cc)
citab


## nngp
xy <- data.frame(st_coordinates(air_svi_merge$geometry_aclima))
coordinates(xy) <- c("X", "Y")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
res <- spTransform(xy, CRS("+proj=utm +zone=10 ellps=WGS84"))
coords<-coordinates(res)/1000



sigma.sq <- 5
tau.sq <- 1
phi <- 3/0.5
starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)
tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)
priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))
cov.model <- "exponential"
theta.alpha <- as.matrix(expand.grid(seq(0.1,1,length.out=15), seq(3/10,3/0.1,length.out=15)))
colnames(theta.alpha) <- c("alpha","phi")


# modeling
start.time <- Sys.time()
model_svi <- spConjNNGP(logno2~race_aa_per+race_asian_per+race_hisp_per+race_other_per+age_median+svi_score,data=air_svi_merge,coords=coords,sigma.sq.IG=c(2, 40), n.neighbors=10,cov.model=cov.model,theta.alpha=theta.alpha,
                        k.fold = 2, score.rule = "crps",fit.rep=TRUE, n.samples=200,n.omp.threads=18,verbose=FALSE)
end.time <- Sys.time()
summary(model_svi)
