## --------------------------------------------------------------------------------------------
library(geepack)



## --------------------------------------------------------------------------------------------
## create census tract id:  1-11 digit of GEOID
data_merge$tract_id<-str_sub(data_merge$geoid,1,11)

## value below detection limit (replace with 5/sqrt(2))
data_merge$value_lod<-with(data_merge,ifelse(value>5,value,5/sqrt(2)))

## log-transformed NO2 
data_merge$log_no2<-log(data_merge$value_lod)

## missing value from ACS data
data_gee1<-na.omit(subset(data_merge,
                      select=c(race_aa_per,
                               race_asian_per,
                               race_hisp_per,
                               race_other_per,
                               age_median,
                               `edu_>high_per`,
                               hhinc_median,  ## use continuous hhinc
                               unemp_per,
                               log_no2,
                              tract_id)))

data_gee2<-na.omit(subset(data_merge,
                      select=c(race_aa_per,
                               race_asian_per,
                               race_hisp_per,
                               race_other_per,
                               age_median,
                               `edu_>high_per`,
                               `hhinc_>100000_per`, ## use categorical hhinc 
                               unemp_per,
                               log_no2,
                              tract_id)))

## incomplete data percentage
1-nrow(data_gee1)/nrow(data_merge) ##2.67%
1-nrow(data_gee2)/nrow(data_merge) ## 1.26%  categorical hhinc has less missing value




## --------------------------------------------------------------------------------------------

 model_lm<-lm(log_no2~race_aa_per+race_asian_per+race_hisp_per+race_other_per+
                          age_median+
                          `edu_>high_per`+
                          `hhinc_>100000_per`+
                          unemp_per,
            data=data_gee2,
                    )

 
 summary(model_lm)




## --------------------------------------------------------------------------------------------

## Choose categorical income variable (with less missing value)
model_full<- geeglm(log_no2~
                          race_aa_per+race_asian_per+race_hisp_per+race_other_per+
                          age_median+
                          `edu_>high_per`+`hhinc_>100000_per`+unemp_per,
                    id=tract_id,
                    data=data_gee2,
                    family=gaussian,
                    corstr="ind")



## --------------------------------------------------------------------------------------------
model_edu<-geeglm(log_no2~
                          race_aa_per+race_asian_per+race_hisp_per+race_other_per+
                          age_median+
                          `edu_>high_per`,
                    id=tract_id,
                    data=data_gee2,
                    family=gaussian,
                    corstr="ind")

model_hinc<-geeglm(log_no2~
                          race_aa_per+race_asian_per+race_hisp_per+race_other_per+
                          age_median+
                          `hhinc_>100000_per`,
                    id=tract_id,
                    data=data_gee2,
                    family=gaussian,
                    corstr="ind")

model_unemp<-geeglm(log_no2~
                          race_aa_per+race_asian_per+race_hisp_per+race_other_per+
                          age_median+
                          unemp_per,
                    id=tract_id,
                    data=data_gee2,
                    family=gaussian,
                    corstr="ind")
 
summary(model_edu)
summary(model_hinc)
summary(model_unemp)

